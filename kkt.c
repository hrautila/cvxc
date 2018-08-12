
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include "convex.h"

// Solution of KKT equations by a dense LDL factorization of the 
// 3 x 3 system.
//    
// Returns a solver that (1) can computes the LDL factorization of
//    
// [ H           A'   GG'*W^{-1} ] 
// [ A           0    0          ],
// [ W^{-T}*GG   0   -I          ] 
//    
// given H, Df, W, where GG' = [Df; G], and (2) compute solution for 
//    
//  [ H     A'   GG'   ]   [ ux ]   [ bx ]
//  [ A     0    0     ] * [ uy ] = [ by ].
//  [ GG    0   -W'*W  ]   [ uz ]   [ bz ]
//    
// H is n x n,  A is p x n, Df is mnl x n, G is N x n where
// N = dims['l'] + sum(dims['q']) + sum( k**2 for k in dims['s'] ).
//

static
int ldl_factor(cvx_kktsolver_t *S, cvx_scaling_t *W, cvx_matrix_t *H, cvx_matrix_t *Df)
{
    cvx_ldlsolver_t *ldl = &S->u.ldl;
    cvx_conelp_problem_t *cp = S->cp;
    cvx_size_t rG, cG;
    cvx_matrix_t Kt, colvec, g0;
    cvx_matgrp_t col_g;
    
    //printf("mnl: %d\n", ldl->mnl);
    cvxm_size(&rG, &cG, cp->G);

    ldl->W = W;
    // TODO: only the LOWER triangular part?
    cvxm_mkconst(&ldl->K, 0.0);
    if (H) {
        // H is symmetric matrix; copy only the lower triangular part
        cvxm_view_map(&Kt, &ldl->K, 0, 0, ldl->n, ldl->n);
        cvxm_copy(&Kt, H, 0);
    }
    // copy A to K
    cvxm_view_map(&Kt, &ldl->K, ldl->n, 0, ldl->p, ldl->n);
    cvxm_copy(&Kt, ldl->A, 0);
    
    // copy scaled [Df; G].T to K ; column by column
    for (int k = 0; k < ldl->n; k++) {
        // g is (mnl + G.Rows(), 1) matrix, Df is (mnl, n), G is (N, n)
        if (ldl->mnl > 0) {
            // copy Df column; Df is mnl-by-n matrix
            cvxm_view_map(&colvec, Df, 0, k, ldl->mnl, 1);
            cvxm_view_map(&g0, &ldl->g, 0, 0, ldl->mnl, 1);
            cvxm_copy(&g0, &colvec, 0);
            // unmap view? no write-back
        }
        cvxm_view_map(&colvec, ldl->G, 0, k, rG, 1);
        cvxm_view_map(&g0, &ldl->g, ldl->mnl, 0, rG, 1);
        cvxm_copy(&g0, &colvec, 0);
        // make it matrix group; like group on z vector
        cvx_mgrp_init(&col_g, &g0, cp->z_g.index);
        // scale column; with W^-T
        cvx_scale(&col_g, W, CVX_INV|CVX_TRANS, &cp->work);
        // copy to K in packed storage;
        cvxm_view_map(&colvec, &ldl->K, ldl->n+ldl->p, k, ldl->ldK-ldl->p-ldl->n, 1);
        cvx_pack(&colvec, &g0, col_g.index);
    }

    // set trailing diagonal to -1
    for (int k = ldl->n + ldl->p; k < ldl->ldK; k++) {
        cvxm_set(&ldl->K, k, k, -1.0);
    }
    if (S->debug > 1) {
        //fprintf(stderr, "pre-factoring K\n");
        //cvxm_printf(stderr, "%6.3f", &ldl->K);
    }
    int err = cvxm_ldlfactor(&ldl->K, ldl->ipiv, CVX_LOWER, &ldl->work);
    if (S->debug > 1) {
        //fprintf(stderr, "LDL(K)\n");
        //cvxm_printf(stderr, "%6.3f", &ldl->K);
    }
    return err;
}

// Solve
//
//     [ H          A'   GG'*W^{-1} ]   [ ux   ]   [ bx        ]
//     [ A          0    0          ] * [ uy   [ = [ by        ]
//     [ W^{-T}*GG  0   -I          ]   [ W*uz ]   [ W^{-T}*bz ]
//
// and return ux, uy, W*uz.
//
// On entry, x, y, z contain bx, by, bz.  On exit, they contain
// the solution ux, uy, W*uz.

static
int ldl_solve(cvx_kktsolver_t *S, cvx_matrix_t *x, cvx_matrix_t *y, cvx_matgrp_t *z_g)
{
    cvx_conelp_problem_t *cp = S->cp;
    cvx_ldlsolver_t *ldl = &S->u.ldl;
    cvx_matrix_t u0;
    int err = 0;

    cvxm_view_map(&u0, &ldl->u, 0, 0, ldl->n, 1);
    cvxm_copy(&u0, x, CVX_ALL);

    cvxm_view_map(&u0, &ldl->u, ldl->n, 0, ldl->p, 1);
    cvxm_copy(&u0, y, CVX_ALL);

    // map z part; scale and copy to temp result vector as packed
    cvx_scale(z_g, ldl->W, CVX_INV|CVX_TRANS, &cp->work);
    cvxm_view_map(&u0, &ldl->u, ldl->n+ldl->p, 0, ldl->ldK-ldl->p-ldl->n, 1);
    cvxm_copy(&u0, z_g->mat, 0);

    //cvx_mat_printf(stdout, "%e", &ldl->u, "ldl solve u");
    err = cvxm_ldlsolve(&ldl->u, &ldl->K, ldl->ipiv, CVX_LOWER, &ldl->work);
    //printf("solution u\n"); cvxm_printf(stdout, "%8.5f", &ldl->u);
    
    cvxm_view_map(&u0, &ldl->u, 0, 0, ldl->n, 1);
    cvxm_copy(x, &u0, CVX_ALL);

    cvxm_view_map(&u0, &ldl->u, ldl->n, 0, ldl->p, 1);
    cvxm_copy(y, &u0,CVX_ALL);
    // unpack z part to result vector
    cvxm_view_map(&u0, &ldl->u, ldl->n+ldl->p, 0, ldl->ldK-ldl->p-ldl->n, 1);
    cvx_unpack(z_g->mat, &u0, z_g->index);
    return err;
}

//int cvx_ldlsolver_init(cvx_kktsolver_t *S, cvx_dims_t *dims, cvx_matrix_t *G, cvx_matrix_t *A, int mnl)
void cvx_ldlsolver_init(cvx_kktsolver_t *S, cvx_conelp_problem_t *cp, cvx_dimset_t *dims, int mnl)
{
    cvx_ldlsolver_t *ldl = &S->u.ldl;
    size_t ar, ac;

    S->cp = cp;
    cvxm_size(&ar, &ac, cp->A);
    ldl->ldK = ar + ac +                                // space for H and A (H can be zeroes)
        cvx_dimset_sum(dims, CVXDIM_NONLINEAR) +        // space for Df (can be zero)
        cvx_dimset_sum(dims, CVXDIM_LINEAR) +           // rest is space for G
        cvx_dimset_sum(dims, CVXDIM_SOCP) +
        cvx_dimset_sum_packed(dims, CVXDIM_SDP);

    cvxm_init(&ldl->K, ldl->ldK, ldl->ldK);
    ldl->ipiv = (int *)calloc(ldl->ldK, sizeof(int));
    
    cvxm_init(&ldl->u, ldl->ldK, 1);
    cvxm_init(&ldl->g, ldl->ldK, 1);

    int nelem = armas_d_bkfactor_work(&ldl->K, (armas_conf_t *)0);
    __mblk_init(&ldl->work, 4*nelem);

    ldl->dims = dims;
    ldl->A = cp->A;
    ldl->G = cp->G;
    ldl->n = ac;
    ldl->p = ar;
    ldl->mnl = mnl;
    S->u.fnc.factor = ldl_factor;
    S->u.fnc.solve = ldl_solve;
    //printf("mnl: %d\n", ldl->mnl);
    return;
}


// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
