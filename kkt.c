
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
// given H, Df, W, where GG = [Df; G], and (2) compute solution for 
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
    cvx_size_t r, c, rA, cA;
    cvx_matrix_t Kt, colvec, g0;
    
    ldl->W = W;
    cvxm_zero(&ldl->K, CVX_LOWER);
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
        if (ldl->mnl > 0) {
            // copy Df column; Df is mnl-by-n matrix
            cvxm_view_map(&colvec, Df, 0, k, ldl->mnl, 1);
            cvxm_view_map(&g0, &ldl->g, 0, 0, ldl->mnl, 1);
            cvxm_copy(&g0, &colvec, 0);
            // unmap view? no write-back
        }
        cvxm_view_map(&colvec, ldl->G, 0, k, n, 1);
        cvxm_view_map(&g0, &ldl->g, ldl->mnl, 0, ldl->n, 1);
        cvxm_copy(&g0, &colvec, 0);
        // scale column; with W^-T
        cvx_scale(&ldl->g, W, CVX_INV|CVX_TRANS);
        // copy to K in packed storage;
        cvx_pack(&ldl->K, &ldl->g, &ldl->dims, 0, k*ldl->ldK + ldl->n + ldl->p, ldl->mnl);
    }

    // set trailing diagonal to -1
    for (int k = ldl->n + ldl->p; k < ldl->ldK; k++) {
        cvxm_set(&ldl->K, k, k, -1.0);
    }

    return cvxm_ldlfactor(&ldl->K, ldl->ipiv, CVX_LOWER);
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
int ldl_solve(cvx_kktsolver_t *S, cvx_matrix_t *x, cvx_matrix_t *y, cvx_matrix_t *z)
{
    cvx_ldlsolver_t *ldl = &S->u.ldl;
    cvx_matrix_t u0;
    int err = 0;

    cvxm_view_map(&u0, u, 0, 0, ldl->n, 1);
    cvxm_copy(&u0, x);

    cvxm_view_map(&u0, u, ldl->n, 0, ldl->p, 1);
    cvxm_copy(&u0, y);

    cvx_scale(z, ldl->W, CVX_INV|CVX_TRANS);

    err = cvxm_ldlsolve(&ldl->u, &ldl->K, ldl->ipiv, CVX_LOWER);
    
    cvxm_view_map(&u0, u, 0, 0, ldl->n, 1);
    cvxm_copy(x, &u0);

    cvxm_view_map(&u0, u, ldl->n, 0, ldl->p, 1);
    cvxm_copy(y, &u0);
    cvx_unpack(z, &ldl->u, ldl->dims, ldl->mnl, ldl->n + ldl->p);
    return err;
}

int cvx_ldlsolver_init(cvx_kktsolver_t *S, cvx_dims_t *dims, cvx_matrix_t *G, cvx_matrix_t *A, int mnl)
{
    cvx_ldlsolver_t *ldl = &S->u.ldl;
    size_t ar, ac;

    cvxm_size(&ar, &ac, A);
    ldl->ldK = ar + ac + mnl +
        cvx_dims_sum(dims, CVXDIM_LINEAR) +
        cvx_dims_sum(dims, CVXDIM_SOCP) +
        cvx_dims_sum_packed(dims, CVXDIM_SDP);
    cvxm_init(&ldl->K, ldl->ldK, ldl->ldK);
    ldl->ipiv = (int *)calloc(ldl->ldK, sizeof(int));
    
    cvxm_init(&ldl->u, ldl->ldK, 1);
    cvxm_init(&ldl->g, ldl->ldK, 1);

    ldl->dims = dims;
    ldl->A = A;
    ldl->G = G;
    ldl->n = ac;
    ldl->p = ar;
    S->u.fnc.factor = ldl_factor;
    S->u.fnc.solve = ldl_solve;
    return 0;
}


// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
