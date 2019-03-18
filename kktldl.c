
// Copyright: Harri Rautila, 2018 <harri.rautila@gmail.com>

#include "convex.h"

// forward declarations
static
int ldl_init(cvx_kktsolver_t *S,
             cvx_problem_t *cp,
             int n,
             int m,
             const cvx_dimset_t *dims);

static
int ldl_factor(cvx_kktsolver_t *S,
               cvx_scaling_t *W,
               cvx_matrix_t *H,
               cvx_matrix_t *Df);
static
int ldl_solve(cvx_kktsolver_t *S,
              cvx_matrix_t *x,
              cvx_matrix_t *y,
              cvx_matgrp_t *z_g);

static
cvx_size_t ldl_bytes(int n, int m, const cvx_dimset_t *dims);

static
cvx_size_t ldl_make(cvx_kktsolver_t *kkt,
                    cvx_problem_t *cp,
                    int n,
                    int m,
                    const cvx_dimset_t *dims,
                    void *mem,
                    cvx_size_t nbytes);

static
cvx_kktsolver_t *ldl_new(cvx_problem_t *cp,
                         int n,
                         int m,
                         const cvx_dimset_t *dims);

static
void ldl_free(cvx_kktsolver_t *S);

// function table
static cvx_kktfuncs_t ldlfunctions = {
    .new    = ldl_new,
    .factor = ldl_factor,
    .solve  = ldl_solve,
    .init   = ldl_init,
    .bytes  = ldl_bytes,
    .make   = ldl_make,
    .free   = ldl_free
};


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
int ldl_factor(cvx_kktsolver_t *S,
               cvx_scaling_t *W,
               cvx_matrix_t *H,
               cvx_matrix_t *Df)
{
    //cvx_ldlsolver_t *ldl = &S->u.ldl;
    cvx_ldlsolver_t *ldl = (cvx_ldlsolver_t *)S;
    cvx_problem_t *cp = ldl->cp;
    cvx_size_t rG, cG;
    cvx_matrix_t Kt, colvec, g0;
    cvx_matgrp_t col_g;

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

    //cvx_mat_printf(stdout, "%4.1f", Df, "Df");
    // copy scaled [Df; G].T to K ; column by column
    for (cvx_size_t k = 0; k < ldl->n; k++) {
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
        cvx_mgrp_init(&col_g, &ldl->g, cp->index_g);
        // scale column; with W^-T
        cvx_scale(&col_g, W, CVX_INV|CVX_TRANS, &cp->work);
        // copy to K in packed storage;
        cvxm_view_map(&colvec, &ldl->K, ldl->n+ldl->p, k, ldl->ldK-ldl->p-ldl->n, 1);
        cvx_pack(&colvec, col_g.mat, col_g.index);
    }

    // set trailing diagonal to -1
    for (cvx_size_t k = ldl->n + ldl->p; k < ldl->ldK; k++) {
        cvxm_set(&ldl->K, k, k, -1.0);
    }
    //cvx_mat_printf(stdout, "%4.1f", &ldl->K, "K");
    //cvx_mat_print_ifenv("LDLKKT_PREFACTOR_K", &ldl->K, "unfactored K");
    int err = cvxm_ldlfactor(&ldl->K, ldl->ipiv, CVX_LOWER, &ldl->work);
    //cvx_mat_print_ifenv("LDLKKT_POSTFACTOR_K", &ldl->K, "factored K");
    //cvx_mat_test_nan("LDLKKT_POSTFACTOR_K", &ldl->K);

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
int ldl_solve(cvx_kktsolver_t *S,
              cvx_matrix_t *x,
              cvx_matrix_t *y,
              cvx_matgrp_t *z_g)
{
    //cvx_ldlsolver_t *ldl = &S->u.ldl;
    cvx_ldlsolver_t *ldl = (cvx_ldlsolver_t *)S;
    cvx_problem_t *cp = ldl->cp;
    cvx_matrix_t u0;
    int err = 0;

    cvxm_view_map(&u0, &ldl->u, 0, 0, ldl->n, 1);
    cvxm_copy(&u0, x, CVX_ALL);

    cvxm_view_map(&u0, &ldl->u, ldl->n, 0, ldl->p, 1);
    cvxm_copy(&u0, y, CVX_ALL);

    // map z part; scale and copy to temp result vector as packed
    cvx_scale(z_g, ldl->W, CVX_INV|CVX_TRANS, &cp->work);
    cvxm_view_map(&u0, &ldl->u, ldl->n+ldl->p, 0, ldl->ldK-ldl->p-ldl->n, 1);
    cvx_pack(&u0, z_g->mat, z_g->index);

    err = cvxm_ldlsolve(&ldl->u, &ldl->K, ldl->ipiv, CVX_LOWER, &ldl->work);

    cvxm_view_map(&u0, &ldl->u, 0, 0, ldl->n, 1);
    cvxm_copy(x, &u0, CVX_ALL);

    cvxm_view_map(&u0, &ldl->u, ldl->n, 0, ldl->p, 1);
    cvxm_copy(y, &u0,CVX_ALL);
    // unpack z part to result vector
    cvxm_view_map(&u0, &ldl->u, ldl->n+ldl->p, 0, ldl->ldK-ldl->p-ldl->n, 1);
    cvx_unpack(z_g->mat, &u0, z_g->index);
    return err;
}


static
int ldl_init(cvx_kktsolver_t *S,
             cvx_problem_t *cp,
             int n,
             int m,
             const cvx_dimset_t *dims)
{
    cvx_ldlsolver_t *ldl = (cvx_ldlsolver_t *)S;
    ldl->fnc = &ldlfunctions;

    ldl->cp = cp;
    ldl->ldK = n + m;                                   // space for H and A (H can be zeroes)
    ldl->ldK += cvx_dimset_sum(dims, CVXDIM_NONLINEAR); // space for Df (can be zero)
    ldl->ldK += cvx_dimset_sum(dims, CVXDIM_LINEAR);    // rest is space for G
    ldl->ldK += cvx_dimset_sum(dims, CVXDIM_SOCP);
    ldl->ldK += cvx_dimset_sum_packed(dims, CVXDIM_SDP);

    cvxm_init(&ldl->K, ldl->ldK, ldl->ldK);
    ldl->ipiv = (int *)calloc(ldl->ldK, sizeof(int));

    cvxm_init(&ldl->u, ldl->ldK, 1);
    cvxm_init(&ldl->g, ldl->ldK, 1);

    // compute size for workspace
#if 0
    armas_wbuf_t wb = ARMAS_WBNULL;
    armas_pivot_t P;

    armas_pivot_make(&P, ldl->ldK, ldl->ipiv);
    armas_d_bkfactor_w(&ldl->K, &P, CVX_LOWER, &wb, (armas_conf_t *)0);
#endif
    cvx_size_t wbytes = cvxm_ldlwork(&ldl->K);
    __mblk_init(&ldl->work, 4*wbytes);

    ldl->dims = dims;
    ldl->A = cp->A;
    ldl->G = cp->G;
    ldl->n = n;
    ldl->p = m;
    ldl->mnl = cvx_dimset_sum(dims, CVXDIM_NONLINEAR);
    return 0;
}


static
cvx_size_t ldl_bytes(int n, int m, const cvx_dimset_t *dims)
{
    cvx_size_t sz = n + m + cvx_dimset_sum_packed(dims, CVXDIM_CONELP);
    cvx_size_t need = sz*sz + 2*sz;  // dense matrix + 2 vectors;
    cvx_size_t ipvlen = sz;
    ipvlen *= sizeof(int);
    ipvlen += __aligned128(ipvlen);
    need   *= sizeof(cvx_float_t);
    need   += __aligned128(need);
    return need + ipvlen;
}

static
cvx_size_t ldl_make(cvx_kktsolver_t *kkt,
                    cvx_problem_t *cp,
                    int n,
                    int m,
                    const cvx_dimset_t *dims,
                    void *mem,
                    cvx_size_t nbytes)
{
    return 0;
}

static
cvx_kktsolver_t *ldl_new(cvx_problem_t *cp,
                         int n,
                         int m,
                         const cvx_dimset_t *dims)
{
    cvx_ldlsolver_t *ldl = (cvx_ldlsolver_t *)calloc(sizeof(cvx_ldlsolver_t), 1);
    ldl_init((cvx_kktsolver_t *)ldl, cp, n, m, dims);
    return  (cvx_kktsolver_t *)ldl;
}

static
void ldl_free(cvx_kktsolver_t *kkt)
{
    if (!kkt)
        return;

    cvx_ldlsolver_t *ldl = (cvx_ldlsolver_t *)kkt;
    cvxm_release(&ldl->K);
    cvxm_release(&ldl->u);
    cvxm_release(&ldl->g);
    __mblk_release(&ldl->work);
    if (ldl->ipiv)
        free(ldl->ipiv);
    free(kkt);
}



cvx_kktfuncs_t *cvx_ldlload(void *ptr)
{
    return &ldlfunctions;
}


// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
