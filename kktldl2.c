
// Copyright: Harri Rautila, 2018 <harri.rautila@gmail.com>

#include "convex.h"

// forward declarations
static
int ldl2_init(cvx_kktsolver_t *S,
             cvx_problem_t *cp,
             int n,
             int m,
             const cvx_dimset_t *dims);

static
int ldl2_factor(cvx_kktsolver_t *S,
               cvx_scaling_t *W,
               cvx_matrix_t *H,
               cvx_matrix_t *Df);
static
int ldl2_solve(cvx_kktsolver_t *S,
              cvx_matrix_t *x,
              cvx_matrix_t *y,
              cvx_matgrp_t *z_g);

static
cvx_size_t ldl2_bytes(int n, int m, const cvx_dimset_t *dims);

static
cvx_size_t ldl2_make(cvx_kktsolver_t *kkt,
                    cvx_problem_t *cp,
                    int n,
                    int m,
                    const cvx_dimset_t *dims,
                    void *mem,
                    cvx_size_t nbytes);

static
cvx_kktsolver_t *ldl2_new(cvx_problem_t *cp,
                         int n,
                         int m,
                         const cvx_dimset_t *dims);

static
void ldl2_free(cvx_kktsolver_t *S);

// function table
static cvx_kktfuncs_t ldl2functions = {
    .new    = ldl2_new,
    .factor = ldl2_factor,
    .solve  = ldl2_solve,
    .init   = ldl2_init,
    .bytes  = ldl2_bytes,
    .make   = ldl2_make,
    .free   = ldl2_free
};

/*
    Solution of KKT equations by a dense LDL factorization of the 2 x 2
    system.

    Returns a function that (1) computes the LDL factorization of

        [ H + GG' * W^{-1} * W^{-T} * GG   A' ]
        [                                     ]
        [ A                                0  ]

    given H, Df, W, where GG = [Df; G], and (2) returns a function for
    solving

        [ H    A'   GG'   ]   [ ux ]   [ bx ]
        [ A    0    0     ] * [ uy ] = [ by ].
        [ GG   0   -W'*W  ]   [ uz ]   [ bz ]

    H is n x n,  A is p x n, Df is mnl x n, G is N x n where
*/

static
int ldl2_factor(cvx_kktsolver_t *S, cvx_scaling_t *W, cvx_matrix_t *H, cvx_matrix_t *Df)
{
    cvx_ldlsolver_t *ldl = (cvx_ldlsolver_t *)S;
    cvx_problem_t *cp = ldl->cp;
    cvx_size_t rG, cG;
    cvx_matrix_t Kt, Dfk, Gk, g0, g1;
    cvx_matgrp_t g_g;
    int mnl = cvx_dimset_sum(cp->dims, CVXDIM_NONLINEAR);
    int err;

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
    if (ldl->p > 0) {
        cvxm_view_map(&Kt, &ldl->K, ldl->n, 0, ldl->p, ldl->n);
        cvxm_copy(&Kt, ldl->A, 0);
    }

    // Build lower triangular part of K column by column
    for (int k = 0; k < ldl->n; k++) {
        // k'th column of K on and below diagonal
        cvxm_view_map(&Kt, &ldl->K, k, k, ldl->n-k, 1);
        if (Df && mnl > 0) {
            // copy Df column; Df is mnl-by-n matrix
            cvxm_view_map(&Dfk, Df, 0, k, mnl, 1);
            cvxm_view_map(&g0, &ldl->g, 0, 0, mnl, 1);
            cvxm_copy(&g0, &Dfk, 0);
        }
        cvxm_view_map(&Gk, ldl->G, 0, k, rG, 1);
        cvxm_view_map(&g1, &ldl->g, mnl, 0, rG, 1);
        cvxm_copy(&g1, &Gk, 0);

        // TODO: works only w/o non-linear part (z_g.index ??)
        cvx_mgrp_init(&g_g, &g1, cp->index_g);
        cvx_scale(&g_g, W, CVX_INV|CVX_TRANS, &cp->work);
        cvx_scale(&g_g, W, CVX_INV, &cp->work);

        if (Df && mnl > 0) {
            cvxm_mvmult(1.0, &Kt, 1.0, Df, &g0, CVX_TRANS);
        }
        // last n-k columns of G
        cvxm_view_map(&Gk, ldl->G, 0, k, rG, cG-k);
        // Kt = Kt + G^T*g
        cvx_sgemv(1.0, &Kt, 1.0, &Gk, &g_g, CVX_TRANS);
    }

    if (ldl->p > 0) {
        err = cvxm_ldlfactor(&ldl->K, ldl->ipiv, CVX_LOWER, &ldl->work);
    } else {
        err = cvxm_cholfactor(&ldl->K, CVX_LOWER);
    }
    return err;
}

/*
   Solve

        [ H + GG' * W^{-1} * W^{-T} * GG    A' ]   [ ux ]
        [                                      ] * [    ]
        [ A                                 0  ]   [ uy ]

                 [ bx + GG' * W^{-1} * W^{-T} * bz ]
             =   [                                 ]
                 [ by                              ]

   and return x, y, W*z = W^{-T} * (GG*x - bz).

 */
static
int ldl2_solve(cvx_kktsolver_t *S, cvx_matrix_t *x, cvx_matrix_t *y, cvx_matgrp_t *z_g)
{
    cvx_ldlsolver_t *ldl = (cvx_ldlsolver_t *)S;
    cvx_problem_t *cp = ldl->cp;
    cvx_matrix_t u0, u1;
    cvx_matgrp_t g_g, x_g;
    cvx_float_t beta = 0.0;
    int mnl = cvx_dimset_sum(cp->dims, CVXDIM_NONLINEAR);
    int err = 0;

    cvxm_view_map(&u0, &ldl->u, 0,      0, ldl->n, 1);
    cvxm_view_map(&u1, &ldl->u, ldl->n, 0, ldl->p, 1);

    // TODO: z_g.index ?? in non-linear case
    cvxm_copy(&ldl->g, z_g->mat, 0);
    cvx_mgrp_init(&g_g, &ldl->g, z_g->index);

    cvx_scale(&g_g, ldl->W, CVX_INV|CVX_TRANS, &cp->work);
    cvx_scale(&g_g, ldl->W, CVX_INV, &cp->work);
    if (mnl > 0) {
        // TODO: non-linear
        beta = 1.0;
    }
    cvx_sgemv(beta, &ldl->u, 1.0, ldl->G, &g_g, CVX_TRANS);
    cvxm_axpy(&u0, 1.0, x);
    cvxm_copy(&u1, y, 0);

    if (ldl->p > 0) {
        err = cvxm_ldlsolve(&ldl->u, &ldl->K, ldl->ipiv, CVX_LOWER, &ldl->work);
    } else {
        err = cvxm_cholsolve(&ldl->u, &ldl->K, CVX_LOWER);
    }

    cvxm_copy(x, &u0, 0);
    cvxm_copy(y, &u1, 0);

    if (mnl > 0) {
    }
    cvx_mgrp_init(&x_g, x, (cvx_index_t *)0);
    // z = G*x - z
    cvx_sgemv(-1.0, z_g->mat, 1.0, cp->G, &x_g, 0);
    cvx_scale(z_g, ldl->W, CVX_TRANS|CVX_INV, &cp->work);

    return err;
}


static
int ldl2_init(cvx_kktsolver_t *S, cvx_problem_t *cp, int n, int m, const cvx_dimset_t *dims)
{
    cvx_ldlsolver_t *ldl = (cvx_ldlsolver_t *)S;
    size_t grows, gcols;


    cvxm_size(&grows, &gcols, cp->G);
    ldl->ldK = n + m;                                // space for H and A (H can be zeroes)

    cvxm_init(&ldl->K, ldl->ldK, ldl->ldK);
    if (m > 0)
        ldl->ipiv = (int *)calloc(ldl->ldK, sizeof(int));
    else
        ldl->ipiv = (int *)0;

    cvxm_init(&ldl->u, ldl->ldK, 1);
    cvxm_init(&ldl->g, grows, 1);

    //int nelem = armas_d_bkfactor_work(&ldl->K, (armas_conf_t *)0);
    int nelem = cvxm_ldlwork(&ldl->K);
    __mblk_init(&ldl->work, 4*nelem);

    ldl->cp = cp;
    ldl->dims = dims;
    ldl->A = cp->A;
    ldl->G = cp->G;
    ldl->n = n;
    ldl->p = m;

    return 0;
}


static
cvx_size_t ldl2_bytes(int n, int m, const cvx_dimset_t *dims)
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
cvx_size_t ldl2_make(cvx_kktsolver_t *kkt,
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
cvx_kktsolver_t *ldl2_new(cvx_problem_t *cp,
                          int n,
                          int m,
                          const cvx_dimset_t *dims)
{
    cvx_ldlsolver_t *ldl = (cvx_ldlsolver_t *)calloc(sizeof(cvx_ldlsolver_t), 1);
    ldl2_init((cvx_kktsolver_t *)ldl, cp, n, m, dims);
    return  (cvx_kktsolver_t *)ldl;
}

static
void ldl2_free(cvx_kktsolver_t *kkt)
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

cvx_kktfuncs_t *cvx_ldl2load()
{
    return &ldl2functions;
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
