/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */

#include "internal.h"

typedef struct cvxc_ldlsolver {
    cvxc_problem_t *cp;
    cvxc_matrix_t K;
    cvxc_matrix_t u;
    cvxc_matrix_t g;
    cvxc_matrix_t Df;
    cvxc_memblk_t work;
    cvxc_scaling_t *W;
    int *ipiv;
    cvxc_size_t ldK;
    cvxc_size_t p;
    cvxc_size_t n;
    cvxc_size_t mnl;
    cvxc_size_t neqn;
} cvxc_ldlsolver_t;


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
int ldl2_factor(cvxc_kktsolver_t *kkt, cvxc_scaling_t *W, cvxc_matrix_t *H, cvxc_matrix_t *Df)
{
    cvxc_ldlsolver_t *ldl = (cvxc_ldlsolver_t *)kkt->private;
    cvxc_problem_t *cp = ldl->cp;
    cvxc_matrix_t Kt, Dfk, Gk, g0, g1;
    cvxc_matgrp_t g_g;
    int err;

    ldl->W = W;
    // Make a mapping of Df to protect against situation where Df itself is mapping.
    // (see function cp_factor in cp.c)
    cvxm_view_map(&ldl->Df, Df, 0, 0, ldl->mnl, ldl->n);

    cvxm_mkconst(&ldl->K, 0.0);
    if (H) {
        // H is symmetric matrix; copy only the lower triangular part
        cvxm_view_map(&Kt, &ldl->K, 0, 0, ldl->n, ldl->n);
        cvxm_copy(&Kt, H, CVXC_LOWER);
    }
    // copy A to K
    if (ldl->p > 0) {
        cvxm_view_map(&Kt, &ldl->K, ldl->n, 0, ldl->p, ldl->n);
        cvxm_copy(&Kt, cp->A, 0);
    }

    // Build lower triangular part of K column by column
    for (int k = 0; k < ldl->n; k++) {
        // k'th column of K on and below diagonal
        cvxm_view_map(&Kt, &ldl->K, k, k, ldl->n-k, 1);
        if (Df && ldl->mnl > 0) {
            // copy Df column; Df is mnl-by-n matrix
            cvxm_view_map(&Dfk, Df, 0, k, ldl->mnl, 1);
            cvxm_view_map(&g0, &ldl->g, 0, 0, ldl->mnl, 1);
            cvxm_copy(&g0, &Dfk, 0);
        }
        if (ldl->neqn - ldl->mnl > 0) {
            cvxm_view_map(&Gk, cp->G, 0, k, ldl->neqn - ldl->mnl, 1);
            cvxm_view_map(&g1, &ldl->g, ldl->mnl, 0, ldl->neqn - ldl->mnl, 1);
            cvxm_copy(&g1, &Gk, 0);
        }

        cvxc_mgrp_init(&g_g, &ldl->g, cp->index_g);
        cvxc_scale(&g_g, W, CVXC_INV|CVXC_TRANS, cp->work);
        cvxc_scale(&g_g, W, CVXC_INV, cp->work);

        if (Df && ldl->mnl > 0) {
            cvxm_view_map(&Dfk, Df, 0, k, ldl->mnl, ldl->n - k);
            cvxm_mvmult(1.0, &Kt, 1.0, &Dfk, &g0, CVXC_TRANS);
        }
        // last n-k columns of G
        if (ldl->neqn - ldl->mnl > 0) {
            cvxm_view_map(&Gk, cp->G, 0, k, ldl->neqn - ldl->mnl, ldl->n - k);
            // Kt = Kt + G^T*g
            cvxc_trisc(&g_g);
            cvxm_mvmult(1.0, &Kt, 1.0, &Gk, &g1, CVXC_TRANS);
            cvxc_triusc(&g_g);
        }
    }

    if (ldl->p > 0) {
        err = cvxm_ldlfactor(&ldl->K, ldl->ipiv, CVXC_LOWER, &ldl->work);
    } else {
        err = cvxm_cholfactor(&ldl->K, CVXC_LOWER);
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
int ldl2_solve(cvxc_kktsolver_t *kkt, cvxc_matrix_t *x, cvxc_matrix_t *y, cvxc_matgrp_t *z_g)
{
    cvxc_ldlsolver_t *ldl = (cvxc_ldlsolver_t *)kkt->private;
    cvxc_problem_t *cp = ldl->cp;
    cvxc_matrix_t u0, u1, g0, g1;
    cvxc_matgrp_t g_g;
    cvxc_float_t beta = 0.0;
    int err = 0;

    cvxm_view_map(&u0, &ldl->u, 0,      0, ldl->n, 1);

    cvxm_copy(&ldl->g, z_g->mat, 0);
    cvxc_mgrp_init(&g_g, &ldl->g, cp->index_g);

    cvxc_scale(&g_g, ldl->W, CVXC_INV|CVXC_TRANS, cp->work);
    cvxc_scale(&g_g, ldl->W, CVXC_INV, cp->work);
    if (ldl->mnl > 0) {
        // u = Df^T * g
        cvxm_view_map(&g0, &ldl->g, 0, 0, ldl->mnl, 1);
        cvxm_mvmult(0.0, &u0, 1.0, &ldl->Df, &g0, CVXC_TRANS);
        beta = 1.0;
    }
    if (ldl->neqn - ldl->mnl > 0) {
        cvxm_view_map(&g1, &ldl->g, ldl->mnl, 0, ldl->neqn - ldl->mnl, 1);

        cvxc_trisc(&g_g);
        cvxm_mvmult(beta, &u0, 1.0, cp->G, &g1, CVXC_TRANS);
        cvxc_triusc(&g_g);
    }

    cvxm_axpy(&u0, 1.0, x);
    if (ldl->p > 0) {
        cvxm_view_map(&u1, &ldl->u, ldl->n, 0, ldl->p, 1);
        cvxm_copy(&u1, y, 0);
        err = cvxm_ldlsolve(&ldl->u, &ldl->K, ldl->ipiv, CVXC_LOWER, &ldl->work);
        cvxm_copy(y, &u1, 0);
    } else {
        err = cvxm_cholsolve(&ldl->u, &ldl->K, CVXC_LOWER);
    }
    cvxm_copy(x, &u0, 0);

    // z = [Df; G]^T*x - z
    if (ldl->mnl > 0) {
        cvxm_view_map(&g0, z_g->mat, 0, 0, ldl->mnl, 1);
        cvxm_mvmult(-1.0, &g0, 1.0, &ldl->Df, x, 0);
    }
    if (ldl->neqn - ldl->mnl > 0) {
        cvxm_view_map(&g1, z_g->mat, ldl->mnl, 0, ldl->neqn - ldl->mnl, 1);
        cvxm_mvmult(-1.0, &g1, 1.0, cp->G, x, 0);
    }
    cvxc_scale(z_g, ldl->W, CVXC_TRANS|CVXC_INV, cp->work);

    return err;
}


static
int ldl2_init(cvxc_kktsolver_t *kkt, cvxc_problem_t *cp, int n, int m, const cvxc_dimset_t *dims)
{
    cvxc_size_t neqn = cvxc_dimset_sum(dims, CVXDIM_NONLINEAR)
        + cvxc_dimset_sum(dims, CVXDIM_LINEAR)
        + cvxc_dimset_sum(dims, CVXDIM_SOCP)
        + cvxc_dimset_sum_squared(dims, CVXDIM_SDP);

    cvxc_size_t sz = n + m;
    cvxc_size_t need = sz*sz + sz + neqn;
    cvxc_size_t ipvlen = m > 0 ? sz : 0;
    cvxc_size_t nwork = m > 0 ? 2 * cvxm_ldl_worksize(sz) : 0;
    cvxc_size_t nbytes =
        (ipvlen + need) * sizeof(cvxc_float_t) + nwork + sizeof(cvxc_ldlsolver_t);

    kkt->private = calloc(nbytes, 1);
    if (!kkt->private)
        return -1;

    cvxc_ldlsolver_t *ldl = (cvxc_ldlsolver_t *)kkt->private;
    cvxc_size_t offset = sizeof(cvxc_ldlsolver_t);
    unsigned char *buf = (unsigned char *)kkt->private;

    // space for H and A (H can be zeroes)
    ldl->ldK = n + m;
    ldl->cp = cp;
    ldl->n = n;
    ldl->p = m;
    ldl->neqn = neqn;
    ldl->mnl = cvxc_dimset_sum(dims, CVXDIM_NONLINEAR);

    offset += cvxm_make(&ldl->K, ldl->ldK, ldl->ldK, &buf[offset], nbytes - offset);
    offset += cvxm_make(&ldl->u, ldl->ldK, 1, &buf[offset], nbytes - offset);
    offset += cvxm_make(&ldl->g, neqn, 1, &buf[offset], nbytes - offset);

    if (m > 0) {
        ldl->ipiv = (int *)&buf[offset];
        offset += sizeof(int) * ldl->ldK;
    } else {
        ldl->ipiv = (int *)0;
    }

    // rest of the space for LDL workspace
    cvxc_mblk_make(&ldl->work, nbytes - offset, &buf[offset], nbytes - offset);

    return 0;
}

static
void ldl2_release(cvxc_kktsolver_t *kkt)
{
    if (!kkt)
        return;

    if (kkt->private)
        free(kkt->private);
}

 // function table
static cvxc_kktfuncs_t ldl2functions = {
    .factor = ldl2_factor,
    .solve  = ldl2_solve,
    .init   = ldl2_init,
    .release   = ldl2_release
};


void cvxc_kktldl2_load(cvxc_kktsolver_t *kkt)
{
    cvxc_kktsolver_init(kkt, &ldl2functions, 0);
}

