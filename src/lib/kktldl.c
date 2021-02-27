/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */
#include "cvxc.h"

// forward declarations
static
int ldl_init(cvxc_kktsolver_t *S,
             cvxc_problem_t *cp,
             int n,
             int m,
             const cvxc_dimset_t *dims);

static
int ldl_factor(cvxc_kktsolver_t *S,
               cvxc_scaling_t *W,
               cvxc_matrix_t *H,
               cvxc_matrix_t *Df);
static
int ldl_solve(cvxc_kktsolver_t *S,
              cvxc_matrix_t *x,
              cvxc_matrix_t *y,
              cvxc_matgrp_t *z_g);

static
cvxc_size_t ldl_bytes(int n, int m, const cvxc_dimset_t *dims);

static
cvxc_size_t ldl_make(cvxc_kktsolver_t *kkt,
                    cvxc_problem_t *cp,
                    int n,
                    int m,
                    const cvxc_dimset_t *dims,
                    void *mem,
                    cvxc_size_t nbytes);

static
cvxc_kktsolver_t *ldl_new(cvxc_problem_t *cp,
                         int n,
                         int m,
                         const cvxc_dimset_t *dims);

static
void ldl_free(cvxc_kktsolver_t *S);

// function table
static cvxc_kktfuncs_t ldlfunctions = {
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
int ldl_factor(cvxc_kktsolver_t *S,
               cvxc_scaling_t *W,
               cvxc_matrix_t *H,
               cvxc_matrix_t *Df)
{
    //cvxc_ldlsolver_t *ldl = &S->u.ldl;
    cvxc_ldlsolver_t *ldl = (cvxc_ldlsolver_t *)S;
    cvxc_problem_t *cp = ldl->cp;
    cvxc_size_t rG, cG;
    cvxc_matrix_t Kt, colvec, g0;
    cvxc_matgrp_t col_g;

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
    for (cvxc_size_t k = 0; k < ldl->n; k++) {
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
        cvxc_mgrp_init(&col_g, &ldl->g, cp->index_g);
        // scale column; with W^-T
        cvxc_scale(&col_g, W, CVXC_INV|CVXC_TRANS, cp->work);
        // copy to K in packed storage;
        cvxm_view_map(&colvec, &ldl->K, ldl->n+ldl->p, k, ldl->ldK-ldl->p-ldl->n, 1);
        cvxc_pack(&colvec, col_g.mat, col_g.index);
    }

    // set trailing diagonal to -1
    for (cvxc_size_t k = ldl->n + ldl->p; k < ldl->ldK; k++) {
        cvxm_set(&ldl->K, k, k, -1.0);
    }
    int err = cvxm_ldlfactor(&ldl->K, ldl->ipiv, CVXC_LOWER, &ldl->work);
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
int ldl_solve(cvxc_kktsolver_t *S,
              cvxc_matrix_t *x,
              cvxc_matrix_t *y,
              cvxc_matgrp_t *z_g)
{
    cvxc_ldlsolver_t *ldl = (cvxc_ldlsolver_t *)S;
    cvxc_problem_t *cp = ldl->cp;
    cvxc_matrix_t u0;
    int err = 0;

    cvxm_view_map(&u0, &ldl->u, 0, 0, ldl->n, 1);
    cvxm_copy(&u0, x, CVXC_ALL);

    cvxm_view_map(&u0, &ldl->u, ldl->n, 0, ldl->p, 1);
    cvxm_copy(&u0, y, CVXC_ALL);

    // map z part; scale and copy to temp result vector as packed
    cvxc_scale(z_g, ldl->W, CVXC_INV|CVXC_TRANS, cp->work);
    cvxm_view_map(&u0, &ldl->u, ldl->n+ldl->p, 0, ldl->ldK-ldl->p-ldl->n, 1);
    cvxc_pack(&u0, z_g->mat, z_g->index);

    err = cvxm_ldlsolve(&ldl->u, &ldl->K, ldl->ipiv, CVXC_LOWER, &ldl->work);

    cvxm_view_map(&u0, &ldl->u, 0, 0, ldl->n, 1);
    cvxm_copy(x, &u0, CVXC_ALL);

    cvxm_view_map(&u0, &ldl->u, ldl->n, 0, ldl->p, 1);
    cvxm_copy(y, &u0,CVXC_ALL);
    // unpack z part to result vector
    cvxm_view_map(&u0, &ldl->u, ldl->n+ldl->p, 0, ldl->ldK-ldl->p-ldl->n, 1);
    cvxc_unpack(z_g->mat, &u0, z_g->index);
    return err;
}


static
int ldl_init(cvxc_kktsolver_t *S,
             cvxc_problem_t *cp,
             int n,
             int m,
             const cvxc_dimset_t *dims)
{
    cvxc_ldlsolver_t *ldl = (cvxc_ldlsolver_t *)S;
    ldl->fnc = &ldlfunctions;

    ldl->cp = cp;
    ldl->ldK = n + m;                                   // space for H and A (H can be zeroes)
    ldl->ldK += cvxc_dimset_sum(dims, CVXDIM_NONLINEAR); // space for Df (can be zero)
    ldl->ldK += cvxc_dimset_sum(dims, CVXDIM_LINEAR);    // rest is space for G
    ldl->ldK += cvxc_dimset_sum(dims, CVXDIM_SOCP);
    ldl->ldK += cvxc_dimset_sum_packed(dims, CVXDIM_SDP);

    cvxm_init(&ldl->K, ldl->ldK, ldl->ldK);
    ldl->ipiv = (int *)calloc(ldl->ldK, sizeof(int));

    cvxm_init(&ldl->u, ldl->ldK, 1);
    cvxm_init(&ldl->g, ldl->ldK, 1);

    // compute size for workspace
    cvxc_size_t wbytes = cvxm_ldlwork(&ldl->K);
    cvxc_mblk_init(&ldl->work, 4*wbytes);

    ldl->dims = dims;
    ldl->A = cp->A;
    ldl->G = cp->G;
    ldl->n = n;
    ldl->p = m;
    ldl->mnl = cvxc_dimset_sum(dims, CVXDIM_NONLINEAR);
    return 0;
}


static
cvxc_size_t ldl_bytes(int n, int m, const cvxc_dimset_t *dims)
{
    cvxc_size_t sz = n + m + cvxc_dimset_sum_packed(dims, CVXDIM_CONELP);
    cvxc_size_t need = sz*sz + 2*sz;  // dense matrix + 2 vectors;
    cvxc_size_t ipvlen = sz;
    ipvlen *= sizeof(int);
    ipvlen += __aligned128(ipvlen);
    need   *= sizeof(cvxc_float_t);
    need   += __aligned128(need);
    return need + ipvlen;
}

static
cvxc_size_t ldl_make(cvxc_kktsolver_t *kkt,
                    cvxc_problem_t *cp,
                    int n,
                    int m,
                    const cvxc_dimset_t *dims,
                    void *mem,
                    cvxc_size_t nbytes)
{
    return 0;
}

static
cvxc_kktsolver_t *ldl_new(cvxc_problem_t *cp,
                         int n,
                         int m,
                         const cvxc_dimset_t *dims)
{
    cvxc_ldlsolver_t *ldl = (cvxc_ldlsolver_t *)calloc(sizeof(cvxc_ldlsolver_t), 1);
    ldl_init((cvxc_kktsolver_t *)ldl, cp, n, m, dims);
    return  (cvxc_kktsolver_t *)ldl;
}

static
void ldl_free(cvxc_kktsolver_t *kkt)
{
    if (!kkt)
        return;

    cvxc_ldlsolver_t *ldl = (cvxc_ldlsolver_t *)kkt;
    cvxm_release(&ldl->K);
    cvxm_release(&ldl->u);
    cvxm_release(&ldl->g);
    cvxc_mblk_release(&ldl->work);
    if (ldl->ipiv)
        free(ldl->ipiv);
    free(kkt);
}



cvxc_kktfuncs_t *cvxc_ldlload(void *ptr)
{
    return &ldlfunctions;
}
