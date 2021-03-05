/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */
#include "internal.h"

typedef struct cvxc_ldlsolver {
    cvxc_kktfuncs_t *fnc;
    cvxc_problem_t *cp;
    cvxc_matrix_t K;
    cvxc_matrix_t u;
    cvxc_matrix_t g;
    cvxc_memblk_t work;
    cvxc_scaling_t *W;
    cvxc_matrix_t *A;
    cvxc_matrix_t *G;
    cvxc_matrix_t *Df;
    const cvxc_dimset_t *dims;
    size_t ldK;
    int *ipiv;
    size_t p;
    size_t n;
    size_t mnl;
} cvxc_ldlsolver_t;


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
    cvxc_ldlsolver_t *ldl = (cvxc_ldlsolver_t *)S->private;
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
    cvxc_ldlsolver_t *ldl = (cvxc_ldlsolver_t *)S->private;
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

extern cvxc_size_t cvxm_ldl_worksize(int);

static
cvxc_size_t ldl_bytes(int n, int m, const cvxc_dimset_t *dims)
{
    cvxc_size_t neqn = cvxc_dimset_sum(dims, CVXDIM_NONLINEAR)
        + cvxc_dimset_sum(dims, CVXDIM_LINEAR)
        + cvxc_dimset_sum(dims, CVXDIM_SOCP)
        + cvxc_dimset_sum_packed(dims, CVXDIM_SDP);

    cvxc_size_t sz = n + m + neqn;
    cvxc_size_t need = sz*sz + 2*sz;  // dense matrix + 2 vectors;
    cvxc_size_t ipvlen = sz;
    cvxc_size_t nwork = 2 * cvxm_ldl_worksize(sz);

    ipvlen *= sizeof(int);
    ipvlen += __aligned128(ipvlen);
    need   *= sizeof(cvxc_float_t);
    need   += __aligned128(need);
    return need + ipvlen + nwork;
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
    cvxc_ldlsolver_t *ldl = kkt->private;
    cvxc_size_t offset = sizeof(cvxc_ldlsolver_t);
    unsigned char *buf = (unsigned char *)mem;

    ldl->cp = cp;
    ldl->ldK = n + m;
    ldl->ldK += cvxc_dimset_sum(dims, CVXDIM_NONLINEAR);
    ldl->ldK += cvxc_dimset_sum(dims, CVXDIM_LINEAR);
    ldl->ldK += cvxc_dimset_sum(dims, CVXDIM_SOCP);
    ldl->ldK += cvxc_dimset_sum_packed(dims, CVXDIM_SDP);

    offset += cvxm_make(&ldl->K, ldl->ldK, ldl->ldK, &buf[offset], nbytes - offset);
    offset += cvxm_make(&ldl->u, ldl->ldK, 1, &buf[offset], nbytes - offset);
    offset += cvxm_make(&ldl->g, ldl->ldK, 1, &buf[offset], nbytes - offset);

    ldl->ipiv = (int *)&buf[offset];
    offset += sizeof(int) * ldl->ldK;

    // rest of the space for LDL workspace
    cvxc_mblk_make(&ldl->work, nbytes - offset, &buf[offset], nbytes - offset);

    ldl->dims = dims;
    ldl->A = cp->A;
    ldl->G = cp->G;
    ldl->n = n;
    ldl->p = m;
    ldl->mnl = cvxc_dimset_sum(dims, CVXDIM_NONLINEAR);
    return 0;
}

static
int ldl_init(cvxc_kktsolver_t *S,
             cvxc_problem_t *cp,
             int n,
             int m,
             const cvxc_dimset_t *dims)
{

    cvxc_size_t nbytes = ldl_bytes(n, m, dims);
    S->private = calloc(nbytes + sizeof(cvxc_ldlsolver_t), 1);
    if (!S->private)
        return -1;

    return ldl_make(S, cp, n, m, dims, S->private, nbytes + sizeof(cvxc_ldlsolver_t));
}

static
void ldl_release(cvxc_kktsolver_t *kkt)
{
    if (!kkt)
        return;

    free(kkt->private);
    kkt->private = 0;
}

// function table
static cvxc_kktfuncs_t ldlfunctions = {
    .factor = ldl_factor,
    .solve  = ldl_solve,
    .init   = ldl_init,
    .bytes  = ldl_bytes,
    .make   = ldl_make,
    .release= ldl_release
};


void cvxc_kktldl_load(cvxc_kktsolver_t *kkt)
{
    kkt->vtable = &ldlfunctions;
    kkt->private = 0;
    kkt->next = 0;
}


cvxc_kktfuncs_t *cvxc_ldlload(void *ptr)
{
    return &ldlfunctions;
}
