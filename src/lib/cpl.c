/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */

#include "internal.h"

#ifndef __ZERO
#define __ZERO 0.0
#endif

#define EFRM "%8.5f"
//#define EFRM "%9.2e"
#define STEP 0.99
#define BETA 0.5
#define ALPHA 0.01
#define EXPON 3
#define MAX_RELAXED_ITERS 8

#define CPL_RESTORE 1
#define CPL_SAVE 0

static inline
cvxc_size_t __WORKBYTES(cvxc_size_t n)
{
    return 8*n*n*sizeof(cvxc_float_t);
}


// Evaluates residuals in Newton equations:
//
//     [ vx ]     [ 0     ]   [ H  A' GG' ] [ ux        ]
//     [ vy ] -=  [ 0     ] + [ A  0  0   ] [ uy        ]
//     [ vz ]     [ W'*us ]   [ GG 0  0   ] [ W^{-1}*uz ]
//
//     vs -= lmbda o (uz + us).
static
int cpl_res(cvxc_problem_t *cp,
            cvxc_matrix_t *ux,
            cvxc_matrix_t *uy,
            cvxc_matgrp_t *uz_g,
            cvxc_matgrp_t *us_g,
            cvxc_matrix_t *vx,
            cvxc_matrix_t *vy,
            cvxc_matgrp_t *vz_g,
            cvxc_matgrp_t *vs_g,
            cvxc_scaling_t *W,
            cvxc_matgrp_t *lmbda_g)
{
    int err = 0;
    cvxc_cpl_internal_t *cpi = cp->u.cpl;

    cvxc_matrix_t *us = us_g->mat, *uz = uz_g->mat;
    cvxc_matrix_t *vs = vs_g->mat, *vz = vz_g->mat;

    cvxc_matgrp_t ux_g = (cvxc_matgrp_t){ .mat = ux, .index = (cvxc_index_t *)0 };

    // vx = vx - H*ux - A^T*uy - G^T*W^-1*uz
    cvxm_mvmult(1.0, vx, -1.0, &cpi->H, ux, 0);
    cvxm_mvmult(1.0, vx, -1.0, cp->A, uy, CVXC_TRANSA);
    cvxm_copy(&cpi->wz3, uz, CVXC_ALL);
    cvxc_scale(&cpi->wz3_g, &cpi->W, CVXC_INV, &cpi->work);
    cvxc_sgemv2(1.0, vx, -1.0, &cpi->Df, cp->G, &cpi->wz3_g, CVXC_TRANS);

    // vy = vy - A*ux
    cvxm_mult(1.0, vy, -1.0, cp->A, ux, 0);

    // vz = vz - W'*us - GG*ux
    cvxm_copy(&cpi->ws3, us, 0);
    cvxc_scale(&cpi->ws3_g, &cpi->W, CVXC_TRANS, &cpi->work);
    cvxm_axpy(vz, -1.0, &cpi->ws3);
    cvxc_sgemv2(1.0, vz, -1.0, &cpi->Df, cp->G, &ux_g, 0);

    // vs = vs - lmbda o (uz + us)
    cvxm_copy(&cpi->ws3, us, 0);
    cvxm_axpy(&cpi->ws3, 1.0, uz);
    cvxc_sprod(&cpi->ws3_g, lmbda_g, CVXC_DIAG, &cpi->work);
    cvxm_axpy(vs, -1.0, &cpi->ws3);

    return err;
}

// Solve
//
//     [ H  A'  GG'  ] [ ux        ]   [ bx                    ]
//     [ A  0   0    ] [ uy        ] = [ by                    ]
//     [ GG 0  -W'*W ] [ W^{-1}*uz ]   [ bz - W'*(lmbda o\ bs) ]
//
//     us = lmbda o\ bs - uz.
static
int f4_no_ir(cvxc_problem_t *cp,
             cvxc_matrix_t *x,
             cvxc_matrix_t *y,
             cvxc_matgrp_t *z_g,
             cvxc_matgrp_t *s_g)
{
    cvxc_matrix_t *s = s_g->mat, *z = z_g->mat;
    cvxc_cpl_internal_t *cpi = cp->u.cpl;

    // s = lmbda o\ s
    cvxc_sinv(s_g, &cpi->lmbda_g, &cpi->work);

    // z := z - W'*s
    cvxm_copy(&cpi->ws3, s, 0);
    cvxc_scale(&cpi->ws3_g, &cpi->W, CVXC_TRANS, &cpi->work);
    cvxm_axpy(z, -1.0, &cpi->ws3);
    // solve
    cvxc_kktsolve(cp->solver, x, y, z_g);

    // s := s - z
    cvxm_axpy(s, -1.0, z);

    return 0;
}


// f6(x, y, z, tau, s, kappa) solves the same system as f6_no_ir,
// but applies iterative refinement.
static
int f4(cvxc_problem_t *cp,
       cvxc_matrix_t *x,
       cvxc_matrix_t *y,
       cvxc_matgrp_t *z_g,
       cvxc_matgrp_t *s_g,
       int refinement)
{
    cvxc_matrix_t *s = s_g->mat, *z = z_g->mat;
    cvxc_cpl_internal_t *cpi = cp->u.cpl;
    int err = 0;

    if (refinement > 0 /*|| (flags | CVXC_DEBUG) != 0*/) {
        cvxm_copy(&cpi->wx, x, 0);
        cvxm_copy(&cpi->wy, y, 0);
        cvxm_copy(&cpi->wz, z, 0);
        cvxm_copy(&cpi->ws, s, 0);
    }
    err = f4_no_ir(cp, x, y, z_g, s_g);
    for (int i = 0; i < refinement; i++) {
        cvxm_copy(&cpi->wx2, &cpi->wx, 0);
        cvxm_copy(&cpi->wy2, &cpi->wy, 0);
        cvxm_copy(&cpi->wz2, &cpi->wz, 0);
        cvxm_copy(&cpi->ws2, &cpi->ws, 0);

        cpl_res(cp, x, y, z_g, s_g, &cpi->wx2, &cpi->wy2, &cpi->wz2_g,
                &cpi->ws2_g, &cpi->W, &cpi->lmbda_g);

        f4_no_ir(cp, &cpi->wx2, &cpi->wy2, &cpi->wz2_g, &cpi->ws2_g);

        cvxm_axpy(x, 1.0, &cpi->wx2);
        cvxm_axpy(y, 1.0, &cpi->wy2);
        cvxm_axpy(s, 1.0, &cpi->ws2);
        cvxm_axpy(z, 1.0, &cpi->wz2);
    }
    return err;
}

static
int cpl_factor(cvxc_kktsolver_t *S, cvxc_scaling_t *W, cvxc_matrix_t *x, cvxc_matrix_t *z)
{
    cvxc_problem_t *cp = (cvxc_problem_t *)S->private;
    cvxc_cpl_internal_t *cpi = cp->u.cpl;
    F2(cp->F, __cvxnil, &cpi->Df, &cpi->H, x, z);
    return cvxc_kktfactor(S->next, W, &cpi->H, &cpi->Df);
}

static
int cpl_solve(cvxc_kktsolver_t *S, cvxc_matrix_t *x, cvxc_matrix_t *y, cvxc_matgrp_t *z_g)
{
    return cvxc_kktsolve(S->next, x, y, z_g);
}

static
void cpl_release(cvxc_kktsolver_t *kkt)
{
    if (!kkt)
        return;

    if (kkt->next && kkt->next->vtable->release)
        (*kkt->next->vtable->release)(kkt->next);
    kkt->private = 0;
}

static cvxc_kktfuncs_t cpl_kktfunctions = {
    .factor = cpl_factor,
    .solve  = cpl_solve,
    .release= cpl_release
};

void cvxc_cpl_solver_init(cvxc_kktsolver_t *cs, cvxc_problem_t *cp, cvxc_kktsolver_t *next)
{
    cs->vtable = &cpl_kktfunctions;
    cs->next = next;
    cs->private = cp;
}

/**
 * @brief Compute memory allocation needed for CONVEXLP problem
 *
 * @param[in] n  Number of variables
 * @param[in] m  Number of equality constrains (rows of A matrix)
 * @param[in] dims Dimensions of inequality constraints, linear, scop and sdp
 * @param[in] nonlinear Non zero if objective function is non-linear
 *
 * @return Number of bytes of memory needed.
 */
cvxc_size_t cvxc_cpl_bytes(int n, int m, const cvxc_dimset_t *dims, int nonlinear)
{
    cvxc_size_t cdim_mnl  = cvxc_dimset_sum_squared(dims, CVXDIM_CONVEXLP);
    cvxc_size_t cdim_diag = cvxc_dimset_sum(dims, CVXDIM_CONELP);
    cvxc_size_t sdim      = cvxc_dimset_sum(dims, CVXDIM_SDP);
    cvxc_size_t maxsdp    = cvxc_dimset_max(dims, CVXDIM_SDP);
    cvxc_size_t mnl       = cvxc_dimset_sum(dims, CVXDIM_NONLINEAR);
    cvxc_size_t total = 0;
    cvxc_size_t nbytes = 0;
    if (nonlinear != 0)
        nonlinear = 1;
    mnl += nonlinear;
    cdim_mnl += nonlinear;

    // size of standard index set (if nonlinear then 2 times)
    nbytes += cvxc_index_bytes(dims, CVXC_INDEX_NORMAL) * (1 + nonlinear);
    // size of packed index set
    nbytes += cvxc_index_bytes(dims, CVXC_INDEX_PACKED);
    // size of diagonal index set
    nbytes += cvxc_index_bytes(dims, CVXC_INDEX_DIAG);
    // size of SDP diagonal index set
    nbytes += cvxc_index_bytes(dims, CVXC_INDEX_SIGS);

    // SIZES
    // cdim_mnl: (23)
    //   z, dz, dz0, dz2, dz20, z0, rz, rz0, newz, newrz, wz, wz2, wz3 (13)
    //   s, ds, ds0, ds2, ds20, s0, news, ws, ws2, ws3  (10)
    // n: (9)
    //   x, dx, x0, rx, rx0, newx, newrx, wx, wx2
    // m: (8)
    //   y, dy, y0, ry, ry0, newy, wy, wy2
    // cdim_diag+mnl: (2)
    //   lmbda, lmbdasq
    // sdim: (2)
    //   sigs, sigz
    //
    // rz -> rznl + rzl; rz0 -> rznl0 + rzl0; wz2 -> wz2nl + wz2l
    total += (nonlinear + 9)*n; // for x, dx, rx, x0, newx, rx0, wx, wx2, newrz
    total += 8*m;               // for y, dy, ry, y0, newy, ry0, wy, wy2
    total += 10*cdim_mnl;       // for s, ds,     s0, news,      ws, ws2, ds0, ds2, ds20, ws3
    total += 13*cdim_mnl;       // for z, dz, rz, z0, newz, rz0, wz, wz2, dz0, dz2, dz20, wz3, newrz0

    total += 2*(cdim_diag+mnl) + 2;   // for lmbda, lmbdasq
    total += 2*sdim;            // for sigs, sigz

    // f, Df, H
    total += mnl*2;   // for f, newf
    total += mnl*n*2; // for Df, newDf
    total += n*n;               // for H

    nbytes += total*sizeof(cvxc_float_t);

    // workspace for SDP constraint scaling; TODO: think about this.
    if (maxsdp > 0)
        nbytes += __WORKBYTES(maxsdp);

    cvxc_size_t isize, szw;
    // calculte space need for scaling matrix (second for state save)
    szw     = cvxc_scaling_bytes(&isize, dims);
    nbytes += 2*szw;

    return nbytes;
}

#define __INIT(a, fn)                                           \
    do {                                                        \
        a = fn;                                                 \
        if ((a) == 0) {                                         \
            fprintf(stderr, "__INIT@%d [%ld,%ld] '%s'\n", __LINE__, nbytes, offset, #fn); \
            abort();                                            \
        }                                                       \
        nbytes -= (a);                                          \
        offset += (a);                                          \
    } while (0)

#define __INITC(a, l, var, fn)                                          \
do {                                                                    \
    if ((l) == 0) {                                                     \
        cvxm_map_data(var, 0, 1, (cvxc_float_t *)0);              \
    }                                                                   \
    else {                                                              \
        a = fn;                                                         \
        if ((a) == 0) {                                                 \
            fprintf(stderr, "__INIT@%d [%ld,%ld] '%s'\n", __LINE__, nbytes, offset, #fn); \
            abort();                                                    \
        }                                                               \
        nbytes -= (a);                                                  \
        offset += (a);                                                  \
    } \
} while (0)

static inline
int cvxm_xvector(int nl, cvxc_matrix_t *x, cvxc_size_t m, cvxc_size_t n, void *data, cvxc_size_t ndata)
{
    if (nl)
        return cvxm_make_epi(x, m, n, data, ndata);
    return cvxm_make(x, m, n, data, ndata);
}

/**
 * @brief Overlay problem variables onto memory block
 *
 * @param[in,out] prob
 *      ConeLP problem structure
 * @param[in] n
 *      Number of variables (length of x -vector)
 * @param[in] m
 *      Number of equality constrains (rows of A matrix)
 * @param[in] 
 *      dims Dimensions of inequality constraints, linear, scop and sdp
 * @param[in]
 *      nl   Nonzero for non-linear objective function
 * @param[in] memory
 *      Pointer to memory block 
 * @param[in] nbytes
 *      Size of memory block in bytes
 *
 * @return Number of bytes of memory used, zero if memory was not large enough.
 *
 */
cvxc_size_t cvxc_cpl_make(cvxc_problem_t *cp,
                        int n,
                        int m,
                        const cvxc_dimset_t *dims,
                        int nl,
                        void *memory,
                        cvxc_size_t msize)
{
    cvxc_cpl_internal_t *cpi = cp->u.cpl;
    cvxc_size_t offset = 0;
    cvxc_size_t used = 0;
    cvxc_size_t nbytes = msize;
    unsigned char *bytes = (unsigned char *)memory;

    if (nl != 0)
        nl = 1;

    cvxc_size_t cdim_mnl  = cvxc_dimset_sum_squared(dims, CVXDIM_CONVEXLP) + nl;

    cp->cdim = cdim_mnl;

    // __INIT macro assumes variables offset and nbytes;
    // overlay index sets
    __INIT(used, cvxc_index_make(&cpi->index_full, dims, CVXC_INDEX_NORMAL, &bytes[offset],  nbytes));
    __INIT(used, cvxc_index_make(&cpi->index_packed, dims, CVXC_INDEX_PACKED, &bytes[offset],  nbytes));
    __INIT(used, cvxc_index_make(&cpi->index_diag, dims, CVXC_INDEX_DIAG, &bytes[offset],  nbytes));
    __INIT(used, cvxc_index_make(&cpi->index_sig, dims, CVXC_INDEX_SIGS, &bytes[offset],  nbytes));

    if (nl) {
        cvxc_dimset_t ldims = *dims;
        ldims.iscpt = 0;
        __INIT(used, cvxc_index_make(&cpi->index_cpt, &ldims, CVXC_INDEX_NORMAL, &bytes[offset],  nbytes));
    }
    //fprintf(stderr, "cpl make: indexes %ld\n", offset);

    // map result matrix;
    __INIT(used, cvxm_xvector(nl, &cpi->x, n, 1, &bytes[offset], nbytes));
    __INITC(used, m, &cpi->y, cvxm_make(&cpi->y, m, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->s, cdim_mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->z, cdim_mnl, 1, &bytes[offset], nbytes));

    if (nl)  {
        __INIT(used, cvxm_make_epi(&cpi->c0, n, 1, &bytes[offset], nbytes));
    }

    // dx, dy, ds, dz
    __INIT(used, cvxm_xvector(nl, &cpi->dx, n, 1, &bytes[offset], nbytes));
    __INITC(used, m, &cpi->dy, cvxm_make(&cpi->dy, m, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->ds, cdim_mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->dz, cdim_mnl, 1, &bytes[offset], nbytes));

    // ds0, dz0, ds2, dz2, ds20, dz20
    __INIT(used, cvxm_make(&cpi->ds0, cdim_mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->dz0, cdim_mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->ds2, cdim_mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->dz2, cdim_mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->ds20, cdim_mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->dz20, cdim_mnl, 1, &bytes[offset], nbytes));

    // rx, ry, rz
    __INIT(used, cvxm_xvector(nl, &cpi->rx, n, 1, &bytes[offset], nbytes));
    __INITC(used, m, &cpi->ry, cvxm_make(&cpi->ry, m, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->rz, cdim_mnl, 1, &bytes[offset], nbytes));

    // x0, y0, z0, s0
    __INIT(used, cvxm_xvector(nl, &cpi->x0, n, 1, &bytes[offset], nbytes));
    __INITC(used, m, &cpi->y0, cvxm_make(&cpi->y0, m, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->s0, cdim_mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->z0, cdim_mnl, 1, &bytes[offset], nbytes));

    // newx, newy, newz, news, newrx
    __INIT(used, cvxm_xvector(nl, &cpi->newx, n, 1, &bytes[offset], nbytes));
    __INITC(used, m, &cpi->newy, cvxm_make(&cpi->newy, m, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->news, cdim_mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->newz, cdim_mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_xvector(nl, &cpi->newrx, n, 1, &bytes[offset], nbytes));

    // rx0, ry0, rz0, newrz0
    __INIT(used, cvxm_xvector(nl, &cpi->rx0, n, 1, &bytes[offset], nbytes));
    __INITC(used, m, &cpi->ry0, cvxm_make(&cpi->ry0, m, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->rz0, cdim_mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->newrz0, cdim_mnl, 1, &bytes[offset], nbytes));

    // wx, wy, ws, wz
    __INIT(used, cvxm_xvector(nl, &cpi->wx, n, 1, &bytes[offset], nbytes));
    __INITC(used, m, &cpi->wy2, cvxm_make(&cpi->wy, m, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->ws, cdim_mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->wz, cdim_mnl, 1, &bytes[offset], nbytes));

    // wx2, wy2, ws2, wz2
    __INIT(used, cvxm_xvector(nl, &cpi->wx2, n, 1, &bytes[offset], nbytes));
    __INITC(used, m, &cpi->wy2, cvxm_make(&cpi->wy2, m, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->ws2, cdim_mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->wz2, cdim_mnl, 1, &bytes[offset], nbytes));

    // ws3, wz3
    __INIT(used, cvxm_make(&cpi->ws3, cdim_mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->wz3, cdim_mnl, 1, &bytes[offset], nbytes));

    cvxc_size_t cdim_diag =
        cvxc_dimset_sum(dims, CVXDIM_CONVEXLP) + nl;

    cp->cdim_diag = cdim_diag;

    // lmbda, lmbdasq
    __INIT(used, cvxm_make(&cpi->lmbda,   cdim_diag, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->lmbdasq, cdim_diag, 1, &bytes[offset], nbytes));

    cvxc_size_t sdim   = cvxc_dimset_sum(dims, CVXDIM_SDP);
    cvxc_size_t maxsdp = cvxc_dimset_max(dims, CVXDIM_SDP);

    // sigs, sigz; space for eigenvalues; zero if no SDP constraints
    __INITC(used, sdim, &cpi->sigs, cvxm_make(&cpi->sigs, sdim, 1, &bytes[offset], nbytes));
    __INITC(used, sdim, &cpi->sigz, cvxm_make(&cpi->sigz, sdim, 1, &bytes[offset], nbytes));

    // f, Df, H matrices
    int mnl = cvxc_dimset_sum(dims, CVXDIM_NONLINEAR) + nl;
    __INIT(used, cvxm_make(&cpi->f, mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->newf, mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->Df, mnl, n, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->newDf, mnl, n, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->H, n, n, &bytes[offset], nbytes));

    // scaling matrix
    __INIT(used, cvxc_scaling_make(&cpi->W,   dims, &bytes[offset], nbytes));
    __INIT(used, cvxc_scaling_make(&cpi->W0, dims, &bytes[offset], nbytes));

    // workspace for SDP contraints handling
    cvxc_mblk_empty(&cpi->work);
    if (maxsdp > 0) {
        __INIT(used, cvxc_mblk_make(&cpi->work, __WORKBYTES(maxsdp), &bytes[offset], nbytes));
    }

    // setup matrix group variables
    cvxc_mgrp_init(&cpi->s_g,   &cpi->s,   &cpi->index_full);
    cvxc_mgrp_init(&cpi->z_g,   &cpi->z,   &cpi->index_full);
    cvxc_mgrp_init(&cpi->z0_g,  &cpi->z0,  &cpi->index_full);
    cvxc_mgrp_init(&cpi->s0_g,  &cpi->s0,  &cpi->index_full);
    cvxc_mgrp_init(&cpi->ds_g,  &cpi->ds,  &cpi->index_full);
    cvxc_mgrp_init(&cpi->dz_g,  &cpi->dz,  &cpi->index_full);
    cvxc_mgrp_init(&cpi->ds0_g, &cpi->ds0, &cpi->index_full);
    cvxc_mgrp_init(&cpi->dz0_g, &cpi->dz0, &cpi->index_full);
    cvxc_mgrp_init(&cpi->ds2_g, &cpi->ds2, &cpi->index_full);
    cvxc_mgrp_init(&cpi->dz2_g, &cpi->dz2, &cpi->index_full);
    cvxc_mgrp_init(&cpi->ds20_g,&cpi->ds20,&cpi->index_full);
    cvxc_mgrp_init(&cpi->dz20_g,&cpi->dz20,&cpi->index_full);
    cvxc_mgrp_init(&cpi->rz_g,  &cpi->rz,  &cpi->index_full);
    cvxc_mgrp_init(&cpi->newz_g,&cpi->newz,&cpi->index_full);
    cvxc_mgrp_init(&cpi->news_g,&cpi->news,&cpi->index_full);
    cvxc_mgrp_init(&cpi->ws_g,  &cpi->ws,  &cpi->index_full);
    cvxc_mgrp_init(&cpi->wz_g,  &cpi->wz,  &cpi->index_full);
    cvxc_mgrp_init(&cpi->ws2_g, &cpi->ws2, &cpi->index_full);
    cvxc_mgrp_init(&cpi->wz2_g, &cpi->wz2, &cpi->index_full);
    cvxc_mgrp_init(&cpi->ws3_g, &cpi->ws3, &cpi->index_full);
    cvxc_mgrp_init(&cpi->wz3_g, &cpi->wz3, &cpi->index_full);

    cvxc_mgrp_init(&cpi->newrz_g,  &cpi->newrz0,&cpi->index_full);

    cvxc_mgrp_init(&cpi->sigs_g,   &cpi->sigs,   &cpi->index_sig);
    cvxc_mgrp_init(&cpi->sigz_g,   &cpi->sigz,   &cpi->index_sig);

    cvxc_mgrp_init(&cpi->lmbda_g,   &cpi->lmbda,   &cpi->index_diag);
    cvxc_mgrp_init(&cpi->lmbdasq_g, &cpi->lmbdasq, &cpi->index_diag);

    return offset;
}

int cvxc_cpl_isok(const cvxc_matrix_t *c,
                 const cvxc_convex_program_t *F,
                 const cvxc_matrix_t *G,
                 const cvxc_matrix_t *h,
                 const cvxc_matrix_t *A,
                 const cvxc_matrix_t *b,
                 const cvxc_dimset_t *dims)
{
    cvxc_size_t mc, nc, mG, nG, mA, nA, mh, nh, mb, nb;

    mc = mG = nG = mA = nA = mh = mb = 0;
    nc = nh = nb = 1;
    if (c)
        cvxm_size(&mc, &nc, c);
    if (G)
        cvxm_size(&mG, &nG, G);
    if (A)
        cvxm_size(&mA, &nA, A);
    if (h)
        cvxm_size(&mh, &nh, h);
    if (b)
        cvxm_size(&mb, &nb, b);

    if (nc > 1 || (c && mc < 1)) {
        return CVXC_ERR_DIMC;
    }

    cvxc_size_t cdim      = cvxc_dimset_sum_squared(dims, CVXDIM_CONELP);
    cvxc_size_t cdim_pckd = cvxc_dimset_sum_packed(dims, CVXDIM_CONELP);

    if (nh > 1 || mh != cdim) {
        return CVXC_ERR_DIMH;
    }

    if (mG != cdim || (c && nG != mc)) {
        return CVXC_ERR_DIMG;
    }
    if (nA != nG || mA != mb) {
        return CVXC_ERR_DIMA;
    }
    if (nb != 1) {
        return CVXC_ERR_DIMB;
    }
    if ( mb > mc || mb + cdim_pckd < mc) {
        return CVXC_ERR_RANK;
    }
    return 0;
}

int cvxc_cpl_setvars(cvxc_problem_t *cp,
                    cvxc_convex_program_t *F,
                    cvxc_size_t n, cvxc_size_t m,
                    cvxc_matrix_t *c,
                    cvxc_matrix_t *G,
                    cvxc_matrix_t *h,
                    cvxc_matrix_t *A,
                    cvxc_matrix_t *b,
                    const cvxc_dimset_t *dims,
                    cvxc_kktsolver_t *kktsolver)
{
    cp->c = c;
    cp->G = G;
    cp->h = h;
    cp->A = A;
    cp->b = b;
    cp->F = F;

    return 0;
}

/**
 * @brief Allocate space for CPL internal variables.
 *
 * @param cp
 *   On exit cp->u.space is pointer to allocated memory block and cp->mlen number of
 *   reserved bytes.
 * @param nl
 *   If non-zero then problem target function is convex otherwise linear.
 * @param n
 *   Number of variables
 * @param m
 *   Number of equality constraints
 * @param dims
 *   Structure and sizes of unequality constraints.
 *
 * @retval 0  Failure
 * @retval >0 Number of bytes used for internal variables.
 */
cvxc_size_t cvxc_cpl_allocate(cvxc_problem_t *cp,
                            int nl,
                            cvxc_size_t n,
                            cvxc_size_t m,
                            cvxc_size_t extra,
                            const cvxc_dimset_t *dims)
{
    cvxc_size_t used, nbytes = cvxc_cpl_bytes(n, m, dims, nl);
    if (extra)
        nbytes += extra;
    nbytes += sizeof(cvxc_cpl_internal_t);

    cp->u.space = calloc(nbytes, 1);
    if (!cp->u.space) {
        cp->error = CVXC_ERR_MEMORY;
        return 0;
    }
    cp->nbytes = nbytes;

    unsigned char *memory = &cp->u.space[sizeof(cvxc_cpl_internal_t)];
    if ((used = cvxc_cpl_make(cp, n, m, dims, nl, memory, nbytes)) == 0) {
        cp->error = CVXC_ERR_MEMORY;
        free(cp->u.space);
        cp->u.space = 0;
        return 0;
    }
    cp->work = &cp->u.cpl->work;
    return used + sizeof(cvxc_cpl_internal_t);
}

/**
 * @brief Setup CPL problem.
 *
 * @retval 0  Failure
 * @retval >0 Number of bytes used for internal variables.
 */
static
cvxc_size_t cvxc_cpl_create(
    cvxc_problem_t *cp,
    cvxc_convex_program_t *F,
    cvxc_matrix_t *c,
    cvxc_umatrix_t *Gf,
    cvxc_matrix_t *h,
    cvxc_umatrix_t *Af,
    cvxc_matrix_t *b,
    const cvxc_dimset_t *dims,
    cvxc_kktsolver_t *kktsolver)
{
    cvxc_size_t mc, nc, mb, nb;
    cvxc_size_t used;

    mc = nc = mb = nb = 0;
    cvxm_size(&mc, &nc, c);
    if (b)
        cvxm_size(&mb, &nb, b);

    if ((used = cvxc_cpl_allocate(cp, 0, mc, mb, 0, dims)) == 0) {
        return 0;
    }
    // cvxc_cpl_setvars(cp, F, mc, mb, c, G, h, A, b, dims, kktsolver);
    cp->c = c;
    cp->F = F;
    cp->Gf = Gf;
    cp->h = h;
    cp->Af = Af;
    cp->b = b;

    cvxc_cpl_internal_t *cpi = cp->u.cpl;
    // provide full index set
    cp->index_g = &cpi->index_full;
    cvxc_mgrp_init(&cpi->h_g, cp->h, cp->index_g);

    // init KKT solver
    if (kktsolver) {
        cp->solver = kktsolver;
    } else {
        cvxc_ldlsolver_init(&cp->__S, cp, mc, mb, dims);
        cp->solver = &cp->__S;
    }
    cvxc_cpl_solver_init(&cpi->cp_solver, cp, cp->solver);
    cp->solver = (cvxc_kktsolver_t *)&cpi->cp_solver;

    return used;
}

cvxc_size_t cvxc_cpl_setup_user(
    cvxc_problem_t *cp,
    cvxc_convex_program_t *F,
    cvxc_matrix_t *c,
    cvxc_umatrix_t *Gf,
    cvxc_matrix_t *h,
    cvxc_umatrix_t *Af,
    cvxc_matrix_t *b,
    const cvxc_dimset_t *dims,
    cvxc_kktsolver_t *kktsolver)
{
    int err;
    if (! cp)
        return 0;

    if ((err = cvxc_cpl_isok(c, F, __cvxnil, h, __cvxnil, b, dims)) != 0) {
        cp->error = err;
        return 0;
    }

    cp->A = __cvxnil;
    cp->G = __cvxnil;
    int stat = cvxc_cpl_create(cp, F, c, &cp->Gu, h, &cp->Au, b, dims, kktsolver);
    return stat;
}

cvxc_size_t cvxc_cpl_setup(
    cvxc_problem_t *cp,
    cvxc_convex_program_t *F,
    cvxc_matrix_t *c,
    cvxc_matrix_t *G,
    cvxc_matrix_t *h,
    cvxc_matrix_t *A,
    cvxc_matrix_t *b,
    const cvxc_dimset_t *dims,
    cvxc_kktsolver_t *kktsolver)
{
    int err;
    if (! cp)
        return 0;

    if ((err = cvxc_cpl_isok(c, F, G, h, A, b, dims)) != 0) {
        cp->error = err;
        return 0;
    }

    cvxc_umat_make(&cp->Au, A);
    cvxc_umat_make(&cp->Gu, G);
    cp->A = A;
    cp->G = G;
    int stat = cvxc_cpl_create(cp, F, c, &cp->Gu, h, &cp->Au, b, dims, kktsolver);
    if (stat <= 0) {
        cvxc_umat_clear(&cp->Au);
        cvxc_umat_clear(&cp->Gu);
    }
    return stat;
}

/**
 * @brief Release resources reserved to program.
 */
void cvxc_cpl_release(cvxc_problem_t *cp)
{
    if (!cp)
        return;
    if (cp->u.space)
        free(cp->u.space);
    cp->u.space = 0;
    cp->nbytes = 0;
    cp->work = 0;

    cvxc_kktrelease(cp->solver);
}

int cvxc_cpl_ready(cvxc_solution_t *sol, cvxc_problem_t *cp, int iter, int stat)
{
    cvxc_cpl_internal_t *cpi = cp->u.cpl;

    if (stat == CVXC_STAT_OPTIMAL && iter == -1) {
        // constructed initial point is feasible and optimal
        cvxc_mksymm(&cpi->s_g);
        cvxc_mksymm(&cpi->z_g);

        // rx = A'*y + G'z + c
        cvxm_copy(&cpi->rx, cp->c, CVXC_ALL);
        cvxm_mvmult(1.0, &cpi->rx, 1.0, cp->A, &cpi->y, CVXC_TRANS);
        cvxm_mvmult(1.0, &cpi->rx, 1.0, cp->G, &cpi->z, CVXC_TRANS);
        cpi->resx = cvxm_nrm2(&cpi->rx);

        // ry = b - A*x  ; TODO - computes -b - A*x ;; check
        cvxm_copy(&cpi->ry, cp->b, CVXC_ALL);
        cvxm_mvmult(-1.0, &cpi->ry, -1.0, cp->A, &cpi->x, CVXC_TRANS);
        cpi->resy = cvxm_nrm2(&cpi->ry);

        // rz = s + G*x - h
        cvxm_copy(&cpi->rz, &cpi->s, CVXC_ALL);
        cvxm_mvmult(1.0, &cpi->rz, 1.0, cp->G, &cpi->x, 0);
        cvxm_axpy(&cpi->rz, 1.0, cp->h);
        cpi->resz = cvxc_snrm2(&cpi->rz_g);

        cpi->pres = MAXF(cpi->resy/cpi->resy0, cpi->resz/cpi->resz0);
        cpi->dres = cpi->resx/cpi->resx0;
        cpi->cx = cvxm_dot(cp->c, &cpi->x);
        cpi->by = cvxm_dot(cp->b, &cpi->y);
        cpi->hz = cvxc_sdot(&cpi->h_g, &cpi->z_g);

        sol->x = &cpi->x;
        sol->s = &cpi->s;
        sol->y = &cpi->y;
        sol->z = &cpi->z;

        sol->status = stat;
        sol->gap = cpi->gap;
        sol->relative_gap = cpi->relgap;
        sol->primal_objective = cpi->cx;
        sol->dual_objective = - (cpi->by + cpi->hz);
        sol->primal_infeasibility = cpi->pres;
        sol->dual_infeasibility = cpi->dres;
        sol->primal_slack = - cpi->ts;
        sol->dual_slack = - cpi->tz;
        sol->primal_residual_cert = __NaN();
        sol->dual_residual_cert = __NaN();
        sol->iterations = 0;
    }
    else if (stat == CVXC_STAT_UNKNOWN || stat == CVXC_STAT_OPTIMAL) {
        cvxm_scale(&cpi->x, 1.0/cpi->tau, CVXC_ALL);
        cvxm_scale(&cpi->y, 1.0/cpi->tau, CVXC_ALL);
        cvxm_scale(&cpi->s, 1.0/cpi->tau, CVXC_ALL);
        cvxm_scale(&cpi->z, 1.0/cpi->tau, CVXC_ALL);
        cvxc_mksymm(&cpi->s_g);
        cvxc_mksymm(&cpi->z_g);
        cpi->ts = cvxc_max_step(&cpi->s_g, __nilgrp, &cpi->work);
        cpi->tz = cvxc_max_step(&cpi->z_g, __nilgrp, &cpi->work);

        cp->error = stat == CVXC_STAT_UNKNOWN ? CVXC_ERR_MAXITER : 0;
        sol->x = &cpi->x;
        sol->s = &cpi->s;
        sol->y = &cpi->y;
        sol->z = &cpi->z;

        sol->status = stat;

        sol->gap = cpi->gap;
        sol->relative_gap = cpi->relgap;
        sol->primal_objective = cpi->cx;
        sol->dual_objective = - (cpi->by + cpi->hz);
        sol->primal_infeasibility = cpi->pres;
        sol->dual_infeasibility = cpi->dres;
        sol->primal_slack = - cpi->ts;
        sol->dual_slack = - cpi->tz;
        if (stat == CVXC_STAT_OPTIMAL) {
            sol->primal_residual_cert = __NaN();
            sol->dual_residual_cert = __NaN();
        } else {
            sol->primal_residual_cert = cpi->pinfres;
            sol->dual_residual_cert = cpi->dinfres;
        }
        sol->iterations = iter;
    }
    return -stat;
}

static
void cvxc_cpl_save_or_restore(cvxc_problem_t *cp, int restore)
{
    cvxc_cpl_internal_t *cpi = cp->u.cpl;
    //cvxc_stats_t *stats = &cp->stats;

    if (restore) {
        cvxc_scaling_copy(&cpi->W, &cpi->W0);
        cvxm_copy(&cpi->x, &cpi->x0, 0);
        cvxm_copy(&cpi->y, &cpi->y0, 0);
        cvxm_copy(&cpi->s, &cpi->s0, 0);
        cvxm_copy(&cpi->z, &cpi->z0, 0);
        cvxm_copy(&cpi->lmbda, &cpi->lmbda0, 0);
        cvxm_copy(&cpi->lmbdasq, &cpi->lmbdasq0, 0);
        cvxm_copy(&cpi->rx, &cpi->rx0, 0);
        cvxm_copy(&cpi->ry, &cpi->ry0, 0);
        cvxm_copy(&cpi->rz, &cpi->rz0, 0);
        cpi->phi   = cpi->phi0;
        cpi->dphi  = cpi->dphi0;
        cpi->sigma = cpi->sigma0;
        cpi->eta   = cpi->eta0;
        cpi->step  = cpi->step0;
        cpi->gap   = cpi->gap0;
        return;
    }
    // save
    cvxc_scaling_copy(&cpi->W0, &cpi->W);
    cvxm_copy(&cpi->x0, &cpi->x, 0);
    cvxm_copy(&cpi->y0, &cpi->y, 0);
    cvxm_copy(&cpi->s0, &cpi->s, 0);
    cvxm_copy(&cpi->z0, &cpi->z, 0);
    cvxm_copy(&cpi->lmbda0, &cpi->lmbda, 0);
    cvxm_copy(&cpi->lmbdasq0, &cpi->lmbdasq, 0);
    cvxm_copy(&cpi->rx0, &cpi->rx, 0);
    cvxm_copy(&cpi->ry0, &cpi->ry, 0);
    cvxm_copy(&cpi->rz0, &cpi->rz, 0);
    cpi->phi0   = cpi->phi;
    cpi->dphi0  = cpi->dphi;
    cpi->sigma0 = cpi->sigma;
    cpi->eta0   = cpi->eta;
    cpi->step0  = cpi->step;
    cpi->gap0   = cpi->gap;
}

static
int cvxc_cpl_linesearch(cvxc_problem_t *cp,
                       int iter,
                       int relaxed_iters)
{
    cvxc_cpl_internal_t *cpi = cp->u.cpl;
    int backtrack;
    cvxc_matrix_t news_nl, newrznl;
    cvxc_float_t newresznl, newgap, newphi, newresx;

    // NOTE: sigma changes!!! parameters??

    // Merit function
    //   phi = theta1 * gap + theta2 * norm(rx) + theta3 * norm(rznl)
    // and its directional derivative dphi.
    cpi->phi = cpi->theta1*cpi->gap + cpi->theta2*cpi->resx + cpi->theta3*cpi->resznl;
    if (iter == 0)
        cpi->dphi = -cpi->phi;
    else
        cpi->dphi =
            - cpi->theta1*(1.0 - cpi->sigma)*cpi->gap
            - cpi->theta2*(1.0 - cpi->eta)*cpi->resx
            - cpi->theta3*(1.0 - cpi->eta)*cpi->resznl;

    for (backtrack = 1; backtrack; backtrack = 0 ) {
        cvxm_copy(&cpi->newx, &cpi->x, 0);
        cvxm_axpy(&cpi->newx, cpi->step, &cpi->dx);
        cvxm_copy(&cpi->newy, &cpi->y, 0);
        cvxm_axpy(&cpi->newy, cpi->step, &cpi->dy);
        cvxm_copy(&cpi->newz, &cpi->z, 0);
        cvxm_axpy(&cpi->newz, cpi->step, &cpi->dz2);
        cvxm_copy(&cpi->news, &cpi->s, 0);
        cvxm_axpy(&cpi->news, cpi->step, &cpi->ds2);

        F1(cp->F, &cpi->newf, &cpi->newDf, &cpi->newx);
        // newrx = c + A'*newy + newDf'*newz[:mnl] + G'*newz[mnl:]
        // ...
        cvxm_copy(&cpi->newrx, cp->c, 0);  // TODO: c == __nil ??
        cvxm_mvmult(1.0, &cpi->newrx, 1.0, cp->A, &cpi->newy, CVXC_TRANS);
        cvxc_sgemv2(1.0, &cpi->newrx, 1.0, &cpi->newDf, cp->G, &cpi->newz_g, CVXC_TRANS);
        newresx = SQRT(cvxm_dot(&cpi->newrx, &cpi->newrx));

        // newrznl = news[:mnl] + newf
        // ...
        cvxc_mgrp_elem(&news_nl, &cpi->news_g, CVXDIM_CONVEX, 0);
        cvxc_mgrp_elem(&newrznl, &cpi->newrz_g, CVXDIM_CONVEX, 0);
        cvxm_copy(&newrznl, &news_nl, 0);
        cvxm_axpy(&newrznl, 1.0, &cpi->newf);
        newresznl = cvxm_nrm2(&newrznl);

        newgap = (1.0 - (1.0 - cpi->sigma)*cpi->step)*cpi->gap +
            cpi->step*cpi->step*cpi->dsdz;
        newphi = cpi->theta1*newgap + cpi->theta2*newresx + cpi->theta3*newresznl;

        if (iter == 0) {
            if (newgap <= (1.0 - ALPHA*cpi->step)*cpi->gap &&
                ((relaxed_iters > 0 && relaxed_iters < MAX_RELAXED_ITERS)
                 || newphi <= cpi->phi+ALPHA*cpi->step*cpi->dphi)) {
                backtrack = 0;
                cpi->sigma = MINF(newgap/cpi->gap, POW((newgap/cpi->gap), EXPON));
                cpi->eta   = 0.0;
            }
            else {
                cpi->step *= BETA;
            }
        }
        else {
            if (relaxed_iters == -1 || (relaxed_iters == 0 && MAX_RELAXED_ITERS == 0)) {
                // standard line search
                if (newphi <= cpi->phi + ALPHA*cpi->step*cpi->dphi) {
                    relaxed_iters = 0;
                    backtrack = 0;
                }
                else {
                    cpi->step *= BETA;
                }
            }
            else if (relaxed_iters == 0 && relaxed_iters < MAX_RELAXED_ITERS) {
                if (newphi <= cpi->phi + ALPHA*cpi->step*cpi->dphi) {
                    relaxed_iters = 0;
                }
                else {
                    // save state
                    cvxc_cpl_save_or_restore(cp, 0);
                    relaxed_iters = 1;
                }
                backtrack = 0;
            }
            else if (relaxed_iters >= 0 && relaxed_iters < MAX_RELAXED_ITERS &&
                     MAX_RELAXED_ITERS > 0) {
                if (newphi <= cpi->phi0 + ALPHA*cpi->step0*cpi->dphi0) {
                    relaxed_iters = 0;
                }
                else {
                    relaxed_iters += 1;
                }
                backtrack = 0;
            }
            else if (relaxed_iters == MAX_RELAXED_ITERS && MAX_RELAXED_ITERS > 0) {
                if (newphi <= cpi->phi0 + ALPHA*cpi->step0*cpi->dphi0) {
                    relaxed_iters = 0;
                    backtrack = 0;
                }
                else if (newphi >= cpi->phi0) {
                    // Resume last saved line search
                    cvxc_cpl_save_or_restore(cp, 1);
                    relaxed_iters = -1;
                }
                else if (newphi <= cpi->phi + ALPHA*cpi->step*cpi->dphi) {
                    backtrack = 0;
                    relaxed_iters = -1;
                }
            }
        }
    }

    return relaxed_iters;
}

int cvxc_cpl_compute_start(cvxc_problem_t *cp)
{
    cvxc_cpl_internal_t *cpi = cp->u.cpl;
    F0(cp->F, &cpi->x0);
    cvxc_mgrp_initial_value(&cpi->z_g, 0);
    cvxc_mgrp_initial_value(&cpi->s_g, 0);
    cvxc_scaling_initial_value(&cpi->W);
    cvxm_copy(&cpi->x,    &cpi->x0, 0);
    cvxm_copy(&cpi->rx,   &cpi->x0, 0);
    cvxm_copy(&cpi->dx,   &cpi->x0, 0);
    cvxm_copy(&cpi->rx0,  &cpi->x0, 0);
    cvxm_copy(&cpi->dx0,  &cpi->x0, 0);
    cvxm_copy(&cpi->newx, &cpi->x0, 0);
    cvxm_copy(&cpi->newrx,&cpi->x0, 0);
    cvxm_copy(&cpi->ry,    cp->b,   0);
    return 0;
}

/**
 * @brief Solve convex problem with linear objective.
 */
int cvxc_cpl_solve(cvxc_solution_t *sol, cvxc_problem_t *cp, cvxc_solopts_t *opts)
{
    cvxc_cpl_internal_t *cpi = cp->u.cpl;

    cvxc_float_t cx, dy, dz, dr;
    cvxc_matrix_t lk, sk, zk, ls, lz, s_l, s_nl, z_l, z_nl, rz_nl, rz_l;

    cvxc_size_t range;
    cvxc_float_t abstol = CVXC_ABSTOL;
    cvxc_float_t reltol = CVXC_RELTOL;
    cvxc_float_t feastol = CVXC_FEASTOL;

    int maxiter = opts->max_iter > 0 ? opts->max_iter : CVXC_MAXITER;
    int refinement = opts->refinement > 0 ? opts->refinement : 0;

    refinement = opts->refinement == 0 &&
        (cvxc_index_count(&cpi->index_full, CVXDIM_SDP) > 0 ||
         cvxc_index_count(&cpi->index_full, CVXDIM_SOCP) > 0);

    cp->error = 0;

    cpi->tau = 1.0;
    cpi->kappa = 1.0;
    cpi->wkappa3 = 0.0;
    cvxm_copy(&cpi->dx,   &cpi->x, 0);
    cvxm_copy(&cpi->rx,   &cpi->x, 0);
    cvxm_copy(&cpi->rx0,  &cpi->x, 0);
    cvxm_copy(&cpi->hrx, cp->c, 0);
    cvxm_copy(&cpi->ry,  cp->b, 0);
    cvxm_copy(&cpi->hry, cp->b, 0);

    cpi->gap = cvxc_sdot(&cpi->s_g, &cpi->z_g);

    range  = cvxc_index_count(&cpi->index_full, CVXDIM_SDP);
    range += cvxc_index_length(&cpi->index_full, CVXDIM_NONLINEAR);
    range += cvxc_index_length(&cpi->index_full, CVXDIM_LINEAR);
    range += cvxc_index_length(&cpi->index_full, CVXDIM_SOCP);
    if (range == 0)
        range = 1;

    // make non-linear and linear mappings
    cvxc_mgrp_elem(&s_nl,  &cpi->s_g,  CVXDIM_NONLINEAR|CVXDIM_NLTARGET, 0);
    cvxc_mgrp_elem(&s_l,   &cpi->s_g,  CVXDIM_CONELP, 0);
    cvxc_mgrp_elem(&z_nl,  &cpi->z_g,  CVXDIM_NONLINEAR|CVXDIM_NLTARGET, 0);
    cvxc_mgrp_elem(&z_l,   &cpi->z_g,  CVXDIM_CONELP, 0);
    cvxc_mgrp_elem(&rz_nl, &cpi->rz_g, CVXDIM_NONLINEAR|CVXDIM_NLTARGET, 0);
    cvxc_mgrp_elem(&rz_l,  &cpi->rz_g, CVXDIM_CONELP, 0);

    // -----------------------------------------------------------------------------

    for (int iter = 0; iter < maxiter; iter++) {

        if (refinement != 0) {
            F2(cp->F, &cpi->f, &cpi->Df, &cpi->H, &cpi->x, &z_nl);
        }
        else  {
            F1(cp->F, &cpi->f, &cpi->Df, &cpi->x);
        }

        cpi->gap = cvxc_sdot(&cpi->s_g, &cpi->z_g);

        // rx = c + A'*y + Df'*z[:mnl] + G'*z[mnl:]
        cvxm_copy(&cpi->rx, cp->c, 0);
        cvxm_mvmult(1.0, &cpi->rx, 1.0, cp->A, &cpi->y, CVXC_TRANS);
        cvxc_sgemv2(1.0, &cpi->rx, 1.0, &cpi->Df, cp->G, &cpi->z_g, CVXC_TRANS);
        cpi->resx = cvxm_nrm2(&cpi->rx);

        // ry = A*x - b
        cvxm_copy(&cpi->ry, cp->b, 0);
        cvxm_mvmult(-1.0, &cpi->ry, 1.0, cp->A, &cpi->x, 0);
        cpi->resy = cvxm_nrm2(&cpi->ry);

        // rz_nl = s_nl + f
        cvxm_copy(&rz_nl, &s_nl, 0);
        cvxm_axpy(&rz_nl, 1.0, &cpi->f);
        cpi->resznl = cvxm_nrm2(&rz_nl);

        // rz_l = s_l + G*x - h
        cvxm_copy(&rz_l, &s_l, 0);
        cvxm_axpy(&rz_l, -1.0, cp->h);
        cvxm_mvmult(1.0, &rz_l, 1.0, cp->G, &cpi->x, 0);
        cpi->reszl = cvxc_snrm2_elem(&cpi->rz_g, CVXDIM_CONELP);

        // Statistics for stopping criteria
        // pcost = c'*x
        // dcost = c'*x + y'*(A*x-b) + znl'*f(x) + zl'*(G*x-h)
        //       = c'*x + y'*(A*x-b) + znl'*(f(x)+snl) + zl'*(G*x-h+sl)
        //         - z'*s
        //       = c'*x + y'*ry + znl'*rznl + zl'*rzl - gap
        cx = cvxm_dot(cp->c, &cpi->x);
        dy = cvxm_dot(&cpi->y, &cpi->ry);
        dz = cvxm_dot(&z_nl, &rz_nl);
        dr = cvxc_sdot_elem(&cpi->z_g, &cpi->rz_g, CVXDIM_CONELP);

        cpi->pcost = cx;
        cpi->dcost = cx + dy + dz +dr - cpi->gap;

        // statistics for stopping
        if (cpi->pcost < 0.0) {
            cpi->relgap = cpi->gap / -cpi->pcost;
        } else if (cpi->dcost > 0.0) {
            cpi->relgap = cpi->gap / cpi->dcost;
        } else {
            cpi->relgap = __NaN();
        }

        cpi->pres = cpi->resy*cpi->resy + cpi->reszl*cpi->reszl + cpi->resznl*cpi->resznl;
        cpi->pres = SQRT(cpi->pres);

        cpi->dres = cpi->resx;
        if (iter == 0) {
            cpi->resx0 = MAXF(1.0, cpi->resx);
            cpi->resznl0 = MAXF(1.0, cpi->resznl);
            cpi->pres0 = MAXF(1.0, cpi->pres);
            cpi->dres0 = MAXF(1.0, cpi->dres);
            cpi->gap0  = cpi->gap;
            cpi->theta1 = 1.0/cpi->gap0;
            cpi->theta2 = 1.0/cpi->resx0;
            cpi->theta3 = 1.0/cpi->resznl0;
        }
        cpi->phi = cpi->theta1*cpi->gap + cpi->theta2*cpi->resx + cpi->theta3*cpi->resznl;
        cpi->pres = cpi->pres / cpi->pres0;
        cpi->dres = cpi->dres / cpi->dres0;

        if (opts->show_progress > 0) {
            if (iter == 0) {
                fprintf(stderr, "%10s %11s %9s %8s %7s\n",
                    "pcost", "dcost", "gap", "pres", "dres");
            }
            fprintf(stderr, "%2d: %11.4e %11.4e %6.0e %7.0e %7.0e\n",
                    iter, cpi->pcost, cpi->dcost, cpi->gap, cpi->pres, cpi->dres);
        }
        // ---------------------------------------------------------------------
        // test for stopping criteria

        if (cpi->pres <= feastol &&
            cpi->dres <= feastol &&
            (cpi->gap <= abstol ||
             (!isnan(cpi->relgap) && cpi->relgap < reltol))) {

            return cvxc_cpl_ready(sol, cp, iter, CVXC_STAT_OPTIMAL);
        }

        // -----------------------------------------------------------------------
        // Compute initial scaling W:
        //
        //     W * z = W^{-T} * s = lambda
        //     dg * tau = 1/dg * kappa = lambdag.
        if (iter == 0) {
            cvxc_compute_scaling(&cpi->W, &cpi->s_g, &cpi->z_g, &cpi->lmbda_g, &cpi->work);
        }

        // lmdasq = lmda o lmbda
        cvxc_ssqr(&cpi->lmbdasq_g, &cpi->lmbda_g);

        // f3(x, y, z) solves
        //
        //     [ H   A'  GG'*W^{-1} ] [ ux ]   [ bx ]
        //     [ A   0   0          ] [ uy ] = [ by ].
        //     [ GG  0  -W'         ] [ uz ]   [ bz ]
        //
        // On entry, x, y, z contain bx, by, bz.
        // On exit, they contain ux, uy, uz.
        int relaxed_iters = 0;
        int err = cvxc_kktfactor(cp->solver, &cpi->W, &cpi->x, &z_nl);

        if (err < 0) {
            int singular_kkt = 0;
            // singular matrix ??
            if (iter == 0) {
                // Rank(A) < p or Rank([H(x); A; Df(x); G] < n]??
                return -2;
            }
            if (relaxed_iters > 0 && relaxed_iters < MAX_RELAXED_ITERS) {
                // restore save point
                cvxc_cpl_save_or_restore(cp, CPL_RESTORE);
                cpi->resznl = cvxm_nrm2(&rz_nl);
                relaxed_iters = -1;
                F2(cp->F, __cvxnil, &cpi->Df, &cpi->H, &cpi->x, &z_nl);
                err = cvxc_kktfactor(cp->solver, &cpi->W, &cpi->H, &cpi->Df);
                if (err)
                    singular_kkt = 1;
            } else {
                singular_kkt = 1;
            }

            if (singular_kkt) {
                // terminate in singular KKT_matrix
                return cvxc_cpl_ready(sol, cp, iter, CVXC_STAT_SINGULAR);
            }
        }

        if (iter == 0) {
            if (refinement > 0 || opts->debug) {
                cvxm_copy(&cpi->wx, cp->c, 0);
                cvxm_copy(&cpi->wy, cp->b, 0);
                cvxm_copy(&cpi->wz, &cpi->z, 0);
                cvxm_copy(&cpi->ws, &cpi->s, 0);
            }
            if (refinement > 0) {
                cvxm_copy(&cpi->wx2, cp->c, 0);
                cvxm_copy(&cpi->wy2, cp->b, 0);
                cvxm_scale(&cpi->wz2, 0.0, 0);
                cvxm_scale(&cpi->ws2, 0.0, 0);
            }
        }

        cpi->sigma = __ZERO;
        cpi->eta   = __ZERO;

        for (int i = 0; i < 2; i++) {
            //cvxc_matrix_t t0, t1;
            cpi->mu = cpi->gap/((cvxc_float_t)range);
            // Solve
            //
            //     [ 0     ]   [ H  A' GG' ] [ dx        ]
            //     [ 0     ] + [ A  0  0   ] [ dy        ] = -(1 - eta)*r
            //     [ W'*ds ]   [ GG 0  0   ] [ W^{-1}*dz ]
            //
            //     lmbda o (dz + ds) = -lmbda o lmbda + sigma*mu*e.
            //
            //  for nl,l,q: ds = -1.0*lmbdasq

            cvxc_mgrp_copy_lambda(&cpi->ds_g, &cpi->lmbdasq_g);
            cvxc_mgrp_scale_sz(&cpi->ds_g, -1.0, 0);
            cvxc_mgrp_update_sz(&cpi->ds_g, cpi->sigma*cpi->mu, 0);

            cvxm_axpby(0.0, &cpi->dx, -1.0+cpi->eta, &cpi->rx);

            cvxm_axpby(0.0, &cpi->dy, -1.0+cpi->eta, &cpi->ry);
            cvxc_mgrp_axpby_sz(0.0, &cpi->dz_g, -1.0+cpi->eta, &cpi->rz_g, CVXDIM_CONVEXPROG);

            // .. computation here
            err = f4(cp, &cpi->dx, &cpi->dy, &cpi->dz_g, &cpi->ds_g, refinement);
            if (err < 0) {
                // terminated ....
                return cvxc_cpl_ready(sol, cp, iter, CVXC_STAT_SINGULAR);
            }
            // line search needs ds'*dz and unscaled steps
            cpi->dsdz = cvxc_sdot(&cpi->ds_g, &cpi->dz_g);
            cvxm_copy(&cpi->dz2, &cpi->dz, 0);
            cvxc_scale(&cpi->dz2_g, &cpi->W, CVXC_INV, &cpi->work);
            cvxm_copy(&cpi->ds2, &cpi->ds, 0);
            cvxc_scale(&cpi->ds2_g, &cpi->W, CVXC_TRANS, &cpi->work);

            // max step to boundary
            cvxc_scale2(&cpi->ds_g, &cpi->lmbda_g, 0, &cpi->work);
            cvxc_scale2(&cpi->dz_g, &cpi->lmbda_g, 0, &cpi->work);

            cpi->ts = cvxc_max_step(&cpi->ds_g, &cpi->sigs_g, &cpi->work);
            cpi->tz = cvxc_max_step(&cpi->dz_g, &cpi->sigz_g, &cpi->work);
            cvxc_float_t t = cvxc_maxvec(3, (cvxc_float_t[]){0.0, cpi->ts, cpi->tz});
            if (t == 0.0)
                cpi->step = 1.0;
            else
                cpi->step = MINF(1.0, STEP/t);

            // backtrack until newx is in domain of f
            for (int backtrack = 1; backtrack; ) {
                cvxm_copy(&cpi->newx, &cpi->x, 0);
                cvxm_axpy(&cpi->newx, cpi->step, &cpi->dx);
                err = F1(cp->F, &cpi->newf, &cpi->newDf, &cpi->newx);
                if (err == 0)
                    backtrack = 0;
                else
                    cpi->step *= BETA;
            }

            // do the line search
            relaxed_iters = cvxc_cpl_linesearch(cp, i, relaxed_iters);
        }

        // update x, y
        cvxm_axpy(&cpi->x, cpi->step, &cpi->dx);
        cvxm_axpy(&cpi->y, cpi->step, &cpi->dy);

        // Replace nonlinear, 'l' and 'q' blocks of ds and dz with the updated
        // variables in the current scaling.
        // Replace 's' blocks of ds and dz with the factors Ls, Lz in a
        // factorization Ls*Ls', Lz*Lz' of the updated variables in the
        // current scaling.
        //
        // ds := e + step*ds for 'l' and 'q' blocks.
        // dz := e + step*dz for 'l' and 'q' blocks.
        cvxc_mgrp_scale_sz(&cpi->ds_g, cpi->step, CVXDIM_CONVEX|CVXDIM_LINEAR|CVXDIM_SOCP);
        cvxc_mgrp_update_sz(&cpi->ds_g, 1.0, CVXDIM_CONVEX|CVXDIM_LINEAR|CVXDIM_SOCP);

        cvxc_mgrp_scale_sz(&cpi->dz_g, cpi->step, CVXDIM_CONVEX|CVXDIM_LINEAR|CVXDIM_SOCP);
        cvxc_mgrp_update_sz(&cpi->dz_g, 1.0, CVXDIM_CONVEX|CVXDIM_LINEAR|CVXDIM_SOCP);

        // ds := H(lambda)^{-1/2} * ds and dz := H(lambda)^{-1/2} * dz.
        //
        // This replaces the 'l' and 'q' components of ds and dz with the
        // updated variables in the current scaling.
        // The 's' components of ds and dz are replaced with
        //
        // diag(lmbda_k)^{1/2} * Qs * diag(lmbda_k)^{1/2}
        // diag(lmbda_k)^{1/2} * Qz * diag(lmbda_k)^{1/2}

        cvxc_scale2(&cpi->ds_g, &cpi->lmbda_g, CVXC_INV, &cpi->work);
        cvxc_scale2(&cpi->dz_g, &cpi->lmbda_g, CVXC_INV, &cpi->work);

        // sigs := ( e + step*sigs ) ./ lambda for 's' blocks.
        // sigz := ( e + step*sigz ) ./ lambda for 's' blocks.
        cvxm_scale(&cpi->sigs, cpi->step, CVXC_ALL);
        cvxm_scale(&cpi->sigz, cpi->step, CVXC_ALL);
        cvxm_add(&cpi->sigs, 1.0, CVXC_ALL);
        cvxm_add(&cpi->sigz, 1.0, CVXC_ALL);

        for (int k = 0; k < cvxc_mgrp_count(&cpi->lmbda_g, CVXDIM_SDP); k++) {
            cvxc_mgrp_elem(&lk, &cpi->lmbda_g, CVXDIM_SDP, k);
            // sigs ./ lmbda
            cvxc_mgrp_elem(&sk, &cpi->sigs_g, CVXDIM_SDP, k);
            cvxm_solve_diag(&sk, 1.0, &lk, CVXC_RIGHT);
            // sigz ./ lmbda
            cvxc_mgrp_elem(&sk, &cpi->sigz_g, CVXDIM_SDP, k);
            cvxm_solve_diag(&sk, 1.0, &lk, CVXC_RIGHT);
        }

        // divide ds, dz by lmbda S blocks;;
        for (int k = 0; k < cvxc_mgrp_count(&cpi->ds_g, CVXDIM_SDP); k++) {
            cvxc_float_t a;
            int m = cvxc_mgrp_elem(&sk, &cpi->ds_g, CVXDIM_SDP, k);
            cvxc_mgrp_elem(&zk, &cpi->dz_g, CVXDIM_SDP, k);
            cvxc_mgrp_elem(&ls, &cpi->sigs_g, CVXDIM_SDP, k);
            cvxc_mgrp_elem(&lz, &cpi->sigz_g, CVXDIM_SDP, k);

            for (int i = 0; i < m; i++) {
                cvxm_view_map(&lk, &sk, 0, i, m, 1);
                a = SQRT(cvxm_get(&ls, i, 0));
                cvxm_scale(&lk, a, 0);

                cvxm_view_map(&lk, &zk, 0, i, m, 1);
                a = SQRT(cvxm_get(&lz, i, 0));
                cvxm_scale(&lk, a, 0);
            }
        }
        cvxc_update_scaling(&cpi->W, &cpi->lmbda_g, &cpi->ds_g, &cpi->dz_g, &cpi->work);

        // Unscale s, z, tau, kappa (unscaled variables are used only to
        // compute feasibility residuals).
        cvxc_mgrp_copy_lambda(&cpi->s_g, &cpi->lmbda_g);
        cvxc_scale(&cpi->s_g, &cpi->W, CVXC_TRANS, &cpi->work);

        cvxc_mgrp_copy_lambda(&cpi->z_g, &cpi->lmbda_g);
        cvxc_scale(&cpi->z_g, &cpi->W, CVXC_INV, &cpi->work);

        cpi->gap = cvxm_dot(&cpi->lmbda, &cpi->lmbda);
    }

    return cvxc_cpl_ready(sol, cp, maxiter, CVXC_STAT_UNKNOWN);
}
