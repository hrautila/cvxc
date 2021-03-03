/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */

#ifndef CVXC_INTERNAL_H
#define CVXC_INTERNAL_H

#include "cvxc.h"

typedef struct cvxc_conelp_internal {
    int primalstart;
    int dualstart;
    cvxc_float_t tau, kappa;
    cvxc_float_t dkappa, dtau;
    cvxc_float_t wkappa3;
    cvxc_float_t nrms, nrmz;
    cvxc_float_t dg, dgi;
    cvxc_float_t ts, tz;                 // step sizes

    // statistics
    cvxc_float_t resx0;
    cvxc_float_t resy0;
    cvxc_float_t resz0;
    cvxc_float_t resx;
    cvxc_float_t resy;
    cvxc_float_t resz;
    cvxc_float_t hresx;
    cvxc_float_t hresy;
    cvxc_float_t hresz;
    cvxc_float_t cx;
    cvxc_float_t by;
    cvxc_float_t hz;
    cvxc_float_t rt;
    cvxc_float_t dres;
    cvxc_float_t pres;
    cvxc_float_t dinfres;
    cvxc_float_t pinfres;
    cvxc_float_t gap;
    cvxc_float_t relgap;
    cvxc_float_t pcost;
    cvxc_float_t dcost;

    // internal result
    cvxc_matrix_t x, y, s, z;

    cvxc_matrix_t dx, dy, ds, dz;
    cvxc_matrix_t x1, y1, z1;
    cvxc_matrix_t rx, ry, rz;
    cvxc_matrix_t hrx, hry, hrz;
    cvxc_matrix_t sigs, sigz;
    cvxc_matrix_t lmbda, lmbdasq;

    cvxc_matrix_t wx, wy, ws, wz;
    cvxc_matrix_t wx2, wy2, ws2, wz2;
    cvxc_matrix_t ws3, wz3;
    cvxc_matrix_t th;

    // matrix groups with correct indexing
    cvxc_matgrp_t h_g;
    cvxc_matgrp_t s_g, z_g;
    cvxc_matgrp_t ds_g, dz_g;
    cvxc_matgrp_t hrz_g;
    cvxc_matgrp_t rz_g;
    cvxc_matgrp_t ws_g, wz_g;
    cvxc_matgrp_t ws2_g, wz2_g;
    cvxc_matgrp_t ws3_g, wz3_g;
    cvxc_matgrp_t th_g;
    cvxc_matgrp_t z1_g;
    cvxc_matgrp_t lmbda_g, lmbdasq_g;
    cvxc_matgrp_t sigs_g, sigz_g;

    cvxc_index_t index_full;             // indexing to full matrix group
    cvxc_index_t index_packed;           // indexing to matrix group with packed 'S' space
    cvxc_index_t index_diag;             // indexing to matrix group with diagonal 'S' space
    cvxc_index_t index_sig;              // indexing to matrix group with only diagonal 'S" space

    cvxc_scaling_t W;                    // scaling matrix group
    cvxc_memblk_t work;                  // workspace

} cvxc_conelp_internal_t;

typedef struct cvxc_gp_params {
    cvxc_matrix_t *F;
    cvxc_matrix_t *g;
    cvxc_gpindex_t *gpi;
    cvxc_matrix_t y;
    cvxc_matrix_t u;
    cvxc_matrix_t Fs;
} cvxc_gp_params_t;

#if 0
typedef struct cvxc_gp_program {
    cvxc_convex_program_t gp;
    cvxc_gp_params_t gp_params;
} cvxc_gp_program_t;
#endif

typedef struct cvxc_cpl_internal {
    cvxc_float_t tau, kappa;
    cvxc_float_t dkappa, dtau;
    cvxc_float_t wkappa3;
    cvxc_float_t nrms, nrmz;
    cvxc_float_t dg, dgi;
    cvxc_float_t ts, tz;                 // step sizes
    cvxc_float_t phi;
    cvxc_float_t dphi;
    cvxc_float_t phi0;
    cvxc_float_t dphi0;
    cvxc_float_t mu;
    cvxc_float_t sigma;
    cvxc_float_t sigma0;
    cvxc_float_t eta;
    cvxc_float_t eta0;
    cvxc_float_t theta1;
    cvxc_float_t theta2;
    cvxc_float_t theta3;
    cvxc_float_t step;
    cvxc_float_t step0;
    cvxc_float_t gap;
    cvxc_float_t gap0;
    cvxc_float_t dsdz;
    cvxc_float_t dsdz0;

    cvxc_float_t relgap;
    cvxc_float_t resznl;
    cvxc_float_t reszl;
    cvxc_float_t resrznl;
    cvxc_float_t resrzl;
    cvxc_float_t relgap0;
    cvxc_float_t resznl0;
    cvxc_float_t reszl0;
    cvxc_float_t resrznl0;
    cvxc_float_t resrzl0;

    cvxc_float_t resx0;
    cvxc_float_t resy0;
    cvxc_float_t resz0;
    cvxc_float_t resx;
    cvxc_float_t resy;
    cvxc_float_t resz;
    cvxc_float_t hresx;
    cvxc_float_t hresy;
    cvxc_float_t hresz;
    cvxc_float_t cx;
    cvxc_float_t by;
    cvxc_float_t hz;
    cvxc_float_t rt;
    cvxc_float_t dres;
    cvxc_float_t pres;
    cvxc_float_t dinfres;
    cvxc_float_t pinfres;
    cvxc_float_t pcost;
    cvxc_float_t dcost;
    cvxc_float_t dres0;
    cvxc_float_t pres0;

    // KKT solver for CPL and CP problems; chaning to proper KKT solver
    cvxc_kktsolver_t cp_solver;

    // solution matrix
    cvxc_matrix_t x, y, s, z;  // ok

    cvxc_matrix_t c0;   // internal vector for non-linear target function
    cvxc_matrix_t z_mnl;

    // f is point in domain of F (mln,1); Df gradient (mnl,n)
    // H is Hessian (n, n);
    cvxc_matrix_t f, Df, H;
    cvxc_matrix_t newf, newDf; //, H;

    cvxc_matrix_t dx, dy, ds, dz;
    cvxc_matrix_t dx0, dy0, ds0, dz0;
    cvxc_matrix_t ds2, dz2, ds20, dz20;

    cvxc_matrix_t x0, y0, s0, z0;
    cvxc_matrix_t x1, y1, s1, z1;

    cvxc_matrix_t rx, ry, rz;
    cvxc_matrix_t rx0, ry0, rz0;
    cvxc_matrix_t newx, newy, newz, news, newrx;
    cvxc_matrix_t newrz0;
    cvxc_matrix_t hrx, hry, hrz;
    cvxc_matrix_t sigs, sigz;
    cvxc_matrix_t lmbda, lmbdasq;
    cvxc_matrix_t lmbda0;
    cvxc_matrix_t lmbdasq0;

    cvxc_matrix_t wx, wy, ws, wz;
    cvxc_matrix_t wx2, wy2, ws2, wz2;
    cvxc_matrix_t ws3, wz3;

    // matrix groups with correct indexing
    cvxc_matgrp_t h_g;
    cvxc_matgrp_t s_g, z_g;
    cvxc_matgrp_t s0_g, z0_g;
    cvxc_matgrp_t ds_g, dz_g;
    cvxc_matgrp_t ds0_g, dz0_g;
    cvxc_matgrp_t ds2_g, dz2_g;
    cvxc_matgrp_t ds20_g, dz20_g;
    cvxc_matgrp_t newz_g;
    cvxc_matgrp_t newrz_g;
    cvxc_matgrp_t news_g;
    cvxc_matgrp_t rz_g;
    cvxc_matgrp_t rzl_g, rznl_g;
    cvxc_matgrp_t ws_g, wz_g;
    cvxc_matgrp_t ws2_g, wz2_g;
    cvxc_matgrp_t ws3_g, wz3_g;
    cvxc_matgrp_t th_g;
    cvxc_matgrp_t z1_g;
    cvxc_matgrp_t lmbda_g, lmbdasq_g;
    cvxc_matgrp_t sigs_g, sigz_g;

    cvxc_index_t index_full;             // indexing to full matrix group
    cvxc_index_t index_packed;           // indexing to matrix group with packed 'S' space
    cvxc_index_t index_diag;             // indexing to matrix group with diagonal 'S' space
    cvxc_index_t index_sig;              // indexing to matrix group with only diagonal 'S" space
    cvxc_index_t index_cpt;              // indexing for G/h matrix with convex target function
    cvxc_scaling_t W;                    // scaling matrix group
    cvxc_scaling_t W0;
    cvxc_memblk_t work;                  // workspace
} cvxc_cpl_internal_t;


static inline
int F0(cvxc_convex_program_t *cp, cvxc_matrix_t *x0)
{
    return cp->F(x0, __cvxnil, __cvxnil, __cvxnil, __cvxnil, cp->user);
}

static inline
int F1(cvxc_convex_program_t *cp, cvxc_matrix_t *f, cvxc_matrix_t *Df, const cvxc_matrix_t *x)
{
    int e = cp->F(f, Df, __cvxnil, x, __cvxnil, cp->user);
    if (e == 0 && cvxm_isepi(x) && f) {
        cvxm_set(f, 0, 0, cvxm_get(f, 0, 0) - x->t);
    }
    return e;
}

static inline
int F2(cvxc_convex_program_t *cp, cvxc_matrix_t *f, cvxc_matrix_t *Df, cvxc_matrix_t *H,
       const cvxc_matrix_t *x, const cvxc_matrix_t *z)
{
    int e = cp->F(f, Df, H, x, z, cp->user);
    if (e == 0 && cvxm_isepi(x) && f) {
        cvxm_set(f, 0, 0, cvxm_get(f, 0, 0) - x->t);
    }
    return e;
}

static inline
cvxc_float_t __NaN()
{
#ifdef NAN
    return NAN;
#else
    return SQRT(-1.0);
#endif
}

static inline
int __isnull(const void *ptr) {
    return ptr == (void *)0;
}

static inline
cvxc_size_t __aligned128(cvxc_size_t n)
{
    return (n & 0xF) != 0 ? 16 - (n & 0xF) : 0;
}

static inline
cvxc_size_t __aligned64(cvxc_size_t n)
{
    return (n & 0x7) != 0 ? 8 - (n & 0x7) : 0;
}

static inline
cvxc_float_t cvxc_maxvec(int n, const cvxc_float_t *vec)
{
    cvxc_float_t r = vec[0];
    for (int k = 1; k < n; k++) {
        r = MAXF(vec[k], r);
    }
    return r;
}

static inline
cvxc_float_t cvxc_minvec(int n, const cvxc_float_t *vec)
{
    cvxc_float_t r = vec[0];
    for (int k = 1; k < n; k++) {
        r = MINF(vec[k], r);
    }
    return r;
}

static inline
void cvxc_mblk_empty(cvxc_memblk_t *m)
{
    if (m) {
        m->memory = (cvxc_float_t *)0;
        m->mlen = 0;
    }
}

/**
 * @brief Allocate space for n float elements
 */
static inline
void cvxc_mblk_init(cvxc_memblk_t *m, cvxc_size_t n)
{
    m->memory = (cvxc_float_t *)0;
    m->mlen = 0;
    cvxc_float_t *space = (cvxc_float_t *)calloc(n, sizeof(cvxc_float_t));
    if (space) {
        m->memory = space;
        m->mlen = n;
        m->__bytes = space;
    }
}

/**
 * @brief Create new block out of provided buffer.
 *
 * @param[out] m
 *    On entry uninitialized block. On exit proper memoery block.
 * @param[in] mlen
 *    Length of the memory block.
 * @param[in] ptr
 *    Pointer to memory buffer.
 * @param[in] nbytes
 *    Size of buffer in bytes.
 */
static inline
cvxc_size_t cvxc_mblk_make(cvxc_memblk_t *m, cvxc_size_t mlen, void *ptr, cvxc_size_t nbytes)
{
    if (!ptr || nbytes < mlen)
        return 0;
    m->memory = (cvxc_float_t *)ptr;
    m->mlen = mlen/sizeof(cvxc_float_t);
    m->__bytes = (void *)0;
    return mlen;
}

/**
 * @brief Release block resources.
 */
static inline
void cvxc_mblk_release(cvxc_memblk_t *m)
{
    if (m->__bytes)
        free(m->__bytes);
    m->memory = (cvxc_float_t *)0;
    m->mlen = 0;
    m->__bytes = (void *)0;
}

/**
 * @brief Get address of item at offset
 */
static inline
cvxc_float_t *cvxc_mblk_offset(cvxc_memblk_t *m, cvxc_size_t off)
{
    if (off >= m->mlen) {
        abort();
        return (cvxc_float_t *)0;
    }
    return &m->memory[off];
}

static inline
void cvxc_mblk_subblk(cvxc_memblk_t *d, cvxc_memblk_t *m, cvxc_size_t off)
{
    if (off >= m->mlen) {
        d->memory = (cvxc_float_t *)0;
        d->mlen = 0;
    } else {
        d->memory = &m->memory[off];
        d->mlen = m->mlen - off;
    }
}

static inline
void cvxc_mblk_clear(cvxc_memblk_t *m)
{
    if (m && m->mlen > 0)
        memset(m->memory, 0, m->mlen*sizeof(cvxc_float_t));
}

#endif /*  CVXC_INTERNAL_H */
