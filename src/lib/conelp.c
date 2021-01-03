/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */
#include <assert.h>
#include "cvxc.h"

#define EFRM "%8.5f"
//#define EFRM "%9.2e"
#define STEP 0.99

#if 0
static
void print_r_rti(cvxc_scaling_t *W, int k, const char *s)
{
    cvxc_matrix_t r, rti;
    cvxc_scaling_elem(&r, W, CVXWS_R, k);
    cvxc_scaling_elem(&rti, W, CVXWS_RTI, k);
    printf("%s r\n", s);
    cvxm_printf(stdout, "%13.6e", &r);
    printf("%s rti\n", s);
    cvxm_printf(stdout, "%13.6e", &rti);
}
#endif

static inline
cvxc_size_t __WORKBYTES(cvxc_size_t n)
{
    return 8*n*n*sizeof(cvxc_float_t);
}

static
int cvxc_res(cvxc_problem_t *cp,
            cvxc_matrix_t *ux,
            cvxc_matrix_t *uy,
            cvxc_matgrp_t *uz_g,
            double utau,
            cvxc_matgrp_t *us_g,
            double ukappa,
            cvxc_matrix_t *vx,
            cvxc_matrix_t *vy,
            cvxc_matgrp_t *vz_g,
            double *vtau,
            cvxc_matgrp_t *vs_g,
            double *vkappa,
            cvxc_scaling_t *W,
            double dg,
            cvxc_matgrp_t *lmbda_g)
{
    cvxc_matrix_t *us = us_g->mat, *uz = uz_g->mat;
    cvxc_matrix_t *vs = vs_g->mat, *vz = vz_g->mat;
    cvxc_matrix_t *lmbda = lmbda_g->mat;

    cvxc_conelp_internal_t *cpi = &cp->u.conelp;

    cvxc_float_t tauplus, lscale, vkplus;
    cvxc_matgrp_t ux_g = (cvxc_matgrp_t){ .mat = ux, .index = (cvxc_index_t *)0 };
    int err = 0;

    // vx = vx - A^T*uy - G^T*W^-1*uz - c*utau/dg
    cvxm_mult(1.0, vx, -1.0, cp->A, uy, CVXC_TRANSA);
    cvxm_copy(&cpi->wz3, uz, CVXC_ALL);
    cvxc_scale(&cpi->wz3_g, &cp->W, CVXC_INV, &cp->work);
    cvxc_sgemv(1.0, vx, -1.0, cp->G, &cpi->wz3_g, CVXC_TRANS);
    cvxm_axpy(vx, -utau/dg, cp->c);

    // vy = vy + A*ux - b*utau/dg
    cvxm_mult(1.0, vy, 1.0, cp->A, ux, 0);
    cvxm_axpy(vy, -utau/dg, cp->b);

    // vz = vz + G*ux - h*utau/dg + W^T*us
    cvxc_sgemv(1.0, vz, 1.0, cp->G, &ux_g, 0);
    cvxm_axpy(vz, -utau/dg, cp->h);
    cvxm_copy(&cpi->ws3, us, 0);
    cvxc_scale(&cpi->ws3_g, &cp->W, CVXC_TRANS, &cp->work);
    cvxm_axpy(vz, 1.0, &cpi->ws3);

    // vtau : vtau + c'*ux + b'uy + h'W^-1*uz + dg*ukappa
    tauplus  = dg*ukappa;
    tauplus += cvxm_dot(cp->c, ux);
    tauplus += cvxm_dot(cp->b, uy);
    tauplus += cvxc_sdot(&cpi->h_g, &cpi->wz3_g);
    *vtau   += tauplus;

    // vs = vs + lmbda o (uz + us)
    cvxm_copy(&cpi->ws3, us, 0);
    cvxm_axpy(&cpi->ws3, 1.0, uz);
    cvxc_sprod(&cpi->ws3_g, &cpi->lmbda_g, CVXC_DIAG, &cp->work);
    cvxm_axpy(vs, 1.0, &cpi->ws3);

    // vkappa += vkappa + lmbdag * (utau + ukappa)
    lscale = cvxm_get(lmbda, cp->cdim_diag, 0);
    vkplus = lscale * (utau + ukappa);
    *vkappa += vkplus;

    return err;
}

// f6_no_ir(x, y, z, tau, s, kappa) solves
//
//    [ 0         ]   [  0   A'  G'  c ] [ ux        ]    [ bx   ]
//    [ 0         ]   [ -A   0   0   b ] [ uy        ]    [ by   ]
//    [ W'*us     ] - [ -G   0   0   h ] [ W^{-1}*uz ] = -[ bz   ]
//    [ dg*ukappa ]   [ -c' -b' -h'  0 ] [ utau/dg   ]    [ btau ]
//
//    lmbda o (uz + us) = -bs
//    lmbdag * (utau + ukappa) = -bkappa.
//
// On entry, x, y, z, tau, s, kappa contain bx, by, bz, btau, 
// bkappa.  On exit, they contain ux, uy, uz, utau, ukappa.
static
int f6_no_ir(cvxc_problem_t *cp,
             cvxc_matrix_t *x,
             cvxc_matrix_t *y,
             cvxc_matgrp_t *z_g,
             cvxc_float_t *tau,
             cvxc_matgrp_t *s_g,
             cvxc_float_t *kappa)
{
    cvxc_conelp_internal_t *cpi = &cp->u.conelp;

    cvxc_matrix_t *s = s_g->mat, *z = z_g->mat;


    //cvxc_matrix_t r;

    // Solve
    //
    // [  0   A'  G'    0   ] [ ux        ]
    // [ -A   0   0     b   ] [ uy        ]
    // [ -G   0   W'*W  h   ] [ W^{-1}*uz ]
    // [ -c' -b' -h'    k/t ] [ utau/dg   ]
    //
    //   [ bx                    ]
    //   [ by                    ]
    // = [ bz - W'*(lmbda o\ bs) ]
    //   [ btau - bkappa/tau     ]
    //
    // us = -lmbda o\ bs - uz
    // ukappa = -bkappa/lmbdag - utau.

    // First solve
    //
    // [ 0  A' G'   ] [ ux        ]   [  bx                    ]
    // [ A  0  0    ] [ uy        ] = [ -by                    ]
    // [ G  0 -W'*W ] [ W^{-1}*uz ]   [ -bz + W'*(lmbda o\ bs) ]

    cvxm_scale(y, -1.0, 0);
    // s = -lmbda o\ s = - lmbda o\ bs
    cvxc_sinv(s_g, &cpi->lmbda_g, &cp->work);
    cvxm_scale(s, -1.0, 0);

    // z = -(z + W'*s) = -bz + W'(lmbda o\ bs)
    cvxm_copy(&cpi->ws3, s, CVXC_ALL);
    cvxc_scale(&cpi->ws3_g, &cp->W, CVXC_TRANS, &cp->work);
    cvxm_axpy(z, 1.0, &cpi->ws3);
    cvxm_scale(z, -1.0, 0);

    cvxc_kktsolve(cp->solver, x, y, z_g);
    // Combine with solution of
    //
    // [ 0   A'  G'    ] [ x1         ]          [ c ]
    // [-A   0   0     ] [ y1         ] = -dgi * [ b ]
    // [-G   0   W'*W  ] [ W^{-1}*dzl ]          [ h ]
    //
    // to satisfy
    //
    // -c'*x - b'*y - h'*W^{-1}*z + dg*tau = btau - bkappa/tau. '
    cvxc_float_t lkappa = - (*kappa) / cvxm_get(&cpi->lmbda, cp->cdim_diag, 0);
    cvxc_float_t ltau  = *tau + lkappa/cp->u.conelp.dgi;
    ltau += cvxm_dot(cp->c, x);
    ltau += cvxm_dot(cp->b, y);
    ltau += cvxc_sdot(&cpi->th_g, z_g);  // TODO: check this
    ltau = cp->u.conelp.dgi * ltau  / (1.0 + cvxc_sdot(&cpi->z1_g, &cpi->z1_g));

    *tau = ltau;
    cvxm_axpy(x, ltau, &cpi->x1);
    cvxm_axpy(y, ltau, &cpi->y1);
    cvxm_axpy(z, ltau, &cpi->z1);

    cvxm_axpy(s, -1.0, z);
    *kappa = lkappa - ltau;

    return 0;
}


// f6(x, y, z, tau, s, kappa) solves the same system as f6_no_ir,
// but applies iterative refinement.
static
int f6(cvxc_problem_t *cp,
       cvxc_matrix_t *x,
       cvxc_matrix_t *y,
       cvxc_matgrp_t *z_g,
       cvxc_float_t *tau,
       cvxc_matgrp_t *s_g,
       cvxc_float_t *kappa,
       int refinement)
{
    cvxc_matrix_t *s = s_g->mat, *z = z_g->mat;
    int /*lags = 0,*/ err = 0;
    double wtau, wkappa, wtau2, wkappa2;

    cvxc_conelp_internal_t *cpi = &cp->u.conelp;

    if (refinement > 0 /*|| (flags | CVXC_DEBUG) != 0*/) {
        cvxm_copy(&cpi->wx, x, 0);
        cvxm_copy(&cpi->wy, y, 0);
        cvxm_copy(&cpi->wz, z, 0);
        cvxm_copy(&cpi->ws, s, 0);
        wtau   = *tau;
        wkappa = *kappa;
    }
    err = f6_no_ir(cp, x, y, z_g, tau, s_g, kappa);
    for (int i = 0; i < refinement; i++) {
        cvxm_copy(&cpi->wx2, &cpi->wx, 0);
        cvxm_copy(&cpi->wy2, &cpi->wy, 0);
        cvxm_copy(&cpi->wz2, &cpi->wz, 0);
        cvxm_copy(&cpi->ws2, &cpi->ws, 0);
        wtau2   = wtau;
        wkappa2 = wkappa;

        cvxc_res(cp, x, y, z_g, *tau, s_g, *kappa,
                &cpi->wx2, &cpi->wy2, &cpi->wz2_g, &wtau2,
                &cpi->ws2_g, &wkappa2, &cp->W, cpi->dg, &cpi->lmbda_g);

        f6_no_ir(cp, &cpi->wx2, &cpi->wy2,
                 &cpi->wz2_g, &wtau2, &cpi->ws2_g, &wkappa2);

        cvxm_axpy(x, 1.0, &cpi->wx2);
        cvxm_axpy(y, 1.0, &cpi->wy2);
        cvxm_axpy(s, 1.0, &cpi->ws2);
        cvxm_axpy(z, 1.0, &cpi->wz2);

        *tau   += wtau2;
        *kappa += wkappa2;
    }
    return err;
}

/**
 * @brief Compute memory allocation needed for CONELP problem
 *
 * @param[in] n  Number of variables
 * @param[in] m  Number of equality constrains (rows of A matrix)
 * @param[in] dims Dimensions of inequality constraints, linear, scop and sdp
 *
 * @return Number of bytes of memory needed.
 */
cvxc_size_t cvxc_conelp_bytes(int n, int m, const cvxc_dimset_t *dims)
{
    cvxc_size_t cdim      = cvxc_dimset_sum_squared(dims, CVXDIM_CONELP);
    cvxc_size_t cdim_diag = cvxc_dimset_sum(dims, CVXDIM_CONELP);
    cvxc_size_t sdim      = cvxc_dimset_sum(dims, CVXDIM_SDP);
    cvxc_size_t maxsdp    = cvxc_dimset_max(dims, CVXDIM_SDP);
    cvxc_size_t total = 0;

    total += 7*n;      // for x, dx, rx, hrx, x1, wx, wx2
    total += 7*m;      // for y, dy, ry, hry, y1, wy, wy2
    total += 5*cdim;   // for s, ds, ws, ws2, ws3
    total += 8*cdim;   // for z, dz, rz, hrz, z1, wz, wz2, wz3
    total += 2*cdim_diag + 2; // for lmbda, lmbdasq
    total += 2*sdim;   // for sigs, sigz
    total += cdim;     // for th

    cvxc_size_t nbytes = total*sizeof(cvxc_float_t);

    // workspace for SDP constraint scaling; TODO: think about this.
    if (maxsdp > 0)
        nbytes += __WORKBYTES(maxsdp);

    cvxc_size_t isize;
    // calculte space need for scaling matrix
    nbytes += cvxc_scaling_bytes(&isize, dims);

    // size of standard index set
    nbytes += cvxc_index_bytes(dims, 0);
    // size of packed index set
    nbytes += cvxc_index_bytes(dims, 1);
    // size of diagonal index set
    nbytes += cvxc_index_bytes(dims, 2);
    // size of SDP diagonal index set
    nbytes += cvxc_index_bytes(dims, 3);

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
 * @param[in] memory
 *      Pointer to memory block
 * @param[in] nbytes
 *      Size of memory block in bytes
 *
 * @return Number of bytes of memory used, zero if memory was not large enough.
 *
 */
cvxc_size_t cvxc_conelp_make(cvxc_problem_t *cp,
                           int n,
                           int m,
                           const cvxc_dimset_t *dims,
                           void *memory,
                           cvxc_size_t nbytes)
{
    cvxc_size_t offset = 0;
    cvxc_size_t used = 0;
    unsigned char *bytes = (unsigned char *)memory;
    cvxc_size_t cdim =
        cvxc_dimset_sum_squared(dims, CVXDIM_CONELP);

    cvxc_conelp_internal_t *cpi = &cp->u.conelp;

    cp->cdim = cdim;

#if 0
    // allocate space for results TODO: always ??
    cp->x = cvxm_new(n, 1);
    cp->y = cvxm_new(m, 1);
    cp->s = cvxm_new(cdim, 1);
    cp->z = cvxm_new(cdim, 1);
#endif
    // __INIT macro assumes variables offset and nbytes;
    // overlay index sets
    __INIT(used, cvxc_index_make(&cpi->index_full, dims, 0, bytes,  nbytes));
    __INIT(used, cvxc_index_make(&cpi->index_packed, dims, 1, &bytes[offset],  nbytes));
    __INIT(used, cvxc_index_make(&cpi->index_diag, dims, 2, &bytes[offset],  nbytes));
    __INIT(used, cvxc_index_make(&cpi->index_sig, dims, 3, &bytes[offset],  nbytes));

    // map result matrix; allocate realy
    __INIT(used, cvxm_make(&cpi->x, n, 1, &bytes[offset], nbytes));
    __INITC(used, m, &cpi->y, cvxm_make(&cpi->y, m, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->s, cdim, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->z, cdim, 1, &bytes[offset], nbytes));

    // dx, dy, ds, dz
    __INIT(used, cvxm_make(&cpi->dx, n, 1, &bytes[offset], nbytes));
    __INITC(used, m, &cpi->dy, cvxm_make(&cpi->dy, m, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->ds, cdim, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->dz, cdim, 1, &bytes[offset], nbytes));

    // rx, ry, rz
    __INIT(used, cvxm_make(&cpi->rx, n, 1, &bytes[offset], nbytes));
    __INITC(used, m, &cpi->ry, cvxm_make(&cpi->ry, m, 1, &bytes[offset], nbytes));
    //__INIT(used, cvxm_make(&cp->rs, cdim, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->rz, cdim, 1, &bytes[offset], nbytes));

    // x1, y1, z1
    __INIT(used, cvxm_make(&cpi->x1, n, 1, &bytes[offset], nbytes));
    __INITC(used, m, &cpi->y1, cvxm_make(&cpi->y1, m, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->z1, cdim, 1, &bytes[offset], nbytes));

    // hrx, hry, hrz
    __INIT(used, cvxm_make(&cpi->hrx, n, 1, &bytes[offset], nbytes));
    __INITC(used, m, &cpi->hry, cvxm_make(&cpi->hry, m, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->hrz, cdim, 1, &bytes[offset], nbytes));

    // wx, wy, ws, wz
    __INIT(used, cvxm_make(&cpi->wx, n, 1, &bytes[offset], nbytes));
    __INITC(used, m, &cpi->wy, cvxm_make(&cpi->wy, m, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->ws, cdim, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->wz, cdim, 1, &bytes[offset], nbytes));

    // wx2, wy2, ws2, wz2
    __INIT(used, cvxm_make(&cpi->wx2, n, 1, &bytes[offset], nbytes));
    __INITC(used, m, &cpi->wy2, cvxm_make(&cpi->wy2, m, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->ws2, cdim, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->wz2, cdim, 1, &bytes[offset], nbytes));

    // ws3, wz3
    __INIT(used, cvxm_make(&cpi->ws3, cdim, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->wz3, cdim, 1, &bytes[offset], nbytes));

    // th
    __INIT(used, cvxm_make(&cpi->th, cdim, 1, &bytes[offset], nbytes));

    cvxc_size_t cdim_diag =
        cvxc_dimset_sum(dims, CVXDIM_CONELP);

    cp->cdim_diag = cdim_diag;

    // lmbda, lmbdasq
    __INIT(used, cvxm_make(&cpi->lmbda,   cdim_diag+1, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->lmbdasq, cdim_diag+1, 1, &bytes[offset], nbytes));

    cvxc_size_t sdim   = cvxc_dimset_sum(dims, CVXDIM_SDP);
    cvxc_size_t maxsdp = cvxc_dimset_max(dims, CVXDIM_SDP);

    // sigs, sigz; space for eigenvalues; zero if no SDP constraints
    __INITC(used, sdim, &cpi->sigs, cvxm_make(&cpi->sigs, sdim, 1, &bytes[offset], nbytes));
    __INITC(used, sdim, &cpi->sigz, cvxm_make(&cpi->sigz, sdim, 1, &bytes[offset], nbytes));

    // scaling matrix
    __INIT(used, cvxc_scaling_make(&cp->W, dims, &bytes[offset], nbytes));

    // workspace for SDP contraints handling
    cvxc_mblk_empty(&cp->work);
    if (maxsdp > 0) {
        __INIT(used, cvxc_mblk_make(&cp->work, __WORKBYTES(maxsdp), &bytes[offset], nbytes));
    }

    // setup matrix group variables
    cvxc_mgrp_init(&cpi->h_g,   cp->h,    &cpi->index_full);
    cvxc_mgrp_init(&cpi->s_g,   &cpi->s,   &cpi->index_full);
    cvxc_mgrp_init(&cpi->z_g,   &cpi->z,   &cpi->index_full);
    cvxc_mgrp_init(&cpi->z1_g,  &cpi->z1,  &cpi->index_full);
    cvxc_mgrp_init(&cpi->th_g,  &cpi->th,  &cpi->index_full);
    cvxc_mgrp_init(&cpi->ds_g,  &cpi->ds,  &cpi->index_full);
    cvxc_mgrp_init(&cpi->dz_g,  &cpi->dz,  &cpi->index_full);
    cvxc_mgrp_init(&cpi->rz_g,  &cpi->rz,  &cpi->index_full);
    cvxc_mgrp_init(&cpi->hrz_g, &cpi->hrz, &cpi->index_full);
    cvxc_mgrp_init(&cpi->ws_g,  &cpi->ws,  &cpi->index_full);
    cvxc_mgrp_init(&cpi->wz_g,  &cpi->wz,  &cpi->index_full);
    cvxc_mgrp_init(&cpi->ws2_g, &cpi->ws2, &cpi->index_full);
    cvxc_mgrp_init(&cpi->wz2_g, &cpi->wz2, &cpi->index_full);
    cvxc_mgrp_init(&cpi->ws3_g, &cpi->ws3, &cpi->index_full);
    cvxc_mgrp_init(&cpi->wz3_g, &cpi->wz3, &cpi->index_full);

    cvxc_mgrp_init(&cpi->sigs_g,   &cpi->sigs,   &cpi->index_sig);
    cvxc_mgrp_init(&cpi->sigz_g,   &cpi->sigz,   &cpi->index_sig);

    cvxc_mgrp_init(&cpi->lmbda_g,   &cpi->lmbda,   &cpi->index_diag);
    cvxc_mgrp_init(&cpi->lmbdasq_g, &cpi->lmbdasq, &cpi->index_diag);


    return offset;
}

int cvxc_conelp_isok(const cvxc_matrix_t *c,
                    const cvxc_matrix_t *G,
                    const cvxc_matrix_t *h,
                    const cvxc_matrix_t *A,
                    const cvxc_matrix_t *b,
                    const cvxc_dimset_t *dims)
{
    cvxc_size_t mc, nc, mG, nG, mA, nA, mh, nh, mb, nb;

    mc = nc = mG = nG = mA = nA = mh = nh = mb = nb = 0;
    if (!c) {
        return CVXC_ERR_NULLCOST;
    }

    cvxm_size(&mc, &nc, c);
    if (G)
        cvxm_size(&mG, &nG, G);
    if (A)
        cvxm_size(&mA, &nA, A);
    if (h)
        cvxm_size(&mh, &nh, h);
    if (b)
        cvxm_size(&mb, &nb, b);

    if (nc > 1 || mc < 1) {
        return CVXC_ERR_DIMC;
    }

    cvxc_size_t cdim      = cvxc_dimset_sum_squared(dims, CVXDIM_CONELP);
    cvxc_size_t cdim_pckd = cvxc_dimset_sum_packed(dims, CVXDIM_CONELP);

    if (nh > 1 || mh != cdim) {
        return CVXC_ERR_DIMH;
    }

    if (mG != cdim || nG != mc) {
        return CVXC_ERR_DIMG;
    }
    if (nA != mc || mA != mb) {
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

cvxc_size_t cvxc_conelp_setup(cvxc_problem_t *cp,
                            cvxc_matrix_t *c,
                            cvxc_matrix_t *G,
                            cvxc_matrix_t *h,
                            cvxc_matrix_t *A,
                            cvxc_matrix_t *b,
                            cvxc_dimset_t *dims,
                            cvxc_kktsolver_t *kktsolver)
{
    if (! cp)
        return 0;

    cvxc_size_t mc, nc, mb, nb;
    int err;
    if ((err = cvxc_conelp_isok(c, G, h, A, b, dims)) != 0) {
        cp->error = err;
        return 0;
    }

    mc = nc = mb = nb = 0;
    cvxm_size(&mc, &nc, c);
    if (b)
        cvxm_size(&mb, &nb, b);

    cvxc_size_t nbytes = cvxc_conelp_bytes(mc, mb, dims);
    void *memory = calloc(nbytes, 1);
    if (!memory) {
        cp->error = CVXC_ERR_MEMORY;
        return 0;
    }

    cp->c = c;
    cp->G = G;
    cp->h = h;
    cp->A = A;
    cp->b = b;

    cp->primal_x = (cvxc_matrix_t *)0;
    cp->primal_s = (cvxc_matrix_t *)0;
    cp->dual_y = (cvxc_matrix_t *)0;
    cp->dual_z = (cvxc_matrix_t *)0;

    if (cvxc_conelp_make(cp, mc, mb, dims, memory, nbytes) == 0) {
        cp->error = CVXC_ERR_MEMORY;
        return 0;
    }
    cp->mlen = nbytes;
    cp->memory = memory;
    // provide full index set
    cp->index_g = &cp->u.conelp.index_full;

    // init KKT solver
    if (kktsolver) {
        cp->solver = kktsolver;
    } else {
        cvxc_ldlsolver_init(&cp->__S, cp, mc, mb, dims);
        cp->solver = &cp->__S;
    }

    return nbytes;
}

void cvxc_conelp_set_start(cvxc_problem_t *prob,
                          cvxc_matrix_t *primal_x,
                          cvxc_matrix_t *primal_s,
                          cvxc_matrix_t *dual_y,
                          cvxc_matrix_t *dual_z)
{
    if (! prob)
        return;
    // TODO: check sizes; and must give primal_x, primal_s or neither (likewise for y,z)
    prob->primal_x = primal_x;
    prob->primal_s = primal_s;
    prob->dual_y = dual_y;
    prob->dual_z = dual_z;
}

void cvxc_conelp_init_scaling(cvxc_scaling_t *W)
{
    cvxc_matrix_t x;
    // D = ones
    cvxc_scaling_elem(&x, W, CVXWS_D, 0);
    cvxm_mkconst(&x, 1.0);
    // DI = ones
    cvxc_scaling_elem(&x, W, CVXWS_DI, 0);
    cvxm_mkconst(&x, 1.0);
    if (W->vcount > 0) {
        // BETA = ones
        cvxc_scaling_elem(&x, W, CVXWS_BETA, 0);
        cvxm_mkconst(&x, 1.0);
        // V is unit vector
        for (int k = 0; k < W->vcount; k++) {
            cvxc_scaling_elem(&x, W, CVXWS_V, k);
            cvxm_mkident(&x);
        }
    }
    // R & RTI are identity
    for (int k = 0; k < W->rcount; k++) {
        cvxc_scaling_elem(&x, W, CVXWS_R, k);
        cvxm_mkident(&x);
        cvxc_scaling_elem(&x, W, CVXWS_RTI, k);
        cvxm_mkident(&x);
    }
}

int cvxc_conelp_compute_start(cvxc_problem_t *cp)
{
    cvxc_conelp_internal_t *cpi = &cp->u.conelp;
    cvxc_stats_t *stats = &cp->stats;

    int primalstart = ! (cp->primal_x && cp->primal_s);
    int dualstart = ! (cp->dual_y && cp->dual_z);


    if (primalstart || dualstart) {
        cvxc_conelp_init_scaling(&cp->W);
        cvxc_kktfactor(cp->solver, &cp->W, __cvxnil, __cvxnil);
    }
    if (primalstart) {
        // minimize    || G * x - h ||^2
        // subject to  A * x = b
        //
        // by solving
        //
        //     [ 0   A'  G' ]   [ x  ]   [ 0 ]
        //     [ A   0   0  ] * [ dy ] = [ b ].
        //     [ G   0  -I  ]   [ -s ]   [ h ]
        cvxm_scale(&cpi->x, 0.0, 0);
        cvxm_copy(&cpi->dy, cp->b, 0);
        cvxm_copy(&cpi->s, cp->h, 0);
        cvxc_kktsolve(cp->solver, &cpi->x, &cpi->dy, &cpi->s_g);
        cvxm_scale(&cpi->s, -1.0, 0);
    } else {
        cvxm_copy(&cpi->x, cp->primal_x, 0);
        cvxm_copy(&cpi->s, cp->primal_s, 0);
    }

    cpi->ts = cvxc_max_step(&cpi->s_g, __nilgrp, &cp->work);

    if (cpi->ts >= 0.0 && ! primalstart) {
        cp->error = CVXC_ERR_NEG_INITIAL_S;
        return -1;
    }
    if (dualstart) {
        // minimize   || z ||^2
        // subject to G'*z + A'*y + c = 0
        //
        // by solving
        //
        //     [ 0   A'  G' ] [ dx ]   [ -c ]
        //     [ A   0   0  ] [ y  ] = [  0 ].
        //     [ G   0  -I  ] [ z  ]   [  0 ]
        cvxm_scale(&cpi->y, 0.0, 0);
        cvxm_copy(&cpi->dx, cp->c, 0);
        cvxm_scale(&cpi->dx, -1.0, 0);
        cvxm_scale(&cpi->z, 0.0, 0);
        cvxc_kktsolve(cp->solver, &cpi->dx, &cpi->y, &cpi->z_g);
    } else {
        cvxm_copy(&cpi->y, cp->dual_y, 0);
        cvxm_copy(&cpi->z, cp->dual_z, 0);
    }

    cpi->tz = cvxc_max_step(&cpi->z_g,  __nilgrp, &cp->work);
    if (cpi->tz >= 0.0 && ! dualstart) {
        cp->error = CVXC_ERR_NEG_INITIAL_Z;
        return -1;
    }

    stats->resx0 = MAXF(1.0, cvxm_nrm2(cp->c));
    stats->resy0 = MAXF(1.0, cvxm_nrm2(cp->b));
    stats->resz0 = MAXF(1.0, cvxc_snrm2(&cpi->h_g));

    cpi->nrms = cvxc_snrm2(&cpi->s_g);
    cpi->nrmz = cvxc_snrm2(&cpi->z_g);

    stats->gap = stats->pcost = stats->dcost = stats->relgap = 0.0;

    return 0;
}



int cvxc_conelp_ready(cvxc_problem_t *cp, cvxc_stats_t *stats, int iter, int stat)
{
    cvxc_conelp_internal_t *cpi = &cp->u.conelp;

    //cvxc_stats_t *stats = &cp->stats;

    if (stat == CVXC_STAT_OPTIMAL && iter == -1) {
        // constructed initial point is feasible and optimal
        cvxc_mksymm(&cpi->s_g);
        cvxc_mksymm(&cpi->z_g);

        // rx = A'*y + G'z + c
        cvxm_copy(&cpi->rx, cp->c, CVXC_ALL);
        cvxm_mvmult(1.0, &cpi->rx, 1.0, cp->A, &cpi->y, CVXC_TRANS);
        cvxm_mvmult(1.0, &cpi->rx, 1.0, cp->G, &cpi->z, CVXC_TRANS);
        stats->resx = cvxm_nrm2(&cpi->rx);

        // ry = b - A*x  ; TODO - computes -b - A*x ;; check
        cvxm_copy(&cpi->ry, cp->b, CVXC_ALL);
        cvxm_mvmult(-1.0, &cpi->ry, -1.0, cp->A, &cpi->x, CVXC_TRANS);
        stats->resy = cvxm_nrm2(&cpi->ry);

        // rz = s + G*x - h
        cvxm_copy(&cpi->rz, &cpi->s, CVXC_ALL);
        cvxm_mvmult(1.0, &cpi->rz, 1.0, cp->G, &cpi->x, 0);
        cvxm_axpy(&cpi->rz, 1.0, cp->h);
        stats->resz = cvxc_snrm2(&cpi->rz_g);

        stats->pres = MAXF(stats->resy/stats->resy0, stats->resz/stats->resz0);
        stats->dres = stats->resx/stats->resx0;
        stats->cx = cvxm_dot(cp->c, &cpi->x);
        stats->by = cvxm_dot(cp->b, &cpi->y);
        stats->hz = cvxc_sdot(&cpi->h_g, &cpi->z_g);

        cp->solution.x = &cpi->x;
        cp->solution.s = &cpi->s;
        cp->solution.y = &cpi->y;
        cp->solution.z = &cpi->z;

        cp->solution.status = stat;
        cp->solution.gap = stats->gap;
        cp->solution.relative_gap = stats->relgap;
        cp->solution.primal_objective = stats->cx;
        cp->solution.dual_objective = - (stats->by + stats->hz);
        cp->solution.primal_infeasibility = stats->pres;
        cp->solution.dual_infeasibility = stats->dres;
        cp->solution.primal_slack = - cpi->ts;
        cp->solution.dual_slack = - cpi->tz;
        cp->solution.primal_residual_cert = __NaN();
        cp->solution.dual_residual_cert = __NaN();
        cp->solution.iterations = 0;
    }
    else if (stat == CVXC_STAT_UNKNOWN || stat == CVXC_STAT_OPTIMAL) {
        cvxm_scale(&cpi->x, 1.0/cpi->tau, CVXC_ALL);
        cvxm_scale(&cpi->y, 1.0/cpi->tau, CVXC_ALL);
        cvxm_scale(&cpi->s, 1.0/cpi->tau, CVXC_ALL);
        cvxm_scale(&cpi->z, 1.0/cpi->tau, CVXC_ALL);
        cvxc_mksymm(&cpi->s_g);
        cvxc_mksymm(&cpi->z_g);
        cpi->ts = cvxc_max_step(&cpi->s_g, __nilgrp, &cp->work);
        cpi->tz = cvxc_max_step(&cpi->s_g, __nilgrp, &cp->work);

        cp->error = stat == CVXC_STAT_UNKNOWN ? CVXC_ERR_MAXITER : 0;
        cp->solution.x = &cpi->x;
        cp->solution.s = &cpi->s;
        cp->solution.y = &cpi->y;
        cp->solution.z = &cpi->z;

        cp->solution.status = stat;

        cp->solution.gap = stats->gap;
        cp->solution.relative_gap = stats->relgap;
        cp->solution.primal_objective = stats->cx;
        cp->solution.dual_objective = - (stats->by + stats->hz);
        cp->solution.primal_infeasibility = stats->pres;
        cp->solution.dual_infeasibility = stats->dres;
        cp->solution.primal_slack = - cpi->ts;
        cp->solution.dual_slack = - cpi->tz;
        if (stat == CVXC_STAT_OPTIMAL) {
            cp->solution.primal_residual_cert = __NaN();
            cp->solution.dual_residual_cert = __NaN();
        } else {
            cp->solution.primal_residual_cert = stats->pinfres;
            cp->solution.dual_residual_cert = stats->dinfres;
        }
        cp->solution.iterations = iter;
    }
    else if (stat == CVXC_STAT_PRIMAL_INFEASIBLE) {
        cp->solution.status = stat;
        cp->error = stat;
        cvxm_scale(&cpi->y, 1.0/(-stats->hz - stats->by), CVXC_ALL);
        cvxm_scale(&cpi->z, 1.0/(-stats->hz - stats->by), CVXC_ALL);
        cvxc_mksymm(&cpi->z_g);

        cp->solution.x = __cvxnil;
        cp->solution.s = __cvxnil;
        cp->solution.y = __cvxnil;
        cp->solution.z = __cvxnil;

        cp->solution.status = stat;

        cp->solution.gap = __NaN();
        cp->solution.relative_gap = __NaN();
        cp->solution.primal_objective = __NaN();
        cp->solution.dual_objective = 1.0;
        cp->solution.primal_infeasibility = __NaN();
        cp->solution.dual_infeasibility = __NaN();
        cp->solution.primal_slack = __NaN();
        cp->solution.dual_slack = - cpi->tz;
        cp->solution.primal_residual_cert = stats->pinfres;
        cp->solution.dual_residual_cert = __NaN();
        cp->solution.iterations = iter;
    }
    else if (stat == CVXC_STAT_DUAL_INFEASIBLE) {
        cp->error = stat;
        cp->solution.status = stat;
        cvxm_scale(&cpi->x, 1.0/-stats->cx, CVXC_ALL);
        cvxm_scale(&cpi->s, 1.0/-stats->cx, CVXC_ALL);
        cvxc_mksymm(&cpi->s_g);

        cp->solution.x = __cvxnil;
        cp->solution.s = __cvxnil;
        cp->solution.y = __cvxnil;
        cp->solution.z = __cvxnil;

        cp->solution.status = stat;

        cp->solution.gap = __NaN();
        cp->solution.relative_gap = __NaN();
        cp->solution.primal_objective = 1.0;
        cp->solution.dual_objective = __NaN();
        cp->solution.primal_infeasibility = __NaN();
        cp->solution.dual_infeasibility = __NaN();
        cp->solution.primal_slack = - cpi->ts;
        cp->solution.dual_slack = __NaN();
        cp->solution.primal_residual_cert = __NaN();
        cp->solution.dual_residual_cert = stats->dinfres;
        cp->solution.iterations = iter;
    }
    return -stat;
}

int cvxc_conelp_solve(cvxc_problem_t *cp, cvxc_solopts_t *opts)
{
    cvxc_conelp_internal_t *cpi = &cp->u.conelp;
    cvxc_stats_t *stats = &cp->stats;

    cvxc_float_t max_nrms, max_nrmz;
    cvxc_matrix_t lk, sk, zk, ls, lz;

    cvxc_float_t abstol = CVXC_ABSTOL;
    cvxc_float_t reltol = CVXC_RELTOL;
    cvxc_float_t feastol = CVXC_FEASTOL;

    int maxiter = opts && opts->max_iter > 0 ? opts->max_iter : CVXC_MAXITER;
    int refinement = opts && opts->refinement > 0 ? opts->refinement : 0;

    if (cvxc_index_count(&cpi->index_full, CVXDIM_SOCP) > 0 ||
        cvxc_index_count(&cpi->index_full, CVXDIM_SDP) > 0)
        refinement = 1;

    int primalstart = ! (cp->primal_x && cp->primal_s);
    int dualstart = ! (cp->dual_y && cp->dual_z);

    cp->error = 0;

    if (primalstart && dualstart) {
        // check for constructed initial point
        stats->gap = cvxc_sdot(&cpi->s_g, &cpi->z_g);
        stats->pcost = cvxm_dot(cp->c, &cpi->x);
        stats->dcost =
            - cvxm_dot(cp->b, &cpi->y) - cvxc_sdot(&cpi->h_g, &cpi->z_g);

        if (stats->pcost < 0.0) {
            stats->relgap = stats->gap / - stats->pcost;
        } else if (stats->dcost > 0.0) {
            stats->relgap = stats->gap / stats->dcost;
        } else {
            stats->relgap = __NaN();
        }

        if (cpi->ts <= 0.0 && cpi->tz < 0 &&
            (stats->gap <= abstol ||
             (!isnan(stats->relgap) && stats->relgap <= reltol))) {

            // initial point is feasible and optimal
            return cvxc_conelp_ready(cp, stats, -1, CVXC_STAT_OPTIMAL);
        }

        max_nrms = MAXF(1.0, cpi->nrms);
        max_nrmz = MAXF(1.0, cpi->nrmz);
        if (cpi->ts >= - 1e-8 * max_nrms) {
            cvxc_mgrp_update_sz(&cpi->s_g, 1.0+cpi->ts, 0);
        }

        if (cpi->tz >= -1e-8 * max_nrmz) {
            cvxc_mgrp_update_sz(&cpi->z_g, 1.0+cpi->tz, 0);
        }
    } else if (primalstart && ! dualstart) {
        max_nrms = MAXF(1.0, cpi->nrms);
        if (cpi->ts >= - 1e-8 * max_nrms) {
            cvxc_mgrp_update_sz(&cpi->s_g, 1.0+cpi->ts, 0);
        }

    } else if (dualstart && ! primalstart) {
        max_nrmz = MAXF(1.0, cpi->nrmz);
        if (cpi->tz >= -1e-8 * max_nrmz) {
            cvxc_mgrp_update_sz(&cpi->z_g, 1.0+cpi->tz, 0);
        }
    }

    cpi->tau = 1.0;
    cpi->kappa = 1.0;
    cpi->wkappa3 = 0.0;
    cvxm_copy(&cpi->rx,  cp->c, 0);
    cvxm_copy(&cpi->hrx, cp->c, 0);
    cvxm_copy(&cpi->ry,  cp->b, 0);
    cvxm_copy(&cpi->hry, cp->b, 0);

    stats->gap = cvxc_sdot(&cpi->s_g, &cpi->z_g);

    // -----------------------------------------------------------------------------

    for (int iter = 0; iter <= maxiter; iter++) {

        // hrx = -A'*y - G'*z
        cvxm_scale(&cpi->hrx, 0.0, CVXC_ALL);
        cvxm_mvmult(0.0, &cpi->hrx, -1.0, cp->A, &cpi->y, CVXC_TRANS);
        cvxc_sgemv(1.0, &cpi->hrx, -1.0, cp->G, &cpi->z_g, CVXC_TRANS);
        stats->hresx = SQRT(cvxm_dot(&cpi->hrx, &cpi->hrx));

        // rx  = hrx - c*tau
        //     = -A'*y - G'*z - c*tau
        cvxm_copy(&cpi->rx, &cpi->hrx, 0);
        cvxm_axpy(&cpi->rx, -cpi->tau, cp->c);
        stats->resx = cvxm_nrm2(&cpi->rx) / cpi->tau;

        // hry = A*x
        cvxm_scale(&cpi->hry, 0.0, CVXC_ALL);
        cvxm_mvmult(0.0, &cpi->hry, 1.0, cp->A, &cpi->x, 0);
        stats->hresy = cvxm_nrm2(&cpi->hry);

        // ry  = hry - b*tau
        //     = A*x - b*tau
        cvxm_copy(&cpi->ry, &cpi->hry, 0);
        cvxm_axpy(&cpi->ry, -cpi->tau, cp->b);
        stats->resy = cvxm_nrm2(&cpi->ry) / cpi->tau;

        // hrz = s + G*x
        cvxm_scale(&cpi->hrz, 0.0, CVXC_ALL);
        cvxm_mvmult(0.0, &cpi->hrz, 1.0, cp->G, &cpi->x, 0);
        cvxm_axpy(&cpi->hrz, 1.0, &cpi->s);
        stats->hresz = cvxc_snrm2(&cpi->hrz_g);

        //  rz = hrz = h*tau
        //     = s + G*x - h*tau
        cvxm_copy(&cpi->rz, &cpi->hrz, 0);
        cvxm_axpy(&cpi->rz, -cpi->tau, cp->h);

        stats->resz = cvxc_snrm2(&cpi->rz_g) / cpi->tau;

        // rt = kappa + c'*x + b'*y + h'*z '
        stats->cx = cvxm_dot(cp->c, &cpi->x);
        stats->by = cvxm_dot(cp->b, &cpi->y);
        stats->hz = cvxc_sdot(&cpi->h_g, &cpi->z_g);
        stats->rt = cpi->kappa + stats->cx + stats->by + stats->hz;

        assert(isfinite(stats->cx));
        assert(isfinite(cpi->tau));

        // statistics for stopping
        stats->pcost = stats->cx / cpi->tau;
        stats->dcost = -(stats->by + stats->hz) / cpi->tau;
        if (stats->pcost < 0.0) {
            stats->relgap = stats->gap / -stats->pcost;
        } else if (stats->dcost > 0.0) {
            stats->relgap = stats->gap / stats->dcost;
        } else {
            stats->relgap = __NaN();
        }

        stats->pres = MAXF((stats->resy/stats->resy0), (stats->resz/stats->resz0));
        stats->dres = stats->resx / stats->resx0;
        stats->pinfres = __NaN();
        if (stats->hz + stats->by < 0.0) {
            stats->pinfres = stats->hresx / stats->resx0 / (-stats->hz - stats->by);
        }
        stats->dinfres = __NaN();
        if (stats->cx < 0.0) {
            stats->dinfres = MAXF((stats->hresy/stats->resy0),  (stats->hresz/stats->resz0)) / (-stats->cx);
        }

        if (opts && opts->show_progress > 0) {
            if (iter == 0) {
                fprintf(stderr, "%10s %12s %9s %8s %7s %5s\n",
                    "pcost", "dcost", "gap", "pres", "dres", "k/t");
            }
            fprintf(stderr, "%2d: %11.4e %11.4e %6.0e %7.0e %7.0e %7.0e\n",
                    iter, stats->pcost, stats->dcost, stats->gap, stats->pres,
                    stats->dres, cpi->kappa/cpi->tau);
        }
        // ---------------------------------------------------------------------
        // test for stopping criteria

        if (stats->pres <= feastol &&
            stats->dres <= feastol &&
            (stats->gap <= abstol ||
             (!isnan(stats->relgap) && stats->relgap < reltol))) {

            return cvxc_conelp_ready(cp, stats, iter, CVXC_STAT_OPTIMAL);
        }
        else if (! isnan(stats->pinfres) && stats->pinfres <= feastol) {
            return cvxc_conelp_ready(cp, stats, iter, CVXC_STAT_PRIMAL_INFEASIBLE);
        }
        else if (! isnan(stats->dinfres) && stats->dinfres <= feastol) {
            return cvxc_conelp_ready(cp, stats, iter, CVXC_STAT_DUAL_INFEASIBLE);
        }

        // -----------------------------------------------------------------------
        // Compute initial scaling W:
        //
        //     W * z = W^{-T} * s = lambda
        //     dg * tau = 1/dg * kappa = lambdag.
        if (iter == 0) {
            cvxc_compute_scaling(&cp->W, &cpi->s_g, &cpi->z_g, &cpi->lmbda_g, &cp->work);
            //     dg = sqrt( kappa / tau )
            //     dgi = sqrt( tau / kappa )
            //     lambda_g = sqrt( tau * kappa )
            //
            // lambda_g is stored in the last position of lmbda.
            cpi->dg  = SQRT(cpi->kappa / cpi->tau);
            cpi->dgi = SQRT(cpi->tau / cpi->kappa);
            cvxm_set(&cpi->lmbda, cp->cdim_diag, 0, SQRT(cpi->tau * cpi->kappa));
        }

        // lmdasq = lmda o lmbda
        cvxc_ssqr(&cpi->lmbdasq_g, &cpi->lmbda_g/*, 0*/);
        cvxc_float_t lambda_g0 = cvxm_get(&cpi->lmbda, cp->cdim_diag, 0);
        cvxm_set(&cpi->lmbdasq, cp->cdim_diag, 0, lambda_g0 * lambda_g0);

        // Solve
        //
        //     [ 0   A'  G'    ] [ x1        ]          [ c ]
        //     [-A   0   0     ]*[ y1        ] = -dgi * [ b ].
        //     [-G   0   W'*W  ] [ W^{-1}*z1 ]          [ h ]

        cvxc_kktfactor(cp->solver, &cp->W, __cvxnil, __cvxnil);

        cvxm_copy(&cpi->x1, cp->c, CVXC_ALL);
        cvxm_scale(&cpi->x1, -1.0, CVXC_ALL);
        cvxm_copy(&cpi->y1, cp->b, CVXC_ALL);
        cvxm_copy(&cpi->z1, cp->h, CVXC_ALL);

        int err = cvxc_kktsolve(cp->solver, &cpi->x1, &cpi->y1, &cpi->z1_g);

        cvxm_scale(&cpi->x1, cpi->dgi, CVXC_ALL);
        cvxm_scale(&cpi->y1, cpi->dgi, CVXC_ALL);
        cvxm_scale(&cpi->z1, cpi->dgi, CVXC_ALL);

        if (err < 0) {
            if (iter == 0 && ! primalstart && ! dualstart) {
                cp->error = CVXC_ERR_RANK;
                return -1;
            } else {
                cp->error = CVXC_ERR_SINGULAR;
                return cvxc_conelp_ready(cp, stats, iter, CVXC_STAT_SINGULAR);
            }
        }

        cvxm_copy(&cpi->th, cp->h, CVXC_ALL);
        cvxc_scale(&cpi->th_g, &cp->W, CVXC_TRANS|CVXC_INV, &cp->work);

        cvxc_float_t nrm, mu, sigma, step, tt, tk;

        // TODO: check this!!!
        nrm = cvxm_nrm2(&cpi->lmbda);
        mu  = (nrm * nrm) / (1.0 + (cvxc_float_t)cp->cdim_diag);
        sigma = 0.0;

        for (int i = 0; i < 2; i++) {
            // Solve
            //
            // [ 0         ]   [  0   A'  G'  c ] [ dx        ]
            // [ 0         ]   [ -A   0   0   b ] [ dy        ]
            // [ W'*ds     ] - [ -G   0   0   h ] [ W^{-1}*dz ]
            // [ dg*dkappa ]   [ -c' -b' -h'  0 ] [ dtau/dg   ]
            //
            //               [ rx   ]
            //               [ ry   ]
            // = - (1-sigma) [ rz   ]
            //               [ rtau ]
            //
            // lmbda o (dz + ds) = -lmbda o lmbda + sigma*mu*e
            // lmbdag * (dtau + dkappa) = - kappa * tau + sigma*mu
            //
            // ds = -lmbdasq                          if i == 0
            //    = -lmbdasq - dsa o dza + sigma*mu*e if i == 1
            // dkappa = -lambdasq[-1]                            if i == 0 
            //        = -lambdasq[-1] - dkappaa*dtaua + sigma*mu if i == 1.

            cvxc_mgrp_copy_lambda(&cpi->ds_g, &cpi->lmbdasq_g);
            cpi->dkappa = cvxm_get(&cpi->lmbdasq, cp->cdim_diag, 0);

            if (i == 1) {
                // ws3 holds ds o dz ; wkappa3 holds dtau*dkappa 
                cvxm_axpy(&cpi->ds, 1.0, &cpi->ws3);

                // update ds with -sigma*mu
                cvxc_mgrp_update_sz(&cpi->ds_g, -sigma*mu, 0);
                cpi->dkappa = cpi->dkappa + cpi->wkappa3 - sigma*mu;
            }
            // (dx, dy, dz, dtau) = (1-sigma)*(rx, ry, rz, rt)
            cvxm_copy(&cpi->dx, &cpi->rx, CVXC_ALL);
            cvxm_scale(&cpi->dx, 1.0 - sigma, CVXC_ALL);
            cvxm_copy(&cpi->dy, &cpi->ry, CVXC_ALL);
            cvxm_scale(&cpi->dy, 1.0 - sigma, CVXC_ALL);
            cvxm_copy(&cpi->dz, &cpi->rz, CVXC_ALL);
            cvxm_scale(&cpi->dz, 1.0 - sigma, CVXC_ALL);
            cpi->dtau = (1.0 - sigma) * stats->rt;

            f6(cp, &cpi->dx, &cpi->dy,
               &cpi->dz_g, &cpi->dtau, &cpi->ds_g, &cpi->dkappa, refinement);

            if (i == 0) {
                // Save ds o dz and dkappa * dtau for Mehrotra correction
                cvxm_copy(&cpi->ws3, &cpi->ds, CVXC_ALL);
                cvxc_sprod(&cpi->ws3_g, &cpi->dz_g, 0, &cp->work);
                cpi->wkappa3 = cpi->dtau * cpi->dkappa;
            }

            // Maximum step to boundary.
            //
            // If i is 1, also compute eigenvalue decomposition of the 's'
            // blocks in ds, dz.  The eigenvectors Qs, Qz are stored in
            // ds_k, dz_k.  The eigenvalues are stored in sigs, sigz.

            cvxc_scale2(&cpi->ds_g, &cpi->lmbda_g, 0, &cp->work);
            cvxc_scale2(&cpi->dz_g, &cpi->lmbda_g, 0, &cp->work);

            if (i == 0) {
                cpi->ts = cvxc_max_step(&cpi->ds_g, __nilgrp, &cp->work);
                cpi->tz = cvxc_max_step(&cpi->dz_g, __nilgrp, &cp->work);
            } else {
                cpi->ts = cvxc_max_step(&cpi->ds_g, &cpi->sigs_g, &cp->work);
                cpi->tz = cvxc_max_step(&cpi->dz_g, &cpi->sigz_g, &cp->work);
            }

            tt = -cpi->dtau / cvxm_get(&cpi->lmbda, cp->cdim_diag, 0);
            tk = -cpi->dkappa / cvxm_get(&cpi->lmbda, cp->cdim_diag, 0);
            cvxc_float_t t = cvxc_maxvec(5, (cvxc_float_t[]){0.0, cpi->ts, cpi->tz, tt, tk});
            if (t == 0.0) {
                step = 1.0;
            } else {
                if (i == 0) {
                    step = MINF(1.0, 1.0/t);
                } else {
                    step = MINF(1.0, STEP/t);
                }
            }
            if (i == 0) {
                // sigma = (1 -step)^3
                sigma = (1.0 - step) * (1.0 - step) * (1.0 - step);
            }
        }
        // update x, y
        cvxm_axpy(&cpi->x, step, &cpi->dx);
        cvxm_axpy(&cpi->y, step, &cpi->dy);

        // Replace 'l' and 'q' blocks of ds and dz with the updated
        // variables in the current scaling.
        // Replace 's' blocks of ds and dz with the factors Ls, Lz in a
        // factorization Ls*Ls', Lz*Lz' of the updated variables in the
        // current scaling.
        //
        // ds := e + step*ds for 'l' and 'q' blocks.
        // dz := e + step*dz for 'l' and 'q' blocks.
        cvxc_mgrp_scale_sz(&cpi->ds_g, step, CVXDIM_LINEAR|CVXDIM_SOCP);
        cvxc_mgrp_update_sz(&cpi->ds_g, 1.0, CVXDIM_LINEAR|CVXDIM_SOCP);

        cvxc_mgrp_scale_sz(&cpi->dz_g, step, CVXDIM_LINEAR|CVXDIM_SOCP);
        cvxc_mgrp_update_sz(&cpi->dz_g, 1.0, CVXDIM_LINEAR|CVXDIM_SOCP);

        // ds := H(lambda)^{-1/2} * ds and dz := H(lambda)^{-1/2} * dz.
        //
        // This replaces the 'l' and 'q' components of ds and dz with the
        // updated variables in the current scaling.
        // The 's' components of ds and dz are replaced with
        //
        // diag(lmbda_k)^{1/2} * Qs * diag(lmbda_k)^{1/2}
        // diag(lmbda_k)^{1/2} * Qz * diag(lmbda_k)^{1/2}

        cvxc_scale2(&cpi->ds_g, &cpi->lmbda_g, CVXC_INV, &cp->work);
        cvxc_scale2(&cpi->dz_g, &cpi->lmbda_g, CVXC_INV, &cp->work);

        // sigs := ( e + step*sigs ) ./ lambda for 's' blocks.
        // sigz := ( e + step*sigz ) ./ lambda for 's' blocks.
        cvxm_scale(&cpi->sigs, step, CVXC_ALL);
        cvxm_scale(&cpi->sigz, step, CVXC_ALL);
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

        // TODO: missing somethings..!!!!
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
        cvxc_update_scaling(&cp->W, &cpi->lmbda_g, &cpi->ds_g, &cpi->dz_g, &cp->work);

        // For kappa, tau block:
        //
        //     dg := sqrt( (kappa + step*dkappa) / (tau + step*dtau) )
        //         = dg * sqrt( (1 - step*tk) / (1 - step*tt) )
        //
        //     lmbda[-1] := sqrt((tau + step*dtau) * (kappa + step*dkappa))
        //                = lmbda[-1] * sqrt(( 1 - step*tt) * (1 - step*tk))
        cpi->dg *= SQRT(1.0 - step*tk)/SQRT(1.0-step*tt);
        cpi->dgi = 1.0/cpi->dg;
        cvxc_float_t a = SQRT(1.0 - step*tk) * SQRT(1.0-step*tt);
        cvxc_float_t g = cvxm_get(&cpi->lmbda, cp->cdim_diag, 0);
        cvxm_set(&cpi->lmbda, cp->cdim_diag, 0, a*g);

        // Unscale s, z, tau, kappa (unscaled variables are used only to
        // compute feasibility residuals).
        cvxc_mgrp_copy_lambda(&cpi->s_g, &cpi->lmbda_g);
        cvxc_scale(&cpi->s_g, &cp->W, CVXC_TRANS, &cp->work);

        cvxc_mgrp_copy_lambda(&cpi->z_g, &cpi->lmbda_g);

        cvxc_scale(&cpi->z_g, &cp->W, CVXC_INV, &cp->work);

        cpi->kappa = cvxm_get(&cpi->lmbda, cp->cdim_diag, 0) / cpi->dgi;
        cpi->tau   = cvxm_get(&cpi->lmbda, cp->cdim_diag, 0) * cpi->dgi;
        cvxm_map_data(&lk, cp->cdim_diag, 1, cvxm_data(&cpi->lmbda, 0));
        g = cvxm_nrm2(&lk) / cpi->tau;
        stats->gap = g*g;
    }

    return cvxc_conelp_ready(cp, stats, maxiter, CVXC_STAT_UNKNOWN);
}
