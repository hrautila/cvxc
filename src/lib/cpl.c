
// Copyright: Harri Rautila, 2018-2019 <harri.rautila@gmail.com>

#include "epi.h"
#include "convex.h"

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
cvx_size_t __WORKBYTES(cvx_size_t n)
{
    return 8*n*n*sizeof(cvx_float_t);
}


// Evaluates residuals in Newton equations:
//
//     [ vx ]     [ 0     ]   [ H  A' GG' ] [ ux        ]
//     [ vy ] -=  [ 0     ] + [ A  0  0   ] [ uy        ]
//     [ vz ]     [ W'*us ]   [ GG 0  0   ] [ W^{-1}*uz ]
//
//     vs -= lmbda o (uz + us).
static
int cpl_res(cvx_problem_t *cp,
            cvx_matrix_t *ux,
            cvx_matrix_t *uy,
            cvx_matgrp_t *uz_g,
            cvx_matgrp_t *us_g,
            cvx_matrix_t *vx,
            cvx_matrix_t *vy,
            cvx_matgrp_t *vz_g,
            cvx_matgrp_t *vs_g,
            cvx_scaling_t *W,
            cvx_matgrp_t *lmbda_g)
{
    int err = 0;
    cvx_cpl_internal_t *cpi = &cp->u.cpl;

    cvx_matrix_t *us = us_g->mat, *uz = uz_g->mat;
    cvx_matrix_t *vs = vs_g->mat, *vz = vz_g->mat;

    cvx_matgrp_t ux_g = (cvx_matgrp_t){ .mat = ux, .index = (cvx_index_t *)0 };

    // vx = vx - H*ux - A^T*uy - G^T*W^-1*uz
    cvxm_mult(1.0, vx, -1.0, &cpi->H, ux, 0);
    cvxm_mult(1.0, vx, -1.0, cp->A, uy, CVX_TRANSA);
    cvxm_copy(&cpi->wz3, uz, CVX_ALL);
    cvx_scale(&cpi->wz3_g, &cp->W, CVX_INV, &cp->work);
    cvx_sgemv2(1.0, vx, -1.0, &cpi->Df, cp->G, &cpi->wz3_g, CVX_TRANS);

    // vy = vy - A*ux
    cvxm_mult(1.0, vy, -1.0, cp->A, ux, 0);

    // vz = vz - W'*us - GG*ux
    cvxm_copy(&cpi->ws3, us, 0);
    cvx_scale(&cpi->ws3_g, &cp->W, CVX_TRANS, &cp->work);
    cvxm_axpy(vz, -1.0, &cpi->ws3);
    cvx_sgemv2(1.0, vz, -1.0, &cpi->Df, cp->G, &ux_g, 0);

    //cvx_mat_printf(stdout, "%e", vz, "res_entry:update(vz)");

    // vs = vs - lmbda o (uz + us)
    cvxm_copy(&cpi->ws3, us, 0);
    cvxm_axpy(&cpi->ws3, 1.0, uz);
    cvx_sprod(&cpi->ws3_g, lmbda_g, CVX_DIAG, &cp->work);
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
int f4_no_ir(cvx_problem_t *cp,
             cvx_matrix_t *x,
             cvx_matrix_t *y,
             cvx_matgrp_t *z_g,
             cvx_matgrp_t *s_g)
{
    cvx_matrix_t *s = s_g->mat, *z = z_g->mat;
    cvx_cpl_internal_t *cpi = &cp->u.cpl;

    // cvx_mat_printf(stdout, "%e", x, "f4_no_ir:in x");

    // cvx_mat_printf(stdout, "%.7f", s, "f4_no_ir:in s");
    // cvx_mat_printf(stdout, "%e", z, "f4_no_ir:in z");
    // s = lmbda o\ s
    // cvx_mat_printf(stdout, "%.7f", cpi->lmbda_g.mat, "f4_no_ir:lmbda)");
    cvx_sinv(s_g, &cpi->lmbda_g, &cp->work);
    // cvx_mat_printf(stdout, "%.7f", s_g->mat, "f4_no_ir.sinv(s)");

    // z := z - W'*s
    cvxm_copy(&cpi->ws3, s, 0);
    cvx_scale(&cpi->ws3_g, &cp->W, CVX_TRANS, &cp->work);
    cvxm_axpy(z, -1.0, &cpi->ws3);
    //cvx_mat_printf(stdout, "%.7f", &cpi->ws3, "pre-solve f4_no_ir ws3");
    // cvx_mat_printf(stdout, "%.7f", z, "pre-solve f4_no_ir z");
    // solve
    cvx_kktsolve(cp->solver, x, y, z_g);
    //cvx_mat_printf(stdout, "%e", z, "post-solve f4_no_ir z");

    // s := s - z
    cvxm_axpy(s, -1.0, z);

    return 0;
}


// f6(x, y, z, tau, s, kappa) solves the same system as f6_no_ir,
// but applies iterative refinement.
static
int f4(cvx_problem_t *cp,
       cvx_matrix_t *x,
       cvx_matrix_t *y,
       cvx_matgrp_t *z_g,
       cvx_matgrp_t *s_g,
       int refinement)
{
    cvx_matrix_t *s = s_g->mat, *z = z_g->mat;
    cvx_cpl_internal_t *cpi = &cp->u.cpl;
    int err = 0;

    //printf("**f6:start tau=%e, kappa=%e, refinement=%d\n", *tau, *kappa, refinement);

    if (refinement > 0 /*|| (flags | CVX_DEBUG) != 0*/) {
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

        //printf("** refinement: %d\n", i);

        //cvx_mat_printf(stdout, "%.7f", &cpi->wx2, "f4:wx2");
        //cvx_mat_printf(stdout, "%.7f", &cpi->wz2, "f4:wz2");
        //cvx_mat_printf(stdout, "%.7f", &cpi->ws2, "f4:ws2");

        cpl_res(cp, x, y, z_g, s_g, &cpi->wx2, &cpi->wy2, &cpi->wz2_g,
                &cpi->ws2_g, &cp->W, &cpi->lmbda_g);

        //cvx_mat_printf(stdout, "%.7f", &cpi->wx2, "f4.res:wx2");
        //cvx_mat_printf(stdout, "%.7f", &cpi->wy2, "f4.res:wy2");
        //cvx_mat_printf(stdout, "%.7f", &cpi->wz2, "f4.res:wz2");
        //cvx_mat_printf(stdout, "%.7f", &cpi->ws2, "f4.res:ws2");

        f4_no_ir(cp, &cpi->wx2, &cpi->wy2, &cpi->wz2_g, &cpi->ws2_g);

        //cvx_mat_printf(stdout, "%.7f", &cpi->wx2, "f4_no_ir:wx2");
        //cvx_mat_printf(stdout, "%.7f", &cpi->wz2, "f4_no_ir:wz2");
        //cvx_mat_printf(stdout, "%.7f", &cpi->ws2, "f4_no_ir:ws2");

        cvxm_axpy(x, 1.0, &cpi->wx2);
        cvxm_axpy(y, 1.0, &cpi->wy2);
        cvxm_axpy(s, 1.0, &cpi->ws2);
        cvxm_axpy(z, 1.0, &cpi->wz2);
    }
    //printf("**f6:end tau=%e, kappa=%e, refinement=%d\n", *tau, *kappa, refinement);
    return err;
}

static
int cpl_factor(cvx_kktsolver_t *S, cvx_scaling_t *W, cvx_matrix_t *x, cvx_matrix_t *z)
{
    cvx_chainsolver_t *cpl_solver = (cvx_chainsolver_t *)S;
    cvx_problem_t *cp = cpl_solver->cp;
    cvx_cpl_internal_t *cpi = &cp->u.cpl;
    F2(cp->F, __cvxnil, &cpi->Df, &cpi->H, x, z);
    return cvx_kktfactor(cpl_solver->next, W, &cpi->H, &cpi->Df);
}

static
int cpl_solve(cvx_kktsolver_t *S, cvx_matrix_t *x, cvx_matrix_t *y, cvx_matgrp_t *z_g)
{
    cvx_chainsolver_t *cpl_solver = (cvx_chainsolver_t *)S;
    return cvx_kktsolve(cpl_solver->next, x, y, z_g);
}

static cvx_kktfuncs_t cpl_kktfunctions = {
    .factor = cpl_factor,
    .solve  = cpl_solve
};

cvx_chainsolver_t *cvx_cpl_solver_init(cvx_chainsolver_t *cs, cvx_problem_t *cp, cvx_kktsolver_t *next)
{
    cs->fnc = &cpl_kktfunctions;
    cs->next = next;
    cs->cp = cp;
    return cs;
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
cvx_size_t cvx_cpl_bytes(int n, int m, const cvx_dimset_t *dims, int nonlinear)
{
    cvx_size_t cdim_mnl  = cvx_dimset_sum_squared(dims, CVXDIM_CONVEXLP);
    cvx_size_t cdim_diag = cvx_dimset_sum(dims, CVXDIM_CONELP);
    cvx_size_t sdim      = cvx_dimset_sum(dims, CVXDIM_SDP);
    cvx_size_t maxsdp    = cvx_dimset_max(dims, CVXDIM_SDP);
    cvx_size_t mnl       = cvx_dimset_sum(dims, CVXDIM_NONLINEAR);
    cvx_size_t total = 0;
    cvx_size_t nbytes = 0;
    if (nonlinear != 0)
        nonlinear = 1;
    mnl += nonlinear;
    cdim_mnl += nonlinear;

    // size of standard index set (if nonlinear then 2 times)
    nbytes += cvx_index_bytes(dims, CVX_INDEX_NORMAL) * (1 + nonlinear);
    // size of packed index set
    nbytes += cvx_index_bytes(dims, CVX_INDEX_PACKED);
    // size of diagonal index set
    nbytes += cvx_index_bytes(dims, CVX_INDEX_DIAG);
    // size of SDP diagonal index set
    nbytes += cvx_index_bytes(dims, CVX_INDEX_SIGS);

    //fprintf(stderr, "cpl size:  indexes %ld\n", nbytes);

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

    //fprintf(stderr, "cpl size: vectors    %ld\n", nbytes+total*sizeof(cvx_float_t));

    total += 2*(cdim_diag+mnl) + 2;   // for lmbda, lmbdasq
    total += 2*sdim;            // for sigs, sigz

    //fprintf(stderr, "cpl size: eigen space %ld\n", nbytes+total*sizeof(cvx_float_t));

    // f, Df, H
    total += mnl*2;   // for f, newf
    total += mnl*n*2; // for Df, newDf
    total += n*n;               // for H

    //fprintf(stderr, "cpl size: f, Df, H %ld  [mnl:%ld, nl:%d, n:%d]\n",
    //      nbytes+total*sizeof(cvx_float_t), mnl, nonlinear, n);

    nbytes += total*sizeof(cvx_float_t);

    // workspace for SDP constraint scaling; TODO: think about this.
    if (maxsdp > 0)
        nbytes += __WORKBYTES(maxsdp);

    cvx_size_t isize, szw;
    // calculte space need for scaling matrix (second for state save)
    szw     = cvx_scaling_bytes(&isize, dims);
    nbytes += 2*szw;

    //fprintf(stderr, "cpl size: scaling %ld [%ld]\n", nbytes, szw);

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
        cvxm_map_data(var, 0, 1, (cvx_float_t *)0);              \
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
int cvxm_xvector(int nl, cvx_matrix_t *x, cvx_size_t m, cvx_size_t n, void *data, cvx_size_t ndata)
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
cvx_size_t cvx_cpl_make(cvx_problem_t *cp,
                        int n,
                        int m,
                        const cvx_dimset_t *dims,
                        int nl,
                        void *memory,
                        cvx_size_t nbytes)
{
    cvx_cpl_internal_t *cpi = &cp->u.cpl;
    cvx_size_t offset = 0;
    cvx_size_t used = 0;
    unsigned char *bytes = (unsigned char *)memory;

    if (nl != 0)
        nl = 1;

    cvx_size_t cdim_mnl  = cvx_dimset_sum_squared(dims, CVXDIM_CONVEXLP) + nl;

    cp->cdim = cdim_mnl;

    // __INIT macro assumes variables offset and nbytes;
    // overlay index sets
    __INIT(used, cvx_index_make(&cpi->index_full, dims, CVX_INDEX_NORMAL, &bytes[offset],  nbytes));
    __INIT(used, cvx_index_make(&cpi->index_packed, dims, CVX_INDEX_PACKED, &bytes[offset],  nbytes));
    __INIT(used, cvx_index_make(&cpi->index_diag, dims, CVX_INDEX_DIAG, &bytes[offset],  nbytes));
    __INIT(used, cvx_index_make(&cpi->index_sig, dims, CVX_INDEX_SIGS, &bytes[offset],  nbytes));

    if (nl) {
        cvx_dimset_t ldims = *dims;
        ldims.iscpt = 0;
        __INIT(used, cvx_index_make(&cpi->index_cpt, &ldims, CVX_INDEX_NORMAL, &bytes[offset],  nbytes));
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

    //fprintf(stderr, "cpl make: vectors %ld\n", offset);

    cvx_size_t cdim_diag =
        cvx_dimset_sum(dims, CVXDIM_CONVEXLP) + nl;

    cp->cdim_diag = cdim_diag;

    // lmbda, lmbdasq
    __INIT(used, cvxm_make(&cpi->lmbda,   cdim_diag, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->lmbdasq, cdim_diag, 1, &bytes[offset], nbytes));

    cvx_size_t sdim   = cvx_dimset_sum(dims, CVXDIM_SDP);
    cvx_size_t maxsdp = cvx_dimset_max(dims, CVXDIM_SDP);

    // sigs, sigz; space for eigenvalues; zero if no SDP constraints
    __INITC(used, sdim, &cpi->sigs, cvxm_make(&cpi->sigs, sdim, 1, &bytes[offset], nbytes));
    __INITC(used, sdim, &cpi->sigz, cvxm_make(&cpi->sigz, sdim, 1, &bytes[offset], nbytes));

    //fprintf(stderr, "cpl make: eigen space %ld\n", offset);

    // f, Df, H matrices
    int mnl = cvx_dimset_sum(dims, CVXDIM_NONLINEAR) + nl;
    __INIT(used, cvxm_make(&cpi->f, mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->newf, mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->Df, mnl, n, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->newDf, mnl, n, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->H, n, n, &bytes[offset], nbytes));

    //fprintf(stderr, "cpl make: f, Df, H %ld [mnl:%d, nl:%d, n:%d]\n", offset, mnl, nl, n);

    // scaling matrix
    __INIT(used, cvx_scaling_make(&cp->W,   dims, &bytes[offset], nbytes));
    __INIT(used, cvx_scaling_make(&cpi->W0, dims, &bytes[offset], nbytes));

    // workspace for SDP contraints handling
    __mblk_empty(&cp->work);
    if (maxsdp > 0) {
        __INIT(used, __mblk_make(&cp->work, __WORKBYTES(maxsdp), &bytes[offset], nbytes));
    }

    // setup matrix group variables
    //cvx_mgrp_init(&cpi->h_g,   cp->h,    &cpi->index_full);
    cvx_mgrp_init(&cpi->s_g,   &cpi->s,   &cpi->index_full);
    cvx_mgrp_init(&cpi->z_g,   &cpi->z,   &cpi->index_full);
    cvx_mgrp_init(&cpi->z0_g,  &cpi->z0,  &cpi->index_full);
    cvx_mgrp_init(&cpi->s0_g,  &cpi->s0,  &cpi->index_full);
    cvx_mgrp_init(&cpi->ds_g,  &cpi->ds,  &cpi->index_full);
    cvx_mgrp_init(&cpi->dz_g,  &cpi->dz,  &cpi->index_full);
    cvx_mgrp_init(&cpi->ds0_g, &cpi->ds0, &cpi->index_full);
    cvx_mgrp_init(&cpi->dz0_g, &cpi->dz0, &cpi->index_full);
    cvx_mgrp_init(&cpi->ds2_g, &cpi->ds2, &cpi->index_full);
    cvx_mgrp_init(&cpi->dz2_g, &cpi->dz2, &cpi->index_full);
    cvx_mgrp_init(&cpi->ds20_g,&cpi->ds20,&cpi->index_full);
    cvx_mgrp_init(&cpi->dz20_g,&cpi->dz20,&cpi->index_full);
    cvx_mgrp_init(&cpi->rz_g,  &cpi->rz,  &cpi->index_full);
    cvx_mgrp_init(&cpi->newz_g,&cpi->newz,&cpi->index_full);
    cvx_mgrp_init(&cpi->news_g,&cpi->news,&cpi->index_full);
    cvx_mgrp_init(&cpi->ws_g,  &cpi->ws,  &cpi->index_full);
    cvx_mgrp_init(&cpi->wz_g,  &cpi->wz,  &cpi->index_full);
    cvx_mgrp_init(&cpi->ws2_g, &cpi->ws2, &cpi->index_full);
    cvx_mgrp_init(&cpi->wz2_g, &cpi->wz2, &cpi->index_full);
    cvx_mgrp_init(&cpi->ws3_g, &cpi->ws3, &cpi->index_full);
    cvx_mgrp_init(&cpi->wz3_g, &cpi->wz3, &cpi->index_full);

    cvx_mgrp_init(&cpi->newrz_g,  &cpi->newrz0,&cpi->index_full);

    cvx_mgrp_init(&cpi->sigs_g,   &cpi->sigs,   &cpi->index_sig);
    cvx_mgrp_init(&cpi->sigz_g,   &cpi->sigz,   &cpi->index_sig);

    cvx_mgrp_init(&cpi->lmbda_g,   &cpi->lmbda,   &cpi->index_diag);
    cvx_mgrp_init(&cpi->lmbdasq_g, &cpi->lmbdasq, &cpi->index_diag);

    cp->f  = &cpi->f;
    cp->Df = &cpi->Df;
    cp->H  = &cpi->H;

    cp->mlen = nbytes;
    cp->memory = memory;

    return offset;
}

int cvx_cpl_isok(const cvx_matrix_t *c,
                 const cvx_convex_program_t *F,
                 const cvx_matrix_t *G,
                 const cvx_matrix_t *h,
                 const cvx_matrix_t *A,
                 const cvx_matrix_t *b,
                 const cvx_dimset_t *dims)
{
    cvx_size_t mc, nc, mG, nG, mA, nA, mh, nh, mb, nb;

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
        return CVX_ERR_DIMC;
    }

    cvx_size_t cdim      = cvx_dimset_sum_squared(dims, CVXDIM_CONELP);
    cvx_size_t cdim_pckd = cvx_dimset_sum_packed(dims, CVXDIM_CONELP);

    if (nh > 1 || mh != cdim) {
        return CVX_ERR_DIMH;
    }

    if (mG != cdim || (c && nG != mc)) {
        return CVX_ERR_DIMG;
    }
    if (nA != nG || mA != mb) {
        return CVX_ERR_DIMA;
    }
    if (nb != 1) {
        return CVX_ERR_DIMB;
    }
    if ( mb > mc || mb + cdim_pckd < mc) {
        return CVX_ERR_RANK;
    }
    return 0;
}

int cvx_cpl_setvars(cvx_problem_t *cp,
                    cvx_convex_program_t *F,
                    cvx_size_t n, cvx_size_t m,
                    cvx_matrix_t *c,
                    cvx_matrix_t *G,
                    cvx_matrix_t *h,
                    cvx_matrix_t *A,
                    cvx_matrix_t *b,
                    const cvx_dimset_t *dims,
                    cvx_kktsolver_t *kktsolver)
{
    cp->c = c;
    cp->G = G;
    cp->h = h;
    cp->A = A;
    cp->b = b;
    cp->F = F;

    cp->primal_x = (cvx_matrix_t *)0;
    cp->primal_s = (cvx_matrix_t *)0;
    cp->dual_y = (cvx_matrix_t *)0;
    cp->dual_z = (cvx_matrix_t *)0;

    return 0;
}

cvx_size_t cvx_cpl_allocate(cvx_problem_t *cp,
                            int nl,
                            cvx_size_t n,
                            cvx_size_t m,
                            cvx_size_t extra,
                            const cvx_dimset_t *dims)
{
    cvx_size_t used, nbytes = cvx_cpl_bytes(n, m, dims, nl);
    if (extra)
        nbytes += extra;
    void *memory = calloc(nbytes, 1);
    if (!memory) {
        cp->error = CVX_ERR_MEMORY;
        return 0;
    }

    if ((used = cvx_cpl_make(cp, n, m, dims, nl, memory, nbytes)) == 0) {
        cp->error = CVX_ERR_MEMORY;
        free(memory);
        return 0;
    }
    return used;
}

cvx_size_t cvx_cpl_setup(cvx_problem_t *cp,
                         cvx_convex_program_t *F,
                         cvx_matrix_t *c,
                         cvx_matrix_t *G,
                         cvx_matrix_t *h,
                         cvx_matrix_t *A,
                         cvx_matrix_t *b,
                         const cvx_dimset_t *dims,
                         cvx_kktsolver_t *kktsolver)
{
    cvx_size_t mc, nc, mb, nb;
    cvx_size_t used;
    int err;

    if (! cp)
        return 0;

    if ((err = cvx_cpl_isok(c, F, G, h, A, b, dims)) != 0) {
        cp->error = err;
        return 0;
    }

    mc = nc = mb = nb = 0;
    cvxm_size(&mc, &nc, c);
    if (b)
        cvxm_size(&mb, &nb, b);

    if ((used = cvx_cpl_allocate(cp, 0, mc, mb, 0, dims)) == 0) {
        return 0;
    }
    cvx_cpl_setvars(cp, F, mc, mb, c, G, h, A, b, dims, kktsolver);

    cvx_cpl_internal_t *cpi = &cp->u.cpl;
    // provide full index set
    cp->index_g = &cpi->index_full;
    cvx_mgrp_init(&cpi->h_g, cp->h, cp->index_g);

    // init KKT solver
    if (kktsolver) {
        cp->solver = kktsolver;
    } else {
        cvx_ldlsolver_init(&cp->__S, cp, mc, mb, dims);
        cp->solver = &cp->__S;
    }
    cvx_cpl_solver_init(&cpi->cp_solver, cp, cp->solver);
    cp->solver = (cvx_kktsolver_t *)&cpi->cp_solver;

    return used;
}

int cvx_cpl_ready(cvx_problem_t *cp,
                  //cvx_stats_t *stats,
                  int iter,
                  int stat)
{
    cvx_cpl_internal_t *cpi = &cp->u.cpl;
    //cvx_stats_t *stats = &prob->stats;

    if (stat == CVX_STAT_OPTIMAL && iter == -1) {
        // constructed initial point is feasible and optimal
        cvx_mksymm(&cpi->s_g);
        cvx_mksymm(&cpi->z_g);

        // rx = A'*y + G'z + c
        cvxm_copy(&cpi->rx, cp->c, CVX_ALL);
        cvxm_mvmult(1.0, &cpi->rx, 1.0, cp->A, &cpi->y, CVX_TRANS);
        cvxm_mvmult(1.0, &cpi->rx, 1.0, cp->G, &cpi->z, CVX_TRANS);
        cpi->resx = cvxm_nrm2(&cpi->rx);

        // ry = b - A*x  ; TODO - computes -b - A*x ;; check
        cvxm_copy(&cpi->ry, cp->b, CVX_ALL);
        cvxm_mvmult(-1.0, &cpi->ry, -1.0, cp->A, &cpi->x, CVX_TRANS);
        cpi->resy = cvxm_nrm2(&cpi->ry);

        // rz = s + G*x - h
        cvxm_copy(&cpi->rz, &cpi->s, CVX_ALL);
        cvxm_mvmult(1.0, &cpi->rz, 1.0, cp->G, &cpi->x, 0);
        cvxm_axpy(&cpi->rz, 1.0, cp->h);
        cpi->resz = cvx_snrm2(&cpi->rz_g);

        cpi->pres = __MAX2(cpi->resy/cpi->resy0, cpi->resz/cpi->resz0);
        cpi->dres = cpi->resx/cpi->resx0;
        cpi->cx = cvxm_dot(cp->c, &cpi->x);
        cpi->by = cvxm_dot(cp->b, &cpi->y);
        cpi->hz = cvx_sdot(&cpi->h_g, &cpi->z_g);

        cp->solution.x = &cpi->x;
        cp->solution.s = &cpi->s;
        cp->solution.y = &cpi->y;
        cp->solution.z = &cpi->z;

        cp->solution.status = stat;
        cp->solution.gap = cpi->gap;
        cp->solution.relative_gap = cpi->relgap;
        cp->solution.primal_objective = cpi->cx;
        cp->solution.dual_objective = - (cpi->by + cpi->hz);
        cp->solution.primal_infeasibility = cpi->pres;
        cp->solution.dual_infeasibility = cpi->dres;
        cp->solution.primal_slack = - cpi->ts;
        cp->solution.dual_slack = - cpi->tz;
        cp->solution.primal_residual_cert = __NaN();
        cp->solution.dual_residual_cert = __NaN();
        cp->solution.iterations = 0;
    }
    else if (stat == CVX_STAT_UNKNOWN || stat == CVX_STAT_OPTIMAL) {
        cvxm_scale(&cpi->x, 1.0/cpi->tau, CVX_ALL);
        cvxm_scale(&cpi->y, 1.0/cpi->tau, CVX_ALL);
        cvxm_scale(&cpi->s, 1.0/cpi->tau, CVX_ALL);
        cvxm_scale(&cpi->z, 1.0/cpi->tau, CVX_ALL);
        cvx_mksymm(&cpi->s_g);
        cvx_mksymm(&cpi->z_g);
        cpi->ts = cvx_max_step(&cpi->s_g, __nilgrp, &cp->work);
        cpi->tz = cvx_max_step(&cpi->s_g, __nilgrp, &cp->work);

        cp->error = stat == CVX_STAT_UNKNOWN ? CVX_ERR_MAXITER : 0;
        cp->solution.x = &cpi->x;
        cp->solution.s = &cpi->s;
        cp->solution.y = &cpi->y;
        cp->solution.z = &cpi->z;

        cp->solution.status = stat;

        cp->solution.gap = cpi->gap;
        cp->solution.relative_gap = cpi->relgap;
        cp->solution.primal_objective = cpi->cx;
        cp->solution.dual_objective = - (cpi->by + cpi->hz);
        cp->solution.primal_infeasibility = cpi->pres;
        cp->solution.dual_infeasibility = cpi->dres;
        cp->solution.primal_slack = - cpi->ts;
        cp->solution.dual_slack = - cpi->tz;
        if (stat == CVX_STAT_OPTIMAL) {
            cp->solution.primal_residual_cert = __NaN();
            cp->solution.dual_residual_cert = __NaN();
        } else {
            cp->solution.primal_residual_cert = cpi->pinfres;
            cp->solution.dual_residual_cert = cpi->dinfres;
        }
        cp->solution.iterations = iter;
    }
    return -stat;
}

static
void cvx_cpl_save_or_restore(cvx_problem_t *cp, int restore)
{
    cvx_cpl_internal_t *cpi = &cp->u.cpl;
    //cvx_stats_t *stats = &cp->stats;

    if (restore) {
        cvx_scaling_copy(&cp->W, &cpi->W0);
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
    cvx_scaling_copy(&cpi->W0, &cp->W);
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
int cvx_cpl_linesearch(cvx_problem_t *cp,
                       int iter,
                       int relaxed_iters)
{
    cvx_cpl_internal_t *cpi = &cp->u.cpl;
    int backtrack;
    cvx_matrix_t news_nl, newrznl;
    cvx_float_t newresznl, newgap, newphi, newresx;

    // NOTE: sigma changes!!! parameters??

    //printf("theta:  %e, %e, %e\n", cpi->theta1, cpi->theta2, cpi->theta3);
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

    //printf("phi  : %e, dphi:  %e\n", cpi->phi, cpi->dphi);

    for (backtrack = 1; backtrack; backtrack = 0 ) {
        cvxm_copy(&cpi->newx, &cpi->x, 0);
        cvxm_axpy(&cpi->newx, cpi->step, &cpi->dx);
        cvxm_copy(&cpi->newy, &cpi->y, 0);
        cvxm_axpy(&cpi->newy, cpi->step, &cpi->dy);
        cvxm_copy(&cpi->newz, &cpi->z, 0);
        cvxm_axpy(&cpi->newz, cpi->step, &cpi->dz2);
        cvxm_copy(&cpi->news, &cpi->s, 0);
        cvxm_axpy(&cpi->news, cpi->step, &cpi->ds2);

        //cvx_mat_printf(stdout, "%.7f", &cpi->newx, "linesearch newx");
        F1(cp->F, &cpi->newf, &cpi->newDf, &cpi->newx);
        //cvx_mat_printf(stdout, "%.7f", &cpi->newf, "linesearch newf");
        // newrx = c + A'*newy + newDf'*newz[:mnl] + G'*newz[mnl:]
        // ...
        cvxm_copy(&cpi->newrx, cp->c, 0);  // TODO: c == __nil ??
        cvxm_mvmult(1.0, &cpi->newrx, 1.0, cp->A, &cpi->newy, CVX_TRANS);
        cvx_sgemv2(1.0, &cpi->newrx, 1.0, &cpi->newDf, cp->G, &cpi->newz_g, CVX_TRANS);
        newresx = SQRT(cvxm_dot(&cpi->newrx, &cpi->newrx));
        //printf("newresx  : %e\n", newresx);

        // newrznl = news[:mnl] + newf
        // ...
        cvx_mgrp_elem(&news_nl, &cpi->news_g, CVXDIM_CONVEX, 0);
        cvx_mgrp_elem(&newrznl, &cpi->newrz_g, CVXDIM_CONVEX, 0);
        cvxm_copy(&newrznl, &news_nl, 0);
        cvxm_axpy(&newrznl, 1.0, &cpi->newf);
        //cvx_mat_printf(stdout, "%.7f", &newrznl, "linsearch newrznl");
        newresznl = cvxm_nrm2(&newrznl);
        //printf("newresznl: %e\n", newresznl);

        //printf("i: %d, gap   : %e, sigma : %e, step: %e, dsdz: %e\n", iter, cpi->gap, cpi->sigma, cpi->step, cpi->dsdz);
        newgap = (1.0 - (1.0 - cpi->sigma)*cpi->step)*cpi->gap +
            cpi->step*cpi->step*cpi->dsdz;
        newphi = cpi->theta1*newgap + cpi->theta2*newresx + cpi->theta3*newresznl;

        //printf("i: %d, newgap: %e, newphi: %e\n", iter, newgap, newphi);

        if (iter == 0) {
            if (newgap <= (1.0 - ALPHA*cpi->step)*cpi->gap &&
                ((relaxed_iters > 0 && relaxed_iters < MAX_RELAXED_ITERS)
                 || newphi <= cpi->phi+ALPHA*cpi->step*cpi->dphi)) {
                backtrack = 0;
                cpi->sigma = __MIN2(newgap/cpi->gap, __POW((newgap/cpi->gap), EXPON));
                //printf("  break 1: sigma=%.7f [%.7f, %.7f]\n", cpi->sigma, newgap, cpi->gap);
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
                    //printf("  break 2 : newphi=%.7f\n", newphi);
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
                    cvx_cpl_save_or_restore(cp, 0);
                    relaxed_iters = 1;
                }
                backtrack = 0;
                //printf("  break 3 : newphi=%.7f\n", newphi);
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
                //printf("break 4 : newphi=%.7f\n", newphi);
            }
            else if (relaxed_iters == MAX_RELAXED_ITERS && MAX_RELAXED_ITERS > 0) {
                if (newphi <= cpi->phi0 + ALPHA*cpi->step0*cpi->dphi0) {
                    relaxed_iters = 0;
                    backtrack = 0;
                    //printf("  break 5 : newphi=%.7f\n", newphi);
                }
                else if (newphi >= cpi->phi0) {
                    // Resume last saved line search
                    cvx_cpl_save_or_restore(cp, 1);
                    relaxed_iters = -1;
                }
                else if (newphi <= cpi->phi + ALPHA*cpi->step*cpi->dphi) {
                    backtrack = 0;
                    relaxed_iters = -1;
                    //printf("  break 6 : newphi=%.7f\n", newphi);
                }
            }
        }
    }

    return relaxed_iters;
}

int cvx_cpl_compute_start(cvx_problem_t *cp)
{
    cvx_cpl_internal_t *cpi = &cp->u.cpl;
    F0(cp->F, &cpi->x0);
    cvx_mgrp_initial_value(&cpi->z_g, 0);
    cvx_mgrp_initial_value(&cpi->s_g, 0);
    cvx_scaling_initial_value(&cp->W);
    cvxm_copy(&cpi->x,    &cpi->x0, 0);
    cvxm_copy(&cpi->rx,   &cpi->x0, 0);
    cvxm_copy(&cpi->dx,   &cpi->x0, 0);
    cvxm_copy(&cpi->rx0,  &cpi->x0, 0);
    cvxm_copy(&cpi->dx0,  &cpi->x0, 0);
    cvxm_copy(&cpi->newx, &cpi->x0, 0);
    cvxm_copy(&cpi->newrx,&cpi->x0, 0);
    cvxm_copy(&cpi->y,     cp->b,   0);
    cvxm_copy(&cpi->dy,   &cpi->y,  0);
    cvxm_copy(&cpi->ry,   &cpi->y,  0);
    cvxm_copy(&cpi->newy, &cpi->y,  0);
    return 0;
}

/*
 * Notes: 
 * phase 1 implementation: CPL (convex program with linear objective)
 *   stardard matrix elements c, x, G, h, A, b, s, z
 *   convex programs CP with F0, F1, F2
 */
int cvx_cpl_solve(cvx_problem_t *cp,
                  cvx_solopts_t *opts)
{
    cvx_cpl_internal_t *cpi = &cp->u.cpl;
    // cvx_stats_t *stats = &cp->stats;

    cvx_float_t cx, dy, dz, dr;
    cvx_matrix_t lk, sk, zk, ls, lz, s_l, s_nl, z_l, z_nl, rz_nl, rz_l;

    cvx_size_t range;
    cvx_float_t abstol = CVX_ABSTOL;
    cvx_float_t reltol = CVX_RELTOL;
    cvx_float_t feastol = CVX_FEASTOL;

    int maxiter = opts->max_iter > 0 ? opts->max_iter : CVX_MAXITER;
    int refinement = opts->refinement > 0 ? opts->refinement : 0;

    refinement = opts->refinement == 0 &&
        (cvx_index_count(&cpi->index_full, CVXDIM_SDP) > 0 ||
         cvx_index_count(&cpi->index_full, CVXDIM_SOCP) > 0);

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

    cpi->gap = cvx_sdot(&cpi->s_g, &cpi->z_g);

    range  = cvx_index_count(&cpi->index_full, CVXDIM_SDP);
    range += cvx_index_length(&cpi->index_full, CVXDIM_NONLINEAR);
    range += cvx_index_length(&cpi->index_full, CVXDIM_LINEAR);
    range += cvx_index_length(&cpi->index_full, CVXDIM_SOCP);

    // make non-linear and linear mappings
    cvx_mgrp_elem(&s_nl,  &cpi->s_g,  CVXDIM_NONLINEAR|CVXDIM_NLTARGET, 0);
    cvx_mgrp_elem(&s_l,   &cpi->s_g,  CVXDIM_CONELP, 0);
    cvx_mgrp_elem(&z_nl,  &cpi->z_g,  CVXDIM_NONLINEAR|CVXDIM_NLTARGET, 0);
    cvx_mgrp_elem(&z_l,   &cpi->z_g,  CVXDIM_CONELP, 0);
    cvx_mgrp_elem(&rz_nl, &cpi->rz_g, CVXDIM_NONLINEAR|CVXDIM_NLTARGET, 0);
    cvx_mgrp_elem(&rz_l,  &cpi->rz_g, CVXDIM_CONELP, 0);

    //cvx_mat_printf(stdout, "%.7f", &cpi->c0, "c0");
    //cvx_mat_printf(stdout, "%.7f", &cpi->x0, "x0");
    //cvx_mat_printf(stdout, "%.7f", &cpi->z, "z");
    //cvx_mat_printf(stdout, "%.7f", &cpi->s, "s");
    //cvx_mat_printf(stdout, "%.7f", cp->h, "h");
    //printf("preloop: resz0=%.4f, resy0=%.4f, resz0=%.4f\n", cpi->resx0, cpi->resy0, cpi->resz0);
    // -----------------------------------------------------------------------------

    for (int iter = 0; iter < maxiter; iter++) {

        //printf("**** start of iteration loop ****\n");

        if (refinement != 0) {
            F2(cp->F, &cpi->f, &cpi->Df, &cpi->H, &cpi->x, &z_nl);
        }
        else  {
            F1(cp->F, &cpi->f, &cpi->Df, &cpi->x);
            //cvx_mat_printf(stdout, "%e", &cpi->f, "f");
            //cvx_mat_printf(stdout, "%e", &cpi->Df, "Df");
        }

        cpi->gap = cvx_sdot(&cpi->s_g, &cpi->z_g);
        //printf( "gap    : %e\n", cpi->gap);

        // rx = c + A'*y + Df'*z[:mnl] + G'*z[mnl:]
        cvxm_copy(&cpi->rx, cp->c, 0);
        cvxm_mvmult(1.0, &cpi->rx, 1.0, cp->A, &cpi->y, CVX_TRANS);
        // cvx_mat_printf(stdout, "%13.6e", &cpi->rx, "rx.0");
        cvx_sgemv2(1.0, &cpi->rx, 1.0, &cpi->Df, cp->G, &cpi->z_g, CVX_TRANS);
        // cvx_mat_printf(stdout, "%13.6e", &cpi->rx, "rx.2");
        cpi->resx = cvxm_nrm2(&cpi->rx);
        //printf( "resx   : %e\n", cpi->resx);

        // ry = A*x - b
        cvxm_copy(&cpi->ry, cp->b, 0);
        cvxm_mvmult(-1.0, &cpi->ry, 1.0, cp->A, &cpi->x, 0);
        cpi->resy = cvxm_nrm2(&cpi->ry);

        // rz_nl = s_nl + f
        cvxm_copy(&rz_nl, &s_nl, 0);
        cvxm_axpy(&rz_nl, 1.0, &cpi->f);
        cpi->resznl = cvxm_nrm2(&rz_nl);
        //printf( "resznl : %e\n", cpi->resznl);

        // rz_l = s_l + G*x - h
        cvxm_copy(&rz_l, &s_l, 0);
        cvxm_axpy(&rz_l, -1.0, cp->h);
        //cvx_mat_printf(stdout, "%.7f", cp->h, "h");
        cvxm_mvmult(1.0, &rz_l, 1.0, cp->G, &cpi->x, 0);
        cpi->reszl = cvx_snrm2_elem(&cpi->rz_g, CVXDIM_CONELP);
        //cvx_mat_printf(stdout, "%.7f", &rz_l, "rz_l");
        //printf( "resrzl : %e\n", cpi->reszl);

        // Statistics for stopping criteria
        // pcost = c'*x
        // dcost = c'*x + y'*(A*x-b) + znl'*f(x) + zl'*(G*x-h)
        //       = c'*x + y'*(A*x-b) + znl'*(f(x)+snl) + zl'*(G*x-h+sl)
        //         - z'*s
        //       = c'*x + y'*ry + znl'*rznl + zl'*rzl - gap
        cx = cvxm_dot(cp->c, &cpi->x);
        dy = cvxm_dot(&cpi->y, &cpi->ry);
        dz = cvxm_dot(&z_nl, &rz_nl);
        dr = cvx_sdot_elem(&cpi->z_g, &cpi->rz_g, CVXDIM_CONELP);

        cpi->pcost = cx;
        cpi->dcost = cx + dy + dz +dr - cpi->gap;
        //printf( "pcost : %e\n", cpi->pcost);
        //printf( "dcost : %e\n", cpi->dcost);

        // statistics for stopping
        if (cpi->pcost < 0.0) {
            cpi->relgap = cpi->gap / -cpi->pcost;
        } else if (cpi->dcost > 0.0) {
            cpi->relgap = cpi->gap / cpi->dcost;
        } else {
            cpi->relgap = __NaN();
        }

        //printf("resy  : %e, reszl : %e,  resznl : %e\n", cpi->resy, cpi->reszl, cpi->resznl);
        cpi->pres = cpi->resy*cpi->resy + cpi->reszl*cpi->reszl + cpi->resznl*cpi->resznl;
        cpi->pres = SQRT(cpi->pres);

        cpi->dres = cpi->resx;
        if (iter == 0) {
            cpi->resx0 = __MAX2(1.0, cpi->resx);
            cpi->resznl0 = __MAX2(1.0, cpi->resznl);
            cpi->pres0 = __MAX2(1.0, cpi->pres);
            cpi->dres0 = __MAX2(1.0, cpi->dres);
            cpi->gap0  = cpi->gap;
            cpi->theta1 = 1.0/cpi->gap0;
            cpi->theta2 = 1.0/cpi->resx0;
            cpi->theta3 = 1.0/cpi->resznl0;
        }
        cpi->phi = cpi->theta1*cpi->gap + cpi->theta2*cpi->resx + cpi->theta3*cpi->resznl;
        //printf( "phi    : %e, pres: %e, pres0: %e\n", cpi->phi, cpi->pres, cpi->pres0);
        cpi->pres = cpi->pres / cpi->pres0;
        cpi->dres = cpi->dres / cpi->dres0;

        if (opts->show_progress > 0) {
            if (iter == 0) {
                fprintf(stderr, "%10s %12s %10s %8s %7s\n",
                    "pcost", "dcost", "gap", "pres", "dres");
            }
            fprintf(stderr, "%2d: %11.4e %11.4e %4.0e %7.0e %7.0e\n",
                    iter, cpi->pcost, cpi->dcost, cpi->gap, cpi->pres, cpi->dres);
        }
        // ---------------------------------------------------------------------
        // test for stopping criteria

        if (iter == maxiter) {
            return cvx_cpl_ready(cp, /*stats,*/ iter, CVX_STAT_UNKNOWN);

        }
        else  if (cpi->pres <= feastol &&
                    cpi->dres <= feastol &&
                    (cpi->gap <= abstol ||
                     (!isnan(cpi->relgap) && cpi->relgap < reltol))) {

            return cvx_cpl_ready(cp, /*stats,*/ iter, CVX_STAT_OPTIMAL);
        }

        // -----------------------------------------------------------------------
        // Compute initial scaling W:
        //
        //     W * z = W^{-T} * s = lambda
        //     dg * tau = 1/dg * kappa = lambdag.
        if (iter == 0) {
            cvx_compute_scaling(&cp->W, &cpi->s_g, &cpi->z_g, &cpi->lmbda_g, &cp->work);
            //cvx_scaling_elem_printf(stdout, "%.7f", &cp->W, CVXWS_DNL, 0, "W.dnl");
        }

        // lmdasq = lmda o lmbda
        cvx_ssqr(&cpi->lmbdasq_g, &cpi->lmbda_g);
        //cvx_mat_printf(stdout, "%.7f", &cpi->lmbdasq, "lmbdasq");
        //cvx_mat_printf(stdout, "%.7f", &cpi->lmbda, "lmbda");

        // f3(x, y, z) solves
        //
        //     [ H   A'  GG'*W^{-1} ] [ ux ]   [ bx ]
        //     [ A   0   0          ] [ uy ] = [ by ].
        //     [ GG  0  -W'         ] [ uz ]   [ bz ]
        //
        // On entry, x, y, z contain bx, by, bz.
        // On exit, they contain ux, uy, uz.
        int relaxed_iters = 0;
        int err = cvx_kktfactor(cp->solver, &cp->W, &cpi->x, &z_nl);

        if (err < 0) {
            int singular_kkt = 0;
            // singular matrix ??
            if (iter == 0) {
                // Rank(A) < p or Rank([H(x); A; Df(x); G] < n]??
                return -2;
            }
            if (relaxed_iters > 0 && relaxed_iters < MAX_RELAXED_ITERS) {
                // restore save point
                cvx_cpl_save_or_restore(cp, CPL_RESTORE);
                cpi->resznl = cvxm_nrm2(&rz_nl);
                relaxed_iters = -1;
                F2(cp->F, __cvxnil, &cpi->Df, &cpi->H, &cpi->x, &z_nl);
                err = cvx_kktfactor(cp->solver, &cp->W, &cpi->H, &cpi->Df);
                if (err)
                    singular_kkt = 1;
            } else {
                singular_kkt = 1;
            }

            if (singular_kkt) {
                // terminate in singular KKT_matrix
                return cvx_cpl_ready(cp, /*stats,*/ iter, CVX_STAT_SINGULAR);
            }
        }

        if (iter == 0) {
            if (refinement > 0 || opts->debug) {
                cvxm_copy(&cpi->wx, cp->c, 0);
                cvxm_copy(&cpi->wy, cp->b, 0);
                cvxm_copy(&cpi->wz, cp->z, 0);
                cvxm_copy(&cpi->ws, cp->s, 0);
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
            //cvx_matrix_t t0, t1;
            cpi->mu = cpi->gap/((cvx_float_t)range);
            //printf( "mu    : %e\n", cpi->mu);
            //printf( "sigma : %e\n", cpi->sigma);
            //printf( "gap   : %e  %e\n", cpi->gap, (cvx_float_t)range);
            // Solve
            //
            //     [ 0     ]   [ H  A' GG' ] [ dx        ]
            //     [ 0     ] + [ A  0  0   ] [ dy        ] = -(1 - eta)*r
            //     [ W'*ds ]   [ GG 0  0   ] [ W^{-1}*dz ]
            //
            //     lmbda o (dz + ds) = -lmbda o lmbda + sigma*mu*e.
            //
            //  for nl,l,q: ds = -1.0*lmbdasq
            //cvxm_view_map(&t0, &cpi->ds, 0, 0, length_no_SDP, 1);
            //cvxm_view_map(&t1, &cpi->lmbdasq, 0, 0, length_no_SDP, 1);
            //cvxm_axpby(0.0, &t0, -1.0, &t1);

            cvx_mgrp_copy_lambda(&cpi->ds_g, &cpi->lmbdasq_g);
            cvx_mgrp_scale_sz(&cpi->ds_g, -1.0, 0);
            cvx_mgrp_update_sz(&cpi->ds_g, cpi->sigma*cpi->mu, 0);
            //cvx_mat_printf(stdout, "%.7f", &cpi->ds, "ds");

            cvxm_axpby(0.0, &cpi->dx, -1.0+cpi->eta, &cpi->rx);
            //cvx_mat_printf(stdout, "%.7f", &cpi->rx, "rx");
            //cvx_mat_printf(stdout, "%.7f", &cpi->dx, "dx");

            cvxm_axpby(0.0, &cpi->dy, -1.0+cpi->eta, &cpi->ry);
            cvx_mgrp_axpby_sz(0.0, &cpi->dz_g, -1.0+cpi->eta, &cpi->rz_g, CVXDIM_CONVEXPROG);
            //cvx_mat_printf(stdout, "%.7f", &cpi->rz, "rz");

            //cvx_mat_printf(stdout, "%.7f", &cpi->ds, "ds_g@pre.f4");
            //cvx_mat_printf(stdout, "%.7f", &cpi->dz, "dz_g@pre.f4");
            //cvx_mat_printf(stdout, "%.7f", &cpi->dx, "dx@pre.f4");
            // .. computation here
            cvx_mat_test_nan("pre solve ds", &cpi->ds);
            err = f4(cp, &cpi->dx, &cpi->dy, &cpi->dz_g, &cpi->ds_g, refinement);
            if (err < 0) {
                // terminated ....
                return cvx_cpl_ready(cp, /*stats,*/ iter, CVX_STAT_SINGULAR);
            }
            cvx_mat_test_nan("post solve ds", &cpi->ds);

            //cvx_mat_printf(stdout, "%.7f", &cpi->dx, "dx@f4");
            //cvx_mat_printf(stdout, "%.7f", &cpi->ds, "ds_g@f4");
            //cvx_mat_printf(stdout, "%.7f", &cpi->dz, "dz_g@f4");
            // line search needs ds'*dz and unscaled steps
            cpi->dsdz = cvx_sdot(&cpi->ds_g, &cpi->dz_g);
            //printf("dsdz   : %e\n", cpi->dsdz);
            cvxm_copy(&cpi->dz2, &cpi->dz, 0);
            cvx_scale(&cpi->dz2_g, &cp->W, CVX_INV, &cp->work);
            cvxm_copy(&cpi->ds2, &cpi->ds, 0);
            cvx_scale(&cpi->ds2_g, &cp->W, CVX_TRANS, &cp->work);

            // max step to boundary
            cvx_mat_test_nan("pre scale2 ds", &cpi->ds);
            cvx_mat_test_nan("pre scale2 dz", &cpi->dz);
            cvx_scale2(&cpi->ds_g, &cpi->lmbda_g, 0, &cp->work);
            cvx_scale2(&cpi->dz_g, &cpi->lmbda_g, 0, &cp->work);
            cvx_mat_test_nan("post scale2 ds", &cpi->ds);
            cvx_mat_test_nan("post scale2 dz", &cpi->dz);
            cpi->ts = cvx_max_step(&cpi->ds_g, &cpi->sigs_g, &cp->work);
            cpi->tz = cvx_max_step(&cpi->dz_g, &cpi->sigz_g, &cp->work);
            //printf( "ts    : %e\n", cpi->ts);
            //printf( "tz    : %e\n", cpi->tz);
            cvx_float_t t = __maxvec(3, (cvx_float_t[]){0.0, cpi->ts, cpi->tz});
            if (t == 0.0)
                cpi->step = 1.0;
            else
                cpi->step = __MIN2(1.0, STEP/t);
            //printf( "step    : %e\n", cpi->step);

            // backtrack until newx is in domain of f
            for (int backtrack = 1; backtrack; ) {
                cvxm_copy(&cpi->newx, &cpi->x, 0);
                cvxm_axpy(&cpi->newx, cpi->step, &cpi->dx);
                err = F1(cp->F, &cpi->newf, &cpi->newDf, &cpi->newx);
                if (err == 0)
                    backtrack = 0;
                else
                    cpi->step *= BETA;
                //printf("newx: backtrack=%d, i=%d\n", backtrack, _i);
            }

            // do the line search
            //printf("--- linesearch\n");
            relaxed_iters = cvx_cpl_linesearch(cp, i, relaxed_iters);
            //printf("--- eol step  : %e\n", cpi->step);
        }

        // update x, y
        cvxm_axpy(&cpi->x, cpi->step, &cpi->dx);
        cvxm_axpy(&cpi->y, cpi->step, &cpi->dy);
        //cvx_mat_printf(stdout, "%e", &cpi->x, "updated x");
        //cvx_mat_printf(stdout, "%e", &cpi->dz, "dz");
        //cvx_mat_printf(stdout, "%e", &cpi->ds, "ds");

        // Replace nonlinear, 'l' and 'q' blocks of ds and dz with the updated
        // variables in the current scaling.
        // Replace 's' blocks of ds and dz with the factors Ls, Lz in a
        // factorization Ls*Ls', Lz*Lz' of the updated variables in the
        // current scaling.
        //
        // ds := e + step*ds for 'l' and 'q' blocks.
        // dz := e + step*dz for 'l' and 'q' blocks.
        cvx_mgrp_scale_sz(&cpi->ds_g, cpi->step, CVXDIM_CONVEX|CVXDIM_LINEAR|CVXDIM_SOCP);
        cvx_mgrp_update_sz(&cpi->ds_g, 1.0, CVXDIM_CONVEX|CVXDIM_LINEAR|CVXDIM_SOCP);

        cvx_mgrp_scale_sz(&cpi->dz_g, cpi->step, CVXDIM_CONVEX|CVXDIM_LINEAR|CVXDIM_SOCP);
        cvx_mgrp_update_sz(&cpi->dz_g, 1.0, CVXDIM_CONVEX|CVXDIM_LINEAR|CVXDIM_SOCP);

        //cvx_mat_printf(stdout, "%e", &cpi->dz, "update dz");
        //cvx_mat_printf(stdout, "%e", &cpi->ds, "update ds");
        // ds := H(lambda)^{-1/2} * ds and dz := H(lambda)^{-1/2} * dz.
        //
        // This replaces the 'l' and 'q' components of ds and dz with the
        // updated variables in the current scaling.
        // The 's' components of ds and dz are replaced with
        //
        // diag(lmbda_k)^{1/2} * Qs * diag(lmbda_k)^{1/2}
        // diag(lmbda_k)^{1/2} * Qz * diag(lmbda_k)^{1/2}
        //cvx_mat_printf(stdout, "%e", &cpi->dz, "pre-scale2 dz");

        cvx_scale2(&cpi->ds_g, &cpi->lmbda_g, CVX_INV, &cp->work);
        cvx_scale2(&cpi->dz_g, &cpi->lmbda_g, CVX_INV, &cp->work);

        //cvx_mat_printf(stdout, "%e", &cpi->dz, "scale2(dz)");
        //cvx_mat_printf(stdout, "%e", &cpi->ds, "scale2(ds)");

        // sigs := ( e + step*sigs ) ./ lambda for 's' blocks.
        // sigz := ( e + step*sigz ) ./ lambda for 's' blocks.
        cvxm_scale(&cpi->sigs, cpi->step, CVX_ALL);
        cvxm_scale(&cpi->sigz, cpi->step, CVX_ALL);
        cvxm_add(&cpi->sigs, 1.0, CVX_ALL);
        cvxm_add(&cpi->sigz, 1.0, CVX_ALL);
        //cvx_mat_printf(stdout, "%e", &cpi->sigs, "scaled sigs");
        //cvx_mat_printf(stdout, "%e", &cpi->sigz, "scaled sigz");

        for (int k = 0; k < cvx_mgrp_count(&cpi->lmbda_g, CVXDIM_SDP); k++) {
            cvx_mgrp_elem(&lk, &cpi->lmbda_g, CVXDIM_SDP, k);
            // sigs ./ lmbda
            cvx_mgrp_elem(&sk, &cpi->sigs_g, CVXDIM_SDP, k);
            cvxm_solve_diag(&sk, 1.0, &lk, CVX_RIGHT);
            // sigz ./ lmbda
            cvx_mgrp_elem(&sk, &cpi->sigz_g, CVXDIM_SDP, k);
            cvxm_solve_diag(&sk, 1.0, &lk, CVX_RIGHT);
        }
        //cvx_mat_printf(stdout, "%e", &cpi->sigs, "scaled sigs");
        //cvx_mat_printf(stdout, "%e", &cpi->sigz, "scaled sigz");

        // TODO: missing somethings..!!!!
        // divide ds, dz by lmbda S blocks;;
        for (int k = 0; k < cvx_mgrp_count(&cpi->ds_g, CVXDIM_SDP); k++) {
            cvx_float_t a;
            int m = cvx_mgrp_elem(&sk, &cpi->ds_g, CVXDIM_SDP, k);
            cvx_mgrp_elem(&zk, &cpi->dz_g, CVXDIM_SDP, k);
            cvx_mgrp_elem(&ls, &cpi->sigs_g, CVXDIM_SDP, k);
            cvx_mgrp_elem(&lz, &cpi->sigz_g, CVXDIM_SDP, k);

            for (int i = 0; i < m; i++) {
                cvxm_view_map(&lk, &sk, 0, i, m, 1);
                a = SQRT(cvxm_get(&ls, i, 0));
                cvxm_scale(&lk, a, 0);

                cvxm_view_map(&lk, &zk, 0, i, m, 1);
                a = SQRT(cvxm_get(&lz, i, 0));
                cvxm_scale(&lk, a, 0);
            }
        }

        //cvx_mat_printf(stdout, "%e", &cpi->lmbda, "lambda");
        //cvx_mat_test_nan("prescaling lambda", &cpi->lmbda);
        //cvx_mat_printf(stdout, "%e", &cpi->lmbda, "pre scaling lambda");
        //cvx_mat_printf(stdout, "%e", &cpi->ds, "pre scaling ds");
        //cvx_mat_printf(stdout, "%e", &cpi->dz, "pre scaling dz");
        cvx_update_scaling(&cp->W, &cpi->lmbda_g, &cpi->ds_g, &cpi->dz_g, &cp->work);
        //cvx_scaling_printf(stdout, "%e", &cpi->W, "W");
        //cvx_mat_printf(stdout, "%e", &cpi->lmbda, "scaled lambda");
        //cvx_mat_test_nan("post scaling lambda", &cpi->lmbda);

        // Unscale s, z, tau, kappa (unscaled variables are used only to
        // compute feasibility residuals).
        cvx_mgrp_copy_lambda(&cpi->s_g, &cpi->lmbda_g);
        cvx_scale(&cpi->s_g, &cp->W, CVX_TRANS, &cp->work);

        cvx_mgrp_copy_lambda(&cpi->z_g, &cpi->lmbda_g);
        cvx_scale(&cpi->z_g, &cp->W, CVX_INV, &cp->work);

        cpi->gap = cvxm_dot(&cpi->lmbda, &cpi->lmbda);

        //printf( "*** end of loop gap: %e ****\n", cpi->gap);
    }

    // we never reach here; TODO: fail if we do
    return 0;
}

