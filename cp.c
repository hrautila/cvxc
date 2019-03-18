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

/**
 * @brief Compute memory allocation needed for CONVEXLP problem
 *
 * @param[in] n  Number of variables
 * @param[in] m  Number of equality constrains (rows of A matrix)
 * @param[in] dims Dimensions of inequality constraints, linear, scop and sdp
 *
 * @return Number of bytes of memory needed.
 */
cvx_size_t cvx_cp_bytes(int n, int m, const cvx_dimset_t *dims)
{
    cvx_size_t cdim_mnl  = cvx_dimset_sum_squared(dims, CVXDIM_CONVEXPROG);
    cvx_size_t cdim_diag = cvx_dimset_sum(dims, CVXDIM_CONELP);
    cvx_size_t sdim      = cvx_dimset_sum(dims, CVXDIM_SDP);
    cvx_size_t maxsdp    = cvx_dimset_max(dims, CVXDIM_SDP);
    cvx_size_t mnl       = cvx_dimset_sum(dims, CVXDIM_NONLINEAR) + 1;
    cvx_size_t total = 0;
    cvx_size_t nbytes = 0;

    // size of standard index set (2 times; 1 with CONVEXLP and 1 for CONVEXPROG)
    nbytes += cvx_index_bytes(dims, CVX_INDEX_NORMAL) * 2;
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
    total += 10*n;              // for x, dx, rx, x0, newx, rx0, wx, wx2, newrz
    total += 8*m;               // for y, dy, ry, y0, newy, ry0, wy, wy2
    total += 10*cdim_mnl;       // for s, ds,     s0, news,      ws, ws2, ds0, ds2, ds20, ws3
    total += 13*cdim_mnl;       // for z, dz, rz, z0, newz, rz0, wz, wz2, dz0, dz2, dz20, wz3, newrz0

    //fprintf(stderr, "cpl size: vectors    %ld\n", nbytes+total*sizeof(cvx_float_t));

    total += 2*(cdim_diag+mnl) + 2;   // for lmbda, lmbdasq
    total += 2*sdim;            // for sigs, sigz

    //fprintf(stderr, "cpl size: eigen space %ld\n", nbytes+total*sizeof(cvx_float_t));

    // f, Df, H
    total += mnl*2;       // for f, newf
    total += mnl*n*2;     // for Df, newDf
    total += n*n;         // for H

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
cvx_size_t cvx_cp_make(cvx_problem_t *cp,
                       int n,
                       int m,
                       const cvx_dimset_t *dims,
                       void *memory,
                       cvx_size_t nbytes)
{
    cvx_cpl_internal_t *cpi = &cp->u.cpl;
    cvx_size_t offset = 0;
    cvx_size_t used = 0;
    unsigned char *bytes = (unsigned char *)memory;

    cvx_size_t cdim_mnl  = cvx_dimset_sum_squared(dims, CVXDIM_CONVEXPROG);

    cp->cdim  = cdim_mnl;
    //printf("cdim_mnl: %ld\n", cdim_mnl);

    // __INIT macro assumes variables offset and nbytes;
    // overlay index sets
    __INIT(used, cvx_index_make(&cpi->index_full, dims, CVX_INDEX_NORMAL, bytes,  nbytes));
    __INIT(used, cvx_index_make(&cpi->index_packed, dims, CVX_INDEX_PACKED, &bytes[offset],  nbytes));
    __INIT(used, cvx_index_make(&cpi->index_diag, dims, CVX_INDEX_DIAG, &bytes[offset],  nbytes));
    __INIT(used, cvx_index_make(&cpi->index_sig, dims, CVX_INDEX_SIGS, &bytes[offset],  nbytes));

    // set convex target function indicator to zero
    cvx_dimset_t *__ldims  = (cvx_dimset_t *)dims;
    __ldims->iscpt = 0;
    __INIT(used, cvx_index_make(&cpi->index_cpt, __ldims, CVX_INDEX_NORMAL, &bytes[offset],  nbytes));
    __ldims->iscpt = 1;

    //fprintf(stderr, "cpl make: indexes %ld\n", offset);

    // map result matrix;
    __INIT(used, cvxm_make_epi(&cpi->x, n, 1, &bytes[offset], nbytes));
    __INITC(used, m, &cpi->y, cvxm_make(&cpi->y, m, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->s, cdim_mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->z, cdim_mnl, 1, &bytes[offset], nbytes));

    __INIT(used, cvxm_make_epi(&cpi->c0, n, 1, &bytes[offset], nbytes));

    // dx, dy, ds, dz
    __INIT(used, cvxm_make_epi(&cpi->dx, n, 1, &bytes[offset], nbytes));
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
    __INIT(used, cvxm_make_epi(&cpi->rx, n, 1, &bytes[offset], nbytes));
    __INITC(used, m, &cpi->ry, cvxm_make(&cpi->ry, m, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->rz, cdim_mnl, 1, &bytes[offset], nbytes));

    // x0, y0, z0, s0
    __INIT(used, cvxm_make_epi(&cpi->x0, n, 1, &bytes[offset], nbytes));
    __INITC(used, m, &cpi->y0, cvxm_make(&cpi->y0, m, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->s0, cdim_mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->z0, cdim_mnl, 1, &bytes[offset], nbytes));

    // newx, newy, newz, news, newrx
    __INIT(used, cvxm_make_epi(&cpi->newx, n, 1, &bytes[offset], nbytes));
    __INITC(used, m, &cpi->newy, cvxm_make(&cpi->newy, m, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->news, cdim_mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->newz, cdim_mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->newrx, n, 1, &bytes[offset], nbytes));

    // rx0, ry0, rz0, newrz0
    __INIT(used, cvxm_make_epi(&cpi->rx0, n, 1, &bytes[offset], nbytes));
    __INITC(used, m, &cpi->ry0, cvxm_make(&cpi->ry0, m, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->rz0, cdim_mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->newrz0, cdim_mnl, 1, &bytes[offset], nbytes));

    // wx, wy, ws, wz
    __INIT(used, cvxm_make_epi(&cpi->wx, n, 1, &bytes[offset], nbytes));
    __INITC(used, m, &cpi->wy2, cvxm_make(&cpi->wy, m, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->ws, cdim_mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->wz, cdim_mnl, 1, &bytes[offset], nbytes));

    // wx2, wy2, ws2, wz2
    __INIT(used, cvxm_make_epi(&cpi->wx2, n, 1, &bytes[offset], nbytes));
    __INITC(used, m, &cpi->wy2, cvxm_make(&cpi->wy2, m, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->ws2, cdim_mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->wz2, cdim_mnl, 1, &bytes[offset], nbytes));

    // ws3, wz3
    __INIT(used, cvxm_make(&cpi->ws3, cdim_mnl, 1, &bytes[offset], nbytes));
    __INIT(used, cvxm_make(&cpi->wz3, cdim_mnl, 1, &bytes[offset], nbytes));

    //fprintf(stderr, "cpl make: vectors %ld\n", offset);

    cvx_size_t cdim_diag =
        cvx_dimset_sum(dims, CVXDIM_CONVEXPROG);

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
    int mnl = cvx_dimset_sum(dims, CVXDIM_NONLINEAR) + 1;
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
    cvx_mgrp_init(&cpi->h_g,    cp->h,    &cpi->index_cpt);
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


    return offset;
}


int cvx_cp_isok(const cvx_convex_program_t *F,
                const cvx_matrix_t *G,
                const cvx_matrix_t *h,
                const cvx_matrix_t *A,
                const cvx_matrix_t *b,
                const cvx_dimset_t *dims)
{
    cvx_size_t mG, nG, mA, nA, mh, nh, mb, nb;

    mG = nG = mA = nA = mh = mb = 0;
    nh = nb = 1;
    if (G)
        cvxm_size(&mG, &nG, G);
    if (A)
        cvxm_size(&mA, &nA, A);
    if (h)
        cvxm_size(&mh, &nh, h);
    if (b)
        cvxm_size(&mb, &nb, b);

    // cdim includes non-linear target; substract one
    cvx_size_t cdim = cvx_dimset_sum_squared(dims, CVXDIM_CONVEXPROG);

    if (nh > 1 || mh != cdim - 1) {
        return CVX_ERR_DIMH;
    }

    if (mG != cdim - 1) {
        return CVX_ERR_DIMG;
    }
    if (nA != nG || mA != mb) {
        return CVX_ERR_DIMA;
    }
    if (nb != 1) {
        return CVX_ERR_DIMB;
    }
    return 0;
}

static
int cp_factor(cvx_kktsolver_t *S, cvx_scaling_t *W, cvx_matrix_t *x, cvx_matrix_t *z)
{
    cvx_chainsolver_t *cp_solver = (cvx_chainsolver_t *)S;
    cvx_matrix_t Dfc;
    cvx_problem_t *cp = cp_solver->cp;
    cvx_cpl_internal_t *cpi = &cp->u.cpl;
    cvx_size_t nr, nc;

    cvxm_size(&nr, &nc, &cpi->Df);
    cvxm_view_map(&Dfc, &cpi->Df, 1, 0, nr-1, nc);
    F2(cp->F, __cvxnil, &cpi->Df, &cpi->H, x, z);
    //cvx_mat_printf(stdout, "%.7f", z, "z");
    //cvx_mat_printf(stdout, "%.7f", &cpi->H, "H");
    //cvx_mat_printf(stdout, "%.7f", x, "x.m()");
    //cvx_mat_printf(stdout, "%.7f", &cpi->Df, "Df");
    //
    return cvx_kktfactor(cp_solver->next, W, &cpi->H, &Dfc);
}

static
int cp_solve(cvx_kktsolver_t *S, cvx_matrix_t *x, cvx_matrix_t *y, cvx_matgrp_t *z_g)
{
    cvx_chainsolver_t *cp_solver = (cvx_chainsolver_t *)S;
    cvx_problem_t *cp = cp_solver->cp;
    cvx_cpl_internal_t *cpi = &cp->u.cpl;
    cvx_matrix_t gradf, dnlt, z_nlt, z_cpl;
    cvx_index_t index_cpl;
    cvx_matgrp_t z_cpl_g;
    cvx_float_t znlt, epi_x, dnlt0;
    cvx_size_t nr, nc;

    //cvx_mat_printf(stdout, "%.7f", z_g->mat, "presolve z");
    cvxm_size(&nr, &nc, &cpi->Df);
    cvxm_view_map(&gradf, &cpi->Df, 0, 0, 1, nc);
    cvx_scaling_elem(&dnlt, &cp->W, CVXWS_DNLT, 0);
    cvx_mgrp_elem(&z_nlt, z_g, CVXDIM_NLTARGET, 0);
    cvx_mgrp_elem(&z_cpl, z_g, CVXDIM_CONVEXLP, 0);

    cvx_subindex(&index_cpl, &cpi->index_full, CVXDIM_CONVEXLP);
    cvx_mgrp_init(&z_cpl_g, &z_cpl, &cpi->index_cpt);

    dnlt0 = cvxm_get(&dnlt, 0, 0);
    epi_x = cvxm_get_epi(x);
    znlt  = cvxm_get(&z_nlt, 0, 0);
    cvxm_axpy(x, epi_x, &gradf);

    //printf("dnlt0 : %.7f\n", dnlt0);
    //printf("znlt  : %.7f\n", znlt);
    //printf("epi(x): %.7f\n", epi_x);
    //cvx_mat_printf(stdout, "%.7f", x, "presolve x");

    int err = cvx_kktsolve(cp_solver->next, x, y, &z_cpl_g);

    //cvx_mat_printf(stdout, "%.7f", x, "postsolve x");

    cvxm_set(&z_nlt, 0, 0, -epi_x*dnlt0);
    //cvx_mat_printf(stdout, "%.7f", &gradf, "gradf0");
    //printf("grd*x    : %.7f\n", cvxm_dot(&gradf, x));
    //printf("d*d*x.t  : %.7f\n", dnlt0*dnlt0*epi_x);
    //printf("z0       : %.7f\n", znlt);
    epi_x = cvxm_dot(&gradf, x) + dnlt0*dnlt0*epi_x - znlt;
    cvxm_set_epi(x, epi_x);
    return err;
}

static cvx_kktfuncs_t cp_kktfunctions = {
    .factor = cp_factor,
    .solve  = cp_solve
};

cvx_chainsolver_t *cvx_cp_solver_init(cvx_chainsolver_t *cs, cvx_problem_t *cp, cvx_kktsolver_t *next)
{
    cs->fnc = &cp_kktfunctions;
    cs->next = next;
    cs->cp = cp;
    return cs;
}


cvx_size_t cvx_cp_setup(cvx_problem_t *cp,
                        cvx_convex_program_t *F,
                        cvx_matrix_t *G,
                        cvx_matrix_t *h,
                        cvx_matrix_t *A,
                        cvx_matrix_t *b,
                        cvx_dimset_t *dims,
                        cvx_kktsolver_t *kktsolver)
{
    cvx_cpl_internal_t *cpi = &cp->u.cpl;
    cvx_size_t mG, nG, mb, nb;
    int err;

    if (! cp)
        return 0;

    if ((err = cvx_cp_isok(F, G, h, A, b, dims)) != 0) {
        cp->error = err;
        return 0;
    }

    mG = nG = mb = nb = 0;
    cvxm_size(&mG, &nG, G);
    if (b)
        cvxm_size(&mb, &nb, b);

    cvx_size_t nbytes = cvx_cp_bytes(nG, mb, dims);
    void *memory = calloc(nbytes, 1);
    if (!memory) {
        cp->error = CVX_ERR_MEMORY;
        return 0;
    }

    cp->c = __cvxnil;
    cp->G = G;
    cp->h = h;
    cp->A = A;
    cp->b = b;
    cp->dims = dims;
    cp->F = F;

    cp->primal_x = (cvx_matrix_t *)0;
    cp->primal_s = (cvx_matrix_t *)0;
    cp->dual_y = (cvx_matrix_t *)0;
    cp->dual_z = (cvx_matrix_t *)0;

    if (cvx_cp_make(cp, nG, mb, dims, memory, nbytes) == 0) {
        cp->error = CVX_ERR_MEMORY;
        return 0;
    }
    cp->mlen = nbytes;
    cp->memory = memory;

    cp->f  = &cp->u.cpl.f;
    cp->Df = &cp->u.cpl.Df;
    cp->H  = &cp->u.cpl.H;

    cp->c = &cpi->c0;

    // provide full index set G/h matrices
    cp->index_g = &cp->u.cpl.index_cpt;

    // init KKT solver
    if (kktsolver) {
        cp->solver = kktsolver;
    } else {
        cvx_ldlsolver_init(&cp->__S, cp, nG, mb, dims);
        cp->solver = &cp->__S;
    }
    cvx_cp_solver_init(&cpi->cp_solver, cp, cp->solver);
    cp->solver = (cvx_kktsolver_t *)&cpi->cp_solver;

    return nbytes;
}

int cvx_cp_compute_start(cvx_problem_t *cp)
{
    cvx_cpl_internal_t *cpi = &cp->u.cpl;
    F0(cp->F, &cpi->x0);
    cvx_mgrp_initial_value(&cpi->z_g, 0);
    cvx_mgrp_initial_value(&cpi->s_g, 0);
    cvx_scaling_initial_value(&cp->W);
    cvxm_set_all(&cpi->c0, 0.0);
    cpi->c0.t = 1.0;

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

int cvx_cp_solve(cvx_problem_t *cp,
                 cvx_solopts_t *opts)
{
    return cvx_cpl_solve(cp, opts);
}
