/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */
#include "cvxc.h"

cvxc_size_t cvxc_gp_program_bytes(cvxc_size_t p, cvxc_size_t m, cvxc_size_t maxf)
{
    return (p + 1)*sizeof(cvxc_size_t) + (m + maxf)*sizeof(cvxc_float_t);
}

cvxc_size_t cvxc_gpi_make(cvxc_gpindex_t *gpi, cvxc_size_t p, void *memory, cvxc_size_t nbytes)
{
    if (nbytes < (p+1)*sizeof(cvxc_size_t))
        return 0;

    gpi->p = p;
    gpi->index = (cvxc_size_t *)memory;
    gpi->index[0] = 0;
    for (cvxc_size_t k = 0; k < p; k++) {
        gpi->index[k+1] = 0;
    }
    return (p+1)*sizeof(cvxc_size_t);
}

/**
 * @brief Get the n'th element in GP matrix F/g.
 *
 * @return Row index of element
 */
cvxc_size_t cvxc_gpi_elem(const cvxc_gpindex_t *gpi, cvxc_matrix_t *e, const cvxc_matrix_t *Fg, cvxc_size_t n)
{
    cvxc_size_t mS, nS;
    if (n < gpi->p) {
        cvxm_size(&mS, &nS, Fg);
        cvxm_view_map(e, Fg, gpi->index[n], 0, gpi->index[n+1]-gpi->index[n], nS);
        return gpi->index[n];
    }
    return 0;
}

/**
 * @brief Get length of n'th GP index element.
 */
cvxc_size_t cvxc_gpi_length(const cvxc_gpindex_t *gpi, cvxc_size_t n)
{
    return n < gpi->p ? gpi->index[n+1] - gpi->index[n] : 0;
}

int cvxc_gpi_setup(cvxc_gpindex_t *gpi, cvxc_size_t *K, cvxc_size_t p)
{
    gpi->index[0] = 0;
    for (cvxc_size_t k = 0; k < p; k++) {
        gpi->index[k+1] = gpi->index[k] + K[k];
    }
    return 0;
}

struct logsum {
    cvxc_float_t ymax;
    cvxc_float_t ysum;
};

cvxc_float_t eval(cvxc_float_t e, void *p)
{
    struct logsum *lgs = (struct logsum *)p;
    e = exp(e - lgs->ymax);
    lgs->ysum += fabs(e);
    return e;
}

int cvxc_gp_f(cvxc_matrix_t *f, cvxc_matrix_t *Df, cvxc_matrix_t *H,
             const cvxc_matrix_t *x, const cvxc_matrix_t *z, void *user)
{
    cvxc_gp_program_t *gp = (cvxc_gp_program_t *)user;
    cvxc_matrix_t iy, iF, iDf, Fsc, *g, *y, *F, *Fs;
    cvxc_gpindex_t *gpi;
    cvxc_float_t ymax;
    cvxc_size_t m, n;
    struct logsum lg;

    if (!x && !z) {
        cvxm_set_all(f, 0.0);
        return 0;
    }

    y  = &gp->gp_params.y;
    Fs = &gp->gp_params.Fs;
    g  = gp->gp_params.g;
    F  = gp->gp_params.F;
    gpi = &gp->gp_params.gpi;

    cvxm_set_all(f, 0.0);
    cvxm_set_all(Df, 0.0);

    // y = g + F*x
    cvxm_copy(y, g, 0);
    cvxm_mvmult(1.0, y, 1.0, F, x, 0);

    cvxm_size(&m, &n, F);
    if (z)
        cvxm_set_all(H, 0.0);

    for (cvxc_size_t k = 0; k < gpi->p; k++) {
        cvxc_gpi_elem(gpi, &iy, y, k);
        cvxc_gpi_elem(gpi, &iF, F, k);
        cvxm_view_map(&iDf, Df, k, 0, 1, n);
        m = cvxc_gpi_length(gpi, k);

        // yi := exp(yi) = exp(Fi*x+gi)
        ymax = cvxm_max(&iy);
        lg = (struct logsum){ .ymax = ymax, .ysum = 0.0 };
        cvxm_apply2(&iy, eval, &lg);

        // f_i = log sum y_i = log sum exp(F_i*x + g_i)
        if (f)
            cvxm_set(f, k, 0, ymax + log(lg.ysum));

        // yi := yi / sum(yi) = exp(Fi*x+gi) / sum(exp(Fi*x+gi))
        cvxm_scale(&iy, 1.0/lg.ysum, 0);

        // gradfi := Fi' * yi
        //         = Fi' * exp(Fi*x+gi) / sum(exp(Fi*x+gi))
        cvxm_mvmult(0.0, &iDf, 1.0, &iF, &iy, CVXC_TRANS);

        if (z) {
            // Hi = Fi' * (diag(yi) - yi*yi') * Fi
            //    = Fisc' * Fisc
            // where
            // Fisc = diag(yi)^1/2 * (I - 1*yi') * Fi
            //      = diag(yi)^1/2 * (Fi - 1*gradfi')
            cvxm_view_map(&Fsc, Fs, 0, 0, m, n);
            cvxm_copy(&Fsc, &iF, 0);

            for (int i = 0; i < m; i++) {
                cvxc_matrix_t Fs_row;
                cvxm_view_map(&Fs_row, &Fsc, i, 0, 1, n);
                cvxm_axpy(&Fs_row, -1.0, &iDf);
                cvxm_scale(&Fs_row, SQRT(cvxm_get(&iy, i, 0)), 0);
            }
            cvxc_float_t zi = cvxm_get(z, k, 0);
            cvxm_update_sym(H, zi, &Fsc, CVXC_TRANS);
        }
    }
    return 0;
}

/**
 * @brief Setup geometric program.
 *
 *   minimize    \f$ \log \sum \exp(F_0*x + g_0) \f$
 *   subject to  \f$ \log \sum \exp(F_i*x + g_i) <= 0, i = 1, ..., m \f$
 *               \f$ G \times x <= h \f$
 *               \f$ A \times x = b \f$
 *
 * @param[in] K
 *   Array of positive integers terminated with zero value [K0, K1, ..., Km, 0]
 * @param[in] F
 *   Matrix with submatrices [F0, F1, ..., Fm] where each submatrix is \$f K_i x n \f$.
 * @param[in] g
 *   Column vector with subvectors [g0, g1, ..., gm] where each subvector is \$f K_i x 1 \f$.
 * @param[in] G
 * @param[in] h
 * @param[in] A
 * @param[in] b
 */
int cvxc_gp_setup(cvxc_problem_t *cp,
                  cvxc_gpindex_t *K,
                  cvxc_matrix_t *F,
                  cvxc_matrix_t *g,
                  cvxc_matrix_t *G,
                  cvxc_matrix_t *h,
                  cvxc_matrix_t *A,
                  cvxc_matrix_t *b,
                  cvxc_kktsolver_t *kktsolver)
{
    cvxc_cpl_internal_t *cpi;
    cvxc_size_t mG, nG, mF, nF, mg, ng, p, mK, mA, nA, maxK, offset, gpbytes;
    cvxc_dimset_t ldims;

    if (!cp)
        return 0;

    cpi = &cp->u.cpl;
    mA = nA = maxK = 0;
    mG = nG = 0;

    cvxm_size(&mF, &nF, F);
    cvxm_size(&mg, &ng, g);
    for (mK = 0, p = 0; p < mF && p < K->p; p++) {
        mK += K->index[p];
        if (K->index[p] > maxK)
            maxK = K->index[p];
    }
    if (mK != mF || mg != mF) {
        return 0;
    }
    if (A)
        cvxm_size(&mA, &nA, A);
    if (G)
        cvxm_size(&mG, &nG, G);

    // p is count of F_i, g_i elements
    gpbytes = cvxc_gp_program_bytes(p, mg, maxK*nF);
    ldims = (cvxc_dimset_t){0};
    ldims.mnl = p - 1;
    ldims.iscpt = 1;
    ldims.ldim = mG;

    offset = cvxc_cpl_allocate(cp, 1, nG, mA, gpbytes, &ldims);
    if (offset == 0)
        return 0;

    cvxc_gp_program_t *gp = &cpi->gp;
    gp->gp_params.F = F;
    gp->gp_params.g = g;
    gp->gp.F = cvxc_gp_f;
    gp->gp.user = gp;

    cvxc_cp_setvars(cp, &gp->gp, nG, mA, G, h, A, b, &ldims, kktsolver);

    cvxc_size_t used, nbytes;

    nbytes = cp->mlen;

    // Workspace for computng gradients and Hessian
    used = cvxm_make(&gp->gp_params.y, mg, 1, &cp->memory[offset], nbytes-offset);
    if (used == 0)
        return 0;
    offset += used;

    used = cvxm_make(&gp->gp_params.Fs, maxK, nG, &cp->memory[offset], nbytes-offset);
    if (used == 0)
        return 0;
    offset += used;

    used = cvxc_gpi_make(&gp->gp_params.gpi, p, &cp->memory[offset], nbytes-offset);
    if (used == 0)
        return 0;
    offset += used;

    cvxc_gpi_setup(&gp->gp_params.gpi, K->index, p);
    return offset;
}


int cvxc_gp_compute_start(cvxc_problem_t *cp)
{
    return cvxc_cp_compute_start(cp);
}

int cvxc_gp_solve(cvxc_problem_t *cp, cvxc_solopts_t *opts)
{
    return cvxc_cp_solve(cp, opts);
}
