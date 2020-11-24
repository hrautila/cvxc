/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */
#include "cvxc.h"

cvxc_size_t cvxc_gp_program_bytes(cvxc_size_t p, cvxc_size_t m)
{
    return (p + 1)*sizeof(cvxc_size_t) + m*sizeof(cvxc_float_t);
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


int cvxc_gp_f(cvxc_matrix_t *f, cvxc_matrix_t *Df, cvxc_matrix_t *H,
             const cvxc_matrix_t *x, const cvxc_matrix_t *z, void *user)
{
    cvxc_gp_program_t *gp = (cvxc_gp_program_t *)user;
    cvxc_matrix_t iy, iF, iDf;
    cvxc_float_t ymax, ysum;
    cvxc_size_t n;

    if (!x && !z) {
        cvxm_set_all(f, 0.0);
        return 0;
    }

    cvxm_set_all(f, 0.0);
    cvxm_set_all(Df, 0.0);

    // y = g + F*x
    cvxm_copy(&gp->gp_params.y, gp->gp_params.g, 0);
    cvxm_mvmult(1.0, &gp->gp_params.y, 1.0, gp->gp_params.F, x, 0);

    for (cvxc_size_t k = 0; k < gp->gp_params.gpi.p; k++) {
        n = cvxc_gpi_elem(&gp->gp_params.gpi, &iy, &gp->gp_params.y, k);
        cvxc_gpi_elem(&gp->gp_params.gpi, &iF, gp->gp_params.F, k);
        cvxc_gpi_elem(&gp->gp_params.gpi, &iDf, Df, k);

        // yi := exp(yi) = exp(Fi*x+gi)
        ymax = cvxm_amax(&iy);
        for (cvxc_size_t i = 0; i < n; i++) {
            cvxm_set(&iy, i, 0, exp(cvxm_get(&iy, i, 0) - ymax));
        }

        // f_i = log sum y_i = log sum exp(F_i*x + g_i)
        ysum = cvxm_asum(&iy);
        cvxm_set(f, k, 0, ymax + log(ysum));

        // yi := yi / sum(yi) = exp(Fi*x+gi) / sum(exp(Fi*x+gi))
        cvxm_scale(&iy, 1.0/ysum, 0);

        // gradfi := Fi' * yi
        //         = Fi' * exp(Fi*x+gi) / sum(exp(Fi*x+gi))
        cvxm_mvmult(0.0, &iDf, 1.0, &iF, &iy, CVXC_TRANS);

        if (z) {
            // Hi = Fi' * (diag(yi) - yi*yi') * Fi
            //    = Fisc' * Fisc
            // where
            // Fisc = diag(yi)^1/2 * (I - 1*yi') * Fi
            //      = diag(yi)^1/2 * (Fi - 1*gradfi')
        }
    }
    return 0;
}



int cvxc_gp_setup(cvxc_problem_t *cp,
                 cvxc_size_t *K,
                 cvxc_matrix_t *F,
                 cvxc_matrix_t *g,
                 cvxc_matrix_t *G,
                 cvxc_matrix_t *h,
                 cvxc_matrix_t *A,
                 cvxc_matrix_t *b,
                 cvxc_dimset_t *dims,
                 cvxc_kktsolver_t *kktsolver)
{
    cvxc_cpl_internal_t *cpi;
    cvxc_size_t mG, nG, mF, nF, mg, ng, p, mK, mA, nA, offset, gpbytes;
    cvxc_dimset_t ldims = *dims;

    if (!cp)
        return 0;

    cpi = &cp->u.cpl;
    mA = nA = 0;
    mG = nG = 0;

    cvxm_size(&mF, &nF, F);
    cvxm_size(&mg, &ng, g);
    for (mK = 0, p = 0; p < mF && K[p] > 0; p++) {
        mK += K[p];
    }
    //
    if (mK != mF || mg != mF) {
        return 0;
    }
    gpbytes = cvxc_gp_program_bytes(p, mg);
    if (A)
        cvxm_size(&mA, &nA, A);
    if (G)
        cvxm_size(&mG, &nG, G);

    ldims = *dims;
    ldims.mnl = p - 1;
    ldims.iscpt = 1;

    offset = cvxc_cpl_allocate(cp, 1, nG, mA, gpbytes, &ldims);
    if (offset == 0)
        return 0;

    cvxc_cpl_setvars(cp, 0, nG, mA, __cvxnil, G, h, A, b, dims, kktsolver);

    cvxc_size_t used, nbytes;

    nbytes = cp->mlen;
    // TODO: not right!!
    cvxc_gp_program_t *gp = &cpi->gp;
    gp->gp_params.F = F;
    gp->gp_params.g = g;
    gp->gp.user = &gp->gp_params;
    cp->F = &gp->gp;

    used = cvxm_make(&gp->gp_params.y, mg, 1, &cp->memory[nbytes], nbytes-offset);
    if (used == 0)
        return 0;

    offset += used;
    used = cvxc_gpi_make(&gp->gp_params.gpi, p, &cp->memory[nbytes], nbytes-offset);
    if (used == 0)
        return 0;
    cvxc_gpi_setup(&gp->gp_params.gpi, K, p);
    return offset;
}


int cvxc_gp_set_start(cvxc_problem_t *cp,
                     cvxc_matrix_t *x0,
                     cvxc_matrix_t *s0,
                     cvxc_matrix_t *y0,
                     cvxc_matrix_t *z0)
{
    return cvxc_cp_compute_start(cp);
}

int cvxc_gp_solve(cvxc_problem_t *cp, cvxc_solopts_t *opts)
{
    return cvxc_cp_solve(cp, opts);
}
