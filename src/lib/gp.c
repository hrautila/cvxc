/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */

#include "internal.h"

cvxc_size_t cvxc_gp_program_bytes(cvxc_size_t p, cvxc_size_t m, cvxc_size_t maxf)
{
    return (p + 1)*sizeof(cvxc_size_t) + (m + maxf)*sizeof(cvxc_float_t);
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
    cvxc_gp_params_t *gp = (cvxc_gp_params_t *)user;
    cvxc_matrix_t iy, iF, iDf, Fsc, *g, *y, *F, *Fs;
    cvxc_gpindex_t *gpi;
    cvxc_float_t ymax;
    cvxc_size_t m, n;
    struct logsum lg;

    if (!x && !z) {
        cvxm_set_all(f, 0.0);
        return 0;
    }

    y  = &gp->y;
    Fs = &gp->Fs;
    gpi = gp->gpi;
    g  = gp->g;
    F  = gp->F;

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

/*
 * GP program memory layout:
 *
 *  +----------------------+
 *  | cvxc_cpl_internal_t  |
 *  +----------------------+
 *  |    CPL program       |
 *  |   internal data      |
 *  |       area           |
 *  +----------------------+
 *  | cvxc_convex_program_t|
 *  +----------------------+
 *  |  cvxc_gp_params_t    |
 *  +----------------------+
 *  |    GP program        |
 *  |   internal data      |
 *  |      area            |
 *  +----------------------+
 */

static
int cvxc_gp_create(
    cvxc_problem_t *cp, cvxc_size_t nvars,  cvxc_gpindex_t *K, cvxc_matrix_t *F,
    cvxc_matrix_t *g, cvxc_umatrix_t *Gf, cvxc_matrix_t *h,
    cvxc_umatrix_t *Af, cvxc_matrix_t *b, cvxc_kktsolver_t *kktsolver)
{
    cvxc_size_t mh, nh, mF, nF, mg, ng, p, mb, nb, maxK, offset, gpbytes;
    cvxc_dimset_t ldims;

    if (!cp)
        return 0;

    mh = nh = mb = nb = maxK = 0;

    cvxm_size(&mF, &nF, F);
    cvxm_size(&mg, &ng, g);
    for (p = 0; p < mF && p < K->p; p++) {
        if (cvxc_gpi_length(K, p) > maxK)
            maxK = cvxc_gpi_length(K, p);
    }
    if (K->index[K->p] != mF || mg != mF) {
        return 0;
    }
    if (b)
        cvxm_size(&mb, &nb, b);
    if (h)
        cvxm_size(&mh, &nh, h);

    gpbytes = cvxc_gp_program_bytes(K->p, mg, maxK*nF);
    gpbytes += sizeof(cvxc_convex_program_t) + sizeof(cvxc_gp_params_t);

    ldims = (cvxc_dimset_t){0};
    ldims.mnl = K->p - 1;
    ldims.iscpt = 1;
    ldims.ldim = mh;

    /* Layout CPL internal memory */
    offset = cvxc_cpl_allocate(cp, 1, nvars, mb, gpbytes, &ldims);
    if (offset == 0)
        return 0;

    cp->F = 0;
    cp->Gf = Gf;
    cp->h = h;
    cp->Af = Af;
    cp->b = b;
    cvxc_cp_finalize_setup(cp, nvars, mb, &ldims, kktsolver);

    /* Append cvxc_convex_program structure. */
    unsigned char *memory = cp->u.space;
    cp->F = (cvxc_convex_program_t *)&memory[offset];
    offset += sizeof(cvxc_convex_program_t);

    /* Append cvxc_gp_params structure. */
    cvxc_gp_params_t *gp = (cvxc_gp_params_t *)&memory[offset];
    offset += sizeof(cvxc_gp_params_t);

    cvxc_convex_program_init(cp->F, cvxc_gp_f, gp);
    gp->F = F;
    gp->g = g;
    gp->gpi = K;

    cvxc_size_t used;
    /* Append GP internal variables. */
    used = cvxm_make(&gp->y, mg, 1, &memory[offset], cp->nbytes-offset);
    if (used == 0)
        return 0;
    offset += used;

    used = cvxm_make(&gp->Fs, maxK, nvars, &memory[offset], cp->nbytes-offset);
    if (used == 0)
        return 0;
    offset += used;

    return offset;
}

int cvxc_gp_setup_user(
    cvxc_problem_t *cp, cvxc_gpindex_t *K, cvxc_matrix_t *F,
    cvxc_matrix_t *g, cvxc_umatrix_t *Gf, cvxc_matrix_t *h,
    cvxc_umatrix_t *Af, cvxc_matrix_t *b, cvxc_kktsolver_t *kktsolver)
{
    if (!cp || !F || !g)
        return 0;

    cvxc_size_t mF = 0, nF = 0;
    cvxm_size(&mF, &nF, F);
    return cvxc_gp_create(cp, nF, K, F, g, Gf, h, Af, b, kktsolver);
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
 *   Submatrix index set.
 * @param[in] F
 *   Matrix with submatrices [F0, F1, ..., Fm] where each submatrix is \$f K_i x n \f$.
 * @param[in] g
 *   Column vector with subvectors [g0, g1, ..., gm] where each subvector is \$f K_i x 1 \f$.
 * @param[in] G
 * @param[in] h
 * @param[in] A
 * @param[in] b
 */
int cvxc_gp_setup(
    cvxc_problem_t *cp, cvxc_gpindex_t *K, cvxc_matrix_t *F,
    cvxc_matrix_t *g, cvxc_matrix_t *G, cvxc_matrix_t *h,
    cvxc_matrix_t *A, cvxc_matrix_t *b, cvxc_kktsolver_t *kktsolver)
{
    if (!cp || !F || !g)
        return 0;

    cvxc_size_t mF = 0, nF = 0;
    cvxm_size(&mF, &nF, F);
    cvxc_umat_make(&cp->Au, A);
    cvxc_umat_make(&cp->Gu, G);
    cp->A = A;
    cp->G = G;

    int stat = cvxc_gp_create(cp, nF, K, F, g, &cp->Gu, h, &cp->Au, b, kktsolver);
    if (stat <= 0) {
        cvxc_umat_clear(&cp->Au);
        cvxc_umat_clear(&cp->Gu);
    }
    return stat;
}

/**
 * @brief Release resources reserved to program.
 */
void cvxc_gp_release(cvxc_problem_t *cp)
{
    if (!cp || !cp->u.cpl)
        return;

    cvxc_kktrelease(cp->solver);

    free(cp->u.space);
    cp->u.space = 0;
    cp->nbytes = 0;
    cp->work = 0;
}

int cvxc_gp_compute_start(cvxc_problem_t *cp)
{
    return cvxc_cp_compute_start(cp);
}

int cvxc_gp_solve(cvxc_solution_t *sol, cvxc_problem_t *cp, cvxc_solopts_t *opts)
{
    return cvxc_cp_solve(sol, cp, opts);
}
