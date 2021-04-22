/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */

#include "internal.h"

/**
 */
int cvxc_cp_isok(
    cvxc_size_t nvars,
    const cvxc_convex_program_t *F,
    const cvxc_matrix_t *G,
    const cvxc_matrix_t *h,
    const cvxc_matrix_t *A,
    const cvxc_matrix_t *b,
    const cvxc_dimset_t *dims)
{
    cvxc_size_t mG, nG, mA, nA, mh, nh, mb, nb;

    mG = mA = mh = mb = 0;
    nG = nA = nvars;
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
    cvxc_size_t cdim = cvxc_dimset_sum_squared(dims, CVXDIM_CONVEXPROG);

    if (nh > 1 || mh != cdim - 1) {
        return CVXC_ERR_DIMH;
    }

    if (mG != cdim - 1) {
        return CVXC_ERR_DIMG;
    }
    if (nA != nG || (A && mA != mb)) {
        return CVXC_ERR_DIMA;
    }
    if (nb != 1) {
        return CVXC_ERR_DIMB;
    }
    return 0;
}

static
int cp_factor(cvxc_cp_kktsolver_t *S, cvxc_scaling_t *W, cvxc_matrix_t *H, cvxc_matrix_t *Df)
{
    cvxc_matrix_t Dfc;
    cvxc_problem_t *cp = (cvxc_problem_t *)S->private;
    cvxc_cpl_internal_t *cpi = cp->u.cpl;
    cvxc_size_t nr, nc;

    cvxm_size(&nr, &nc, &cpi->Df);
    cvxm_view_map(&Dfc, &cpi->Df, 1, 0, nr-1, nc);
    // F2(cp->F, &cpi->f, &cpi->Df, &cpi->H, x, z);

    return cvxc_kktfactor(S->next, W, &cpi->H, &Dfc);
}

static
int cp_solve(cvxc_cp_kktsolver_t *S, cvxc_epigraph_t *x, cvxc_matrix_t *y, cvxc_matgrp_t *z_g)
{
    cvxc_problem_t *cp = (cvxc_problem_t *)S->private;
    cvxc_cpl_internal_t *cpi = cp->u.cpl;
    cvxc_matrix_t gradf, dnlt, z_nlt, z_cpl;
    cvxc_index_t index_cpl;
    cvxc_matgrp_t z_cpl_g;
    cvxc_float_t znlt, epi_x, dnlt0;
    cvxc_size_t nr, nc;

    cvxm_size(&nr, &nc, &cpi->Df);
    cvxm_view_map(&gradf, &cpi->Df, 0, 0, 1, nc);
    cvxc_scaling_elem(&dnlt, &cpi->W, CVXWS_DNLT, 0);
    cvxc_mgrp_elem(&z_nlt, z_g, CVXDIM_NLTARGET, 0);
    cvxc_mgrp_elem(&z_cpl, z_g, CVXDIM_CONVEXLP, 0);

    cvxc_subindex(&index_cpl, &cpi->index_full, CVXDIM_CONVEXLP);
    cvxc_mgrp_init(&z_cpl_g, &z_cpl, &cpi->index_cpt);

    dnlt0 = cvxm_get(&dnlt, 0, 0);
    epi_x = cvxm_get_epival(x);
    znlt  = cvxm_get(&z_nlt, 0, 0);
    cvxm_axpy(x->m, epi_x, &gradf);

    int err = cvxc_kktsolve(S->next, x->m, y, &z_cpl_g);

    cvxm_set(&z_nlt, 0, 0, -epi_x*dnlt0);
    epi_x = cvxm_dot(&gradf, x->m) + dnlt0*dnlt0*epi_x - znlt;
    cvxm_set_epival(x, epi_x);
    return err;
}

static
void cp_release(cvxc_cp_kktsolver_t *kkt)
{
    if (!kkt)
        return;

    if (kkt->next && kkt->next->vtable->release)
        (*kkt->next->vtable->release)(kkt->next);
    kkt->private = 0;
}

static cvxc_cp_kktfuncs_t cp_kktfunctions = {
    .factor = cp_factor,
    .solve  = cp_solve,
    .release= cp_release
};

static
void cvxc_cp_solver_init(cvxc_cp_kktsolver_t *cs, cvxc_problem_t *cp, cvxc_kktsolver_t *next)
{
    cs->vtable = &cp_kktfunctions;
    cs->next = next;
    cs->private = cp;
}

/**
 * @brief Initialize internal variables for convex program.
 */
int cvxc_cp_finalize_setup(
    cvxc_problem_t *cp,
    cvxc_size_t n,
    cvxc_size_t m,
    const cvxc_dimset_t *dims,
    cvxc_kktsolver_t *kktsolver)
{
    // fix index set G/h matrices
    cvxc_cpl_internal_t *cpi = cp->u.cpl;
    cp->index_g = &cpi->index_cpt;
    cp->c = &cpi->c0;
    cvxc_mgrp_init(&cpi->h_g, cp->h, &cpi->index_cpt);

    // init KKT solver
    if (kktsolver) {
        cp->solver = kktsolver;
    } else {
        cvxc_ldlsolver_init(&cp->__S, cp, n, m, dims);
        cp->solver = &cp->__S;
    }
    cvxc_cp_solver_init(&cpi->cp_solver, cp, cp->solver);
    cp->solver = (cvxc_kktsolver_t *)&cpi->cp_solver;

    return 0;
}

/*
 * @brief Allocate space and initialize variables for convex program.
 */
static
cvxc_size_t cvxc_cp_create(
    cvxc_problem_t *cp,
    cvxc_convex_program_t *F,
    cvxc_size_t nvars,
    cvxc_umatrix_t *Gf,
    cvxc_matrix_t *h,
    cvxc_umatrix_t *Af,
    cvxc_matrix_t *b,
    const cvxc_dimset_t *dims,
    cvxc_kktsolver_t *kktsolver)
{
    cvxc_size_t mb, nb;
    cvxc_size_t used;

    mb = nb = 0;
    if (b)
        cvxm_size(&mb, &nb, b);

    if ((used = cvxc_cpl_allocate(cp, 1, nvars, mb, 0, dims)) == 0) {
        return 0;
    }
    // cvxc_cp_setvars(cp, F, nG, mb, G, h, A, b, dims, kktsolver);
    cp->F = F;
    cp->Gf = Gf;
    cp->h = h;
    cp->Af = Af;
    cp->b = b;

    cvxc_cp_finalize_setup(cp, nvars, mb, dims, kktsolver);
    /* // fix index set G/h matrices */
    /* cvxc_cpl_internal_t *cpi = cp->u.cpl; */
    /* cp->index_g = &cpi->index_cpt; */
    /* cp->c = &cpi->c0; */
    /* cvxc_mgrp_init(&cpi->h_g, cp->h, &cpi->index_cpt); */

    /* // init KKT solver */
    /* if (kktsolver) { */
    /*     cp->solver = kktsolver; */
    /* } else { */
    /*     cvxc_ldlsolver_init(&cp->__S, cp, nvars, mb, dims); */
    /*     cp->solver = &cp->__S; */
    /* } */
    /* cvxc_cp_solver_init(&cpi->cp_solver, cp, cp->solver); */
    /* cp->solver = (cvxc_kktsolver_t *)&cpi->cp_solver; */

    return used;
}

cvxc_size_t cvxc_cp_setup_user(
    cvxc_problem_t *cp,
    cvxc_convex_program_t *F,
    cvxc_size_t nvars,
    cvxc_umatrix_t *Gf,
    cvxc_matrix_t *h,
    cvxc_umatrix_t *Af,
    cvxc_matrix_t *b,
    const cvxc_dimset_t *dims,
    cvxc_kktsolver_t *kktsolver)
{
    if (! cp)
        return 0;

    int err;
    if ((err = cvxc_cp_isok(nvars, F, __cvxnil, h, __cvxnil, b, dims)) != 0) {
        cp->error = err;
        return 0;
    }

    return cvxc_cp_create(cp, F, nvars, &cp->Gu, h, &cp->Au, b, dims, kktsolver);
}

cvxc_size_t cvxc_cp_setup(
    cvxc_problem_t *cp,
    cvxc_convex_program_t *F,
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

    cvxc_size_t mG = 0, nG = 0;
    cvxm_size(&mG, &nG, G);
    if ((err = cvxc_cp_isok(nG, F, G, h, A, b, dims)) != 0) {
        cp->error = err;
        return 0;
    }

    cvxc_umat_make(&cp->Au, A);
    cvxc_umat_make(&cp->Gu, G);
    cp->A = A;
    cp->G = G;

    int stat = cvxc_cp_create(cp, F, nG, &cp->Gu, h, &cp->Au, b, dims, kktsolver);
    if (stat <= 0) {
        cvxc_umat_clear(&cp->Au);
        cvxc_umat_clear(&cp->Gu);
    }
    return stat;
}

/**
 * @brief Release resources reserved to program.
 */
void cvxc_cp_release(cvxc_problem_t *cp)
{
    if (!cp || !cp->u.cpl)
        return;

    cvxc_kktrelease(cp->solver);

    free(cp->u.space);
    cp->u.space = 0;
    cp->nbytes = 0;
    cp->work = 0;
}

/**
 * @brief Compute starting point for solver.
 */
int cvxc_cp_compute_start(cvxc_problem_t *cp)
{
    cvxc_cpl_internal_t *cpi = cp->u.cpl;
    F0(cp->F, &cpi->x0_e);
    cvxc_mgrp_initial_value(&cpi->z_g, 0);
    cvxc_mgrp_initial_value(&cpi->s_g, 0);
    cvxc_scaling_initial_value(&cpi->W);
    cvxm_set_all(&cpi->c0, 0.0);
    cvxm_set_epival(&cpi->c_e, 1.0);

    cvxm_epi_copy(&cpi->x_e,    &cpi->x0_e, 0);
    cvxm_epi_copy(&cpi->rx_e,   &cpi->x0_e, 0);
    cvxm_epi_copy(&cpi->dx_e,   &cpi->x0_e, 0);
    cvxm_epi_copy(&cpi->rx0_e,  &cpi->x0_e, 0);
    // cvxm_epi_copy(&cpi->dx0_e,  &cpi->x0_e, 0);
    cvxm_epi_copy(&cpi->newx_e, &cpi->x0_e, 0);
    cvxm_epi_copy(&cpi->newrx_e,&cpi->x0_e, 0);
    cvxm_copy(&cpi->ry,    cp->b,   0);
    return 0;
}

/**
 * @brief Solve convex program with convex target function.
 */
int cvxc_cp_solve(cvxc_solution_t *sol, cvxc_problem_t *cp, cvxc_solopts_t *opts)
{
    return cvxc_cpl_solve(sol, cp, opts);
}
