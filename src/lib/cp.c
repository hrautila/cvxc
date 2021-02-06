/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */
#include "cvxc.h"

#ifndef __ZERO
#define __ZERO 0.0
#endif


int cvxc_cp_isok(const cvxc_convex_program_t *F,
                const cvxc_matrix_t *G,
                const cvxc_matrix_t *h,
                const cvxc_matrix_t *A,
                const cvxc_matrix_t *b,
                const cvxc_dimset_t *dims)
{
    cvxc_size_t mG, nG, mA, nA, mh, nh, mb, nb;

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
    cvxc_size_t cdim = cvxc_dimset_sum_squared(dims, CVXDIM_CONVEXPROG);

    if (nh > 1 || mh != cdim - 1) {
        return CVXC_ERR_DIMH;
    }

    if (mG != cdim - 1) {
        return CVXC_ERR_DIMG;
    }
    if (nA != nG || mA != mb) {
        return CVXC_ERR_DIMA;
    }
    if (nb != 1) {
        return CVXC_ERR_DIMB;
    }
    return 0;
}

static
int cp_factor(cvxc_kktsolver_t *S, cvxc_scaling_t *W, cvxc_matrix_t *x, cvxc_matrix_t *z)
{
    cvxc_chainsolver_t *cp_solver = (cvxc_chainsolver_t *)S;
    cvxc_matrix_t Dfc;
    cvxc_problem_t *cp = cp_solver->cp;
    cvxc_cpl_internal_t *cpi = &cp->u.cpl;
    cvxc_size_t nr, nc;

    cvxm_size(&nr, &nc, &cpi->Df);
    cvxm_view_map(&Dfc, &cpi->Df, 1, 0, nr-1, nc);
    F2(cp->F, &cpi->f, &cpi->Df, &cpi->H, x, z);

    return cvxc_kktfactor(cp_solver->next, W, &cpi->H, &Dfc);
}

static
int cp_solve(cvxc_kktsolver_t *S, cvxc_matrix_t *x, cvxc_matrix_t *y, cvxc_matgrp_t *z_g)
{
    cvxc_chainsolver_t *cp_solver = (cvxc_chainsolver_t *)S;
    cvxc_problem_t *cp = cp_solver->cp;
    cvxc_cpl_internal_t *cpi = &cp->u.cpl;
    cvxc_matrix_t gradf, dnlt, z_nlt, z_cpl;
    cvxc_index_t index_cpl;
    cvxc_matgrp_t z_cpl_g;
    cvxc_float_t znlt, epi_x, dnlt0;
    cvxc_size_t nr, nc;

    cvxm_size(&nr, &nc, &cpi->Df);
    cvxm_view_map(&gradf, &cpi->Df, 0, 0, 1, nc);
    cvxc_scaling_elem(&dnlt, &cp->W, CVXWS_DNLT, 0);
    cvxc_mgrp_elem(&z_nlt, z_g, CVXDIM_NLTARGET, 0);
    cvxc_mgrp_elem(&z_cpl, z_g, CVXDIM_CONVEXLP, 0);

    cvxc_subindex(&index_cpl, &cpi->index_full, CVXDIM_CONVEXLP);
    cvxc_mgrp_init(&z_cpl_g, &z_cpl, &cpi->index_cpt);

    dnlt0 = cvxm_get(&dnlt, 0, 0);
    epi_x = cvxm_get_epi(x);
    znlt  = cvxm_get(&z_nlt, 0, 0);
    cvxm_axpy(x, epi_x, &gradf);

    int err = cvxc_kktsolve(cp_solver->next, x, y, &z_cpl_g);

    cvxm_set(&z_nlt, 0, 0, -epi_x*dnlt0);
    epi_x = cvxm_dot(&gradf, x) + dnlt0*dnlt0*epi_x - znlt;
    cvxm_set_epi(x, epi_x);
    return err;
}

static cvxc_kktfuncs_t cp_kktfunctions = {
    .factor = cp_factor,
    .solve  = cp_solve
};

cvxc_chainsolver_t *cvxc_cp_solver_init(cvxc_chainsolver_t *cs, cvxc_problem_t *cp, cvxc_kktsolver_t *next)
{
    cs->fnc = &cp_kktfunctions;
    cs->next = next;
    cs->cp = cp;
    return cs;
}

int cvxc_cp_setvars(cvxc_problem_t *cp,
                   cvxc_convex_program_t *F,
                   cvxc_size_t n, cvxc_size_t m,
                   cvxc_matrix_t *G,
                   cvxc_matrix_t *h,
                   cvxc_matrix_t *A,
                   cvxc_matrix_t *b,
                   const cvxc_dimset_t *dims,
                   cvxc_kktsolver_t *kktsolver)
{
    cvxc_cpl_setvars(cp, F, n, m, __cvxnil, G, h, A, b, dims, kktsolver);

    // fix index set G/h matrices
    cvxc_cpl_internal_t *cpi = &cp->u.cpl;
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

cvxc_size_t cvxc_cp_setup(cvxc_problem_t *cp,
                        cvxc_convex_program_t *F,
                        cvxc_matrix_t *G,
                        cvxc_matrix_t *h,
                        cvxc_matrix_t *A,
                        cvxc_matrix_t *b,
                        const cvxc_dimset_t *dims,
                        cvxc_kktsolver_t *kktsolver)
{
    cvxc_size_t mG, nG, mb, nb;
    cvxc_size_t used;
    int err;

    if (! cp)
        return 0;

    if ((err = cvxc_cp_isok(F, G, h, A, b, dims)) != 0) {
        cp->error = err;
        return 0;
    }

    mG = nG = mb = nb = 0;
    cvxm_size(&mG, &nG, G);
    if (b)
        cvxm_size(&mb, &nb, b);

    if ((used = cvxc_cpl_allocate(cp, 1, nG, mb, 0, dims)) == 0) {
        return 0;
    }
    cvxc_cp_setvars(cp, F, nG, mb, G, h, A, b, dims, kktsolver);

    return used;
}

int cvxc_cp_compute_start(cvxc_problem_t *cp)
{
    cvxc_cpl_internal_t *cpi = &cp->u.cpl;
    F0(cp->F, &cpi->x0);
    cvxc_mgrp_initial_value(&cpi->z_g, 0);
    cvxc_mgrp_initial_value(&cpi->s_g, 0);
    cvxc_scaling_initial_value(&cp->W);
    cvxm_set_all(&cpi->c0, 0.0);
    cpi->c0.t = 1.0;

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

int cvxc_cp_solve(cvxc_problem_t *cp,
                 cvxc_solopts_t *opts)
{
    return cvxc_cpl_solve(cp, opts);
}
