
// Copyright: Harri Rautila, 2018-2019 <harri.rautila@gmail.com>

#include "cvxc.h"

#ifndef __ZERO
#define __ZERO 0.0
#endif


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

    int err = cvx_kktsolve(cp_solver->next, x, y, &z_cpl_g);

    cvxm_set(&z_nlt, 0, 0, -epi_x*dnlt0);
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

int cvx_cp_setvars(cvx_problem_t *cp,
                   cvx_convex_program_t *F,
                   cvx_size_t n, cvx_size_t m,
                   cvx_matrix_t *G,
                   cvx_matrix_t *h,
                   cvx_matrix_t *A,
                   cvx_matrix_t *b,
                   const cvx_dimset_t *dims,
                   cvx_kktsolver_t *kktsolver)
{
    cvx_cpl_setvars(cp, F, n, m, __cvxnil, G, h, A, b, dims, kktsolver);

    // fix index set G/h matrices
    cvx_cpl_internal_t *cpi = &cp->u.cpl;
    cp->index_g = &cpi->index_cpt;
    cp->c = &cpi->c0;
    cvx_mgrp_init(&cpi->h_g, cp->h, &cpi->index_cpt);

    // init KKT solver
    if (kktsolver) {
        cp->solver = kktsolver;
    } else {
        cvx_ldlsolver_init(&cp->__S, cp, n, m, dims);
        cp->solver = &cp->__S;
    }
    cvx_cp_solver_init(&cpi->cp_solver, cp, cp->solver);
    cp->solver = (cvx_kktsolver_t *)&cpi->cp_solver;
    return 0;
}

cvx_size_t cvx_cp_setup(cvx_problem_t *cp,
                        cvx_convex_program_t *F,
                        cvx_matrix_t *G,
                        cvx_matrix_t *h,
                        cvx_matrix_t *A,
                        cvx_matrix_t *b,
                        const cvx_dimset_t *dims,
                        cvx_kktsolver_t *kktsolver)
{
    cvx_size_t mG, nG, mb, nb;
    cvx_size_t used;
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

    if ((used = cvx_cpl_allocate(cp, 1, nG, mb, 0, dims)) == 0) {
        return 0;
    }
    cvx_cp_setvars(cp, F, nG, mb, G, h, A, b, dims, kktsolver);

    return used;
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
