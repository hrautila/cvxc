/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */
#include <getopt.h>
#include <string.h>
#include <stdio.h>
#include <dlfcn.h>
#include "solver.h"

int solve_conelp(cvxc_params_t *params, cvxc_solopts_t *opts, struct solver_args *args)
{
    cvxc_size_t nbytes;
    cvxc_problem_t cp;
    cvxc_solution_t solution = {0};

    nbytes = cvxc_conelp_setup(&cp, params->c, params->G, params->h,
                               params->A, params->b, params->dims, (cvxc_kktsolver_t *)0);
    if (nbytes == 0) {
        fprintf(stderr, "error: conelp problem setup failed\n");
        return -1;
    }
    cvxc_conelp_compute_start(&cp);
    cvxc_conelp_solve(&solution, &cp, params->opts ? params->opts : opts);

    solver_write_solution(args->output, &solution);

    cvxc_conelp_release(&cp);
    return 0;
}

int solver_cp_init(void *handle, cvxc_convex_program_t *F, const char *sym, void *ptr)
{
    cvxc_cp_initfunc init_fn = (cvxc_cp_initfunc)dlsym(handle, sym);
    if (!init_fn) {
        perror(sym);
        return -1;
    }
    return (*init_fn)(F, ptr);
}

int solve_cpl(cvxc_params_t *params, cvxc_solopts_t *opts, struct solver_args *args)
{
    cvxc_size_t nbytes;
    cvxc_convex_program_t F;
    cvxc_solution_t solution = {0};

    void *handle = solver_load_shared(params->module);
    if (!handle)
        return -1;

    if (solver_cp_init(handle, &F, "cp_load", params->args) < 0) {
        dlclose(handle);
        return -1;
    }

    cvxc_problem_t cpl;
    nbytes = cvxc_cpl_setup(&cpl, &F, params->c, params->G, params->h,
                            params->A, params->b, params->dims, (cvxc_kktsolver_t *)0);
    if (nbytes == 0) {
        fprintf(stderr, "error: cpl problem setup failed\n");
        dlclose(handle);
        return -1;
    }
    cvxc_cpl_compute_start(&cpl);

    if (params->opts)
        opts = params->opts;
    opts->show_progress = args->verbose;
    opts->max_iter = args->maxiter;

    cvxc_cpl_solve(&solution, &cpl, opts);
    solver_write_solution(args->output, &solution);

    dlclose(handle);
    cvxc_cpl_release(&cp);
    return 0;
}

int solve_cp(cvxc_params_t *params, cvxc_solopts_t *opts, struct solver_args *args)
{
    cvxc_size_t nbytes;
    cvxc_convex_program_t F;
    cvxc_solution_t solution = {0};

    void *handle = solver_load_shared(params->module);
    if (!handle)
        return -1;

    if (solver_cp_init(handle, &F, "cp_load", params->args) < 0) {
        dlclose(handle);
        return -1;
    }

    cvxc_problem_t cp;
    params->dims->iscpt = 1;
    nbytes = cvxc_cp_setup(&cp, &F, params->G, params->h,
                           params->A, params->b, params->dims, (cvxc_kktsolver_t *)0);
    if (nbytes == 0) {
        fprintf(stderr, "error: cp problem setup failed\n");
        dlclose(handle);
        return -1;
    }
    cvxc_cp_compute_start(&cp);
    cvxc_cp_solve(&solution, &cp, params->opts ? params->opts : opts);
    solver_write_solution(args->output, &solution);

    dlclose(handle);
    cvxc_cp_release(&cp);

    return 0;
}

int solve_gp(cvxc_params_t *params, cvxc_solopts_t *opts, struct solver_args *args)
{
    cvxc_size_t nbytes;
    cvxc_problem_t cp;
    cvxc_solution_t solution = {0};

    nbytes = cvxc_gp_setup(&cp, params->K, params->F, params->g, params->G, params->h,
                           params->A, params->b, (cvxc_kktsolver_t *)0);
    if (nbytes == 0) {
        fprintf(stderr, "error: conelp problem setup failed\n");
        return -1;
    }
    cvxc_gp_compute_start(&cp);
    cvxc_gp_solve(&solution, &cp, params->opts ? params->opts : opts);

    solver_write_solution(args->output, &solution);
    cvxc_gp_release(&cp);
    return 0;
}
