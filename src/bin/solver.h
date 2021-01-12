/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */
#ifndef CVXC_SOLVER_H
#define CVXC_SOLVER_H 1

#include "cvxc.h"

#define PATH_SEPARATOR '/'
#define SO_EXT "so"

struct solver_args {
    char *name_G;
    char *name_h;
    char *name_c;
    char *name_A;
    char *name_b;
    char *name_dims;
    char *name_prob;
    char *name_cp;
    char *args_cp;
    char *type;
    char *output;
    cvxc_size_t L;
    cvxc_size_t nL;
    cvxc_size_t nQ;
    cvxc_size_t nS;
    cvxc_size_t *Q;
    cvxc_size_t *S;
    cvxc_size_t maxiter;
    int verbose;
    int refinement;
    int json;
};


extern cvxc_size_t *parse_intlist(cvxc_size_t *n, char *arg);
extern int solver_read_matrix(cvxc_matrix_t **A, const char *path);
extern int solver_read_params(cvxc_params_t **pptr, const char *path);
extern int solver_write_solution(const char *path, const cvxc_solution_t *sol);
extern int solver_init_dims(cvxc_dimset_t **dims, struct solver_args *args, cvxc_size_t hrows);
extern void print_params(cvxc_params_t *p);
extern int solver_have_file_args(struct solver_args *args);
extern int solver_read_args(struct solver_args *args, cvxc_params_t **params, cvxc_solopts_t *opts);
extern int solver_is_shared_name(const char *name);
extern char *solver_type_from_name(char *name);
extern char *solver_find_shared(const char *name, const char *var);
extern void *solver_load_shared(const char *name);
extern int solver_cp_init(void *handle, cvxc_convex_program_t *F, const char *sym, void *ptr);

extern int solve_conelp(cvxc_params_t *params, cvxc_solopts_t *opts, struct solver_args *args);
extern int solve_cpl(cvxc_params_t *params, cvxc_solopts_t *opts, struct solver_args *args);
extern int solve_cp(cvxc_params_t *params, cvxc_solopts_t *opts, struct solver_args *args);

#endif // CVXC_SOLVER_H
