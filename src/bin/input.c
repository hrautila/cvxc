/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */

#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include "solver.h"

int solver_read_matrix(cvxc_matrix_t **A, const char *path)
{
    int err;
    cvxc_stream_t ios;
    if (!path) {
        return 0;
    }
    if (*A) {
        cvxm_init(*A, 0, 0);
    }
    int json = strstr(path, ".json") != 0;

    FILE *fp = fopen(path, "r");
    if (!fp) {
        perror(path);
        return -1;
    }
    cvxc_file_stream(&ios, fp);
    if (json) {
        err = cvxc_json_matrix_read(A, &ios);
    } else {
        err = 0; // cvxm_mmload(A, fp);
    }
    fclose(fp);
    return err;
}

int solver_read_params(cvxc_params_t **pptr, const char *path)
{
    FILE *fp;
    cvxc_stream_t ios;

    if (!path || *path == '-') {
        fp = stdin;
    } else {
        fp = fopen(path, "r");
        if (!fp) {
            perror(path);
            return -1;
        }
    }
    cvxc_file_stream(&ios, fp);
    int err = cvxc_json_read_params(pptr, &ios);
    if (fp != stdin)
        fclose(fp);
    return err;
}

int solver_write_solution(const char *path, const cvxc_solution_t *sol)
{
    FILE *fp;
    cvxc_stream_t ios;

    if (!path || *path == '-') {
        fp = stdout;
    } else {
        fp = fopen(path, "w");
        if (!fp) {
            perror(path);
            return -1;
        }
    }
    cvxc_file_stream(&ios, fp);
    int err = cvxc_json_write_solution(&ios, sol);
    if (fp != stdout)
        fclose(fp);
    return err;
}

int solver_init_dims(cvxc_dimset_t **dims, struct solver_args *args, cvxc_size_t hrows)
{
    if (args->name_dims) {
        cvxc_stream_t ios;
        FILE *fp = fopen(args->name_dims, "r");
        if (fp) {
            cvxc_file_stream(&ios, fp);
        } else {
            cvxc_str_stream(&ios, args->name_dims, strlen(args->name_dims));
        }
        int err = cvxc_json_dimset_read(dims, &ios);
        if (fp)
            fclose(fp);
        return err;
    }
    if (!*dims) {
        *dims = cvxc_dimset_new(0, 0, 0, 0);
    }
    if (args->L + args->nL + args->nQ + args->nS > 0) {
        if (!cvxc_dimset_alloc(*dims, args->L, args->Q, args->S))
            return -1;
        if (args->nL > 0)
            (*dims)->mnl = args->nL;
    } else {
        /* No dimensions defined, assume linear with hrows contraints. */
        cvxc_dimset_create(*dims, 0, hrows, 0, 0);
    }
    return 0;
}

void print_params(cvxc_params_t *p)
{
    cvxc_stream_t ios;
    cvxc_file_stream(&ios, stdout);
    cvxc_json_write_params(&ios, p);
}

int solver_have_file_args(struct solver_args *args)
{
    return args->name_G
        || args->name_h
        || args->name_c
        || args->name_A
        || args->name_b
        || args->name_prob;
}

int solver_read_args(struct solver_args *args, cvxc_params_t **params, cvxc_solopts_t *opts)
{
    cvxc_params_t *pars;

    opts->refinement = args->refinement;
    opts->show_progress = args->verbose;
    opts->max_iter = args->maxiter;

    if (!solver_have_file_args(args) || args->name_prob) {
        if (solver_read_params(params, args->name_prob) < 0)
            return -1;
    }

    pars = *params;
    if (!pars) {
        *params = malloc(sizeof(cvxc_params_t));
        if (!(*params))
            return -1;
        pars = *params;
    }

    if (!pars->G && solver_read_matrix(&pars->G, args->name_G) < 0)
        return -1;
    if (!pars->h && solver_read_matrix(&pars->h, args->name_h) < 0)
        return -1;
    if (!pars->c && solver_read_matrix(&pars->c, args->name_c) < 0)
        return -1;
    if (!pars->A && solver_read_matrix(&pars->A, args->name_A) < 0)
        return -1;
    if (!pars->b && solver_read_matrix(&pars->b, args->name_b) < 0)
        return -1;

    cvxc_size_t rows, cols;
    cvxm_size(&rows, &cols, pars->h);
    if (!pars->dims && solver_init_dims(&pars->dims, args, rows) < 0)
        return -1;

    if (!pars->module) {
        pars->module = args->name_cp;
        pars->args = args->args_cp;
    }
    return 0;
}
