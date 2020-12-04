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
#include "cvxc.h"

enum long_opt_enum {
    OPT_SPDA = 260,
    OPT_JSON = 261
};


static struct option long_options[] = {
    {"matrix-G",   required_argument, 0, 'G'},
    {"matrix-A",   required_argument, 0, 'A'},
    {"matrix-b",   required_argument, 0, 'b'},
    {"matrix-c",   required_argument, 0, 'c'},
    {"matrix-h",   required_argument, 0, 'h'},
    {"dims",       required_argument, 0, 'd'},
    {"program",    required_argument, 0, 'P'},
    {"mm-mask",    required_argument, 0, 'm'},
    {"max-iter",   required_argument, 0, 'i'},
    {"linear",     required_argument, 0, 'L'},
    {"socp",       required_argument, 0, 'Q'},
    {"sdp" ,       required_argument, 0, 'S'},
    {"json",       no_argument,       0, 'j'},
    {"verbose",    no_argument,       0, 'v'},
    {"refinement", required_argument, 0, 'r'},
    {0,            0,                 0, 0}
};

#define SHORT_OPTS "A:G:c:h:b:d:P:L:Q:S:m:i:jv"

struct solver_args {
    char *name_G;
    char *name_h;
    char *name_c;
    char *name_A;
    char *name_b;
    char *name_dims;
    char *name_prob;
    char *name_cp;
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

cvxc_size_t *parse_intlist(cvxc_size_t *n, char *arg)
{
    char *tok, *cp;
    int k, count = 0;
    if (strlen(arg) == 0)
        return (cvxc_size_t *)0;

    for (cp = strchr(arg, ','); cp; cp = strchr(cp+1, ',')) {
        count += 1;
    }
    // we have count+1 values, need count+2 space to zero fill last entry
    cvxc_size_t *p = (cvxc_size_t *)calloc(count+2, sizeof(cvxc_size_t));
    if (!p)
        return p;;
    if (n)
        *n = count + 1;

    if (count == 1) {
        p[0] = atoi(arg);
        p[1] = 0;
        return p;
    }

    cp = arg; k = 0;
    for (tok = strsep(&cp, ","); *tok; tok = strsep(&cp, ",")) {
        p[k++] = atoi(tok);
        if (k == count+1)
            break;
    }
    return p;
}

/*
 * read matrix market formatted file
 */
int solver_read_matrix(cvxc_matrix_t **A, const char *path)
{
    int err;
    cvxc_stream_t ios;
    if (*A) {
        cvxm_init(*A, 0, 0);
    }
    if (!path) {
        return 0;
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
    cvxc_json_write_params(&ios, p->opts, p->dims,
                           p->c, p->G, p->h, p->A, p->b);
}

static
int solver_have_file_args(struct solver_args *args)
{
    return args->name_G
        || args->name_h
        || args->name_c
        || args->name_A
        || args->name_b
        || args->name_prob;
}

int solve_conelp(struct solver_args *args)
{
    cvxc_solopts_t opts = {0};
    cvxc_params_t *pptr, params = {0};
    pptr = &params;

    opts.refinement = args->refinement;
    opts.show_progress = args->verbose;
    opts.max_iter = args->maxiter;

    if (!solver_have_file_args(args) || args->name_prob) {
        if (solver_read_params(&pptr, args->name_prob) < 0)
            return -1;
    }
    if (solver_read_matrix(&params.G, args->name_G) < 0)
        return -1;
    if (solver_read_matrix(&params.h, args->name_h) < 0)
        return -1;
    if (solver_read_matrix(&params.c, args->name_c) < 0)
        return -1;
    if (solver_read_matrix(&params.A, args->name_A) < 0)
        return -1;
    if (solver_read_matrix(&params.b, args->name_b) < 0)
        return -1;

    cvxc_size_t rows, cols;
    cvxm_size(&rows, &cols, params.h);
    if (solver_init_dims(&params.dims, args, rows) < 0)
        return -1;

    if (!params.G || !params.h || !params.c) {
        fprintf(stderr, "matrices c, G, h need to defined\n");
        return -1;
    }
    if ((!params.A && params.b) || (params.A && !params.b)) {
        fprintf(stderr, "must have both  A and b\n");
        return -1;
    }

    if (!params.A && !params.b) {
        cvxc_size_t rows, cols;
        cvxm_size(&rows, &cols, params.G);
        params.A = cvxm_new(0, cols);
        params.b = cvxm_new(0, 1);
    }
    // print_params(&params);

    cvxc_problem_t cp;
    cvxc_size_t nbytes =
        cvxc_conelp_setup(&cp, params.c, params.G, params.h,
                          params.A, params.b, params.dims, (cvxc_kktsolver_t *)0);
    if (nbytes == 0) {
        fprintf(stderr, "Conelp problem setup failed\n");
        return -1;
    }
    cvxc_conelp_compute_start(&cp);
    cvxc_conelp_solve(&cp, params.opts ? params.opts : &opts);

    solver_write_solution("-", &cp.solution);

    return 0;
}

int main(int argc, char **argv)
{
    int opt;
    int option_index;

    struct solver_args args = {0};

    while ((opt = getopt_long(argc, argv, SHORT_OPTS, long_options, &option_index)) != -1) {
        switch (opt) {
        case 'A':
            if (optarg)
                args.name_A = optarg;
            break;
        case 'G':
            if (optarg)
                args.name_G = optarg;
            break;
        case 'h':
            if (optarg)
                args.name_h = optarg;
            break;
        case 'b':
            if (optarg)
                args.name_b = optarg;
            break;
        case 'c':
            if (optarg)
                args.name_c = optarg;
            break;
        case 'd':
            if (optarg)
                args.name_dims = optarg;
        case 'm':
            // mm_mask = optarg;
            break;
        case 'i':
            args.maxiter = atoi(optarg);
            break;
        case 'L':
            args.L = atoi(optarg);
            break;
        case 'Q':
            args.Q = parse_intlist(&args.nQ, optarg);
            break;
        case 'S':
            args.S = parse_intlist(&args.nS, optarg);
            break;
        case 'j':
            args.json = 1;
            break;
        case 'r':
            args.refinement = atoi(optarg);
            break;
        case 'v':
            args.verbose++;
            break;
        default:
            fprintf(stderr, "unknown option: %c\n", opt);
            exit(1);
        }
    }
    int err = solve_conelp(&args);
    return err;
}
