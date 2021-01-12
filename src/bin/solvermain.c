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
#include <unistd.h>
#include "solver.h"

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
    {"convex",     required_argument, 0, 'F'},
    {"args",       required_argument, 0, 'a'},
    {"output",     required_argument, 0, 'o'},
    {"mm-mask",    required_argument, 0, 'm'},
    {"max-iter",   required_argument, 0, 'i'},
    {"linear",     required_argument, 0, 'L'},
    {"socp",       required_argument, 0, 'Q'},
    {"sdp" ,       required_argument, 0, 'S'},
    {"json",       no_argument,       0, 'j'},
    {"verbose",    no_argument,       0, 'v'},
    {"version",    no_argument,       0, 'V'},
    {"refinement", required_argument, 0, 'r'},
    {"type",       required_argument, 0, 'T'},
    {"shell-template", no_argument,   0,   0},
    {0,            0,                 0, 0}
};

#define SHORT_OPTS "A:G:c:h:b:d:F:P:a:o:L:Q:S:m:i:T:jvV"


static char *shell_template[] = {
    "#!/bin/bash",
    "",
    "images=$(docker images | awk '$1 ~ /cvxc/ {printf \"%s:%s\\n\", $1, $2}')",
    "",
    "if [ \"$images\" = \"\" ]; then",
    "    echo \"No cvxc images. Pull one!\" >&2",
    "    exit 1",
    "fi",
    "",
    "docker_image=$(echo \"$images\" | tr ' ' '\\n' | awk -F: '$2 ~ /latest/ {print $0}')",
    "if [ \"$docker_image\" = \"\" ]; then",
    "    docker_image=$(echo \"$images\" | tr ' ' '\\n' | head -1)",
    "fi",
    "",
    "istty=\"-t\"",
    "if [ ! -t 0 ]; then",
    "    istty=",
    "fi",
    "",
    "docker run --rm -i $istty --volume $PWD:/var/run $docker_image $*",
    "",
    0
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
        if (k == count + 1)
            break;
    }
    return p;
}

int solver_is_convex(const char *type)
{
    return strncmp(type, "cpl", 3) == 0 || strncmp(type, "cp", 2) == 0;
}

int handle_long_only_option(const char *name, char *optarg, struct solver_args *args)
{
    if (strncmp(name, "shell-template", 14) == 0) {
        for (int i = 0; shell_template[i]; i++) {
            fprintf(stdout, "%s\n", shell_template[i]);
        }
        exit(0);
    }
    return 0;
}

int main(int argc, char **argv)
{
    int opt;
    int option_index;
    cvxc_size_t rows, cols;

    struct solver_args args = {0};
    cvxc_solopts_t opts = {0};
    cvxc_params_t *pptr, params = {0};
    pptr = &params;

    while ((opt = getopt_long(argc, argv, SHORT_OPTS, long_options, &option_index)) != -1) {
        switch (opt) {
        case 0:
            handle_long_only_option(long_options[option_index].name, optarg, &args);
            break;
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
            break;
        case 'P':
            if (optarg)
                args.name_prob = optarg;
            break;
        case 'F':
            if (optarg)
                args.name_cp = optarg;
            break;
        case 'a':
            if (optarg)
                args.args_cp = optarg;
            break;
        case 'o':
            if (optarg)
                args.output = optarg;
            break;
        case 'T':
            if (optarg)
                args.type = optarg;
            break;
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
        case 'V':
            // print version and exit
            break;
        default:
            fprintf(stderr, "unknown option: %c\n", opt);
            exit(1);
        }
    }
    int err;

    if (!args.type) {
        args.type = solver_type_from_name(argv[0]);
    }

    if (solver_read_args(&args, &pptr, &opts) < 0)
        return -1;

    if (!params.G || !params.h) {
        fprintf(stderr, "error: matrices G, h must be defined\n");
        return -1;
    }
    if ((!params.A && params.b) || (params.A && !params.b)) {
        fprintf(stderr, "error: must have both  A and b\n");
        return -1;
    }

    if (!params.A && !params.b) {
        cvxm_size(&rows, &cols, params.G);
        params.A = cvxm_new(0, cols);
        params.b = cvxm_new(0, 1);
    }

    if (params.module) {
        /* This is either CP or CPL */
        cvxm_size(&rows, &cols, params.c);
        if (rows * cols == 0) {
            err = solve_cp(&params, &opts, &args);
        } else {
            err = solve_cpl(&params, &opts, &args);
        }
    } else {
        /* This is CONELP problem. */
        err = solve_conelp(&params, &opts, &args);
    }
    // solver_release_params(&params);

    return err;
}
