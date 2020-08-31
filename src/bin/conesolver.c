
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include <getopt.h>
#include <string.h>

#include "cvxc.h"

enum long_opt_enum {
    OPT_SPDA = 260,
    OPT_JSON = 261
};


static struct option long_options[] = {
#if 0
    {"matrixG",   required_argument, 0, 'G'},
    {"matrixA",   required_argument, 0, 'A'},
    {"matrixb",   required_argument, 0, 'b'},
    {"matrixc",   required_argument, 0, 'c'},
    {"matrixh",   required_argument, 0, 'h'},
#endif
    {"mm-mask",   required_argument, 0, 'M'},
    {"max-iter",  required_argument, 0, 'I'},
    {"linear",    required_argument, 0, 'L'},
    {"socp",      required_argument, 0, 'Q'},
    {"sdp" ,      required_argument, 0, 'S'},
    {"json",      required_argument, 0, 'J'},
};

//#define SHORT_OPTS "A:G:c:h:b:L:Q:S:"
#define SHORT_OPTS "L:Q:S:M:I:"

int *parse_intlist(int *n, char *arg)
{
    char *tok, *cp;
    int k, count = 0;
    if (strlen(arg) == 0)
        return (int *)0;

    for (cp = strchr(arg, ','); cp; cp = strchr(cp+1, ',')) {
        count += 1;
    }
    // we have count+1 values, need count+2 space to zero fill last entry
    int *p = (int *)calloc(count+2, sizeof(int));
    if (!p)
        return (int *)0;
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

// read matrix market formatted file
int mmread_file(cvxc_matrix_t *A, const char *path)
{
    cvxm_init(A, 0, 0);
    if (!path) {
        return 0;
    }
    FILE *fp = fopen(path, "r");
    if (!fp) {
        perror(path);
        return -1;
    }
    int err = cvxm_mmload(A, fp);
    fclose(fp);
    return err;
}

#define MM_MASK_DEFAULT "matrix_%c.mtx"

void print_solution(cvxc_solution_t *sol)
{
    printf("status      : %2d [%s]\n", sol->status, solution_name[sol->status]);
    printf("primal obj  : %13.6e\n", sol->primal_objective);
    printf("dual obj    : %13.6e\n", sol->dual_objective);
    printf("primal inf  : %13.6e\n", sol->primal_infeasibility);
    printf("dual int    : %13.6e\n", sol->dual_infeasibility);
    printf("primal slack: %13.6e\n", sol->primal_slack);
    printf("dual slack  : %13.6e\n", sol->dual_slack);
    printf("primal cert : %13.6e\n", sol->primal_residual_cert);
    printf("dual cert   : %13.6e\n", sol->dual_residual_cert);
    printf("gap         : %13.6e\n", sol->gap);
    printf("relative gap: %13.6e\n", sol->relative_gap);
    printf("iterations  : %d\n", sol->iterations);
    if (sol->status != CVXC_STAT_OPTIMAL) 
        return;
    cvxc_mat_printf(stdout, "%13.6e", sol->x, "x");
    cvxc_mat_printf(stdout, "%13.6e", sol->s, "s");
    cvxc_mat_printf(stdout, "%13.6e", sol->y, "y");
    cvxc_mat_printf(stdout, "%13.6e", sol->z, "z");
}


int main(int argc, char **argv)
{
    char *mm_mask = MM_MASK_DEFAULT;
    char name[512];

    int *dim_Q = (int *)0;
    int *dim_S = (int *)0;
    int nQ, nS;
    int opt;
    int dim_L;
    int option_index;
    int max_iter = 100;

    cvxc_matrix_t A, G, c, b, h;
    cvxc_dimset_t dims;

    while ((opt = getopt_long(argc, argv, SHORT_OPTS, long_options, &option_index)) != -1) {
        switch (opt) {
        case 'M':
            mm_mask = optarg;
            break;
        case 'I':
            max_iter = atoi(optarg);
            break;
        case 'L':
            dim_L = atoi(optarg);
            break;
        case 'Q':
            dim_Q = parse_intlist(&nQ, optarg);
            break;
        case 'S':
            dim_S = parse_intlist(&nS, optarg);
            break;
        }
    }

    cvxc_dimset_alloc(&dims, dim_L, dim_Q, dim_S);

    snprintf(name, sizeof(name), mm_mask, 'c');
    mmread_file(&c, name);

    snprintf(name, sizeof(name), mm_mask, 'G');
    mmread_file(&G, name);

    snprintf(name, sizeof(name), mm_mask, 'h');
    mmread_file(&h, name);

    snprintf(name, sizeof(name), mm_mask, 'A');
    mmread_file(&A, name);

    snprintf(name, sizeof(name), mm_mask, 'b');
    mmread_file(&b, name);

    cvxc_conelp_problem_t cp;
    cvxc_solopts_t opts = (cvxc_solopts_t){
        .abstol = 0.0,
        .reltol = 0.0,
        .feastol = 0.0,
        .max_iter = max_iter,
        .debug = 0,
        .refinement = 0,
        .kkt_solver_name = 0,
        .show_progress = 1
    };

    cvxc_size_t nbytes =
        cvxc_conelp_setup(&cp, &c, &G, &h, &A, &b, &dims, (cvxc_kktsolver_t *)0);
    if (nbytes == 0) {
        fprintf(stderr, "Conelp problem setup failed\n");
        exit(1);
    }
    cvxc_conelp_compute_start(&cp);
    cvxc_conelp_solve(&cp, &opts);
    print_solution(&cp.solution);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
