
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include <unistd.h>
#include <math.h>
#include "cvxc.h"

extern int print_solution(cvxc_solution_t *sol);

int main(int argc, char **argv)
{
    cvxc_matrix_t G, h, A, b, F, g;
    cvxc_problem_t cp;
    cvxc_gpindex_t gpi;
    int opt;

    cvxc_float_t aflr = 1000.0;
    cvxc_float_t awall = 100.0;
    cvxc_float_t alpha = 0.5;
    cvxc_float_t beta = 2.0;
    cvxc_float_t gamma = 0.5;
    cvxc_float_t delta = 2.0;

    // 8x3
    cvxc_float_t fdata[24] = {
        -1.0, 1.0, 1.0, 0.0, -1.0, 1.0, 0.0, 0.0,
        -1.0, 1.0, 0.0, 1.0, 1.0, -1.0, 1.0, -1.0,
        -1.0, 0.0, 1.0, 1.0, 0.0, 0.0, -1.0, 1.0
    };

    cvxc_float_t gdata[8] = {
        1.0, 2.0/awall, 2.0/awall, 1.0/aflr, alpha, 1.0/beta, gamma, 1.0/delta
    };
    cvxc_size_t K[8] = { 1, 2, 1, 1, 1, 1, 1, 0 };

    cvxc_solopts_t opts = (cvxc_solopts_t){
        .abstol = 0.0,
        .reltol = 0.0,
        .feastol = 0.0,
        .max_iter = 30,
        .debug = 0,
        .refinement = 0,
        .kkt_solver_name = 0,
        .show_progress = 1
    };

    while ((opt = getopt(argc, argv, "N:")) != -1) {
        switch (opt) {
        case 'N':
            opts.max_iter = atoi(optarg);
            break;
        default:
            break;
        }
    }

    cvxc_gpi_init(&gpi, K, 7);

    cvxm_map_data(&F, 8, 3, fdata);
    cvxm_map_data(&g, 8, 1, gdata);
    cvxm_map_data(&h, 0, 1, (cvxc_float_t *)0);
    cvxm_map_data(&G, 0, 3, (cvxc_float_t *)0);
    cvxm_map_data(&A, 0, 3, (cvxc_float_t *)0);
    cvxm_map_data(&b, 0, 1,  (cvxc_float_t *)0);

    // g = log(g)
    cvxm_apply(&g, log, 0);

    if (opts.max_iter == 0)
        return 0;

    cvxc_solution_t solution = {0};

    cvxc_gp_setup(&cp, &gpi, &F, &g, &G, &h, &A, &b, (cvxc_kktsolver_t *)0);
    cvxc_gp_compute_start(&cp);
    cvxc_gp_solve(&solution, &cp, &opts);

    // x = exp(x)
    cvxm_apply(solution.x, exp, 0);

    print_solution(&solution);
    cvxc_gpi_release(&gpi);

    int ok = solution.status == CVXC_STAT_OPTIMAL && solution.iterations == 13;
    printf("test_gp: %s\n", ok ? "OK" : "FAILED");

    exit(1 - ok);
}
