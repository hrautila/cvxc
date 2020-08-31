
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include <unistd.h>
#include "convex.h"

int main(int argc, char **argv)
{
    cvxc_matrix_t G, h, A, b, F, g;
    cvxc_matrix_t x, z;
    cvxc_problem_t cp;
    cvxc_dimset_t dims;
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

    cvxc_float_t xdata[3] = {1.0, 1.0, 1.0};
    cvxc_float_t zdata[7] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    cvxc_solopts_t opts = (cvxc_solopts_t){
        .abstol = 0.0,
        .reltol = 0.0,
        .feastol = 0.0,
        .max_iter = 1,
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

    cvxc_dimset_alloc(&dims, 0, (cvxc_size_t *)0, (cvxc_size_t *)0);

    cvxm_map_data(&F, 8, 3, fdata);
    cvxm_map_data(&g, 8, 1, gdata);
    cvxm_map_data(&h, 0, 1, (cvxc_float_t *)0);
    cvxm_map_data(&G, 0, 3, (cvxc_float_t *)0);
    cvxm_map_data(&A, 0, 3, (cvxc_float_t *)0);
    cvxm_map_data(&b, 0, 1,  (cvxc_float_t *)0);

    if (opts.max_iter == 0)
        return 0;

    cvxc_gp_setup(&cp, K, &F, &g, &G, &h, &A, &b, &dims, (cvxc_kktsolver_t *)0);

    cvxc_size_t m, n;
    cvxm_map_data(&x, 3, 1, xdata);
    cvxm_map_data(&z, 7, 1, zdata);

    cvxm_size(&m, &n, &cp.u.cpl.f);
    printf("f : [%ld,%ld]\n", m, n);

    cvxm_size(&m, &n, &cp.u.cpl.Df);
    printf("Df: [%ld,%ld]\n", m, n);

    cvxm_size(&m, &n, &cp.u.cpl.H);
    printf("H : [%ld,%ld]\n", m, n);

    //cvxc_gp_f(&cp.u.cpl.f, &cp.u.cpl.Df, &cp.u.cpl.H, &x, &z);
    //cvxc_gp_set_start(&cp, __cvxnil, __cvxnil, __cvxnil, __cvxnil);
    //cvxc_gp_solve(&cp, &opts);
    return 0;
}

