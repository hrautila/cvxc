
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include "cvxc.h"


cvxc_float_t gdata[] = {
    2.0, 1.0, -1.0, 0.0,
    1.0, 2.0, 0.0, -1.0
};
cvxc_float_t cdata[] = {-4.0, -5.0};
cvxc_float_t hdata[] = {3.0, 3.0, 0.0, 0.0};

extern int print_solution(cvxc_solution_t *sol);

int main(int argc, char **argv)
{
    cvxc_matrix_t c, G, h, A, b;
    cvxc_problem_t cp;
    cvxc_dimset_t dims;

    cvxc_solopts_t opts = (cvxc_solopts_t){
        .abstol = 0.0,
        .reltol = 0.0,
        .feastol = 0.0,
        .max_iter = 7,
        .debug = 0,
        .refinement = 0,
        .kkt_solver_name = 0,
        .show_progress = 1
    };

    //
    cvxm_map_data(&c, 2, 1, cdata);
    // inequality constraints; G*x <= h
    cvxm_map_data(&G, 4, 2, gdata);
    cvxm_map_data(&h, 4, 1, hdata);
    // equality constraints: Ax = b  (empty matrices if missing)
    cvxm_map_data(&A, 0, 2, (cvxc_float_t *)0);
    cvxm_map_data(&b, 0, 1, (cvxc_float_t *)0);

    cvxc_dimset_alloc(&dims, 4, (cvxc_size_t *)0, (cvxc_size_t *)0);

    cvxc_conelp_setup(&cp, &c, &G, &h, &A, &b, &dims, (cvxc_kktsolver_t *)0);

    cvxc_solution_t solution = {0};

    cvxc_conelp_compute_start(&cp);
    cvxc_conelp_solve(&solution, &cp, &opts);
    print_solution(&solution);

    int ok = solution.status == CVXC_STAT_OPTIMAL && solution.iterations == 4;
    printf("test_lp: %s\n", ok ? "OK" : "FAILED");

    exit(1 - ok);
}
