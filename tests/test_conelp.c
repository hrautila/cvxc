
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include <unistd.h>
#include "cvxc.h"


// one linear [2,1], two socp of size [4,1], one sdp of size [3,3]
cvxc_float_t gdata[] = {
    // 1st 
    16., 7., /**/ 24., -8., 8., -1., /**/ 0., -1., 0., 0.,/**/ 7., -5., 1., -5., 1., -7., 1., -7., -4.,
    // 2nd 
    -14., 2.,/**/ 7., -13., -18., 3.,/**/ 0., 0., -1., 0.,/**/ 3., 13., -6., 13., 12., -10., -6., -10., -28.,
    // 3rd
    5., 0., /**/-15., 12., -6., 17., /**/ 0., 0., 0., -1.,/**/ 9.,  6., -6., 6., -7., -7., -6., -7., -11.
};
cvxc_float_t cdata[] = {-6., -4., -5.};
cvxc_float_t hdata[] = {
    -3., 5.,/**/ 12., -2., -14., -13.,/**/ 10., 0., 0., 0.,/**/ 68., -30., -19., -30., 99., 23., -19., 23., 10.};

extern int print_solution(cvxc_solution_t *sol);

int main(int argc, char **argv)
{
    cvxc_matrix_t c, G, h, A, b;
    cvxc_problem_t cp;
    cvxc_dimset_t dims;
    int opt;

    cvxc_solopts_t opts = (cvxc_solopts_t){
        .abstol = 0.0,
        .reltol = 0.0,
        .feastol = 0.0,
        .max_iter = 15,
        .refinement = 0,
        .bits = CVXC_OPROGRESS
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
    //
    cvxm_map_data(&c, 3, 1, cdata);
    // inequality constraints; G*x <= h
    cvxm_map_data(&G, 19, 3, gdata);
    cvxm_map_data(&h, 19, 1, hdata);
    // equality constraints: Ax = b  (empty matrices if missing)
    cvxm_map_data(&A, 0, 3, (cvxc_float_t *)0);
    cvxm_map_data(&b, 0, 1, (cvxc_float_t *)0);

    cvxc_dimset_alloc(&dims, 2, (cvxc_size_t[]){4, 4, 0}, (cvxc_size_t[]){3, 0});

    cvxc_size_t nbytes = cvxc_conelp_setup(&cp, &c, &G, &h, &A, &b, &dims, (cvxc_kktsolver_t *)0);
    if (nbytes == 0)
        exit(1);

    cvxc_solution_t solution = {0};

    cvxc_conelp_compute_start(&cp);
    cvxc_conelp_solve(&solution, &cp, &opts);
    print_solution(&solution);

    int ok = solution.status == CVXC_STAT_OPTIMAL && solution.iterations == 12;
    printf("test_conelp: %s\n", ok ? "OK" : "FAILED");

    exit(1 - ok);
}
