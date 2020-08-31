
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include <unistd.h>
#include "cvxc.h"

char *solution_name[] = {
    "Optimal",
    "Unknown",
    "Primal Infeasible",
    "Dual Infeasible",
    "Singular"
};

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
    cvxc_matrix_t c, G, h, A, b;
    cvxc_problem_t cp;
    cvxc_dimset_t dims;
    int opt;

    cvxc_solopts_t opts = (cvxc_solopts_t){
        .abstol = 0.0,
        .reltol = 0.0,
        .feastol = 0.0,
        .max_iter = 15,
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
    printf("conelp setup: %ld bytes\n", nbytes);
    if (nbytes == 0)
        exit(1);

    //cp.solver->debug = 2;

    cvxc_conelp_compute_start(&cp);
    if (cvxc_conelp_solve(&cp, &opts) == 0)
        print_solution(&cp.solution);
}
