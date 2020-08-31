
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include "cvxm.h"
#include "cvxc.h"

char *solution_name[] = {
    "Optimal",
    "Unknown",
    "Primal Infeasible",
    "Dual Infeasible",
    "Singular"
};

cvxc_float_t gdata[] = {
    2.0, 1.0, -1.0, 0.0,
    1.0, 2.0, 0.0, -1.0
};
cvxc_float_t cdata[] = {-4.0, -5.0};
cvxc_float_t hdata[] = {3.0, 3.0, 0.0, 0.0};

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
    //cp.solver->debug = 2;

    cvxc_conelp_compute_start(&cp);
    cvxc_conelp_solve(&cp, &opts);
    print_solution(&cp.solution);
}
