
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include <stdio.h>
#include <unistd.h>
#include "convex.h"

extern void cvx_conelp_init_scaling(cvx_scaling_t *);

cvx_float_t gdata[] = {
    2.0, 1.0, -1.0, 0.0,
    1.0, 2.0, 0.0, -1.0
};
cvx_float_t cdata[] = {-4.0, -5.0};
cvx_float_t hdata[] = {3.0, 3.0, 0.0, 0.0};

int test_lp_init(int verbose)
{
    cvx_matrix_t c, G, h, A, b;
    cvx_conelp_problem_t cp;
    cvx_dimset_t dims;

    cvx_solopts_t opts = (cvx_solopts_t){
        .abstol = 0.0,
        .reltol = 0.0,
        .feastol = 0.0,
        .max_iter = 2,
        .debug = 0,
        .refinement = 1,
        .kkt_solver_name = 0
    };
    
    // 
    cvxm_map_data(&c, 2, 1, cdata);
    // inequality constraints; G*x <= h
    cvxm_map_data(&G, 4, 2, gdata);
    cvxm_map_data(&h, 4, 1, hdata);
    // equality constraints: Ax = b  (empty matrices if missing)
    cvxm_map_data(&A, 0, 2, (cvx_float_t *)0);
    cvxm_map_data(&b, 0, 1, (cvx_float_t *)0);

    cvx_dimset_alloc(&dims, 4, (int *)0, (int *)0);

    cvx_conelp_setup(&cp, &c, &G, &h, &A, &b, &dims, (cvx_kktsolver_t *)0);

    cvx_conelp_init_scaling(&cp.W);

    cp.__S.debug = verbose;
    int err = cvx_kktfactor(&cp.__S, &cp.W, __cvxnil, __cvxnil);
    if (err < 0) 
        printf("KKT factor error\n");
    if (verbose > 1) {
        fprintf(stderr, "post-factoring K\n");
        cvxm_printf(stderr, "%6.3f", &cp.__S.u.ldl.K);
    }
    
    return 0;
}

int main(int argc, char **argv)
{
    int verbose = 0;
    int opt;

    while ((opt = getopt(argc, argv, "v")) != -1) {
        switch (opt) {
        case 'v':
            verbose += 1;
            break;
        default:
            fprintf(stderr, "usage: tstmul [-v -f path] \n");
            exit(1);
        }
    }
    
    test_lp_init(verbose);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
