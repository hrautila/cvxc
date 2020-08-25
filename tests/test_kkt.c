
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include <unistd.h>
#include "cvxm.h"
#include "convex.h"

char *solution_name[] = {
    "Optimal",
    "Unknown",
    "Primal Infeasible",
    "Dual Infeasible",
    "Singular"
};

cvx_float_t W_d[] = {
    7.330014e+00,
    5.394851e+00
};
cvx_float_t W_di[] = {
    1.364254e-01,
    1.853619e-01
};
cvx_float_t W_beta[] = {
    1.073844e+01,
    6.125768e+00
};
cvx_float_t W_v[] = {
    1.198764e+00,
    2.946933e-01,
    4.979940e-02,
    -5.896700e-01,
    /**/
    1.000156e+00,
    -8.571120e-03,
    1.046440e-02,
    -1.134278e-02
};
cvx_float_t W_r[] = {
    2.919851e+00,
    -3.036526e+00,
    -1.083566e+00, 
    // 2nd column
    -2.263271e+00, 
    -2.680673e+00,
    1.738179e-01, 
    // 3rd column
    5.926005e-01,
    -1.344990e-02,
    2.402415e+00
};
cvx_float_t W_rti[] = {
    // 1st
    1.722480e-01,
    -1.482364e-01,
    -4.331807e-02, 
    // 2nd 
    -1.955742e-01, 
    -2.048651e-01,
    4.709509e-02, 
    // 3rd
    9.183940e-02,
    -5.203710e-02,
    3.933025e-01
};

// initial set
cvx_float_t K_x[] = {
    6.000000e+00,
    4.000000e+00,
    5.000000e+00
};
cvx_float_t K_y[] = {
};
cvx_float_t K_z[] = {
    -3.000000e+00,
    5.000000e+00,
    1.200000e+01,
    -2.000000e+00,
    -1.400000e+01,
    -1.300000e+01,
    1.000000e+01,
    0.000000e+00,
    0.000000e+00,
    0.000000e+00,
    6.800000e+01,
    -3.000000e+01,
    -1.900000e+01,
    -3.000000e+01,
    9.900000e+01,
    2.300000e+01,
    -1.900000e+01,
    2.300000e+01,
    1.000000e+01
};
// result set
cvx_float_t R_x[] = {
    2.183512e-01,
    1.632292e-01,
    3.860374e-01
};
cvx_float_t R_y[] = {
};
cvx_float_t R_z[] = {
    8.374605e-01,
    -5.829787e-01,
    3.218721e-01,
    4.393832e-01,
    1.007852e+00,
    1.482787e+00,
    -1.634949e+00,
    -6.364560e-02,
    7.539681e-03,
    -1.000742e-01,   // SDP constraint lower triangular
    -6.389662e+00,
    -8.770923e-01,
    4.591148e-01,
    0.0, //-3.000000e+01,
    -3.675809e+00,
    6.494153e-01,
    0.0, //-1.900000e+01,
    0.0, //2.300000e+01,
    -1.829548e+00
};

// one linear [2,1], two socp of size [4,1], one sdp of size [3,3]
cvx_float_t gdata[] = {
    // 1st 
    16., 7., /**/ 24., -8., 8., -1., /**/ 0., -1., 0., 0.,/**/ 7., -5., 1., -5., 1., -7., 1., -7., -4.,
    // 2nd 
    -14., 2.,/**/ 7., -13., -18., 3.,/**/ 0., 0., -1., 0.,/**/ 3., 13., -6., 13., 12., -10., -6., -10., -28.,
    // 3rd
    5., 0., /**/-15., 12., -6., 17., /**/ 0., 0., 0., -1.,/**/ 9.,  6., -6., 6., -7., -7., -6., -7., -11.
};
cvx_float_t cdata[] = {-6., -4., -5.};
cvx_float_t hdata[] = {
    -3., 5.,/**/ 12., -2., -14., -13.,/**/ 10., 0., 0., 0.,/**/ 68., -30., -19., -30., 99., 23., -19., 23., 10.
};

void make_scaling(cvx_scaling_t *W, cvx_dimset_t *dims)
{
    cvx_matrix_t s, m;
    cvx_scaling_init(W, dims);

    cvxm_map_data(&s, 2, 1, W_d);
    cvx_scaling_elem(&m, W, CVXWS_D, 0);
    cvxm_copy(&m, &s, 0);

    cvxm_map_data(&s, 2, 1, W_di);
    cvx_scaling_elem(&m, W, CVXWS_DI, 0);
    cvxm_copy(&m, &s, 0);

    cvxm_map_data(&s, 2, 1, W_beta);
    cvx_scaling_elem(&m, W, CVXWS_BETA, 0);
    cvxm_copy(&m, &s, 0);

    cvxm_map_data(&s, 4, 1, &W_v[0]);
    cvx_scaling_elem(&m, W, CVXWS_V, 0);
    cvxm_copy(&m, &s, 0);

    cvxm_map_data(&s, 4, 1, &W_v[4]);
    cvx_scaling_elem(&m, W, CVXWS_V, 1);
    cvxm_copy(&m, &s, 0);

    cvxm_map_data(&s, 3, 3, W_r);
    cvx_scaling_elem(&m, W, CVXWS_R, 0);
    cvxm_copy(&m, &s, 0);

    cvxm_map_data(&s, 3, 3, W_rti);
    cvx_scaling_elem(&m, W, CVXWS_RTI, 0);
    cvxm_copy(&m, &s, 0);

    //cvx_scaling_printf(stdout, "%13.6e", W, "result");
}

void make_data(cvx_matrix_t *x, cvx_matrix_t *y, cvx_matrix_t *z)
{
    cvx_matrix_t m;
    cvxm_map_data(&m, 3, 1, K_x);
    cvxm_copy(x, &m, 0);
    
    cvxm_map_data(&m, 0, 1, K_y);
    cvxm_copy(y, &m, 0);

    cvxm_map_data(&m, 19, 1, K_z);
    cvxm_copy(z, &m, 0);
}

void compare_data(cvx_matrix_t *x, cvx_matrix_t *y, cvx_matrix_t *z)
{
    cvx_matrix_t m;
    cvxm_map_data(&m, 3, 1, R_x);
    cvxm_axpy(x, -1.0, &m);
    
    cvxm_map_data(&m, 0, 1, R_y);
    cvxm_axpy(y, -1.0, &m);

    cvxm_map_data(&m, 19, 1, R_z);
    cvxm_axpy(z, -1.0, &m);
}

extern cvx_kktfuncs_t *cvx_ldl2load();

int main(int argc, char **argv)
{
    cvx_matrix_t c, G, h, A, b;
    cvx_matrix_t x, y, z, m;
    cvx_matgrp_t z_g;
    cvx_dimset_t dims;
    cvx_problem_t cp;
    char *solver = "ldl";
    int opt;

    while ((opt = getopt(argc, argv, "N:S:")) != -1) {
        switch (opt) {
        case 'N':
            break;
        case 'S':
            solver = optarg;
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
    cvxm_map_data(&A, 0, 3, (cvx_float_t *)0);
    cvxm_map_data(&b, 0, 1, (cvx_float_t *)0);
    cvx_dimset_alloc(&dims, 2, (int[]){4, 4, 0}, (int[]){3, 0});
    
    make_scaling(&cp.W, &dims);
    cp.A = &A;
    cp.b = &b;
    cp.G = &G;
    cp.h = &h;
    cp.c = &c;
    cp.dims = &dims;
    cvx_index_init(&cp.index_full, &dims, CVX_INDEX_NORMAL);
    __mblk_init(&cp.work, 100);

    cvxm_init(&x, 3, 1);
    cvxm_init(&y, 0, 1);
    cvxm_init(&z, 19, 1);

    make_data(&x, &y, &z);
    cvx_mgrp_init(&z_g, &z, &cp.index_full);
    cvx_mgrp_init(&cp.z_g, &z, &cp.index_full);
    

    cvx_kktfuncs_t *kkt;
    printf("KKT solver: %s\n", solver);
    if (strcmp(solver, "ldl2") == 0) {
        kkt = cvx_ldl2load();
    } else {
        kkt = cvx_ldlload((void *)0);
    }
    
    cvx_kktsolver_t *S = kkt->new(&cp, 3, 0, &dims);
    kkt->factor(S, &cp.W, __cvxnil, __cvxnil);

    kkt->solve(S, &x, &y, &z_g);
    cvx_mgrp_elem(&m, &z_g, CVXDIM_SDP, 0);
    cvxm_make_trm(&m, CVX_LOWER);

    cvx_mat_printf(stdout, "%13.6e", &x, "x");
    cvx_mat_printf(stdout, "%13.6e", &z, "z");
    
    compare_data(&x, &y, &z);
    printf("||x||_2: %e\n", cvxm_nrm2(&x));
    printf("||y||_2: %e\n", cvxm_nrm2(&y));
    printf("||z||_2: %e\n", cvxm_nrm2(&z));
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
