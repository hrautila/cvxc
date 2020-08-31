
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include <unistd.h>
#include "cvxc.h"

extern int print_solution(cvxc_solution_t *sol);

/*
 * The analytic centering with cone constraints example in section 9.1 of CVXBook.
 */
typedef struct acenter {
    int rows;
    int cols;
} acenter_t;

int acenter_F(cvxc_matrix_t *f,
              cvxc_matrix_t *Df,
              cvxc_matrix_t *H,
              const cvxc_matrix_t *x,
              const cvxc_matrix_t *z,
              void *user)
{
    cvxc_size_t r, c;
    acenter_t *p = (acenter_t *)user;
    if (!x && !z) {
        cvxm_size(&r, &c, f);
        if ((int)r != p->rows || (int)c != p->cols)
            return -1;
        cvxm_set_all(f, 0.0);
        return 0;
    }
    if (cvxm_amax(x) >= 1.0)
        return -1;

    // compute u := 1 - x**2 ; f = u
    cvxc_float_t z0, xv, u, usum = 0.0;
    if (z) {
        z0 = cvxm_get(z, 0, 0);
        cvxm_set_all(H, 0.0);
        //printf("acenter z0: %.7f\n", z0);
    }

    //  f = [1, 1]; Df = [1, 3]; H = [3, 3]
    //
    //  f = - sum( log(1 - x_i^2) )
    // Df = 2*x*u = 2*x/(1 - x**2)
    for (int i = 0; i < p->rows; i++) {
        xv = cvxm_get(x, i, 0);
        u = 1.0 - xv*xv;
        cvxm_set(Df, 0, i, 2.0*xv/u);
        if (z) {
            // H[i,i] =  2*z0* (1 + u**2)/u**2 == 2*z0*(1 + 1/u**2)
            cvxm_set(H, i, i, 2*z0*((1.0 + u*u)/(u*u)));
        }
        usum += log(u);
    }
    if (f)
        cvxm_set(f, 0, 0, -usum);

    return 0;
 }



int main(int argc, char **argv)
{
    cvxc_matrix_t G, h, A, b;
    cvxc_problem_t cp;
    cvxc_dimset_t dims;
    int opt;

    // 39x3
    cvxc_float_t gdata[39] = {
        0., -1.,  0.,  0., -21., -11.,   0., -11.,  10.,   8.,   0.,   8., 5.,
        0.,  0., -1.,  0.,   0.,  10.,  16.,  10., -10., -10.,  16., -10., 3.,
        0.,  0.,  0., -1.,  -5.,   2., -17.,   2.,  -6.,   8., -17.,  -7., 6.
    };
    cvxc_float_t hdata[13] = {
        1.0, 0.0, 0.0, 0.0, 20., 10., 40., 10., 80., 10., 40., 10., 15.
    };

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
    cvxc_convex_program_t F;
    acenter_t Acenter = (acenter_t){ .rows = 3, .cols = 1};

    while ((opt = getopt(argc, argv, "N:")) != -1) {
        switch (opt) {
        case 'N':
            opts.max_iter = atoi(optarg);
            break;
        default:
            break;
        }
    }

    cvxc_dimset_alloc(&dims, 0, (cvxc_size_t []){4, 0}, (cvxc_size_t []){3, 0});
    dims.iscpt = 1;
    cvxc_convex_program_init(&F, acenter_F, &Acenter);

    cvxm_map_data(&h, 13, 1, hdata);
    cvxm_map_data(&G, 13, 3, gdata);;
    cvxm_map_data(&A, 0, 3, (cvxc_float_t *)0);
    cvxm_map_data(&b, 0, 1,  (cvxc_float_t *)0);

    if (opts.max_iter == 0)
        return 0;

    cvxc_cp_setup(&cp, &F, &G, &h, &A, &b, &dims, (cvxc_kktsolver_t *)0);
    cvxc_cp_compute_start(&cp);
    cvxc_cp_solve(&cp, &opts);
    return print_solution(&cp.solution);
}
