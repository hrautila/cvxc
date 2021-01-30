
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include <unistd.h>
#include "cvxc.h"

extern int print_solution(cvxc_solution_t *sol);

/*
 * The analytic centering with cone constraints example in section 9.1 of CVXBook.
 */

static
cvxc_float_t inv(cvxc_float_t v)
{
    return 1.0/v;
}

static
cvxc_float_t inv_square(cvxc_float_t v)
{
    return (1.0/v)*(1.0/v);
}

static
int log_accumulate(cvxc_float_t e, void *p)
{
    cvxc_float_t *accum = (cvxc_float_t *)p;
    if (e <= 0.0)
        return -1;
    *accum += (cvxc_float_t)log((double)e);
    return 0;
}

static int cvxm_logsum(cvxc_float_t *result, const cvxc_matrix_t *x)
{
    cvxc_float_t sum = 0.0;
    if (armas_iterate(&x->data, log_accumulate, &sum, 0) < 0)
        return -1;
    *result = sum;
    return 0;
}

static
int acenter_F(cvxc_matrix_t *f,
              cvxc_matrix_t *Df,
              cvxc_matrix_t *H,
              const cvxc_matrix_t *x,
              const cvxc_matrix_t *z,
              void *user)
{
    cvxc_float_t logsum;
    cvxc_matrix_t d;

    if (!x) {
        cvxm_set_all(f, 1.0);
        return 0;
    }

    logsum = 0.0;
    if (cvxm_logsum(&logsum, x) < 0)
        return -1;
    if (f)
        cvxm_set(f, 0, 0, -logsum);

    // Df = 1 / x
    cvxm_copy(Df, x, 0);
    cvxm_apply(Df, inv, 0);
    cvxm_scale(Df, -1.0, 0);
    // f and Df

    if (!z)
        return 0;

    // Hessian: diag(z[0] * x^-2)
    cvxm_view_diag(&d, H, 0);
    cvxm_copy(&d, x, 0);
    cvxm_apply(&d, inv_square, 0);
    cvxm_scale(&d, cvxm_get(z, 0, 0), 0);

    return 0;
}

int read_mat(cvxc_matrix_t *m, const char *path)
{
    cvxc_stream_t ios;
    FILE *fp = fopen(path, "r");
    if (!fp) {
        perror(path);
        return -1;
    }
    cvxc_file_stream(&ios, fp);
    int err = cvxc_json_matrix_read(&m, &ios);
    fclose(fp);
    return err;
}

int main(int argc, char **argv)
{
    cvxc_matrix_t G, h, A, b, y, r, x;
    cvxc_problem_t cp;
    cvxc_dimset_t dims;
    int opt;
    cvxc_size_t m, n;
    char *Aname = 0, *bname = 0;
    m = 15;
    n = 30;

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

    while ((opt = getopt(argc, argv, "N:A:b:")) != -1) {
        switch (opt) {
        case 'N':
            opts.max_iter = atoi(optarg);
            break;
        case 'A':
            Aname = optarg;
            break;
        case 'b':
            bname = optarg;
            break;
        default:
            break;
        }
    }

    if (argc > 1)
        m = atoi(argv[1]);
    if (argc > 2)
        n = atoi(argv[2]);

    if (Aname && bname) {
        read_mat(&A, Aname);
        read_mat(&b, bname);
        cvxm_size(&m, &n, &A);
    } else {
        cvxm_init(&b, m, 1);
        cvxm_init(&A, m, n);
        cvxm_init(&y, m, 1);
        cvxm_init(&r, n, 1);
        cvxm_init(&x, n, 1);

        cvxm_set_from(&y, armas_normal);
        cvxm_set_from(&A, armas_normal);
        cvxm_set_from(&r, armas_uniform);

        // r = r - A^T * y
        // A = A - (1/y^Ty) * y * r^T
        cvxm_mvmult(1.0, &r, -1.0, &A, &y, CVXC_TRANS);
        cvxc_float_t doty = cvxm_dot(&y, &y);
        cvxm_mvupdate(&A, -(1.0 / doty), &y, &r);

        cvxm_set_from(&x, armas_uniform);
        cvxm_mvmult(0.0, &b, 1.0, &A, &x, 0);
    }

    dims = (cvxc_dimset_t){0};
    dims.iscpt = 1;
    cvxc_convex_program_init(&F, acenter_F, (void *)0);

    cvxm_map_data(&h, 0, 1, (cvxc_float_t *)0);
    cvxm_map_data(&G, 0, n, (cvxc_float_t *)0);

    if (opts.max_iter == 0)
        return 0;

    cvxc_cp_setup(&cp, &F, &G, &h, &A, &b, &dims, (cvxc_kktsolver_t *)0);
    cvxc_cp_compute_start(&cp);
    cvxc_cp_solve(&cp, &opts);
    return print_solution(&cp.solution);
}
