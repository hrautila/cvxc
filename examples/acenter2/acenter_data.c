
// Copyright: Harri Rautila, 2021 <harri.rautila@gmail.com>

#include <stdio.h>
#include <unistd.h>
#include "cvxc.h"


int main(int argc, char **argv)
{
    cvxc_matrix_t G, h, A, b, y, r, x;
    cvxc_params_t params;
    int opt;
    cvxc_size_t m, n;
    char *name = 0;
    m = 15;
    n = 30;

    while ((opt = getopt(argc, argv, "m:n:")) != -1) {
        switch (opt) {
        case 'm':
            m = atoi(optarg);
            break;
        case 'n':
            n = atoi(optarg);
            break;
        default:
            break;
        }
    }

    if (argc > optind)
        name = argv[optind];

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

    cvxm_map_data(&h, 0, 1, (cvxc_float_t *)0);
    cvxm_map_data(&G, 0, n, (cvxc_float_t *)0);

    params = (cvxc_params_t){
        .A = &A,
        .b = &b,
        .dims = 0,
        .G = 0,
        .h = 0,
        .c = 0,
        .module = "acenter_eq",
        .args = 0,
        .opts = 0
    };

    FILE *fp;
    if (name)
        fp = fopen(name, "w");
    else
        fp = stdout;

    if (fp) {
        cvxc_stream_t ios;
        cvxc_file_stream(&ios, fp);
        cvxc_json_write_params(&ios, &params);
        fclose(fp);
    } else {
        perror(name);
    }

    cvxm_release(&b);
    cvxm_release(&A);
    cvxm_release(&y);
    cvxm_release(&r);
    cvxm_release(&x);
    return 0;
}
