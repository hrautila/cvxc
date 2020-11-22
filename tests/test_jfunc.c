
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include <stdio.h>
#include "convex.h"



int test_jdot(int N, int verbose)
{
    cvxc_matrix_t x, y;
    cvxc_float_t jval = 0.0, exp;
    int ok, k, dot = 0;
    
    cvxm_init(&x, N, 1);
    cvxm_init(&y, N, 1);
    for (k = 1; k < N; k++) {
        cvxm_set(&x, k, 0, N-k);
        cvxm_set(&y, k, 0, N-k);
        dot += (N-k)*(N-k);
    }
    exp = (cvxc_float_t)dot;
    cvxm_set(&x, 0, 0, SQRT(exp));
    cvxm_set(&y, 0, 0, SQRT(exp));

    jval = cvxc_jdot(&x, &y);
    ok = abs(jval) < DBL_EPSILON;
    printf("jdot :  %e == %e: %s\n", jval, 0.0, ok ? "OK" : "NOT OK");
    return ok;
}


/*
 * Compute sqrt(x^T * J * x), where J = [1, 0; 0, -I]
 *
 *   (x0 X1^T) (1  0) (x0)  = (x0 -X1^T) (x0)  = sqrt(x0*x0 - X1^T*X1) = 
 *             (0 -I) (X1)               (X1)
 *
 *   sqrt( x0*x0 - |X1|*|X1| ) = sqrt(x0 - |X1|)*sqrt(x0 + |X1|)
 */
int test_jnrm2(int N, int verbose)
{
    cvxc_matrix_t x, x1;
    cvxc_float_t exp, jval = 0.0, abs_x1;
    int ok, k;
    
    cvxm_init(&x, N, 1);
    for (k = 0; k < N; k++) {
        cvxm_set(&x, k, 0, N-k);
    }
    cvxm_view_map(&x1, &x, 1, 0, N-1, 1);
    abs_x1 = cvxm_nrm2(&x1);
    // set x[0] = |X1|
    cvxm_set(&x, 0, 0, abs_x1);
    jval = cvxc_jnrm2(&x);
    exp = 0.0;
    ok = abs(jval - 0.0) < DBL_EPSILON;
    printf("jnrm2:  %e == %e: %s\n", jval, exp, ok ? "OK" : "NOT OK");

    return ok;
}

#define sqrt2    1.41421356237309504880
#define invsq2 1.0/1.41421356237309504880

/*
 */
int test_pack(int N, int verbose)
{
    cvxc_float_t data[] = {
        1.0,    invsq2, invsq2, invsq2, invsq2,
        invsq2, 1.0,    invsq2, invsq2, invsq2,
        invsq2, invsq2, 1.0,    invsq2, invsq2,
        invsq2, invsq2, invsq2, 1.0,    invsq2,
        invsq2, invsq2, invsq2, invsq2, 1.0
    };
    cvxc_matrix_t x, xc, y, x0, y0;
    cvxc_dimset_t dims;
    cvxc_index_t index;
    cvxc_float_t res;

    // one 5x5 symmetric matrix
    cvxc_dimset_alloc(&dims, 0, (int *)0, (int[]){5, 0});
    cvxc_index_init(&index, &dims, CVXC_INDEX_NORMAL);

    cvxm_map_data(&x, 5, 5, data);
    cvxm_make_trm(&x, CVXC_LOWER);
    cvxm_init(&y,  5, 5);
    cvxm_init(&xc, 5*3, 1);   // packed size is N*(N+1)/2

    //printf("x\n"); cvxm_printf(stdout, "%6.2f", &x);
    cvxc_pack(&xc, &x, &index);
    //printf("xc\n"); cvxm_printf(stdout, "%6.2f", &xc);

    // result: packed storage of ones
    res = cvxm_asum(&xc);
    int ok = abs(res - 1.0*15) < DBL_EPSILON;
    printf("pack:  %s\n", ok ? "OK" : "NOT OK");
    cvxc_unpack(&y, &xc, &index);
    //printf("y\n"); cvxm_printf(stdout, "%6.2f", &y);

    // compute ||y -x||; map matrices to vectors; substract
    cvxm_map_data(&x0, 25, 1, cvxm_data(&x, 0));
    cvxm_map_data(&y0, 25, 1, cvxm_data(&y, 0));
    cvxm_axpy(&y0, -1.0, &x0);
    res = cvxm_nrm2(&y0);
    printf("  ||unpack(pack(x)) - x||: %e\n", res);
    ok = abs(res) < DBL_EPSILON;
    printf("unpack:  %s\n", ok ? "OK" : "NOT OK");

    cvxc_dimset_release(&dims);
    cvxc_index_release(&index);
    return ok;
}

extern int cvxm_read_sbuffer(cvxc_matrix_t *, const char *);

void test_read()
{
    const char *s = "{2 2 [1.0 2.0 3.0 4.0]}";
    cvxc_matrix_t m;
    cvxm_read_sbuffer(&m, s);
    printf("m\n");
    cvxm_printf(stdout, "%6.2f", &m);
}

int main(int argc, char **argv)
{
    int verbose = 1;

    test_jdot(7, verbose);
    test_jnrm2(7, verbose);
    test_read();
    test_pack(7, verbose);
}
