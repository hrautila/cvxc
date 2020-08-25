
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include <stdio.h>
#include "convex.h"


#define sqrt2    1.41421356237309504880
#define invsq2 1.0/1.41421356237309504880


int test_sdot(int N, int verbose)
{
    cvx_matrix_t x, y, xk, yk;
    cvx_matgrp_t x_g, y_g;
    cvx_index_t index;
    cvx_float_t jval = 0.0, exp;
    cvx_dimset_t dims;
    int ok, k;
    
    cvx_dimset_alloc(&dims, N, (int[]){N, 0}, (int[]){N, 0});

    // create matrix and index set from dimensions
    cvx_index_create(&x, &index, &dims, CVX_INDEX_NORMAL);
    cvxm_init(&y, x.rows, x.cols);

    cvx_mgrp_init(&x_g, &x, &index);
    cvx_mgrp_init(&y_g, &y, &index);
    
    // set all x_k = sqrt(2), y_k = 1/sqrt(2)
    for (k = 0; k < x.rows; k++) {
        cvxm_set(&x, k, 0, sqrt2);
        cvxm_set(&y, k, 0, invsq2);
    }

    exp = (cvx_float_t)(2*N + N*N);
    jval = cvx_sdot(&x_g, &y_g);
    ok = abs(jval - exp) < DBL_EPSILON;
    printf("sdot           :  %e == %e: %s\n", jval, exp, ok ? "OK" : "NOT OK");

    // make 'S' matrices lower triangular
    for (k = 0; k < cvx_mgrp_count(&x_g, CVXDIM_SDP); k++) {
        cvx_mgrp_elem(&xk, &x_g, CVXDIM_SDP, k);
        cvx_mgrp_elem(&yk, &y_g, CVXDIM_SDP, k);
        cvxm_make_trm(&xk, CVX_LOWER);
        cvxm_make_trm(&yk, CVX_LOWER);
    }

    jval = cvx_sdot(&x_g, &y_g);
    ok = abs(jval - exp) < DBL_EPSILON;
    printf("sdot tril('S') :  %e == %e: %s\n", jval, exp, ok ? "OK" : "NOT OK");

    cvx_index_release(&index);
    cvxm_release(&x);
    cvxm_release(&y);
    return ok;
}

int test_snrm2(int N, int verbose)
{
    cvx_matrix_t x, xk;
    cvx_matgrp_t x_g;
    cvx_index_t index;
    cvx_float_t jval = 0.0, exp;
    cvx_dimset_t dims;
    int ok, k;
    
    cvx_dimset_alloc(&dims, N, (int[]){N, 0}, (int[]){N, 0});

    // create matrix and index set from dimensions
    cvx_index_create(&x, &index, &dims, CVX_INDEX_NORMAL);
    cvx_mgrp_init(&x_g, &x, &index);
    
    // set all x_k = sqrt(2), y_k = 1/sqrt(2)
    for (k = 0; k < x.rows; k++) {
        cvxm_set(&x, k, 0, sqrt2);
    }

    exp  = (cvx_float_t)sqrt(2*N + N*N)*sqrt2;
    jval = cvx_snrm2(&x_g);
    ok   = abs(jval - exp) < DBL_EPSILON;
    printf("snrm2           :  %e == %e: %s\n", jval, exp, ok ? "OK" : "NOT OK");

    // make 'S' matrices lower triangular
    for (k = 0; k < cvx_mgrp_count(&x_g, CVXDIM_SDP); k++) {
        cvx_mgrp_elem(&xk, &x_g, CVXDIM_SDP, k);
        cvxm_make_trm(&xk, CVX_LOWER);
    }

    jval = cvx_snrm2(&x_g);
    ok   = abs(jval - exp) < DBL_EPSILON;
    printf("snrm2 tril('S') :  %e == %e: %s\n", jval, exp, ok ? "OK" : "NOT OK");

    cvx_index_release(&index);
    cvxm_release(&x);
    return ok;
    //exp = sqrt(sqrt2*sqrt2*slen);
}

int test_trisc(int N, int verbose)
{
    cvx_matrix_t x;
    cvx_matgrp_t x_g;
    cvx_index_t index;
    cvx_float_t asum;
    cvx_dimset_t dims;
    int ok, k;
    
    cvx_dimset_alloc(&dims, N, (int[]){N, 0}, (int[]){N, 0});

    // create matrix and index set from dimensions
    cvx_index_create(&x, &index, &dims, CVX_INDEX_NORMAL);
    cvx_mgrp_init(&x_g, &x, &index);
    
    // set all x_k = sqrt(2), y_k = 1/sqrt(2)
    for (k = 0; k < x.rows; k++) {
        cvxm_set(&x, k, 0, 2.0);
    }


    // expected sum of entries is slen*2.0
    asum = cvxm_asum(&x);
    ok = (2.0*x.rows - asum) < DBL_EPSILON;
    printf("sum(|x_i|)      : exp %f == %f : %s\n", 2.0*x.rows, asum,
           ok ? "OK" : "NOT OK");
    
    cvx_trisc(&x_g);
    // strictly upper part of 'S' matrices are zero; expected sum slen*2.0
    asum = cvxm_asum(&x);
    ok = (2.0*x.rows - asum) < DBL_EPSILON;
    printf("sum(|trisc(x)|) : exp %f == %f : %s\n", 2.0*x.rows, asum,
           ok ? "OK" : "NOT OK");

    // unscale strictly lower part; upper part still zeros 
    cvx_triusc(&x_g);
    // size of strictly upper (lower) part of triangular matrix is N*(N-1)/2
    // expected sum = (slen - N*(N-1)/2)*2.0
    asum = cvxm_asum(&x);
    ok = (2.0*(x.rows-(N*(N-1)/2)) - asum) < DBL_EPSILON;
    printf("sum(|triusc(x)|): exp %f == %f : %s\n", 2.0*(x.rows-(N*(N-1)/2)), asum,
           ok ? "OK" : "NOT OK");
    
    return ok;
}

void mgrp_print(cvx_matgrp_t *g, char *s)
{
    cvx_matrix_t xk;
    
    cvx_mgrp_elem(&xk, g, CVXDIM_LINEAR, 0);
    printf("L(%s)\n", s); cvxm_printf(stdout, "%.2e", &xk);

    cvx_mgrp_elem(&xk, g, CVXDIM_SOCP, 0);
    printf("Q(%s)\n", s); cvxm_printf(stdout, "%.2e", &xk);

    cvx_mgrp_elem(&xk, g, CVXDIM_SDP, 0);    
    printf("S(%s)\n", s); cvxm_printf(stdout, "%.2e", &xk); 
}
/*
 */
int test_sprod(int N, int verbose)
{
    cvx_matrix_t x, y, xk, yk;
    cvx_matgrp_t x_g, y_g;
    cvx_memblk_t wrk;
    cvx_index_t index;
    cvx_float_t exp_sum_l, exp_sum_q, exp_sum_s, exp, asum;
    cvx_dimset_t dims;
    int ok, k;
    
    cvx_dimset_alloc(&dims, N, (int[]){N, 0}, (int[]){N, 0});
    __mblk_init(&wrk, 64);
    
    // create matrix and index set from dimensions
    cvx_index_create(&x, &index, &dims, CVX_INDEX_NORMAL);
    cvxm_init(&y, x.rows, x.cols);

    cvx_mgrp_init(&x_g, &x, &index);
    cvx_mgrp_init(&y_g, &y, &index);
    
    // set all x_k = sqrt(2), y_k = 1/sqrt(2)
    for (k = 0; k < x.rows; k++) {
        cvxm_set(&x, k, 0, sqrt2);
        cvxm_set(&y, k, 0, sqrt2);
    }

    exp_sum_l = N*2.0;
    exp_sum_q = N*2.0 + (N-1)*2.0 + (N-1)*2.0;
    // full matrix; 
    exp_sum_s = N*2.0*(N*(N+1)/2);
    // exp_sum_sd = 2.0*(N*(N+1)/2);
    exp =  exp_sum_l + exp_sum_q + exp_sum_s;

    cvx_sprod(&x_g, &y_g, 0, &wrk);

    asum = cvxm_asum(&x);
    ok = abs(asum - exp) < DBL_EPSILON;
    printf("sprod          :  %e == %e: %s\n", asum, exp, ok ? "OK" : "NOT OK");
    mgrp_print(&x_g, "x");

    // reset X
    for (k = 0; k < x.rows; k++) {
        cvxm_set(&x, k, 0, sqrt2);
    }
    // make 'S' matrices lower triangular
    for (k = 0; k < cvx_mgrp_count(&x_g, CVXDIM_SDP); k++) {
        cvx_mgrp_elem(&xk, &x_g, CVXDIM_SDP, k);
        cvx_mgrp_elem(&yk, &y_g, CVXDIM_SDP, k);
        cvxm_make_trm(&xk, CVX_LOWER);
        cvxm_make_trm(&yk, CVX_LOWER);
    }

    cvx_sprod(&x_g, &y_g, 0, &wrk);
    asum = cvxm_asum(&x);
    ok = abs(asum - exp) < DBL_EPSILON;
    printf("sdot tril('S') :  %e == %e: %s\n", asum, exp, ok ? "OK" : "NOT OK");
    mgrp_print(&x_g, "x");

    __mblk_release(&wrk);
    cvx_index_release(&index);
    cvxm_release(&x);
    cvxm_release(&y);
    return ok;
}

int test_sprod_diag(int N, int verbose)
{
    cvx_matrix_t x, y, xk;
    cvx_matgrp_t x_g, y_g;
    cvx_memblk_t wrk;
    cvx_index_t index, index_d;
    cvx_float_t exp_sum_l, exp_sum_q, exp_sum_s, exp, asum;
    cvx_dimset_t dims;
    int ok, k;
    
    cvx_dimset_alloc(&dims, N, (int[]){N, 0}, (int[]){N, 0});
    __mblk_init(&wrk, 64);
    
    // create matrix and index set from dimensions
    cvx_index_create(&x, &index, &dims, CVX_INDEX_NORMAL);
    cvx_index_create(&y, &index_d, &dims, CVX_INDEX_DIAG);

    cvx_mgrp_init(&x_g, &x, &index);
    cvx_mgrp_init(&y_g, &y, &index_d);
    
    // set all x_k = sqrt(2), y_k = 1/sqrt(2)
    for (k = 0; k < x.rows; k++) {
        cvxm_set(&x, k, 0, sqrt2);
        if (k < y.rows)
            cvxm_set(&y, k, 0, sqrt2);
    }

    exp_sum_l = N*2.0;
    exp_sum_q = N*2.0 + (N-1)*2.0 + (N-1)*2.0;
    // full matrix; 
    exp_sum_s = 2.0*(N*(N+1)/2);
    // exp_sum_sd = 2.0*(N*(N+1)/2);
    exp =  exp_sum_l + exp_sum_q + exp_sum_s;

    cvx_sprod(&x_g, &y_g, CVX_DIAG, &wrk);

    asum = cvxm_asum(&x);
    ok = abs(asum - exp) < DBL_EPSILON;
    printf("sprod          :  %e == %e: %s\n", asum, exp, ok ? "OK" : "NOT OK");
    mgrp_print(&x_g, "x");

    // reset X
    for (k = 0; k < x.rows; k++) {
        cvxm_set(&x, k, 0, sqrt2);
    }
    // make 'S' matrices lower triangular
    for (k = 0; k < cvx_mgrp_count(&x_g, CVXDIM_SDP); k++) {
        cvx_mgrp_elem(&xk, &x_g, CVXDIM_SDP, k);
        cvxm_make_trm(&xk, CVX_LOWER);
    }

    cvx_sprod(&x_g, &y_g, CVX_DIAG, &wrk);
    asum = cvxm_asum(&x);
    ok = abs(asum - exp) < DBL_EPSILON;
    printf("sprod tril('S'):  %e == %e: %s\n", asum, exp, ok ? "OK" : "NOT OK");
    mgrp_print(&x_g, "x");

    __mblk_release(&wrk);
    cvx_index_release(&index);
    cvxm_release(&x);
    cvxm_release(&y);

    return ok;
}


int test_sinv_diag(int N, int verbose)
{
    cvx_matrix_t x, y;
    cvx_matgrp_t x_g, y_g;
    cvx_memblk_t wrk;
    cvx_index_t index, index_d;
    cvx_float_t exp, asum;
    cvx_dimset_t dims;
    int ok, k;
    
    cvx_dimset_alloc(&dims, N, (int[]){N, 0}, (int[]){N, 0});
    __mblk_init(&wrk, 64);
    
    // create matrix and index set from dimensions
    cvx_index_create(&x, &index, &dims, CVX_INDEX_NORMAL);
    cvx_index_create(&x, &index_d, &dims, CVX_INDEX_DIAG);

    cvx_mgrp_init(&x_g, &x, &index);
    cvx_mgrp_init(&y_g, &y, &index_d);
    
    // set all x_k = sqrt(2), y_k = 1/sqrt(2)
    for (k = 0; k < y.rows; k++) {
        cvxm_set(&x, k, 0, sqrt2);
        if (k < y.rows)
            cvxm_set(&y, k, 0, sqrt2);
    }

    // expect original ie. all entries sqrt2
    exp =  (N + N + N*(N+1)/2)*sqrt2;

    cvx_sprod(&x_g, &y_g, CVX_DIAG, &wrk);
    cvx_sinv(&x_g, &y_g, &wrk);
    
    asum = cvxm_asum(&x);
    ok = abs(asum - exp) < DBL_EPSILON;
    printf("sinv            :  %e == %e: %s\n", asum, exp, ok ? "OK" : "NOT OK");
    mgrp_print(&x_g, "x");

    __mblk_release(&wrk);
    cvx_index_release(&index);
    cvxm_release(&x);
    cvxm_release(&y);

    return ok;
}


int test_ssqr_diag(int N, int verbose)
{
    cvx_matrix_t x, y;
    cvx_matgrp_t x_g, y_g;
    cvx_index_t index;
    cvx_float_t exp, asum;
    cvx_dimset_t dims;
    int ok, k;
    
    cvx_dimset_alloc(&dims, N, (int[]){N, 0}, (int[]){N, 0});
    
    // create matrix and index set from dimensions
    cvx_index_create(&x, &index, &dims, CVX_INDEX_DIAG);
    cvxm_init(&y, x.rows, 1);

    cvx_mgrp_init(&x_g, &x, &index);
    cvx_mgrp_init(&y_g, &y, &index);
    
    // set all x_k = sqrt(2), y_k = 1/sqrt(2)
    for (k = 0; k < y.rows; k++) {
        cvxm_set(&y, k, 0, sqrt2);
    }
    mgrp_print(&y_g, "y");

    // expect original ie. all entries sqrt2
    exp =  (N + N + N)*2.0;

    cvx_ssqr(&x_g, &y_g);
    
    asum = cvxm_asum(&x);
    ok = abs(asum - exp) < DBL_EPSILON;
    printf("ssqr            :  %e == %e: %s\n", asum, exp, ok ? "OK" : "NOT OK");
    mgrp_print(&x_g, "x");

    cvx_index_release(&index);
    cvxm_release(&x);
    cvxm_release(&y);

    return ok;
}

int main(int argc, char **argv)
{
    int verbose = 1;
    
    //test_sdot(4, verbose);
    //test_snrm2(4, verbose);
    //test_trisc(4, verbose);
    //test_sprod(4, verbose);
    //test_sprod_diag(4, verbose);
    //test_sinv_diag(4, verbose);
    test_ssqr_diag(4, verbose);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
