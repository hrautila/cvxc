
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include <stdio.h>
#include <assert.h>
#include "convex.h"

#define SQRT2    1.41421356237309504880
#define INVSQRT2 1.0/1.41421356237309504880

static inline
int nEPS(cvxc_float_t x, int n)
{
    return x < n*DBL_EPSILON;
}

cvxc_float_t sqrt2(cvxc_float_t v)
{
    return SQRT2;
}

cvxc_float_t invsqrt2(cvxc_float_t v)
{
    return INVSQRT2;
}

cvxc_float_t two(cvxc_float_t v)
{
    return 2.0;
}

void mgrp_print(cvxc_matgrp_t *g, char *s)
{
    cvxc_matrix_t xk;
    
    if (cvxc_mgrp_elem(__cvxnil, g, CVXDIM_LINEAR, 0) > 0) {
        cvxc_mgrp_elem(&xk, g, CVXDIM_LINEAR, 0);
        printf("L(%s)\n", s); cvxm_printf(stdout, "%9.2e", &xk);
    }

    if (cvxc_mgrp_count(g, CVXDIM_SOCP) > 0) {
        cvxc_mgrp_elem(&xk, g, CVXDIM_SOCP, 0);
        printf("Q(%s)\n", s); cvxm_printf(stdout, "%9.2e", &xk);
    }

    if (cvxc_mgrp_count(g, CVXDIM_SDP) > 0) {
        cvxc_mgrp_elem(&xk, g, CVXDIM_SDP, 0);    
        printf("S(%s)\n", s); cvxm_printf(stdout, "%9.2e", &xk);
    }
}

void test_compute_linear(int N, int verbose)
{
    cvxc_matrix_t s, z, lmbda, d, di;
    cvxc_matgrp_t s_g, z_g, lmbda_g;
    cvxc_scaling_t W;
    cvxc_index_t index;
    cvxc_dimset_t dims;

    
    cvxc_dimset_alloc(&dims, N, (int *)0, (int *)0);
    cvxc_index_create(&s, &index, &dims, CVXC_INDEX_DIAG);
    cvxc_scaling_alloc(&W, &dims);
    
    cvxm_init(&z, s.rows, 1);
    cvxm_init(&lmbda, s.rows, 1);

    cvxc_mgrp_init(&s_g, &s, &index);
    cvxc_mgrp_init(&z_g, &z, &index);
    cvxc_mgrp_init(&lmbda_g, &lmbda, &index);

    cvxm_apply(&s, sqrt2, CVXC_ALL);
    cvxm_apply(&z, invsqrt2, CVXC_ALL);
    cvxm_apply(&lmbda, sqrt2, CVXC_ALL);

    mgrp_print(&s_g, "s");
    mgrp_print(&z_g, "z");
    mgrp_print(&lmbda_g, "lmbda");

    assert(cvxm_get(&s, 0, 0) == SQRT2);
    assert(cvxm_get(&z, 0, 0) == INVSQRT2);

    cvxc_compute_scaling(&W, &s_g, &z_g, &lmbda_g, (cvxc_memblk_t *)0);
    
    cvxc_scaling_elem(&d, &W, CVXWS_D, 0);
    cvxc_scaling_elem(&di, &W, CVXWS_DI, 0);

    // d = sqrt(s ./ z) == sqrt(2)
    printf("D: "); cvxm_printf(stdout, "%6.2f", &d);
    // di = 1/d == 1.0/sqrt(2)
    printf("DI:"); cvxm_printf(stdout, "%6.2f", &di);
    assert(cvxm_get(&d, 0, 0) == SQRT2);
    assert(cvxm_get(&di, 0, 0) == INVSQRT2);

    printf("l: "); cvxm_printf(stdout, "%6.2f", &lmbda);
    assert(cvxm_get(&lmbda, 0, 0) == 1.0);
    
    // test scaling
    cvxm_apply(&s, two, CVXC_ALL);
    cvxm_copy(&z, &s, CVXC_ALL);
    
    cvxc_memblk_t work;
    __mblk_init(&work, 4*N);
    
    cvxc_scale(&s_g, &W, 0, &work);
    printf("W*s\n"); cvxm_printf(stdout, "%6.2f", &s);
    
    cvxc_scale(&s_g, &W, CVXC_INV, &work);
    printf("W^-1*s\n"); cvxm_printf(stdout, "%6.2f", &s);

    cvxc_scale(&s_g, &W, CVXC_TRANS, &work);
    printf("W^T*s\n"); cvxm_printf(stdout, "%6.2f", &s);
    
    cvxc_scale(&s_g, &W, CVXC_TRANS|CVXC_INV, &work);
    printf("W^-T*s\n"); cvxm_printf(stdout, "%6.2f", &s);
}

void test_compute_socp(int N, int verbose)
{
    cvxc_matrix_t s, z, ds, dz, lmbda, v, beta;
    cvxc_matgrp_t s_g, ds_g, z_g, dz_g, lmbda_g;
    cvxc_scaling_t W;
    cvxc_index_t index;
    cvxc_dimset_t dims;

    
    cvxc_dimset_alloc(&dims, 0, (int []){N, 0}, (int *)0);
    cvxc_index_create(&s, &index, &dims, CVXC_INDEX_DIAG);
    cvxm_init(&ds, s.rows, 1);
    cvxm_init(&z, s.rows, 1);
    cvxm_init(&dz, s.rows, 1);
    cvxc_scaling_alloc(&W, &dims);
    
    cvxm_init(&lmbda, s.rows, 1);

    cvxc_mgrp_init(&s_g, &s, &index);
    cvxc_mgrp_init(&z_g, &z, &index);
    cvxc_mgrp_init(&ds_g, &ds, &index);
    cvxc_mgrp_init(&dz_g, &dz, &index);
    cvxc_mgrp_init(&lmbda_g, &lmbda, &index);

    cvxm_apply(&s, sqrt2, CVXC_ALL);
    cvxm_apply(&z, invsqrt2, CVXC_ALL);
    cvxm_apply(&lmbda, sqrt2, CVXC_ALL);

    cvxm_set(&s, 0, 0, N*SQRT2); 
    //cvxm_set(&s, 1, 0, 2*N*SQRT2); 
    //cvxm_set(&s, 2, 0, 1*N*SQRT2); 
    cvxm_set(&z, 0, 0, N*INVSQRT2); 
    //cvxm_set(&z, 1, 0, 2*N*INVSQRT2); 
    //cvxm_set(&z, 2, 0, 1*N*INVSQRT2); 

    // ds = s/10; dz = z/10
    cvxm_copy(&ds, &s, CVXC_ALL);
    cvxm_scale(&ds, 1/10.0, CVXC_ALL);
    cvxm_copy(&dz, &z, CVXC_ALL);
    cvxm_scale(&dz, 1/10.0, CVXC_ALL);

    cvxc_compute_scaling(&W, &s_g, &z_g, &lmbda_g, (cvxc_memblk_t *)0);
    
    cvxc_scaling_elem(&beta, &W, CVXWS_BETA, 0);
    cvxc_scaling_elem(&v, &W, CVXWS_V, 0);


    printf("V   : "); cvxm_printf(stdout, "%6.2f", &v);
    printf("BETA: "); cvxm_printf(stdout, "%6.2f", &beta);
    mgrp_print(&s_g, "s");
    mgrp_print(&z_g, "z");
    mgrp_print(&lmbda_g, "lmbda");

    // sqrt(v'*J*v) == 1, beta > 0.0
    //printf("|v|_j: %f\n", cvxc_jnrm2(&v));

    // beta * H * z = lambda; beta * H \ s = lambda
    // H = 2*v*v' - J;
    //  (2*v*v' - J)*z = 2*v*v'*z - J*z = 2*w*[v0; v1] - [z0, -z1]
    //  

#if 0
    cvxm_init(&x, s.rows, 1);
    cvxc_float_t w = 2*cvxm_dot(&v, &z);
    cvxm_copy(&x, &z, CVXC_ALL);
    // x = - [z0 , -z1]
    cvxm_set(&x, 0, 0, - cvxm_get(&z, 0, 0));
    // x = - [z0, -z1] + 2*w*v
    cvxm_axpy(&x, &v, w);
    // x = beta*x == lmbda
    cvxm_scale(&x, cvxm_get(&beta, 0, 0), CVXC_ALL);       
    printf("x: "); cvxm_printf(stdout, "%6.2f", &x);
    cvxm_axpy(&x, &lmbda, -1.0);
    printf("||beta*H*z - lmbda||: %f\n", cvxm_nrm2(&x));
#endif
    //
    //cvxm_apply(&s, two, CVXC_ALL);
    //cvxm_set(&s, 0, 0, N*2.0);
    //cvxm_copy(&z, &s, CVXC_ALL);
    
    cvxc_memblk_t work;
    __mblk_init(&work, 4*N);
#if 0
    printf("s\n"); cvxm_printf(stdout, "%6.2f", &s);
    
    cvxc_scale(&s_g, &W, 0, &work);
    printf("W*s\n"); cvxm_printf(stdout, "%6.2f", &s);
    
    cvxc_scale(&s_g, &W, CVXC_INV, &work);
    printf("W^-1*s\n"); cvxm_printf(stdout, "%6.2f", &s);

    cvxc_scale(&s_g, &W, CVXC_TRANS, &work);
    printf("W^T*s\n"); cvxm_printf(stdout, "%6.2f", &s);
    
    cvxc_scale(&s_g, &W, CVXC_TRANS|CVXC_INV, &work);
    printf("W^-T*s\n"); cvxm_printf(stdout, "%6.2f", &s);
    
    cvxc_scale2(&s_g, &lmbda_g, 0, &work);
    printf("H(lmbda^{1/2})*sk\n"); cvxm_printf(stdout, "%6.2f", &s);

    cvxc_scale2(&s_g, &lmbda_g, CVXC_INV, &work);
    printf("H(lmbda^{-1/2})*sk\n"); cvxm_printf(stdout, "%6.2f", &s);
#endif

    printf("scale s, z ..\n");
    cvxc_scale(&s_g, &W, CVXC_INV|CVXC_TRANS, &work);
    cvxc_scale(&z_g, &W, 0, &work);
    mgrp_print(&s_g, "s");
    mgrp_print(&z_g, "z");

    //cc=1.067490, vs=1.067490, vz=1.067490, vq=1.000000, vu=0.000000, wk0=1.000000, dd=0.000000
    printf("update scaling ..\n");
    cvxc_update_scaling(&W,  &lmbda_g, &s_g, &z_g,&work);
    printf("V    : "); cvxm_printf(stdout, "%6.2f", &v);
    printf("BETA : "); cvxm_printf(stdout, "%6.2f", &beta);
    printf("lmbda: "); cvxm_printf(stdout, "%6.2f", &lmbda);
    
#if 0
    // beta * H *(z + dz) == lmbda
    //cvxm_axpy(&z, &dz, 1.0);
    w = 2*cvxm_dot(&v, &dz);
    cvxm_copy(&x, &dz, CVXC_ALL);
    // x = - [z0 , -z1]
    cvxm_set(&x, 0, 0, - cvxm_get(&dz, 0, 0));
    // x = - [z0, -z1] + 2*w*v
    cvxm_axpy(&x, &v, w);
    // x = beta*x == lmbda
    cvxm_scale(&x, cvxm_get(&beta, 0, 0), CVXC_ALL);       
    printf("x: "); cvxm_printf(stdout, "%6.2f", &x);
    cvxm_axpy(&x, &lmbda, -1.0);
    printf("||beta*H*z - lmbda||: %f\n", cvxm_nrm2(&x));
#endif
}

int check_scaling(cvxc_scaling_t *W, cvxc_matgrp_t *s_g, cvxc_memblk_t *wrk)
{
    cvxc_matrix_t st, *s;
    cvxc_float_t nrm;
    s = s_g->mat;

    cvxm_init(&st, s->rows, s->cols);
    cvxm_copy(&st, s, CVXC_ALL);
    
    cvxc_scale(s_g, W, 0, wrk);
    cvxc_scale(s_g, W, CVXC_INV, wrk);
    cvxm_axpy(s, -1.0, &st);
    cvxm_norm(&nrm, s, 0);
    cvxm_copy(s, &st, CVXC_ALL);
    printf("           ||W^-1*W*s - s|| : [%8.2e]\n", nrm);

    cvxc_scale(s_g, W, CVXC_TRANS, wrk);
    cvxc_scale(s_g, W, CVXC_INV|CVXC_TRANS, wrk);
    cvxm_axpy(s, -1.0, &st);
    cvxm_norm(&nrm, s, 0);
    cvxm_copy(s, &st, CVXC_ALL);
    printf("         ||W^-T*W^T*s - s|| : [%8.2e]\n", nrm);
    
    cvxm_release(&st);
    return 1;
}

//
// check
// ||rti*r - I|| ~ 0
// ||W^-T*s - lmbda|| ~ 0
// ||W*z - lmbda|| ~ 0

int check_s(cvxc_scaling_t *W, cvxc_matgrp_t *s_g, cvxc_matgrp_t *z_g, cvxc_matgrp_t *lmbda_g, int k, int verbose)
{
    cvxc_matrix_t s, z, l, r, rti, lk, T, sD;
    cvxc_size_t m;
    cvxc_float_t nrm;
    cvxc_memblk_t wrk;
    int ok;
    
    m = cvxc_mgrp_elem(&s, s_g, CVXDIM_SDP, k);
    cvxc_mgrp_elem(&z, z_g, CVXDIM_SDP, k);
    cvxc_mgrp_elem(&l, lmbda_g, CVXDIM_SDP, k);

    __mblk_init(&wrk, m*m);
    cvxm_init(&lk, m, 1);

    cvxc_scale(s_g, W, CVXC_INV|CVXC_TRANS, &wrk);
    cvxm_make_trm(&s, CVXC_LOWER);
    cvxm_view_diag(&sD, &s, 0);
    cvxm_copy(&lk, &sD, CVXC_ALL);
    cvxm_axpy(&sD, -1.0, &l);
    nrm = cvxm_nrm2(&sD);
    ok = nEPS(nrm, 100);
    printf("S%d ||diag(W^-T*s) - lmbda|| : [%8.2e] %s\n", k, nrm, (ok ? "OK" : "NOT OK"));
    // ||W.^T*s - lmbda||_inf
    cvxm_norm(&nrm, &s, -1);
    printf("S%d       ||W^-T*s - lmbda|| : [%8.2e] \n", k, nrm);
    // restore diagonal
    cvxm_copy(&sD, &lk, CVXC_ALL);
    
    cvxc_scale(z_g, W, 0, &wrk);
    cvxm_make_trm(&z, CVXC_LOWER);
    cvxm_view_diag(&sD, &z, 0);
    cvxm_copy(&lk, &sD, CVXC_ALL);
    cvxm_axpy(&sD, -1.0, &l);
    nrm = cvxm_nrm2(&sD);
    ok = nEPS(nrm, 120);
    printf("S%d    ||diag(W*z) - lmbda|| : [%8.2e] %s\n", k, nrm, (ok ? "OK" : "NOT OK"));
    cvxm_norm(&nrm, &z, -1);
    printf("S%d          ||W*z - lmbda|| : [%8.2e] \n", k, nrm);
    // restore diagonal
    cvxm_copy(&sD, &lk, CVXC_ALL);

    // compute ||r^-1*r - I||_inf
    cvxc_scaling_elem(&r, W, CVXWS_R, 0);
    cvxc_scaling_elem(&rti, W, CVXWS_RTI, 0);
    cvxm_init(&T, m, m);
    cvxm_mult(0.0, &T, 1.0, &rti, &r, CVXC_TRANS);
    if (verbose) {
        printf("rti^T*r\n"); cvxm_printf(stdout, "%9.2e", &T);
    }
    cvxm_view_diag(&sD, &T, 0);
    cvxm_add(&sD, -1.0, CVXC_ALL);
    cvxm_norm(&nrm, &T, -1);
    ok = nEPS(nrm, 100);
    printf("S%d          ||rti^T*r - I|| : [%8.2e] %s\n", k, nrm, (ok ? "OK" : "NOT OK"));
    
    __mblk_release(&wrk);
    cvxm_release(&lk);
    cvxm_release(&T);
    return 1;
}

void test_compute_s(int N, int verbose)
{
    cvxc_matrix_t s, z, lmbda, lmbda2, sk, zk;
    cvxc_matgrp_t s_g, z_g, lmbda_g;
    cvxc_scaling_t W;
    cvxc_index_t index, index_d;
    cvxc_dimset_t dims;
    cvxc_matrix_t L1, L2, sD, wrk;
    cvxc_memblk_t work;
    cvxc_size_t wlen;
    cvxc_float_t d;
    int ok;
    
    cvxc_float_t L1data[7*7] = {
        1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
        0.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0,
        0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0
    };
    cvxc_float_t L2data[7*7] = {
        1.0, -2.0,  2.0, -2.0,  2.0, -2.0,  2.0,
        0.0,  1.0, -2.0,  2.0, -2.0,  2.0, -2.0,
        0.0,  0.0,  1.0, -2.0,  2.0, -2.0,  2.0,
        0.0,  0.0,  0.0,  1.0, -2.0,  2.0, -2.0,
        0.0,  0.0,  0.0,  0.0,  1.0, -2.0,  2.0,
        0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -2.0,
        0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0
    };
    
    if (N != 7)
        N = 7;
    
    cvxc_dimset_alloc(&dims, 0, (int *)0, (int[]){N, 0});
    cvxc_scaling_alloc(&W, &dims);

    cvxc_index_create(&s, &index, &dims, CVXC_INDEX_NORMAL);
    cvxc_index_create(&lmbda, &index_d, &dims, CVXC_INDEX_DIAG);

    cvxm_init(&z, s.rows, s.cols);
    cvxm_map_data(&L1, 7, 7, L1data);
    cvxm_map_data(&L2, 7, 7, L2data);

    cvxc_mgrp_init(&s_g, &s, &index);
    cvxc_mgrp_init(&z_g, &z, &index);
    cvxc_mgrp_init(&lmbda_g, &lmbda, &index_d);

    __mblk_init(&work, 4*N*N);
    
    // lambda = sqrt(2)
    cvxm_apply(&lmbda, sqrt2, CVXC_ALL);
    
    // s = L1*L1.T => chol(s) == L1
    cvxc_mgrp_elem(&sk, &s_g, CVXDIM_SDP, 0);
    cvxm_mult(0.0, &sk, 1.0, &L1, &L1, CVXC_TRANSB);
    // z = L2*L2.T => chol(z) == L2
    cvxc_mgrp_elem(&zk, &z_g, CVXDIM_SDP, 0);
    cvxm_mult(0.0, &zk, 1.0, &L2, &L2, CVXC_TRANSB);

    printf("compute scaling for s=L1*L1^T and z=L2*L2^T ...\n");
    mgrp_print(&s_g, "s");
    cvxc_compute_scaling(&W, &s_g, &z_g, &lmbda_g, &work);
    printf("lmbda:\n"); cvxm_printf(stdout, "%9.2e", &lmbda);
    check_s(&W, &s_g, &z_g, &lmbda_g, 0, 0);
    check_scaling(&W, &s_g, &work);
    mgrp_print(&s_g, "s");
    //cvxm_cholfactor(&sk, CVXC_LOWER);
    //cvxm_cholfactor(&zk, CVXC_LOWER);
    //mgrp_print(&s_g, "chol(W^-T*s)");

    //printf("lmbda\n"); cvxm_printf(stdout, "%6.2f", &lmbda);
    
    printf("\nupdate scaling for s=L1, z=L2 ...\n");

    cvxm_copy(&sk, &L1, CVXC_ALL);
    cvxm_copy(&zk, &L2, CVXC_ALL);
    cvxm_scale(&sk, SQRT2, 0);
    cvxm_scale(&zk, SQRT2, 0);
    
    cvxm_init(&lmbda2, lmbda.rows, lmbda.cols);
    cvxm_copy(&lmbda2, &lmbda, CVXC_ALL);
    cvxm_make_trm(&sk, CVXC_LOWER);
    cvxm_make_trm(&zk, CVXC_LOWER);

    mgrp_print(&s_g, "1.05*s");
    
    cvxc_update_scaling(&W, &lmbda_g, &s_g, &z_g, &work);
    printf("lmbda.new:\n"); cvxm_printf(stdout, "%9.2e", &lmbda);
    mgrp_print(&s_g, "s");
    //mgrp_print(&z_g, "z");
    check_s(&W, &s_g, &z_g, &lmbda_g, 0, 1);
    check_scaling(&W, &s_g, &work);

#if 0
    cvxm_axpy(&lmbda2, -1.0, &lmbda);
    d = cvxm_nrm2(&lmbda2);
    ok = nEPS(d, 100);
    printf("expect ||lmbda.new - lmbda|| ~ 0 : [%8.2e]\n", d);
    //printf("lmbda: "); cvxm_printf(stdout, "%6.2f", &lmbda);
#endif
}

int main(int argc, char **argv)
{
    //test_compute_linear(7, 0);
    //test_compute_socp(7, 0);
    test_compute_s(7, 0);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
