/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */
#if ! defined(__STDC_VERSION__)
#  if defined(__CVXC_INLINE)
#    undef __CVXC_INLINE
#  endif
#  define __CVXC_INLINE
#endif

#include <ctype.h>
#include <errno.h>

#include "cvxc.h"

#if __STDC_VERSION__
// inline function reference
void cvxm_init(cvxc_matrix_t *m, cvxc_size_t r, cvxc_size_t c);
void cvxm_release(cvxc_matrix_t *X);
void cvxm_free(cvxc_matrix_t *X);
size_t cvxm_sizeof();

int cvxm_isepi(const cvxc_matrix_t *m);
cvxc_float_t cvxm_get_epi(const cvxc_matrix_t *m);
void cvxm_set_epi(cvxc_matrix_t *m, cvxc_float_t v);

cvxc_matrix_t *cvxm_new(cvxc_size_t r, cvxc_size_t c);
cvxc_matrix_t *cvxm_newcopy(const cvxc_matrix_t *A);
cvxc_float_t *cvxm_data(const cvxc_matrix_t *A, cvxc_size_t k);

cvxc_size_t cvxm_make(cvxc_matrix_t *A, cvxc_size_t rows, cvxc_size_t cols,  void *mem, cvxc_size_t nbytes);
cvxc_size_t cvxm_make_epi(cvxc_matrix_t *A, cvxc_size_t rows, cvxc_size_t cols,  void *mem, cvxc_size_t nbytes);
cvxc_matrix_t *cvxm_map_data(cvxc_matrix_t *A, cvxc_size_t rows, cvxc_size_t cols,  cvxc_float_t *data);

cvxc_matrix_t *cvxm_view_map(cvxc_matrix_t *A, const cvxc_matrix_t *B, cvxc_size_t r, cvxc_size_t c, cvxc_size_t nr, cvxc_size_t nc);
void cvxm_view_unmap(cvxc_matrix_t *A, cvxc_matrix_t *B);
cvxc_matrix_t *cvxm_view_diag(cvxc_matrix_t *A, const cvxc_matrix_t *B, int nd);
void cvxm_size(size_t *r, size_t *c, const cvxc_matrix_t *A);
cvxc_float_t cvxm_get(const cvxc_matrix_t *A, cvxc_size_t r, cvxc_size_t c);
void cvxm_set(cvxc_matrix_t *A, cvxc_size_t r, cvxc_size_t c, cvxc_float_t val);
void cvxm_apply(cvxc_matrix_t *A, cvxm_oper_t f, int flags);

cvxc_float_t cvxm_dot(const cvxc_matrix_t *X, const cvxc_matrix_t *Y);
cvxc_float_t cvxm_nrm2(const cvxc_matrix_t *X);
cvxc_float_t cvxm_amax(const cvxc_matrix_t *X);
cvxc_float_t cvxm_asum(const cvxc_matrix_t *X);
int cvxm_scale(cvxc_matrix_t *X, cvxc_float_t alpha, int flags);
int cvxm_add(cvxc_matrix_t *X, cvxc_float_t alpha, int flags);
void cvxm_make_trm(cvxc_matrix_t *X, int flags);
void cvxm_copy(cvxc_matrix_t *X, const cvxc_matrix_t *Y, int flags);

int cvxm_axpy(cvxc_matrix_t *Y, cvxc_float_t alpha, const cvxc_matrix_t *X);
int cvxm_axpby(cvxc_float_t beta, cvxc_matrix_t *Y, cvxc_float_t alpha, const cvxc_matrix_t *X);
cvxc_float_t cvxc_get(cvxc_matrix_t *A, cvxc_size_t r, cvxc_size_t c);
cvxc_float_t cvxc_set(cvxc_matrix_t *A, cvxc_size_t r, cvxc_size_t c, cvxc_float_t val);
int cvxm_mvsolve_trm(cvxc_matrix_t *X, cvxc_float_t alpha, const cvxc_matrix_t *A, int flags);
int cvxm_mvupdate(cvxc_matrix_t *A, cvxc_float_t alpha, const cvxc_matrix_t *X, const cvxc_matrix_t *Y);
int cvxm_solve_diag(cvxc_matrix_t *X, cvxc_float_t alpha, const cvxc_matrix_t *A, int flags);
int cvxm_mult_diag(cvxc_matrix_t *X, cvxc_float_t alpha, const cvxc_matrix_t *A, int flags);
int cvxm_solve_trm(cvxc_matrix_t *X, cvxc_float_t alpha, const cvxc_matrix_t *A, int flags);
int cvxm_mult_trm(cvxc_matrix_t *X, cvxc_float_t alpha, const cvxc_matrix_t *A, int flags);
int cvxm_mult(cvxc_float_t beta, cvxc_matrix_t *C, cvxc_float_t alpha, const cvxc_matrix_t *A,
              const cvxc_matrix_t *B, int flags);
int cvxm_update_sym(cvxc_matrix_t *C, cvxc_float_t alpha, const cvxc_matrix_t *A, int flags);
int cvxm_update2_sym(cvxc_float_t beta, cvxc_matrix_t *C, cvxc_float_t alpha, const cvxc_matrix_t *A,
                     const cvxc_matrix_t *B, int flags);
int cvxm_mplus(cvxc_float_t alpha, cvxc_matrix_t *A, cvxc_float_t beta, const cvxc_matrix_t *B);


#endif

cvxc_float_t cvxm_min(const cvxc_matrix_t *x)
{
    const cvxc_float_t *vec =  x->data.elems;
    int n = x->data.rows * x->data.cols;
    int inc = x->data.rows == 1 ? x->data.step : 1;
    cvxc_float_t min = FLOAT_BIG;

    for (int i = 0; i < n; i++) {
        min = FMIN(min, vec[i * inc]);
    }
    return min;
}

cvxc_float_t cvxm_max(const cvxc_matrix_t *x)
{
    const cvxc_float_t *vec =  x->data.elems;
    int n = x->data.rows * x->data.cols;
    int inc = x->data.rows == 1 ? x->data.step : 1;
    cvxc_float_t max = -FLOAT_BIG;

    for (int i = 0; i < n; i++) {
        max = FMAX(max, vec[i * inc]);
    }
    return max;
}

//  \brief matrix-vector multiply; Y = beta*Y + alpha*A*X
int cvxm_mvmult(cvxc_float_t beta, cvxc_matrix_t *Y, cvxc_float_t alpha, const cvxc_matrix_t *A,
                const cvxc_matrix_t *X, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    int err = armas_mvmult(beta, &Y->data, alpha, &A->data, &X->data, flags, &cf);
    if (cvxm_isepi(Y) && cvxm_isepi(X)) {
        if ((flags & CVXC_TRANS) != 0) {
            Y->t *= beta;
        }
    }
    return err;
}

//  \brief symmetric matrix-vector multiply; Y = beta*Y + alpha*A*X
int cvxm_mvmult_sym(cvxc_float_t beta, cvxc_matrix_t *Y, cvxc_float_t alpha, const cvxc_matrix_t *A,
                    const cvxc_matrix_t *X, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    int err = armas_mvmult_sym(beta, &Y->data, alpha, &A->data, &X->data, flags, &cf);
    if (cvxm_isepi(Y) && cvxm_isepi(X)) {
        Y->t +=  beta * Y->t;
    }
    return err;
}

// get size of the work space needed to factor matrix A
cvxc_size_t cvxm_ldlwork(const cvxc_matrix_t *A)
{
    armas_conf_t cf = *armas_conf_default();
    armas_wbuf_t wb = (armas_wbuf_t){
        .buf = (char *)0,
        .offset = 0,
        .bytes = 0
    };
    armas_bkfactor_w((armas_dense_t *)&A->data, (armas_pivot_t *)0, ARMAS_LOWER, &wb, &cf);
    return wb.bytes;
}

cvxc_size_t cvxm_ldl_worksize(int N)
{
    armas_conf_t cf = *armas_conf_default();
    armas_wbuf_t wb = {0};
    armas_dense_t K = {0};
    K.rows = K.cols = N;
    armas_bkfactor_w(&K, (armas_pivot_t *)0, ARMAS_LOWER, &wb, &cf);
    return wb.bytes;
}

int cvxm_ldlfactor(cvxc_matrix_t *A, int *ipiv, int flags, cvxc_memblk_t *wrk)
{
    armas_pivot_t pv;
    size_t rows, cols;
    armas_conf_t cf = *armas_conf_default();
    armas_wbuf_t wb = (armas_wbuf_t){
        .buf = (char *)wrk->memory,
        .offset = 0,
        .bytes = wrk->mlen*sizeof(cvxc_float_t)
    };
    cvxm_size(&rows, &cols, A);
    armas_pivot_make(&pv, rows, ipiv);
    // min workspace is 2*N
    if (armas_bkfactor_w(&A->data, &pv, flags, &wb, &cf) < 0) {
        printf("bkfactor error: %d\n", cf.error);
        return -1;
    }
    return 0;
}
int cvxm_ldlsolve(cvxc_matrix_t *B, const cvxc_matrix_t *A, const int *ipiv, int flags, cvxc_memblk_t *wrk)
{
    armas_pivot_t pv;
    size_t rows, cols;
    armas_conf_t cf = *armas_conf_default();
    armas_wbuf_t wb = (armas_wbuf_t){
        .buf = (char *)wrk->memory,
        .offset = 0,
        .bytes = wrk->mlen*sizeof(cvxc_float_t)
    };

    cvxm_size(&rows, &cols, A);
    armas_pivot_make(&pv, rows, (int *)ipiv);
    // min workspace is 2*N
    if (armas_bksolve_w(&B->data, &A->data, &pv, flags, &wb, &cf) < 0) {
        printf("bksolve error: %d\n", cf.error);
    }
    return 0;
}

// \brief Cholesky factorization
int cvxm_cholfactor(cvxc_matrix_t *A, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_cholfactor(&A->data, ARMAS_NOPIVOT, flags, &cf);
}

int cvxm_cholsolve(cvxc_matrix_t *x, cvxc_matrix_t *A, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_cholsolve(&x->data, &A->data, ARMAS_NOPIVOT, flags, &cf);
}

// \brief SVD factorization
int cvxm_svd(cvxc_matrix_t *D, cvxc_matrix_t *U, cvxc_matrix_t *V, cvxc_matrix_t *A, int flags, cvxc_memblk_t *wrk)
{
    armas_conf_t cf = *armas_conf_default();
    armas_wbuf_t wb = (armas_wbuf_t){
        .buf = (char *)wrk->memory,
        .offset = 0,
        .bytes = wrk->mlen*sizeof(cvxc_float_t)
    };
    return armas_svd_w(&D->data, &U->data, &V->data, &A->data, flags, &wb, &cf);
}

// \brief Compute eigenvalues of symmetric matrix.
int cvxm_evd_sym(cvxc_matrix_t *D, cvxc_matrix_t *S, int flags, cvxc_memblk_t *work)
{
    armas_conf_t cf = *armas_conf_default();
    armas_wbuf_t wb = (armas_wbuf_t){
        .buf = (char *)work->memory,
        .offset = 0,
        .bytes = work->mlen*sizeof(cvxc_float_t)
    };
    if (armas_eigen_sym_w(&D->data, &S->data, flags, &wb, &cf) < 0) {
        printf("evd_sym error: %d\n", cf.error);
        return -1;
    }
    return 0;
}

// \brief Compute selected eigenvalues of symmetric matrix.
//  Parameter ival[2] hold  lower and upper limits of speficied interval.
int cvxm_evd_sym_selected(cvxc_matrix_t *D, cvxc_matrix_t *S, int *index, int flags, cvxc_memblk_t *wrk)
{
    armas_conf_t cf = *armas_conf_default();
    armas_wbuf_t wb = (armas_wbuf_t){
        .buf = (char *)wrk->memory,
        .offset = 0,
        .bytes = wrk->mlen*sizeof(cvxc_float_t)
    };
    return armas_eigen_sym_selected_w(&D->data, &S->data, ARMAS_EIGEN_INT(index[0], index[1]), flags, &wb, &cf);
}

int cvxm_lqfactor(cvxc_matrix_t *A, cvxc_matrix_t *tau, cvxc_memblk_t *wrk)
{
    armas_conf_t cf = *armas_conf_default();
    armas_wbuf_t wb = (armas_wbuf_t){
        .buf = (char *)wrk->memory,
        .offset = 0,
        .bytes = wrk->mlen*sizeof(cvxc_float_t)
    };
    return armas_lqfactor_w(&A->data, &tau->data, &wb, &cf);
}

int cvxm_lqmult(cvxc_matrix_t *C,
                const cvxc_matrix_t *A,
                const cvxc_matrix_t *tau,
                int flags,
                cvxc_memblk_t *wrk)
{
    armas_conf_t cf = *armas_conf_default();
    armas_wbuf_t wb = (armas_wbuf_t){
        .buf = (char *)wrk->memory,
        .offset = 0,
        .bytes = wrk->mlen*sizeof(cvxc_float_t)
    };
    return armas_lqmult_w(&C->data, &A->data, &tau->data, flags, &wb, &cf);
}

cvxc_matrix_t *cvxm_mkconst(cvxc_matrix_t *x, cvxc_float_t val)
{
    for (int j = 0; j < x->data.cols; j++) {
        for (int k = 0; k < x->data.rows; k++) {
            x->data.elems[k + j*x->data.step] = val;
        }
    }
    return x;
}

// \brief Make identity matrix or unit column vector 
cvxc_matrix_t *cvxm_mkident(cvxc_matrix_t *x)
{
    cvxm_mkconst(x, 0.0);
    for (int k = 0; k < x->data.cols; k++) {
        x->data.elems[k + k*x->data.step] = 1.0;
    }
    return x;
}


void cvxm_set_all(cvxc_matrix_t *A, cvxc_float_t val)
{
    cvxc_size_t n, m, i, j;
    cvxm_size(&m, &n, A);
    for (j = 0; j < n; j++) {
        for (i = 0; i < m; i++) {
            cvxm_set(A, i, j, val);
        }
    }
}

void cvxm_set_from(cvxc_matrix_t *A, cvxm_generator_t gen)
{
    cvxc_size_t n, m, i, j;
    cvxm_size(&m, &n, A);
    for (j = 0; j < n; j++) {
        for (i = 0; i < m; i++) {
            cvxm_set(A, i, j, gen());
        }
    }
}

int cvxm_apply2(cvxc_matrix_t *A, cvxm_operator2_t func, void *arg)
{
    return armas_apply2(&A->data, func, arg, 0);
}

void cvxm_unit_vector(cvxc_matrix_t *A)
{
    cvxc_size_t n, m, i;
    cvxm_size(&m, &n, A);

    cvxm_set(A, 0, 0, 1.0);
    for (i = 1; i < m; i++) {
        cvxm_set(A, i, 0, 0.0);
    }
}

// \brief Create new identity matrix
cvxc_matrix_t *cvxm_identity(cvxc_size_t n)
{
    cvxc_size_t k;
    cvxc_matrix_t *m = cvxm_new(n, n);
    for (k = 0; k < n; k++) {
        cvxm_set(m, k, k, 1.0);
    }
    return m;
}

// \brief print matrix to file stream
void cvxm_printf(FILE *fp, const char *fmt, const cvxc_matrix_t *x)
{
    armas_printf(fp, fmt, &x->data);
    if (cvxm_isepi(x)) {
        fprintf(fp, "/");
        fprintf(fp, fmt, x->t);
        fprintf(fp, "/\n");
    }
}

// \brief Matrix or vector norm
int cvxm_norm(cvxc_float_t *nrm, const cvxc_matrix_t *A, int norm)
{
    cvxc_size_t rows, cols;

    *nrm = 0;
    cvxm_size(&rows, &cols, A);

    if (rows != 1 && cols != 1 && norm == 2) {
        return -1;
    }

    switch (norm) {
    case 0:
        norm = ARMAS_NORM_FRB;
        break;
    case 1:
        norm = ARMAS_NORM_ONE;
        break;
    case 2:
        norm = ARMAS_NORM_TWO;
        break;
    default:
        norm = ARMAS_NORM_INF;
        break;
    }
    *nrm = armas_mnorm(&A->data, norm, NULLCONF);
    return 0;
}

void cvxm_mksymm(cvxc_matrix_t *x, int n)
{
    cvxc_matrix_t xcol, xrow;
    for (int i = 0; i < n-1; i++) {
        // i'th column and i'th row
        cvxm_view_map(&xcol, x, i+1, i, n-1-i, 1);
        cvxm_view_map(&xrow, x, i, i+1, 1, n-1-i);
        cvxm_copy(&xrow, &xcol, CVXC_ALL);
    }
}

