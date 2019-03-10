
// Copyright: Harri Rautila, 2018 <harri.rautila@gmail.com>

#if ! defined(__STDC_VERSION__)
#  if defined(__CVX_INLINE)
#    undef __CVX_INLINE
#  endif
#  define __CVX_INLINE
#endif

#include <ctype.h>
#include <errno.h>

#include "convex.h"
#include "cvxm.h"

#if __STDC_VERSION__
// inline function reference
void cvxm_init(cvx_matrix_t *m, cvx_size_t r, cvx_size_t c);
void cvxm_release(cvx_matrix_t *X);
void cvxm_free(cvx_matrix_t *X);
size_t cvxm_sizeof();

int cvxm_isepi(const cvx_matrix_t *m);
cvx_float_t cvxm_get_epi(const cvx_matrix_t *m);
void cvxm_set_epi(cvx_matrix_t *m, cvx_float_t v);

cvx_matrix_t *cvxm_new(cvx_size_t r, cvx_size_t c);
cvx_matrix_t *cvxm_newcopy(const cvx_matrix_t *A);
cvx_float_t *cvxm_data(const cvx_matrix_t *A, cvx_size_t k);

cvx_size_t cvxm_make(cvx_matrix_t *A, cvx_size_t rows, cvx_size_t cols,  void *mem, cvx_size_t nbytes);
cvx_size_t cvxm_make_epi(cvx_matrix_t *A, cvx_size_t rows, cvx_size_t cols,  void *mem, cvx_size_t nbytes);
cvx_matrix_t *cvxm_map_data(cvx_matrix_t *A, cvx_size_t rows, cvx_size_t cols,  cvx_float_t *data);

cvx_matrix_t *cvxm_view_map(cvx_matrix_t *A, const cvx_matrix_t *B, cvx_size_t r, cvx_size_t c, cvx_size_t nr, cvx_size_t nc);
void cvxm_view_unmap(cvx_matrix_t *A, cvx_matrix_t *B);
cvx_matrix_t *cvxm_view_diag(cvx_matrix_t *A, const cvx_matrix_t *B, int nd);
void cvxm_size(size_t *r, size_t *c, const cvx_matrix_t *A);
cvx_float_t cvxm_get(const cvx_matrix_t *A, cvx_size_t r, cvx_size_t c);
void cvxm_set(cvx_matrix_t *A, cvx_size_t r, cvx_size_t c, cvx_float_t val);
void cvxm_apply(cvx_matrix_t *A, cvx_oper_t f, int flags);

cvx_float_t cvxm_dot(const cvx_matrix_t *X, const cvx_matrix_t *Y);
cvx_float_t cvxm_nrm2(const cvx_matrix_t *X);
cvx_float_t cvxm_amax(const cvx_matrix_t *X);
int cvxm_scale(cvx_matrix_t *X, cvx_float_t alpha, int flags);
int cvxm_add(cvx_matrix_t *X, cvx_float_t alpha, int flags);
void cvxm_make_trm(cvx_matrix_t *X, int flags);
void cvxm_copy(cvx_matrix_t *X, const cvx_matrix_t *Y, int flags);

int cvxm_axpy(cvx_matrix_t *Y, cvx_float_t alpha, const cvx_matrix_t *X);
int cvxm_axpby(cvx_float_t beta, cvx_matrix_t *Y, cvx_float_t alpha, const cvx_matrix_t *X);
cvx_float_t cvx_get(cvx_matrix_t *A, cvx_size_t r, cvx_size_t c);
cvx_float_t cvx_set(cvx_matrix_t *A, cvx_size_t r, cvx_size_t c, cvx_float_t val);
int cvxm_mvsolve_trm(cvx_matrix_t *X, cvx_float_t alpha, const cvx_matrix_t *A, int flags);
int cvxm_mvupdate(cvx_matrix_t *A, cvx_float_t alpha, const cvx_matrix_t *X, const cvx_matrix_t *Y);
int cvxm_solve_diag(cvx_matrix_t *X, cvx_float_t alpha, const cvx_matrix_t *A, int flags);
int cvxm_mult_diag(cvx_matrix_t *X, cvx_float_t alpha, const cvx_matrix_t *A, int flags);
int cvxm_solve_trm(cvx_matrix_t *X, cvx_float_t alpha, const cvx_matrix_t *A, int flags);
int cvxm_mult_trm(cvx_matrix_t *X, cvx_float_t alpha, const cvx_matrix_t *A, int flags);
int cvxm_mult(cvx_float_t beta, cvx_matrix_t *C, cvx_float_t alpha, const cvx_matrix_t *A, 
              const cvx_matrix_t *B, int flags);
int cvxm_update_sym(cvx_matrix_t *C, cvx_float_t alpha, const cvx_matrix_t *A, int flags);
int cvxm_update2_sym(cvx_float_t beta, cvx_matrix_t *C, cvx_float_t alpha, const cvx_matrix_t *A, 
                     const cvx_matrix_t *B, int flags);


#endif

//  \brief matrix-vector multiply; Y = beta*Y + alpha*A*X
int cvxm_mvmult(cvx_float_t beta, cvx_matrix_t *Y, cvx_float_t alpha, const cvx_matrix_t *A,
                const cvx_matrix_t *X, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    int err = armas_d_mvmult(beta, &Y->data, alpha, &A->data, &X->data, flags, &cf);
    if (cvxm_isepi(Y) && cvxm_isepi(X)) {
        if ((flags & CVX_TRANS) != 0) {
            Y->t *= beta;
        }
    }
    return err;
}

//  \brief symmetric matrix-vector multiply; Y = beta*Y + alpha*A*X
int cvxm_mvmult_sym(cvx_float_t beta, cvx_matrix_t *Y, cvx_float_t alpha, const cvx_matrix_t *A,
                    const cvx_matrix_t *X, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    int err = armas_d_mvmult_sym(beta, &Y->data, alpha, &A->data, &X->data, flags, &cf);
    if (cvxm_isepi(Y) && cvxm_isepi(X)) {
        Y->t +=  beta * Y->t;
    }
    return err;
}

// get size of the work space needed to factor matrix A
cvx_size_t cvxm_ldlwork(const cvx_matrix_t *A)
{
    armas_conf_t cf = *armas_conf_default();
    armas_wbuf_t wb = (armas_wbuf_t){
        .buf = (char *)0,
        .offset = 0,
        .bytes = 0
    };
    armas_d_bkfactor_w((armas_d_dense_t *)&A->data, (armas_pivot_t *)0, ARMAS_LOWER, &wb, &cf);
    return wb.bytes;
}

int cvxm_ldlfactor(cvx_matrix_t *A, int *ipiv, int flags, cvx_memblk_t *wrk)
{
    armas_pivot_t pv;
    size_t rows, cols;
    //armas_d_dense_t W;
    armas_conf_t cf = *armas_conf_default();
    armas_wbuf_t wb = (armas_wbuf_t){
        .buf = (char *)wrk->memory,
        .offset = 0,
        .bytes = wrk->mlen*sizeof(cvx_float_t)
    };
    cvxm_size(&rows, &cols, A);
    armas_pivot_make(&pv, rows, ipiv);
    // min workspace is 2*N
    //armas_d_make(&W, wrk->mlen, 1, wrk->mlen, (cvx_float_t*)wrk->memory);
    if (armas_d_bkfactor_w(&A->data, &pv, flags, &wb, &cf) < 0) {
        printf("bkfactor error: %d\n", cf.error);
    }
    return 0;
}
int cvxm_ldlsolve(cvx_matrix_t *B, const cvx_matrix_t *A, const int *ipiv, int flags, cvx_memblk_t *wrk)
{
    armas_pivot_t pv;
    size_t rows, cols;
    //armas_d_dense_t W;
    armas_conf_t cf = *armas_conf_default();
    armas_wbuf_t wb = (armas_wbuf_t){
        .buf = (char *)wrk->memory,
        .offset = 0,
        .bytes = wrk->mlen*sizeof(cvx_float_t)
    };

    cvxm_size(&rows, &cols, A);
    armas_pivot_make(&pv, rows, (int *)ipiv);
    // min workspace is 2*N
    //armas_d_make(&W, wrk->mlen, 1, wrk->mlen, (cvx_float_t*)wrk->memory);
    if (armas_d_bksolve_w(&B->data, &A->data, &pv, flags, &wb, &cf) < 0) {
        printf("bksolve error: %d\n", cf.error);
    }
    return 0;
}

// \brief Cholesky factorization
int cvxm_cholfactor(cvx_matrix_t *A, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_d_cholfactor(&A->data, ARMAS_NOPIVOT, flags, &cf);
}

int cvxm_cholsolve(cvx_matrix_t *x, cvx_matrix_t *A, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_d_cholsolve(&x->data, &A->data, ARMAS_NOPIVOT, flags, &cf);
}

// \brief SVD factorization
int cvxm_svd(cvx_matrix_t *D, cvx_matrix_t *U, cvx_matrix_t *V, cvx_matrix_t *A, int flags, cvx_memblk_t *wrk)
{
    //cvx_matrix_t w;
    armas_conf_t cf = *armas_conf_default();
    armas_wbuf_t wb = (armas_wbuf_t){
        .buf = (char *)wrk->memory,
        .offset = 0,
        .bytes = wrk->mlen*sizeof(cvx_float_t)
    };
    //armas_d_make(&w, wrk->mlen, 1, wrk->mlen, wrk->memory);
    return armas_d_svd_w(&D->data, &U->data, &V->data, &A->data, flags, &wb, &cf);
}

// \brief Compute eigenvalues of symmetric matrix.
int cvxm_evd_sym(cvx_matrix_t *D, cvx_matrix_t *S, int flags, cvx_memblk_t *work)
{
    //cvx_matrix_t wrk;
    //armas_d_make(&wrk, work->mlen, 1, work->mlen, work->memory);
    armas_conf_t cf = *armas_conf_default();
    armas_wbuf_t wb = (armas_wbuf_t){
        .buf = (char *)work->memory,
        .offset = 0,
        .bytes = work->mlen*sizeof(cvx_float_t)
    };
    if (armas_d_eigen_sym_w(&D->data, &S->data, flags, &wb, &cf) < 0) {
        printf("evd_sym error: %d\n", cf.error);
        return -1;
    }
    return 0;
}

// \brief Compute selected eigenvalues of symmetric matrix.
//  Parameter ival[2] hold  lower and upper limits of speficied interval.
int cvxm_evd_sym_selected(cvx_matrix_t *D, cvx_matrix_t *S, int *index, int flags, cvx_memblk_t *wrk)
{
    //cvx_matrix_t wrk;
    //armas_d_make(&wrk, work->mlen, 1, work->mlen, work->memory);
    armas_conf_t cf = *armas_conf_default();
    armas_wbuf_t wb = (armas_wbuf_t){
        .buf = (char *)wrk->memory,
        .offset = 0,
        .bytes = wrk->mlen*sizeof(cvx_float_t)
    };
    return armas_d_eigen_sym_selected_w(&D->data, &S->data, ARMAS_EIGEN_INT(index[0], index[1]), flags, &wb, &cf);
}

int cvxm_lqfactor(cvx_matrix_t *A, cvx_matrix_t *tau, cvx_memblk_t *wrk)
{
    //cvx_matrix_t W;
    //armas_d_make(&W, wrk->mlen, 1, wrk->mlen, wrk->memory);
    armas_conf_t cf = *armas_conf_default();
    armas_wbuf_t wb = (armas_wbuf_t){
        .buf = (char *)wrk->memory,
        .offset = 0,
        .bytes = wrk->mlen*sizeof(cvx_float_t)
    };
    return armas_d_lqfactor_w(&A->data, &tau->data, &wb, &cf);
}

int cvxm_lqmult(cvx_matrix_t *C,
                const cvx_matrix_t *A,
                const cvx_matrix_t *tau,
                int flags,
                cvx_memblk_t *wrk)
{
    //armas_conf_t *cf = armas_conf_default();
    //cvx_matrix_t W;
    //armas_d_make(&W, wrk->mlen, 1, wrk->mlen, wrk->memory);
    armas_conf_t cf = *armas_conf_default();
    armas_wbuf_t wb = (armas_wbuf_t){
        .buf = (char *)wrk->memory,
        .offset = 0,
        .bytes = wrk->mlen*sizeof(cvx_float_t)
    };
    return armas_d_lqmult_w(&C->data, &A->data, &tau->data, flags, &wb, &cf);
}

cvx_matrix_t *cvxm_mkconst(cvx_matrix_t *x, cvx_float_t val)
{
    for (int j = 0; j < x->data.cols; j++) {
        for (int k = 0; k < x->data.rows; k++) {
            x->data.elems[k + j*x->data.step] = val;
        }
    }
    return x;
}

// \brief Make identity matrix or unit column vector 
cvx_matrix_t *cvxm_mkident(cvx_matrix_t *x)
{
    // this does not work for non-continues submatrix or row vector
    //memset(x->elems, 0, x->rows*x->cols*sizeof(x->elems[0]));
    cvxm_mkconst(x, 0.0);
    for (int k = 0; k < x->data.cols; k++) {
        x->data.elems[k + k*x->data.step] = 1.0;
    }
    return x;
}


void cvxm_set_all(cvx_matrix_t *A, cvx_float_t val)
{
    cvx_size_t n, m, i, j;
    cvxm_size(&m, &n, A);
    for (j = 0; j < n; j++) {
        for (i = 0; i < m; i++) {
            cvxm_set(A, i, j, val);
        }
    }
}

void cvxm_unit_vector(cvx_matrix_t *A)
{
    cvx_size_t n, m, i;
    cvxm_size(&m, &n, A);

    cvxm_set(A, 0, 0, 1.0);
    for (i = 1; i < m; i++) {
        cvxm_set(A, i, 0, 0.0);
    }
}

// \brief Create new identity matrix
cvx_matrix_t *cvxm_identity(cvx_size_t n)
{
    cvx_size_t k;
    cvx_matrix_t *m = cvxm_new(n, n);
    for (k = 0; k < n; k++) {
        cvxm_set(m, k, k, 1.0);
    }
    return m;
}

// \brief Create new column vector set to `val`
cvx_matrix_t *cvxm_new_vector(cvx_size_t n, cvx_float_t val)
{
    cvx_size_t k;
    cvx_matrix_t *m = cvxm_new(n, 1);
    if (val != 0.0) {
        for (k = 0; k < n; k++) {
            cvxm_set(m, k, 0, val);
        }
    }
    return m;
}

// \brief Create new column unit vector with `val`
cvx_matrix_t *cvxm_new_unit_vector(cvx_size_t n, cvx_float_t val)
{
    cvx_matrix_t *m = cvxm_new(n, 1);
    cvxm_set(m, 0, 0, val);
    return m;
}

// \brief print matrix to file stream
void cvxm_printf(FILE *fp, const char *fmt, const cvx_matrix_t *x)
{
#if 0
    if (armas_d_isvector(x) && x->cols == 1) {
        armas_d_dense_t row;
        armas_d_col_as_row(&row, x);
        armas_d_printf(fp, fmt, &row);
        return;
    }
#endif
    armas_d_printf(fp, fmt, &x->data);
    if (cvxm_isepi(x)) {
        fprintf(fp, "/");
        fprintf(fp, fmt, x->t);
        fprintf(fp, "/\n");
    }
}

// \brief Matrix or vector norm
int cvxm_norm(cvx_float_t *nrm, const cvx_matrix_t *A, int norm)
{
    cvx_size_t rows, cols;

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
    *nrm = armas_d_norm(&A->data, norm, 0, NULLCONF);
    return 0;
}

void cvxm_mksymm(cvx_matrix_t *x, int n)
{
    cvx_matrix_t xcol, xrow;
    for (int i = 0; i < n-1; i++) {
        // i'th column and i'th row
        cvxm_view_map(&xcol, x, i+1, i, n-1-i, 1);
        cvxm_view_map(&xrow, x, i, i+1, 1, n-1-i);
        cvxm_copy(&xrow, &xcol, CVX_ALL);
    }
}


// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
