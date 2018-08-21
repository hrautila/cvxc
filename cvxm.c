
// Copyright: Harri Rautila, 2016

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
size_t cvxm_sizeof();
cvx_matrix_t *cvxm_new(cvx_size_t r, cvx_size_t c);
cvx_matrix_t *cvxm_newcopy(const cvx_matrix_t *A);
cvx_float_t *cvxm_data(const cvx_matrix_t *A, cvx_size_t k);

cvx_size_t cvxm_make(cvx_matrix_t *A, cvx_size_t rows, cvx_size_t cols,  void *mem, cvx_size_t nbytes);
cvx_matrix_t *cvxm_map_data(cvx_matrix_t *A, cvx_size_t rows, cvx_size_t cols,  cvx_float_t *data);

cvx_matrix_t *cvxm_view_map(cvx_matrix_t *A, const cvx_matrix_t *B, cvx_size_t r, cvx_size_t c, cvx_size_t nr, cvx_size_t nc);
void cvxm_view_unmap(cvx_matrix_t *A, cvx_matrix_t *B);
cvx_matrix_t *cvxm_view_diag(cvx_matrix_t *A, const cvx_matrix_t *B, int nd);
void cvxm_size(size_t *r, size_t *c, const cvx_matrix_t *A);
cvx_float_t cvxm_get(cvx_matrix_t *A, cvx_size_t r, cvx_size_t c);
void cvxm_set(cvx_matrix_t *A, cvx_size_t r, cvx_size_t c, cvx_float_t val);
void cvxm_apply(cvx_matrix_t *A, cvx_oper_t f, int flags);

cvx_float_t cvxm_dot(const cvx_matrix_t *X, const cvx_matrix_t *Y);
cvx_float_t cvxm_nrm2(const cvx_matrix_t *X);
cvx_float_t cvxm_asum(const cvx_matrix_t *X);
int cvxm_scale(cvx_matrix_t *X, cvx_float_t alpha, int flags);
int cvxm_add(cvx_matrix_t *X, cvx_float_t alpha, int flags);
void cvxm_make_trm(cvx_matrix_t *X, int flags);
void cvxm_copy(cvx_matrix_t *X, const cvx_matrix_t *Y, int flags);

int cvxm_axpy(cvx_matrix_t *Y, cvx_float_t alpha, const cvx_matrix_t *X);
cvx_float_t cvx_get(cvx_matrix_t *A, cvx_size_t r, cvx_size_t c);
cvx_float_t cvx_set(cvx_matrix_t *A, cvx_size_t r, cvx_size_t c, cvx_float_t val);
int cvxm_mvsolve_trm(cvx_matrix_t *X, cvx_float_t alpha, const cvx_matrix_t *A, int flags);
int cvxm_mvmult(cvx_float_t beta, cvx_matrix_t *Y, cvx_float_t alpha, const cvx_matrix_t *A, const cvx_matrix_t *X, int flags);
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

int cvxm_ldlfactor(cvx_matrix_t *A, int *ipiv, int flags, cvx_memblk_t *wrk)
{
    armas_pivot_t pv;
    size_t rows, cols;
    armas_d_dense_t W;
    armas_conf_t cf = *armas_conf_default();
    cvxm_size(&rows, &cols, A);
    armas_pivot_make(&pv, rows, ipiv);
    // min workspace is 2*N
    armas_d_make(&W, wrk->mlen, 1, wrk->mlen, (cvx_float_t*)wrk->memory);
    if (armas_d_bkfactor(A, &W, &pv, flags, &cf) < 0) {
        printf("bkfactor error: %d\n", cf.error);
    }
    return 0;
}
int cvxm_ldlsolve(cvx_matrix_t *B, const cvx_matrix_t *A, const int *ipiv, int flags, cvx_memblk_t *wrk)
{
    armas_pivot_t pv;
    size_t rows, cols;
    armas_d_dense_t W;
    armas_conf_t cf = *armas_conf_default();
    cvxm_size(&rows, &cols, A);
    armas_pivot_make(&pv, rows, (int *)ipiv);
    // min workspace is 2*N
    armas_d_make(&W, wrk->mlen, 1, wrk->mlen, (cvx_float_t*)wrk->memory);
    if (armas_d_bksolve(B, (cvx_matrix_t *)A, &W, &pv, flags, &cf) < 0) {
        printf("bksolve error: %d\n", cf.error);
    }
    return 0;
}

// \brief Cholesky factorization
int cvxm_cholfactor(cvx_matrix_t *A, int flags)
{
    return armas_d_cholfactor(A, (armas_d_dense_t *)0, ARMAS_NOPIVOT, flags, NULLCONF);
}

// \brief SVD factorization
int cvxm_svd(cvx_matrix_t *D, cvx_matrix_t *U, cvx_matrix_t *V, cvx_matrix_t *A, int flags, cvx_memblk_t *wrk)
{
    cvx_matrix_t w;
    armas_d_make(&w, wrk->mlen, 1, wrk->mlen, wrk->memory);
    return armas_d_svd(D, U, V, A, &w, flags, NULLCONF);
}

// \brief Compute eigenvalues of symmetric matrix.
int cvxm_evd_sym(cvx_matrix_t *D, cvx_matrix_t *S, int flags, cvx_memblk_t *work)
{
    cvx_matrix_t wrk;
    armas_d_make(&wrk, work->mlen, 1, work->mlen, work->memory);
    return armas_d_eigen_sym(D, S, &wrk, flags, NULLCONF);
}

// \brief Compute selected eigenvalues of symmetric matrix.
//  Parameter ival[2] hold  lower and upper limits of speficied interval.
int cvxm_evd_sym_selected(cvx_matrix_t *D, cvx_matrix_t *S, int *index, int flags, cvx_memblk_t *work)
{
    cvx_matrix_t wrk;
    armas_d_make(&wrk, work->mlen, 1, work->mlen, work->memory);
    return armas_d_eigen_sym_selected(D, S, &wrk, ARMAS_EIGEN_INT(index[0], index[1]), flags, NULLCONF);
}

cvx_matrix_t *cvxm_mkconst(cvx_matrix_t *x, cvx_float_t val)
{
    for (int j = 0; j < x->cols; j++) {
        for (int k = 0; k < x->rows; k++) {
            x->elems[k + j*x->step] = val;
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
    for (int k = 0; k < x->cols; k++) {
        x->elems[k + k*x->step] = 1.0;
    }
    return x;
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
    armas_d_printf(fp, fmt, x);
}

// \brief Compute workspace needed for SVD factorization of matrix of size [r, c]
cvx_size_t cvxm_svd_workspace(cvx_size_t r, cvx_size_t c)
{
    cvx_matrix_t Tmp;
    armas_d_make(&Tmp, r, c,  r, (cvx_float_t *)0);
    return armas_d_svd_work(&Tmp, ARMAS_WANTU, NULLCONF);
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
    *nrm = armas_d_norm(A, norm, 0, NULLCONF);
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
