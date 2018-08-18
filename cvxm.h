
// Copyright: Harri Rautila, 2016

#ifndef __CVX_CVXM_H
#define __CVX_CVXM_H 1

#include <stdio.h>
#include <stddef.h>

#if ! defined(__CVX_INLINE)
#  if __GNUC__
#    if ! __STDC_VERSION__
#      define __CVX_INLINE extern inline
#    else
#      define __CVX_INLINE inline
#      define __ARMAS_INLINE inline
#    endif
#  else
#    define __CVX_INLINE inline
#  endif
#endif

#include <armas/dmatrix.h>

typedef uint64_t cvx_size_t;
typedef double cvx_float_t;

typedef armas_d_dense_t cvx_matrix_t;
typedef armas_d_operator_t cvx_oper_t;

typedef struct cvx_memblk cvx_memblk_t;

#define NULLCONF (armas_conf_t *)0
#define __cvxnil (cvx_matrix_t *)0

typedef enum cvx_flags {
    CVX_ALL    = 0,
    CVX_TRANS  = ARMAS_TRANS,
    CVX_TRANSA = ARMAS_TRANSA,
    CVX_TRANSB = ARMAS_TRANSB,
    CVX_LOWER  = ARMAS_LOWER, 
    CVX_UPPER  = ARMAS_UPPER,
    CVX_UNIT   = ARMAS_UNIT,
    CVX_DIAG   = ARMAS_UNIT,
    CVX_INV    = ARMAS_UNIT,
    CVX_LEFT   = ARMAS_LEFT,
    CVX_RIGHT  = ARMAS_RIGHT,
    CVX_WANTU  = ARMAS_WANTU,
    CVX_WANTV  = ARMAS_WANTV
} cvx_flags_t;

__CVX_INLINE
void cvxm_init(cvx_matrix_t *m, cvx_size_t r, cvx_size_t c)
{
    armas_d_init(m, (int)r, (int)c);
}

__CVX_INLINE
size_t cvxm_sizeof()
{
    return sizeof(armas_d_dense_t);
}

// new matrix; call with (0,0) allowed, returns empty matrix
__CVX_INLINE
cvx_matrix_t *cvxm_new(cvx_size_t r, cvx_size_t c)
{
    return armas_d_alloc((int)r, (int)c);
}

__CVX_INLINE
void cvxm_release(cvx_matrix_t *X)
{
    armas_d_release(X);
}

// create a copy of matrix 
__CVX_INLINE
cvx_matrix_t *cvxm_newcopy(const cvx_matrix_t *A)
{
    return armas_d_newcopy((cvx_matrix_t *)A);
}

// \brief Get pointer to underlying storage of k'th element in a matrix.
__CVX_INLINE
cvx_float_t *cvxm_data(const cvx_matrix_t *A, cvx_size_t k)
{
    return &armas_d_data(A)[k];
}

// \brief Create matrix view of elements pointed by data
__CVX_INLINE
cvx_matrix_t *cvxm_map_data(cvx_matrix_t *A, cvx_size_t rows, cvx_size_t cols,  cvx_float_t *data)
{
    armas_d_make(A, (int)rows, (int)cols, (int)rows, data);
    return A;
}

// create submatrix (matrix view)
__CVX_INLINE
cvx_matrix_t *cvxm_view_map(cvx_matrix_t *A, const cvx_matrix_t *B,
                            cvx_size_t r, cvx_size_t c, cvx_size_t nr, cvx_size_t nc)
{
    armas_d_submatrix(A, (cvx_matrix_t *)B, (int)r, (int)c, (int)nr, (int)nc);
    return A;
}

// release submatrix view
__CVX_INLINE
void cvxm_view_unmap(cvx_matrix_t *A, cvx_matrix_t *B)
{
    // empty
}

// create diagonal view
__CVX_INLINE
cvx_matrix_t *cvxm_view_diag(cvx_matrix_t *A, const cvx_matrix_t *B, int nd)
{
    armas_d_diag(A, B, nd);
    return A;
}

// get matrix size
__CVX_INLINE
void cvxm_size(size_t *r, size_t *c, const cvx_matrix_t *A)
{
    *r = A->rows;
    *c = A->cols;
}

__CVX_INLINE
cvx_float_t cvxm_get(cvx_matrix_t *A, cvx_size_t r, cvx_size_t c)
{
    return armas_d_get_unsafe(A, r, c);
}

__CVX_INLINE
void cvxm_set(cvx_matrix_t *A, cvx_size_t r, cvx_size_t c, cvx_float_t val)
{
    armas_d_set_unsafe(A, r, c, val);
}

// elementwise
__CVX_INLINE
void cvxm_apply(cvx_matrix_t *A, cvx_oper_t f, int flags)
{
    armas_d_apply(A, f, flags);
}

// vector-vector dot product
__CVX_INLINE
cvx_float_t cvxm_dot(const cvx_matrix_t *X, const cvx_matrix_t *Y)
{
    return armas_d_dot(X, Y, NULLCONF);
}

// \brief sum of absolute values
__CVX_INLINE
cvx_float_t cvxm_asum(const cvx_matrix_t *X)
{
    return armas_d_asum(X, NULLCONF);
}

// \brief vector nrm2
__CVX_INLINE
cvx_float_t cvxm_nrm2(const cvx_matrix_t *X)
{
    return armas_d_nrm2(X, NULLCONF);
}

// \brief element-wise scale with constant
__CVX_INLINE
int cvxm_scale(cvx_matrix_t *X, cvx_float_t alpha, int flags)
{
    return armas_d_mscale(X, alpha, flags);
}

// \brief Element-wise add constant
__CVX_INLINE
int cvxm_add(cvx_matrix_t *X, cvx_float_t alpha, int flags)
{
    return armas_d_madd(X, alpha, flags);
}

__CVX_INLINE
void cvxm_make_trm(cvx_matrix_t *X, int flags)
{
    armas_d_make_trm(X, flags);
}

// \brief copy 
__CVX_INLINE
void cvxm_copy(cvx_matrix_t *X, const cvx_matrix_t *Y, int flags)
{
    armas_d_mcopy(X, (cvx_matrix_t *)Y);
}

// \brief axpy; y = y + alpha*x
__CVX_INLINE
int cvxm_axpy(cvx_matrix_t *Y, cvx_float_t alpha, const cvx_matrix_t *X)
{
    return armas_d_axpy(Y, alpha, X, NULLCONF);
}

// \brief matrix-vector solve; X = alpha*A.-1*X
__CVX_INLINE
int cvxm_mvsolve_trm(cvx_matrix_t *X, cvx_float_t alpha, const cvx_matrix_t *A, int flags)
{
    return armas_d_mvsolve_trm(X, alpha, A, flags, NULLCONF);
}

//  \brief matrix-vector multiply; Y = alpha*A*X
__CVX_INLINE
int cvxm_mvmult(cvx_float_t beta, cvx_matrix_t *Y, cvx_float_t alpha, const cvx_matrix_t *A,
                const cvx_matrix_t *X, int flags)
{
    return armas_d_mvmult(1.0, Y, alpha, A, X, flags, NULLCONF);
}

//  \brief matrix-vector rank update; A = A + alpha*X*Y
__CVX_INLINE
int cvxm_mvupdate(cvx_matrix_t *A, cvx_float_t alpha, const cvx_matrix_t *X, const cvx_matrix_t *Y)
{
    return armas_d_mvupdate(A, alpha, X, Y, NULLCONF);
}

//  \brief diag solve; X = alpha*diag(A).-1*X
__CVX_INLINE
int cvxm_solve_diag(cvx_matrix_t *X, cvx_float_t alpha, const cvx_matrix_t *A, int flags)
{
    return armas_d_solve_diag(X, alpha, A, flags, NULLCONF);
}

//  \brief diag mult; X = alpha*diag(A)*X
__CVX_INLINE
int cvxm_mult_diag(cvx_matrix_t *X, cvx_float_t alpha, const cvx_matrix_t *A, int flags)
{
    return armas_d_mult_diag(X, alpha, A, flags, NULLCONF);
}

//  \brief matrix-matrix solve; X = alpha*A.-1*X
__CVX_INLINE
int cvxm_solve_trm(cvx_matrix_t *X, cvx_float_t alpha, const cvx_matrix_t *A, int flags)
{
    return armas_d_solve_trm(X, alpha, A, flags, NULLCONF);
}

//  \brief triangiar matrix multiply; X = alpha*A*X
__CVX_INLINE
int cvxm_mult_trm(cvx_matrix_t *X, cvx_float_t alpha, const cvx_matrix_t *A, int flags)
{
    return armas_d_mult_trm(X, alpha, A, flags, NULLCONF);
}

//  \brief matrix-matrix multiply; C = beta*C + alpha*A*B
__CVX_INLINE
int cvxm_mult(cvx_float_t beta, cvx_matrix_t *C, cvx_float_t alpha, const cvx_matrix_t *A, 
              const cvx_matrix_t *B, int flags)
{
    return armas_d_mult(beta, C, alpha, A, B, flags, NULLCONF);
}

//  \brief symmetric rank-k update; C = C + alpha*A*A^T
__CVX_INLINE
int cvxm_update_sym(cvx_matrix_t *C, cvx_float_t alpha, const cvx_matrix_t *A, int flags)
{
    return armas_d_update_sym(1.0, C, alpha, A, flags, NULLCONF);
}

//  \brief symmetric rank 2k update; C = C + alpha*B*B^T
__CVX_INLINE
int cvxm_update2_sym(cvx_float_t beta, cvx_matrix_t *C, cvx_float_t alpha, const cvx_matrix_t *A, 
                     const cvx_matrix_t *B, int flags)
{
    return armas_d_update2_sym(beta, C, alpha, A, B, flags, NULLCONF);
}


extern int cvxm_write_file(FILE *fp, const cvx_matrix_t *m);
extern int cvxm_json_write_file(FILE *fp, const cvx_matrix_t *m);
extern int cvxm_read_file(cvx_matrix_t *m, FILE *fp);

// other functions
extern cvx_matrix_t *cvxm_mkident(cvx_matrix_t *x);
extern cvx_matrix_t *cvxm_mkconst(cvx_matrix_t *x, cvx_float_t val);

extern cvx_matrix_t *cvxm_identity(cvx_size_t n);
extern cvx_matrix_t *cvxm_new_vector(cvx_size_t n, cvx_float_t val);
extern cvx_matrix_t *cvxm_new_unit_vector(cvx_size_t n, cvx_float_t val);
extern void cvxm_printf(FILE *, const char *fmt, const cvx_matrix_t *);

extern void cvxm_mksymm(cvx_matrix_t *A, int n);
extern int cvxm_norm(cvx_float_t *nrm, const cvx_matrix_t *A, int norm);
extern void cvxm_zero(cvx_matrix_t *A, int flags);

// \brief LDL factorization
extern int cvxm_ldlfactor(cvx_matrix_t *A, int *ipiv, int flags, cvx_memblk_t *work);
// \brief Solve equations using computed LDL factorization
extern int cvxm_ldlsolve(cvx_matrix_t *B, const cvx_matrix_t *A, const int *ipiv, int flags, cvx_memblk_t *wrk);
// \brief Eigenvalues of symmetric matrix
extern int cvxm_evd_sym(cvx_matrix_t *D, cvx_matrix_t *S, int flags, cvx_memblk_t *work);
// \brief Selected eigenvalues of symmetric matrix
extern int cvxm_evd_sym_selected(cvx_matrix_t *D, cvx_matrix_t *S, int *ival, int flags, cvx_memblk_t *work);
// \brief Cholesky factorization
extern int cvxm_cholfactor(cvx_matrix_t *A, int flags);
// \brief SVD factorization
extern int cvxm_svd(cvx_matrix_t *D, cvx_matrix_t *U, cvx_matrix_t *V, cvx_matrix_t *A, int flags, cvx_memblk_t *wrk);

// compute workspace needed for SVD factorization of matrix [r, c]
extern cvx_size_t cvxm_svd_workspace(cvx_size_t r, cvx_size_t c);

#endif // __CVX_CVXM_H

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
