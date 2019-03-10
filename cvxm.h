
// Copyright: Harri Rautila, 2018, <harri.rautila@gmail.com>

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
typedef int64_t cvx_int_t;

//typedef armas_d_dense_t cvx_matrix_t;
typedef struct cvx_matrix {
    armas_d_dense_t data;
    int bits;
    cvx_float_t t;
} cvx_matrix_t;


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
    CVX_WANTV  = ARMAS_WANTV,
    CVX_DEBUG  = 0x8000000
} cvx_flags_t;

enum cvxm_types {
    CVXM_SPARSE = 0x1,
    CVXM_EPIGRAPH = 0x4
};

__CVX_INLINE
int cvxm_isepi(const cvx_matrix_t *m) {
    return (m->bits & CVXM_EPIGRAPH ) != 0;
}

__CVX_INLINE
void cvxm_init(cvx_matrix_t *m, cvx_size_t r, cvx_size_t c)
{
    armas_d_init(&m->data, (int)r, (int)c);
    m->bits = 0;
    m->t  = 0.0;
}

__CVX_INLINE
size_t cvxm_sizeof()
{
    return sizeof(cvx_matrix_t);
}

// new matrix; call with (0,0) allowed, returns empty matrix
__CVX_INLINE
cvx_matrix_t *cvxm_new(cvx_size_t r, cvx_size_t c)
{
    cvx_matrix_t *m = (cvx_matrix_t *)calloc(1, sizeof(cvx_matrix_t));
    if (m) {
        cvxm_init(m, r, c);
    }
    return m;
}

__CVX_INLINE
void cvxm_release(cvx_matrix_t *X)
{
    armas_d_release(&X->data);
}

__CVX_INLINE
void cvxm_free(cvx_matrix_t *X)
{
    if (X) {
        armas_d_release(&X->data);
        free(X);
    }
}

// create a copy of matrix 
//__CVX_INLINE
//cvx_matrix_t *cvxm_newcopy(const cvx_matrix_t *A)
//{
//    return armas_d_newcopy((cvx_matrix_t *)A);
//}

// \brief Get pointer to underlying storage of k'th element in a matrix.
__CVX_INLINE
cvx_float_t *cvxm_data(const cvx_matrix_t *A, cvx_size_t k)
{
    return &armas_d_data(&A->data)[k];
}

__CVX_INLINE
cvx_size_t cvxm_make(cvx_matrix_t *A, cvx_size_t rows, cvx_size_t cols,  void *data, cvx_size_t nbytes)
{
    cvx_size_t nb = rows * cols * sizeof(cvx_float_t);
    if (nbytes < nb)
        return 0;
    armas_d_make(&A->data, (int)rows, (int)cols, (int)rows, (cvx_float_t *)data);
    A->bits = 0;
    A->t = 0.0;
    return nb;
}

__CVX_INLINE
cvx_size_t cvxm_make_epi(cvx_matrix_t *A, cvx_size_t rows, cvx_size_t cols,  void *data, cvx_size_t nbytes)
{
    cvx_size_t nb = rows * cols * sizeof(cvx_float_t);
    if (nbytes < nb)
        return 0;
    armas_d_make(&A->data, (int)rows, (int)cols, (int)rows, (cvx_float_t *)data);
    A->bits = CVXM_EPIGRAPH;
    A->t = 0.0;
    return nb;
}

// \brief Create matrix view of elements pointed by data
__CVX_INLINE
cvx_matrix_t *cvxm_map_data(cvx_matrix_t *A, cvx_size_t rows, cvx_size_t cols,  cvx_float_t *data)
{
    armas_d_make(&A->data, (int)rows, (int)cols, (int)rows, data);
    A->bits = 0;
    A->t = 0.0;
    return A;
}

__CVX_INLINE
cvx_matrix_t *cvxm_map_data_epi(cvx_matrix_t *A, cvx_size_t rows, cvx_size_t cols,  cvx_float_t *data)
{
    armas_d_make(&A->data, (int)rows, (int)cols, (int)rows, data);
    A->bits = CVXM_EPIGRAPH;
    A->t = 0.0;
    return A;
}

// create submatrix (matrix view)
__CVX_INLINE
cvx_matrix_t *cvxm_view_map(cvx_matrix_t *A, const cvx_matrix_t *B,
                            cvx_size_t r, cvx_size_t c, cvx_size_t nr, cvx_size_t nc)
{
    armas_d_submatrix(&A->data, (armas_d_dense_t *)&B->data, (int)r, (int)c, (int)nr, (int)nc);
    A->bits = 0;
    A->t = 0.0;
    return A;
}

// release submatrix view
//__CVX_INLINE
//void cvxm_view_unmap(cvx_matrix_t *A, cvx_matrix_t *B)
//{
// empty
//}

// create diagonal view
__CVX_INLINE
cvx_matrix_t *cvxm_view_diag(cvx_matrix_t *A, const cvx_matrix_t *B, int nd)
{
    armas_d_diag(&A->data, &B->data, nd);
    A->bits = 0;
    A->t = 0.0;
    return A;
}

// get matrix size
__CVX_INLINE
void cvxm_size(size_t *r, size_t *c, const cvx_matrix_t *A)
{
    *r = A->data.rows;
    *c = A->data.cols;
}

__CVX_INLINE
cvx_float_t cvxm_get(const cvx_matrix_t *A, cvx_size_t r, cvx_size_t c)
{
    return armas_d_get_unsafe(&A->data, r, c);
}

__CVX_INLINE
cvx_float_t cvxm_get_epi(const cvx_matrix_t *A)
{
    return A->t;
}

__CVX_INLINE
void cvxm_set(cvx_matrix_t *A, cvx_size_t r, cvx_size_t c, cvx_float_t val)
{
    armas_d_set_unsafe(&A->data, r, c, val);
}

__CVX_INLINE
void cvxm_set_epi(cvx_matrix_t *A, cvx_float_t v)
{
    A->t = v;
}

// elementwise
__CVX_INLINE
void cvxm_apply(cvx_matrix_t *A, cvx_oper_t f, int flags)
{
    armas_d_apply(&A->data, f, flags);
}

// vector-vector dot product
__CVX_INLINE
cvx_float_t cvxm_dot(const cvx_matrix_t *X, const cvx_matrix_t *Y)
{
    armas_conf_t cf = *armas_conf_default();
    cvx_float_t val = armas_d_dot(&X->data, &Y->data, &cf);
    if (cvxm_isepi(X) && cvxm_isepi(Y)) {
        val += X->t * Y->t;
    }
    return  val;
}

// \brief sum of absolute values
__CVX_INLINE
cvx_float_t cvxm_amax(const cvx_matrix_t *X)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_d_amax(&X->data, &cf);
}

// \brief vector nrm2
__CVX_INLINE
cvx_float_t cvxm_nrm2(const cvx_matrix_t *X)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_d_nrm2(&X->data, &cf);
}

// \brief element-wise scale with constant
__CVX_INLINE
int cvxm_scale(cvx_matrix_t *X, cvx_float_t alpha, int flags)
{
    int err = armas_d_mscale(&X->data, alpha, flags);
    if (cvxm_isepi(X)) {
        X->t *= alpha;
    }
    return err;
}

// \brief Element-wise add constant
__CVX_INLINE
int cvxm_add(cvx_matrix_t *X, cvx_float_t alpha, int flags)
{
    return armas_d_madd(&X->data, alpha, flags);
}

__CVX_INLINE
void cvxm_make_trm(cvx_matrix_t *X, int flags)
{
    armas_d_make_trm(&X->data, flags);
}

// \brief copy 
__CVX_INLINE
void cvxm_copy(cvx_matrix_t *X, const cvx_matrix_t *Y, int flags)
{
    armas_d_mcopy(&X->data, (armas_d_dense_t *)&Y->data);
    X->t = Y->t;
}

// \brief axpy; y = y + alpha*x
__CVX_INLINE
int cvxm_axpy(cvx_matrix_t *Y, cvx_float_t alpha, const cvx_matrix_t *X)
{
    armas_conf_t cf = *armas_conf_default();
    int err = armas_d_axpy(&Y->data, alpha, &X->data, &cf);
    if (cvxm_isepi(Y) && cvxm_isepi(X)) {
        Y-> t += X->t * alpha;
    }
    return err;
}

__CVX_INLINE
int cvxm_axpby(cvx_float_t beta, cvx_matrix_t *Y, cvx_float_t alpha, const cvx_matrix_t *X)
{
    armas_conf_t cf = *armas_conf_default();
    int err = armas_d_axpby(beta, &Y->data, alpha, &X->data, &cf);
    if (cvxm_isepi(X) && cvxm_isepi(Y)) {
        Y->t = beta*Y->t + alpha*X->t;
    }
    return err;
}

extern int cvxm_mvmult(cvx_float_t beta, cvx_matrix_t *Y, cvx_float_t alpha,
                       const cvx_matrix_t *A, const cvx_matrix_t *X, int flags);
extern int cvxm_mvmult_sym(cvx_float_t beta, cvx_matrix_t *Y, cvx_float_t alpha,
                           const cvx_matrix_t *A, const cvx_matrix_t *X, int flags);

// \brief matrix-vector solve; X = alpha*A.-1*X
__CVX_INLINE
int cvxm_mvsolve_trm(cvx_matrix_t *X, cvx_float_t alpha, const cvx_matrix_t *A, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_d_mvsolve_trm(&X->data, alpha, &A->data, flags, &cf);
}


//  \brief matrix-vector rank update; A = A + alpha*X*Y
__CVX_INLINE
int cvxm_mvupdate(cvx_matrix_t *A, cvx_float_t alpha, const cvx_matrix_t *X, const cvx_matrix_t *Y)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_d_mvupdate(&A->data, alpha, &X->data, &Y->data, &cf);
}

//  \brief diag solve; X = alpha*diag(A).-1*X
__CVX_INLINE
int cvxm_solve_diag(cvx_matrix_t *X, cvx_float_t alpha, const cvx_matrix_t *A, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_d_solve_diag(&X->data, alpha, &A->data, flags, &cf);
}

//  \brief diag mult; X = alpha*diag(A)*X
__CVX_INLINE
int cvxm_mult_diag(cvx_matrix_t *X, cvx_float_t alpha, const cvx_matrix_t *A, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_d_mult_diag(&X->data, alpha, &A->data, flags, &cf);
}

//  \brief matrix-matrix solve; X = alpha*A.-1*X
__CVX_INLINE
int cvxm_solve_trm(cvx_matrix_t *X, cvx_float_t alpha, const cvx_matrix_t *A, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_d_solve_trm(&X->data, alpha, &A->data, flags, &cf);
}

//  \brief triangiar matrix multiply; X = alpha*A*X
__CVX_INLINE
int cvxm_mult_trm(cvx_matrix_t *X, cvx_float_t alpha, const cvx_matrix_t *A, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_d_mult_trm(&X->data, alpha, &A->data, flags, &cf);
}

//  \brief matrix-matrix multiply; C = beta*C + alpha*A*B
__CVX_INLINE
int cvxm_mult(cvx_float_t beta, cvx_matrix_t *C, cvx_float_t alpha, const cvx_matrix_t *A, 
              const cvx_matrix_t *B, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_d_mult(beta, &C->data, alpha, &A->data, &B->data, flags, &cf);
}

//  \brief symmetric rank-k update; C = C + alpha*A*A^T
__CVX_INLINE
int cvxm_update_sym(cvx_matrix_t *C, cvx_float_t alpha, const cvx_matrix_t *A, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_d_update_sym(1.0, &C->data, alpha, &A->data, flags, &cf);
}

//  \brief symmetric rank 2k update; C = C + alpha*B*B^T
__CVX_INLINE
int cvxm_update2_sym(cvx_float_t beta, cvx_matrix_t *C, cvx_float_t alpha, const cvx_matrix_t *A, 
                     const cvx_matrix_t *B, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_d_update2_sym(beta, &C->data, alpha, &A->data, &B->data, flags, &cf);
}


extern int cvxm_mmload(cvx_matrix_t *A, FILE *fp);
extern int cvxm_write_file(FILE *fp, const cvx_matrix_t *m);
extern int cvxm_json_write_file(FILE *fp, const cvx_matrix_t *m);
extern int cvxm_read_file(cvx_matrix_t *m, FILE *fp);

extern void cvxm_set_all(cvx_matrix_t *A, cvx_float_t val);
extern void cvxm_unit_vector(cvx_matrix_t *A);

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

extern cvx_size_t cvxm_ldlwork(const cvx_matrix_t *A);
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
extern int cvxm_cholsolve(cvx_matrix_t *x, cvx_matrix_t *A, int flags);
// \brief SVD factorization
extern int cvxm_svd(cvx_matrix_t *D, cvx_matrix_t *U, cvx_matrix_t *V, cvx_matrix_t *A, int flags, cvx_memblk_t *wrk);
extern int cvxm_lqfactor(cvx_matrix_t *A, cvx_matrix_t *tau, cvx_memblk_t *wrk);
extern int cvxm_lqmult(cvx_matrix_t *C, const cvx_matrix_t *A, const cvx_matrix_t *tau, int flags, cvx_memblk_t *wrk);
extern int cvxm_qrfactor(cvx_matrix_t *A, cvx_matrix_t *tau, cvx_memblk_t *wrk);
extern int cvxm_qrmult(cvx_matrix_t *C, const cvx_matrix_t *A, const cvx_matrix_t *tau, int flags, cvx_memblk_t *wrk);

// compute workspace needed for SVD factorization of matrix [r, c]
extern cvx_size_t cvxm_svd_workspace(cvx_size_t r, cvx_size_t c);

#endif // __CVX_CVXM_H

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
