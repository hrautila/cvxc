/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */

#ifndef __CVXC_CVXM_H
#define __CVXC_CVXM_H 1

#include <float.h>
#include <stdio.h>
#include <stddef.h>

#if ENABLE_FLOAT32

#include <armas/sdense.h>
typedef uint32_t cvxc_size_t;
typedef float cvxc_float_t;
typedef int32_t cvxc_int_t;

#define FMIN fminf
#define FMAX fmaxf
#define FLOAT_BIG FLT_MAX

#else

#include <armas/ddense.h>
typedef uint64_t cvxc_size_t;
typedef double cvxc_float_t;
typedef int64_t cvxc_int_t;

#define FMIN fmin
#define FMAX fmax
#define FLOAT_BIG DBL_MAX

#endif /* ENABLE_FLOAT32 */

//typedef armas_dense_t cvxc_matrix_t;
typedef struct cvxc_matrix {
    armas_dense_t data;
    int bits;
    cvxc_float_t t;
} cvxc_matrix_t;

typedef struct cvxc_epigraph {
    cvxc_matrix_t *m;
    int bits;
    cvxc_float_t t;
} cvxc_epigraph_t;

typedef armas_operator_t cvxm_oper_t;
typedef armas_generator_t cvxm_generator_t;
typedef armas_operator2_t cvxm_operator2_t;

typedef struct cvxc_memblk cvxc_memblk_t;

#define NULLCONF (armas_conf_t *)0
#define __cvxnil (cvxc_matrix_t *)0

typedef armas_iostream_t cvxc_stream_t;


typedef enum cvxc_flags {
    CVXC_ALL    = 0,
    CVXC_TRANS  = ARMAS_TRANS,
    CVXC_TRANSA = ARMAS_TRANSA,
    CVXC_TRANSB = ARMAS_TRANSB,
    CVXC_LOWER  = ARMAS_LOWER, 
    CVXC_UPPER  = ARMAS_UPPER,
    CVXC_UNIT   = ARMAS_UNIT,
    CVXC_DIAG   = ARMAS_UNIT,
    CVXC_INV    = ARMAS_UNIT,
    CVXC_LEFT   = ARMAS_LEFT,
    CVXC_RIGHT  = ARMAS_RIGHT,
    CVXC_WANTU  = ARMAS_WANTU,
    CVXC_WANTV  = ARMAS_WANTV,
    CVXC_DEBUG  = 0x8000000
} cvxc_flags_t;

enum cvxm_types {
    CVXM_SPARSE = 0x1,
    CVXM_EPIGRAPH = 0x4
};

__CVXC_INLINE
int cvxm_isepi(const cvxc_matrix_t *m)
{
    return (m->bits & CVXM_EPIGRAPH ) != 0;
}

__CVXC_INLINE
int cvxm_is_epigraph(const cvxc_epigraph_t *m)
{
    return (m->bits & CVXM_EPIGRAPH ) != 0;
}

__CVXC_INLINE
void cvxm_init(cvxc_matrix_t *m, cvxc_size_t r, cvxc_size_t c)
{
    armas_init(&m->data, (int)r, (int)c);
    m->bits = 0;
    m->t  = 0.0;
}

__CVXC_INLINE
void cvxm_epi_make(cvxc_epigraph_t *x_e, cvxc_matrix_t *x, int epi)
{
    x_e->m = x;
    x_e->t = 0.0;
    x_e->bits = epi ? CVXM_EPIGRAPH : 0;
}

// __CVXC_INLINE

__CVXC_INLINE
size_t cvxm_sizeof()
{
    return sizeof(cvxc_matrix_t);
}

// new matrix; call with (0,0) allowed, returns empty matrix
__CVXC_INLINE
cvxc_matrix_t *cvxm_new(cvxc_size_t r, cvxc_size_t c)
{
    cvxc_matrix_t *m = (cvxc_matrix_t *)calloc(1, sizeof(cvxc_matrix_t));
    if (m) {
        cvxm_init(m, r, c);
    }
    return m;
}

__CVXC_INLINE
void cvxm_release(cvxc_matrix_t *X)
{
    armas_release(&X->data);
}

__CVXC_INLINE
void cvxm_free(cvxc_matrix_t *X)
{
    if (X) {
        armas_release(&X->data);
        free(X);
    }
}

// create a copy of matrix 
//__CVXC_INLINE
//cvxc_matrix_t *cvxm_newcopy(const cvxc_matrix_t *A)
//{
//    return armas_newcopy((cvxc_matrix_t *)A);
//}

// \brief Get pointer to underlying storage of k'th element in a matrix.
__CVXC_INLINE
cvxc_float_t *cvxm_data(const cvxc_matrix_t *A, cvxc_size_t k)
{
    return &armas_data(&A->data)[k];
}

__CVXC_INLINE
cvxc_size_t cvxm_make(cvxc_matrix_t *A, cvxc_size_t rows, cvxc_size_t cols,  void *data, cvxc_size_t nbytes)
{
    cvxc_size_t nb = rows * cols * sizeof(cvxc_float_t);
    if (nbytes < nb)
        return 0;
    armas_make(&A->data, (int)rows, (int)cols, (int)rows, (cvxc_float_t *)data);
    A->bits = 0;
    A->t = 0.0;
    return nb;
}

__CVXC_INLINE
cvxc_size_t cvxm_make_epi(cvxc_matrix_t *A, cvxc_size_t rows, cvxc_size_t cols,  void *data, cvxc_size_t nbytes)
{
    cvxc_size_t nb = rows * cols * sizeof(cvxc_float_t);
    if (nbytes < nb)
        return 0;
    armas_make(&A->data, (int)rows, (int)cols, (int)rows, (cvxc_float_t *)data);
    A->bits = CVXM_EPIGRAPH;
    A->t = 0.0;
    return nb;
}

// \brief Create matrix view of elements pointed by data
__CVXC_INLINE
cvxc_matrix_t *cvxm_map_data(cvxc_matrix_t *A, cvxc_size_t rows, cvxc_size_t cols,  cvxc_float_t *data)
{
    armas_make(&A->data, (int)rows, (int)cols, (int)rows, data);
    A->bits = 0;
    A->t = 0.0;
    return A;
}

__CVXC_INLINE
cvxc_matrix_t *cvxm_map_data_epi(cvxc_matrix_t *A, cvxc_size_t rows, cvxc_size_t cols,  cvxc_float_t *data)
{
    armas_make(&A->data, (int)rows, (int)cols, (int)rows, data);
    A->bits = CVXM_EPIGRAPH;
    A->t = 0.0;
    return A;
}

// create submatrix (matrix view)
__CVXC_INLINE
cvxc_matrix_t *cvxm_view_map(cvxc_matrix_t *A, const cvxc_matrix_t *B,
                            cvxc_size_t r, cvxc_size_t c, cvxc_size_t nr, cvxc_size_t nc)
{
    armas_submatrix(&A->data, (armas_dense_t *)&B->data, (int)r, (int)c, (int)nr, (int)nc);
    A->bits = 0;
    A->t = 0.0;
    return A;
}

// release submatrix view
//__CVXC_INLINE
//void cvxm_view_unmap(cvxc_matrix_t *A, cvxc_matrix_t *B)
//{
// empty
//}

// create diagonal view
__CVXC_INLINE
cvxc_matrix_t *cvxm_view_diag(cvxc_matrix_t *A, const cvxc_matrix_t *B, int nd)
{
    armas_diag(&A->data, &B->data, nd);
    A->bits = 0;
    A->t = 0.0;
    return A;
}

// get matrix size
__CVXC_INLINE
void cvxm_size(size_t *r, size_t *c, const cvxc_matrix_t *A)
{
    *r = A ? A->data.rows : 0;
    *c = A ? A->data.cols : 0;
}

__CVXC_INLINE
cvxc_float_t cvxm_get(const cvxc_matrix_t *A, cvxc_size_t r, cvxc_size_t c)
{
    return armas_get_unsafe(&A->data, r, c);
}

__CVXC_INLINE
cvxc_float_t cvxm_get_epi(const cvxc_matrix_t *A)
{
    return A->t;
}

__CVXC_INLINE
cvxc_float_t cvxm_get_epival(const cvxc_epigraph_t *y)
{
    return y->t;
}

__CVXC_INLINE
void cvxm_set(cvxc_matrix_t *A, cvxc_size_t r, cvxc_size_t c, cvxc_float_t val)
{
    armas_set_unsafe(&A->data, r, c, val);
}

__CVXC_INLINE
void cvxm_set_epi(cvxc_matrix_t *A, cvxc_float_t v)
{
    A->t = v;
}

__CVXC_INLINE
void cvxm_set_epival(cvxc_epigraph_t *y, cvxc_float_t v)
{
    y->t = v;
}

// elementwise
__CVXC_INLINE
void cvxm_apply(cvxc_matrix_t *A, cvxm_oper_t f, int flags)
{
    armas_apply(&A->data, f, flags);
}

// vector-vector dot product
__CVXC_INLINE
cvxc_float_t cvxm_dot(const cvxc_matrix_t *X, const cvxc_matrix_t *Y)
{
    if (!X || !Y)
        return 0;

    armas_conf_t cf = *armas_conf_default();
    cvxc_float_t val = armas_dot(&X->data, &Y->data, &cf);
    if (cvxm_isepi(X) && cvxm_isepi(Y)) {
        val += X->t * Y->t;
    }
    return  val;
}

// vector-vector dot product
__CVXC_INLINE
cvxc_float_t cvxm_epi_dot(const cvxc_epigraph_t *x, const cvxc_epigraph_t *y)
{
    if (!x || !y)
        return 0;

    armas_conf_t cf = *armas_conf_default();
    cvxc_float_t val = armas_dot(&x->m->data, &y->m->data, &cf);
    if (cvxm_is_epigraph(x) && cvxm_is_epigraph(y)) {
        val += x->t * y->t;
    }
    return  val;
}

// \brief absolute maximum element
__CVXC_INLINE
cvxc_float_t cvxm_amax(const cvxc_matrix_t *X)
{
    armas_conf_t cf = *armas_conf_default();
    if (X)
        return armas_amax(&X->data, &cf);
    return 0.0;
}

// \brief sum of absolute values
__CVXC_INLINE
cvxc_float_t cvxm_asum(const cvxc_matrix_t *X)
{
    armas_conf_t cf = *armas_conf_default();
    if (X)
        return armas_asum(&X->data, &cf);
    return 0.0;
}

// \brief vector nrm2
__CVXC_INLINE
cvxc_float_t cvxm_nrm2(const cvxc_matrix_t *X)
{
    armas_conf_t cf = *armas_conf_default();
    if (X)
        return armas_nrm2(&X->data, &cf);
    return 0.0;
}

// \brief element-wise scale with constant
__CVXC_INLINE
int cvxm_scale(cvxc_matrix_t *X, cvxc_float_t alpha, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    int err = armas_mscale(&X->data, alpha, flags, &cf);
    if (cvxm_isepi(X)) {
        X->t *= alpha;
    }
    return err;
}

__CVXC_INLINE
int cvxm_epi_scale(cvxc_epigraph_t *X, cvxc_float_t alpha, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    int err = armas_mscale(&X->m->data, alpha, flags, &cf);
    if (cvxm_is_epigraph(X)) {
        X->t *= alpha;
    }
    return err;
}

// \brief Element-wise add constant
__CVXC_INLINE
int cvxm_add(cvxc_matrix_t *X, cvxc_float_t alpha, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_madd(&X->data, alpha, flags, &cf);
}

__CVXC_INLINE
void cvxm_make_trm(cvxc_matrix_t *X, int flags)
{
    armas_make_trm(&X->data, flags);
}

// \brief copy 
__CVXC_INLINE
void cvxm_copy(cvxc_matrix_t *X, const cvxc_matrix_t *Y, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    if (X && Y) {
        armas_mcopy(&X->data, (armas_dense_t *)&Y->data, flags, &cf);
        X->t = Y->t;
    }
}

// \brief copy
__CVXC_INLINE
void cvxm_epi_copy(cvxc_epigraph_t *X, const cvxc_epigraph_t *Y, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    if (X && Y) {
        armas_mcopy(&X->m->data, (armas_dense_t *)&Y->m->data, flags, &cf);
        X->t = Y->t;
    }
}

// \brief axpy; y = y + alpha*x
__CVXC_INLINE
int cvxm_axpy(cvxc_matrix_t *Y, cvxc_float_t alpha, const cvxc_matrix_t *X)
{
    int err = 0;
    if (X && Y) {
        armas_conf_t cf = *armas_conf_default();
        err = armas_axpy(&Y->data, alpha, &X->data, &cf);
        if (cvxm_isepi(Y) && cvxm_isepi(X)) {
            Y-> t += X->t * alpha;
        }
    }
    return err;
}

// \brief axpy; y = y + alpha*x
__CVXC_INLINE
int cvxm_epi_axpy(cvxc_epigraph_t *y, cvxc_float_t alpha, const cvxc_epigraph_t *x)
{
    int err = 0;
    if (x && y) {
        armas_conf_t cf = *armas_conf_default();
        err = armas_axpy(&y->m->data, alpha, &x->m->data, &cf);
        if (cvxm_is_epigraph(y) && cvxm_is_epigraph(x)) {
            y->t += x->t * alpha;
        }
    }
    return err;
}

__CVXC_INLINE
int cvxm_axpby(cvxc_float_t beta, cvxc_matrix_t *Y, cvxc_float_t alpha, const cvxc_matrix_t *X)
{
    int err = 0;
    if (X && Y) {
        armas_conf_t cf = *armas_conf_default();
        err = armas_axpby(beta, &Y->data, alpha, &X->data, &cf);
        if (cvxm_isepi(X) && cvxm_isepi(Y)) {
            Y->t = beta*Y->t + alpha*X->t;
        }
    }
    return err;
}

__CVXC_INLINE
int cvxm_epi_axpby(cvxc_float_t beta, cvxc_epigraph_t *y, cvxc_float_t alpha, const cvxc_epigraph_t *x)
{
    int err = 0;
    if (x && y) {
        armas_conf_t cf = *armas_conf_default();
        err = armas_axpby(beta, &y->m->data, alpha, &x->m->data, &cf);
        if (cvxm_is_epigraph(x) && cvxm_is_epigraph(y)) {
            y->t = beta*y->t + alpha*x->t;
        }
    }
    return err;
}

extern int cvxm_epi_mvmult(cvxc_float_t beta, cvxc_epigraph_t *Y, cvxc_float_t alpha,
                       const cvxc_matrix_t *A, const cvxc_epigraph_t *X, int flags);
extern int cvxm_mvmult(cvxc_float_t beta, cvxc_matrix_t *Y, cvxc_float_t alpha,
                       const cvxc_matrix_t *A, const cvxc_matrix_t *X, int flags);
extern int cvxm_mvmult_sym(cvxc_float_t beta, cvxc_matrix_t *Y, cvxc_float_t alpha,
                           const cvxc_matrix_t *A, const cvxc_matrix_t *X, int flags);

// \brief matrix-vector solve; X = alpha*A.-1*X
__CVXC_INLINE
int cvxm_mvsolve_trm(cvxc_matrix_t *X, cvxc_float_t alpha, const cvxc_matrix_t *A, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_mvsolve_trm(&X->data, alpha, &A->data, flags, &cf);
}


//  \brief matrix-vector rank update; A = A + alpha*X*Y
__CVXC_INLINE
int cvxm_mvupdate(cvxc_matrix_t *A, cvxc_float_t alpha, const cvxc_matrix_t *X, const cvxc_matrix_t *Y)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_mvupdate(1.0, &A->data, alpha, &X->data, &Y->data, &cf);
}

//  \brief diag solve; X = alpha*diag(A).-1*X
__CVXC_INLINE
int cvxm_solve_diag(cvxc_matrix_t *X, cvxc_float_t alpha, const cvxc_matrix_t *A, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_solve_diag(&X->data, alpha, &A->data, flags, &cf);
}

//  \brief diag mult; X = alpha*diag(A)*X
__CVXC_INLINE
int cvxm_mult_diag(cvxc_matrix_t *X, cvxc_float_t alpha, const cvxc_matrix_t *A, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_mult_diag(&X->data, alpha, &A->data, flags, &cf);
}

//  \brief matrix-matrix solve; X = alpha*A.-1*X
__CVXC_INLINE
int cvxm_solve_trm(cvxc_matrix_t *X, cvxc_float_t alpha, const cvxc_matrix_t *A, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_solve_trm(&X->data, alpha, &A->data, flags, &cf);
}

//  \brief triangiar matrix multiply; X = alpha*A*X
__CVXC_INLINE
int cvxm_mult_trm(cvxc_matrix_t *X, cvxc_float_t alpha, const cvxc_matrix_t *A, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_mult_trm(&X->data, alpha, &A->data, flags, &cf);
}

//  \brief matrix-matrix multiply; C = beta*C + alpha*A*B
__CVXC_INLINE
int cvxm_mult(cvxc_float_t beta, cvxc_matrix_t *C, cvxc_float_t alpha, const cvxc_matrix_t *A, 
              const cvxc_matrix_t *B, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_mult(beta, &C->data, alpha, &A->data, &B->data, flags, &cf);
}

//  \brief symmetric rank-k update; C = C + alpha*A*A^T
__CVXC_INLINE
int cvxm_update_sym(cvxc_matrix_t *C, cvxc_float_t alpha, const cvxc_matrix_t *A, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_update_sym(1.0, &C->data, alpha, &A->data, flags, &cf);
}

//  \brief symmetric rank 2k update; C = C + alpha*B*B^T
__CVXC_INLINE
int cvxm_update2_sym(cvxc_float_t beta, cvxc_matrix_t *C, cvxc_float_t alpha, const cvxc_matrix_t *A, 
                     const cvxc_matrix_t *B, int flags)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_update2_sym(beta, &C->data, alpha, &A->data, &B->data, flags, &cf);
}

__CVXC_INLINE
int cvxm_mplus(cvxc_float_t alpha, cvxc_matrix_t *A, cvxc_float_t beta, const cvxc_matrix_t *B)
{
    armas_conf_t cf = *armas_conf_default();
    return armas_mplus(alpha, &A->data, beta, &B->data, 0, &cf);
}


extern int cvxm_mm_read_file(cvxc_matrix_t *A, FILE *fp);
extern int cvxm_write_file(FILE *fp, const cvxc_matrix_t *m);
extern int cvxm_json_write_file(FILE *fp, const cvxc_matrix_t *m);
extern int cvxm_json_read_file(cvxc_matrix_t *m, FILE *fp);
extern int cvxm_read_file(cvxc_matrix_t *m, FILE *fp);

extern void cvxm_set_all(cvxc_matrix_t *A, cvxc_float_t val);
extern void cvxm_unit_vector(cvxc_matrix_t *A);
extern cvxc_float_t cvxm_min(const cvxc_matrix_t *x);
extern cvxc_float_t cvxm_max(const cvxc_matrix_t *x);
extern void cvxm_set_from(cvxc_matrix_t *A, cvxm_generator_t func);
extern int cvxm_apply2(cvxc_matrix_t *A, cvxm_operator2_t func, void *p);

// other functions
extern cvxc_matrix_t *cvxm_mkident(cvxc_matrix_t *x);
extern cvxc_matrix_t *cvxm_mkconst(cvxc_matrix_t *x, cvxc_float_t val);
extern cvxc_matrix_t *cvxm_identity(cvxc_size_t n);
extern void cvxm_printf(FILE *, const char *fmt, const cvxc_matrix_t *);

extern void cvxm_mksymm(cvxc_matrix_t *A, int n);
extern int cvxm_norm(cvxc_float_t *nrm, const cvxc_matrix_t *A, int norm);
extern void cvxm_zero(cvxc_matrix_t *A, int flags);

extern cvxc_size_t cvxm_ldlwork(const cvxc_matrix_t *A);
extern cvxc_size_t cvxm_ldl_worksize(int N);
// \brief LDL factorization
extern int cvxm_ldlfactor(cvxc_matrix_t *A, int *ipiv, int flags, cvxc_memblk_t *work);
// \brief Solve equations using computed LDL factorization
extern int cvxm_ldlsolve(cvxc_matrix_t *B, const cvxc_matrix_t *A, const int *ipiv, int flags, cvxc_memblk_t *wrk);
// \brief Eigenvalues of symmetric matrix
extern int cvxm_evd_sym(cvxc_matrix_t *D, cvxc_matrix_t *S, int flags, cvxc_memblk_t *work);
// \brief Selected eigenvalues of symmetric matrix
extern int cvxm_evd_sym_selected(cvxc_matrix_t *D, cvxc_matrix_t *S, int *ival, int flags, cvxc_memblk_t *work);
// \brief Cholesky factorization
extern int cvxm_cholfactor(cvxc_matrix_t *A, int flags);
extern int cvxm_cholsolve(cvxc_matrix_t *x, cvxc_matrix_t *A, int flags);
// \brief SVD factorization
extern int cvxm_svd(cvxc_matrix_t *D, cvxc_matrix_t *U, cvxc_matrix_t *V, cvxc_matrix_t *A, int flags, cvxc_memblk_t *wrk);
extern int cvxm_lqfactor(cvxc_matrix_t *A, cvxc_matrix_t *tau, cvxc_memblk_t *wrk);
extern int cvxm_lqmult(cvxc_matrix_t *C, const cvxc_matrix_t *A, const cvxc_matrix_t *tau, int flags, cvxc_memblk_t *wrk);
extern int cvxm_qrfactor(cvxc_matrix_t *A, cvxc_matrix_t *tau, cvxc_memblk_t *wrk);
extern int cvxm_qrmult(cvxc_matrix_t *C, const cvxc_matrix_t *A, const cvxc_matrix_t *tau, int flags, cvxc_memblk_t *wrk);

// compute workspace needed for SVD factorization of matrix [r, c]
extern cvxc_size_t cvxm_svd_workspace(cvxc_size_t r, cvxc_size_t c);

extern void cvxm_libstart();
extern void cvxm_libstop();

#endif // __CVXC_CVXM_H
