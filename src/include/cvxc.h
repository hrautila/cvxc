
#ifndef __CVX_CONVEX_H
#define __CVX_CONVEX_H 1

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <math.h>

/**
 * LibTool version numbering.
 */
#define CVXC_ABI_CURRENT  0
#define CVXC_ABI_REVISION  0
#define CVXC_ABI_AGE  0

#define SQRT sqrt
#define __POW  pow

#define SQRT2           1.41421356237309504880
#define CVX_MAXITER     100
#define CVX_ABSTOL      1e-7
#define CVX_RELTOL      1e-6
#define CVX_FEASTOL     1e-7
#define CVX_STEP        0.99

#define CVX_STAT_OPTIMAL                0
#define CVX_STAT_UNKNOWN                1
#define CVX_STAT_PRIMAL_INFEASIBLE      2
#define CVX_STAT_DUAL_INFEASIBLE        3
#define CVX_STAT_SINGULAR               4
#define CVX_STAT_MAX                    4

#define CVX_ERR_NULLCOST        1
#define CVX_ERR_DIMC            2
#define CVX_ERR_DIMH            3
#define CVX_ERR_DIMG            4
#define CVX_ERR_DIMA            5
#define CVX_ERR_DIMB            6
#define CVX_ERR_RANK            7
#define CVX_ERR_MEMORY          8
#define CVX_ERR_NEG_INITIAL_S   9
#define CVX_ERR_NEG_INITIAL_Z   10
#define CVX_ERR_MAXITER         11
#define CVX_ERR_SINGULAR        12

#include "cvxm.h"

static inline
cvx_float_t __NaN()
{
#ifdef NAN
    return NAN;
#else
    return sqrt(-1.0);
#endif
}

static inline
int __isnull(const void *ptr) {
    return ptr == (void *)0;
}

static inline
cvx_size_t __aligned128(cvx_size_t n)
{
    return (n & 0xF) != 0 ? 16 - (n & 0xF) : 0;
}

static inline
cvx_size_t __aligned64(cvx_size_t n)
{
    return (n & 0x7) != 0 ? 8 - (n & 0x7) : 0;
}

static inline
cvx_float_t __maxvec(int n, const cvx_float_t *vec)
{
    cvx_float_t r = vec[0];
    for (int k = 1; k < n; k++) {
        if (vec[k] > r)
            r = vec[k];
    }
    return r;
}

static inline
cvx_float_t __minvec(int n, const cvx_float_t *vec)
{
    cvx_float_t r = vec[0];
    for (int k = 1; k < n; k++) {
        if (vec[k] < r)
            r = vec[k];
    }
    return r;
}

#define __MAX2(a, b)  __maxvec(2, (cvx_float_t []){(a), (b)})
#define __MIN2(a, b)  __minvec(2, (cvx_float_t []){(a), (b)})

// -------------------------------------------------------------------------------------
// problem dimension data structures
typedef enum {
    L = 0x1,
    S = 0x2,
    Q = 0x4,
    NL = 0x8,
    CVXDIM_LINEAR = 0x1,
    CVXDIM_SDP = 0x2,
    CVXDIM_SOCP = 0x4,
    CVXDIM_NONLINEAR = 0x8,
    CVXDIM_NLTARGET = 0x10,
    CVXDIM_CONVEX = 0x18,
    CVXDIM_CONELP = 0x7,
    CVXDIM_CONVEXLP = 0xF,
    CVXDIM_CONVEXPROG = 0x1F
} cvx_dim_enum;

typedef enum {
    D = 0,
    DI = 1,
    BETA = 2,
    V = 3,
    R = 4,
    RTI = 5,
    DNL = 6,
    DNLI = 7,
    DNLT = 8,
    DNLTI = 9,
    CVXWS_D = 0,
    CVXWS_DI = 1,
    CVXWS_BETA = 2,
    CVXWS_V = 3,
    CVXWS_R = 4,
    CVXWS_RTI = 5,
    CVXWS_DNL = 6,
    CVXWS_DNLI = 7,
    CVXWS_DNLT = 8,
    CVXWS_DNLTI = 9
} cvx_mset_enum;

typedef enum  {
    CVX_KKTSOL_LDL = 1,
    CVX_KKTSOL_LDL2 = 2,
    CVX_KKTSOL_QR = 3,
    CVX_KKTSOL_CHOL = 4,
    CVX_KKTSOL_CHOL2 = 5
} cvx_kktsolver_name_t;

typedef enum  {
    CVX_INDEX_NORMAL = 0,
    CVX_INDEX_PACKED = 1,
    CVX_INDEX_DIAG = 2,
    CVX_INDEX_SIGS = 3
} cvx_index_type;

// -------------------------------------------------------------------------------------

// very simple
typedef struct cvx_memblk {
    cvx_float_t *memory;
    cvx_size_t mlen;
    void *__bytes;
} cvx_memblk_t;

static inline
void __mblk_empty(cvx_memblk_t *m)
{
    if (m) {
        m->memory = (cvx_float_t *)0;
        m->mlen = 0;
    }
}
static inline
void __mblk_init(cvx_memblk_t *m, cvx_size_t n)
{
    m->memory = (cvx_float_t *)0;
    m->mlen = 0;
    cvx_float_t *space = (cvx_float_t *)calloc(n, sizeof(cvx_float_t));
    if (space) {
        m->memory = space;
        m->mlen = n;
        m->__bytes = space;
    }
}

static inline
cvx_size_t __mblk_make(cvx_memblk_t *m, cvx_size_t mlen, void *ptr, cvx_size_t nbytes)
{
    if (!ptr || nbytes < mlen)
        return 0;
    m->memory = (cvx_float_t *)ptr;
    m->mlen = mlen/sizeof(cvx_float_t);
    m->__bytes = (void *)0;
    return mlen;
}

static inline
void __mblk_release(cvx_memblk_t *m)
{
    if (m->__bytes)
        free(m->__bytes);
    m->memory = (cvx_float_t *)0;
    m->mlen = 0;
    m->__bytes = (void *)0;
}

static inline
cvx_float_t *__mblk_offset(cvx_memblk_t *m, cvx_size_t off)
{
    if (off >= m->mlen) {
        abort();
        return (cvx_float_t *)0;
    }
    return &m->memory[off];
}

static inline
void __mblk_subblk(cvx_memblk_t *d, cvx_memblk_t *m, cvx_size_t off)
{
    if (off >= m->mlen) {
        d->memory = (cvx_float_t *)0;
        d->mlen = 0;
    } else {
        d->memory = &m->memory[off];
        d->mlen = m->mlen - off;
    }
}

static inline
void __mblk_clear(cvx_memblk_t *m)
{
    if (m && m->mlen > 0)
        memset(m->memory, 0, m->mlen*sizeof(cvx_float_t));
}

//----------------------------------------------------------------------------------------

typedef struct cvx_dimset {
    cvx_size_t *qdims;         ///< Second Order Cone constraint dimensions
    cvx_size_t *sdims;         ///< SDP constraint dimensions
    cvx_size_t qlen;           ///< Number of SOCP dimensions
    cvx_size_t slen;           ///< Number of SDP dimensions
    cvx_size_t ldim;           ///< Linear constraints
    cvx_size_t mnl;            ///< Non-linear constraints
    cvx_size_t iscpt;          ///< Target function is non-linear convex (CP)
} cvx_dimset_t;

extern cvx_dimset_t *cvx_dimset_new(cvx_size_t nonlinear, cvx_size_t linear, cvx_size_t socp, cvx_size_t sdp);
extern void cvx_dimset_free(cvx_dimset_t *dims);
extern cvx_size_t cvx_dimset_count(const cvx_dimset_t *dims, cvx_dim_enum name);
extern cvx_dimset_t *cvx_dimset_create(cvx_dimset_t *dims, cvx_size_t nonlinear, cvx_size_t linear, cvx_size_t socp, cvx_size_t sdp);
extern cvx_dimset_t *cvx_dimset_alloc(cvx_dimset_t *dims, cvx_size_t linear, const cvx_size_t *socp, const cvx_size_t *sdp);
extern void cvx_dimset_release(cvx_dimset_t *dims);
extern cvx_size_t cvx_dimset_max(const cvx_dimset_t *dims, cvx_dim_enum name);
extern cvx_size_t cvx_dimset_sum(const cvx_dimset_t *dims, cvx_dim_enum name);
extern cvx_size_t cvx_dimset_sum_squared(const cvx_dimset_t *dims, cvx_dim_enum name);
extern cvx_size_t cvx_dimset_sum_packed(const cvx_dimset_t *dims, cvx_dim_enum name);

typedef struct cvx_index {
    cvx_size_t *index;          ///< Index data
    cvx_size_t *indq;           ///< SOCP vector
    cvx_size_t *inds;           ///< Indexes of SDP constraints
    cvx_size_t *indl;           ///< Index to linear
    cvx_size_t *indnl;          ///< Index to non-linear
    cvx_size_t *indnlt;         ///< Index to non-linear target
    cvx_size_t indlen;          ///< Total number of  index elements
    cvx_size_t qlen;            ///< Number of SOCP dimensions
    cvx_size_t slen;            ///< Number of SDP dimensions
    cvx_index_type type;        ///< Index type
    void *__bytes;
} cvx_index_t;

extern cvx_size_t cvx_index_bytes(const cvx_dimset_t *dims, int kind);
extern cvx_size_t cvx_index_count(const cvx_index_t *index, cvx_dim_enum name);
extern cvx_size_t cvx_index_make(cvx_index_t *ind, const cvx_dimset_t *dims, int kind, void *buf, cvx_size_t nbytes);
extern cvx_index_t *cvx_index_init(cvx_index_t *index, const cvx_dimset_t *dims, int packed);
extern void cvx_index_release(cvx_index_t *dims);
extern cvx_size_t cvx_index_elem(cvx_matrix_t *x, const cvx_matrix_t *y, const cvx_index_t *ind, cvx_dim_enum name, int k);
extern cvx_size_t cvx_index_length(const cvx_index_t *ind, cvx_dim_enum name);

extern void cvx_index_create(cvx_matrix_t *x, cvx_index_t *index, const cvx_dimset_t *dims, cvx_index_type kind);
extern void cvx_subindex(cvx_index_t *ind, const cvx_index_t *src, int parts);

typedef struct cvx_scaling {
    cvx_float_t *data;          // Scaling matrix data space, one big block
    cvx_size_t *indexes;        // Data space for V, R, RTI sizes and offsets
    cvx_float_t *dnlt;          // Pointer to DNLT vector
    cvx_float_t *dnlti;         // Pointer to DNLTI vector
    int dnltsz;                 // Length of DNLT/DNLTI vectors (0 or 1)
    cvx_float_t *dnl;           // Pointer to DNL vector
    cvx_float_t *dnli;          // Pointer to DNLI vector
    int dnlsz;                  // Length of DNL/DNLI vectors; equal to dims->mnl
    cvx_float_t *d;             // Pointer to D vector
    cvx_float_t *di;            // Pointer to DI vector
    int dsz;                    // Length of D/DI vectors; equal to dims->ldim
    cvx_float_t *beta;          // Pointer to BETA vector
    cvx_size_t *indv;           // Offsets to V vector;
    int vcount;                 // # of V vectors; equal to dims->qlen
    cvx_size_t *indr;           // Offsets to R matrices
    cvx_size_t *indrti;         // Offsets to RTI matrices
    int rcount;                 // # of R/RTI matrices, equal to dims->slen
    unsigned char *__bytes;
    cvx_size_t nbytes;
} cvx_scaling_t;

extern cvx_size_t cvx_scaling_bytes(cvx_size_t *isize, const cvx_dimset_t *dims);
extern cvx_size_t cvx_scaling_make(cvx_scaling_t *W, const cvx_dimset_t *dims, void *mem, size_t nbytes);

extern cvx_size_t cvx_scaling_elem(cvx_matrix_t *A, const cvx_scaling_t *W, cvx_mset_enum name, int ind);
extern cvx_scaling_t *cvx_scaling_init(cvx_scaling_t *W, const cvx_dimset_t *dims);
extern void cvx_scaling_initial_value(cvx_scaling_t *W);
extern void cvx_scaling_release(cvx_scaling_t *W);
extern void cvx_scaling_printf(FILE *f, const char *format, cvx_scaling_t *W, const char *s);
extern void cvx_scaling_elem_printf(FILE *f, const char *format,
                                    const cvx_scaling_t *W, cvx_mset_enum name, int index, const char *s);

extern int cvx_scaling_copy(cvx_scaling_t *W, const cvx_scaling_t *Ws);


/*
 * Matrix group stored in a continuous memory matrix.
 */
typedef struct cvx_matgrp {
    cvx_matrix_t *mat;
    cvx_index_t *index;
} cvx_matgrp_t;

#define __nilgrp (cvx_matgrp_t *)0

static inline
void cvx_mgrp_init(cvx_matgrp_t *g, cvx_matrix_t *m, cvx_index_t *ind)
{
    g->mat = m; g->index = ind;
}

static inline
cvx_size_t cvx_mgrp_elem(cvx_matrix_t *x, cvx_matgrp_t *g, cvx_dim_enum name, int k)
{
    return cvx_index_elem(x, g->mat, g->index, name, k);
}

static inline
int cvx_mgrp_count(cvx_matgrp_t *g, cvx_dim_enum name)
{
    return cvx_index_count(g->index, name);
}

extern void cvx_mgrp_copy(cvx_matgrp_t *x_g, cvx_matgrp_t *y_g, cvx_dim_enum name);
extern void cvx_mgrp_copy_lambda(cvx_matgrp_t *ds_g, cvx_matgrp_t *lmbda);
extern void cvx_mgrp_scale_sz(cvx_matgrp_t *ds_g, cvx_float_t v, int flags);
extern void cvx_mgrp_update_sz(cvx_matgrp_t *ds_g, cvx_float_t v, int flags);
extern void cvx_mgrp_axpby_sz(cvx_float_t beta, cvx_matgrp_t *y_g, cvx_float_t alpha, cvx_matgrp_t *x_g, int flags);
extern void cvx_mgrp_initial_value(cvx_matgrp_t *x_g, int flags);
extern void cvx_mgrp_printf(FILE *f, const char *format, cvx_matgrp_t *g, const char *s);
extern void cvx_mat_printf(FILE *f, const char *format, cvx_matrix_t *g, const char *s);
extern void cvx_mat_test_nan(const char *name, const cvx_matrix_t *A);
extern void cvx_mat_print_ifenv(const char *name, const cvx_matrix_t *A, const char *s);

// -------------------------------------------------------------------------------------

// \brief hyperbolic ||x||_2
extern cvx_float_t cvx_jnrm2(cvx_matrix_t *x);
// \brief hyperbolic x'*y
extern cvx_float_t cvx_jdot(cvx_matrix_t *x, cvx_matrix_t *y);
extern cvx_float_t cvx_snrm2(cvx_matgrp_t *x);
extern cvx_float_t cvx_sdot(cvx_matgrp_t *x, cvx_matgrp_t *y);
extern cvx_float_t cvx_sdot_elem(cvx_matgrp_t *x, cvx_matgrp_t *y, cvx_dim_enum name);
extern cvx_float_t cvx_snrm2_elem(cvx_matgrp_t *x, cvx_dim_enum name);

extern void cvx_pack(cvx_matrix_t *y, cvx_matrix_t *x, const cvx_index_t *index);
extern void cvx_unpack(cvx_matrix_t *y, cvx_matrix_t *x, const cvx_index_t *index);
extern void cvx_pack2(cvx_matrix_t *x, const cvx_index_t *index, cvx_memblk_t *work);

extern int cvx_scale(cvx_matgrp_t *x, cvx_scaling_t *W, int flags, cvx_memblk_t *work);
extern int cvx_scale_part(cvx_matgrp_t *x, cvx_scaling_t *W, int flags, int parts, cvx_memblk_t *work);
extern int cvx_scale2(cvx_matgrp_t *x_g, cvx_matgrp_t *lmbda_g, int flags, cvx_memblk_t *work);
extern int cvx_sprod(cvx_matgrp_t *x, cvx_matgrp_t *y, int flags, cvx_memblk_t *wrk);
extern int cvx_ssqr(cvx_matgrp_t *x, cvx_matgrp_t *y);
extern int cvx_sinv(cvx_matgrp_t *x, cvx_matgrp_t *y, cvx_memblk_t *work);

extern int cvx_sgemv(cvx_float_t beta, cvx_matrix_t *y,
                     cvx_float_t alpha, const cvx_matrix_t *A, cvx_matgrp_t *x, int flags);
extern int cvx_sgemv2(cvx_float_t beta, cvx_matrix_t *y, cvx_float_t alpha,
                      const cvx_matrix_t *A, const cvx_matrix_t *B, cvx_matgrp_t *x, int flags);
extern int cvx_triusc(cvx_matgrp_t *x);
extern int cvx_trisc(cvx_matgrp_t *x);
extern int cvx_mksymm(cvx_matgrp_t *x);

extern int
cvx_compute_scaling(cvx_scaling_t *W, cvx_matgrp_t *s_g, cvx_matgrp_t *z_g, cvx_matgrp_t *l_g, cvx_memblk_t *wrk);
extern int
cvx_update_scaling(cvx_scaling_t *W, cvx_matgrp_t *s_g, cvx_matgrp_t *z_g, cvx_matgrp_t *l_g, cvx_memblk_t *wrk);

// ------------------------------------------------------------------------------------------
// solver things

typedef int (*cvx_convex_f0)(cvx_matrix_t *x0, void *user);
typedef int (*cvx_convex_f1)(cvx_matrix_t *f, cvx_matrix_t *Df, const cvx_matrix_t *x, void *user);
typedef int (*cvx_convex_f2)(cvx_matrix_t *f, cvx_matrix_t *Df, cvx_matrix_t *H,
                             const cvx_matrix_t *x, const cvx_matrix_t *z, void *user);

typedef struct cvx_convex_program {
    cvx_convex_f2 F;
    void *user;
}  cvx_convex_program_t;

static inline
void cvx_convex_program_init(cvx_convex_program_t *cp, cvx_convex_f2 F, void *user)
{
    cp->F = F;
    cp->user = user;
}

static inline
int F0(cvx_convex_program_t *cp, cvx_matrix_t *x0)
{
    return cp->F(x0, __cvxnil, __cvxnil, __cvxnil, __cvxnil, cp->user);
}

static inline
int F1(cvx_convex_program_t *cp, cvx_matrix_t *f, cvx_matrix_t *Df, const cvx_matrix_t *x)
{
    int e = cp->F(f, Df, __cvxnil, x, __cvxnil, cp->user);
    if (e == 0 && cvxm_isepi(x)) {
        cvxm_set(f, 0, 0, cvxm_get(f, 0, 0) - x->t);
    }
    return e;
}

static inline
int F2(cvx_convex_program_t *cp, cvx_matrix_t *f, cvx_matrix_t *Df, cvx_matrix_t *H,
       const cvx_matrix_t *x, const cvx_matrix_t *z)
{
    int e = cp->F(f, Df, H, x, z, cp->user);
    if (e == 0 && cvxm_isepi(x) && f) {
        cvxm_set(f, 0, 0, cvxm_get(f, 0, 0) - x->t);
    }
    return e;
}


/* Forward declarations */
typedef struct cvx_problem cvx_problem_t;
typedef struct cvx_kktsolver cvx_kktsolver_t;

typedef struct cvx_kktfuncs {
    cvx_kktsolver_t *(*new)(cvx_problem_t *cp, int n, int m, const cvx_dimset_t *dims);
    int (*init)(cvx_kktsolver_t *S, cvx_problem_t *cp, int n, int m, const cvx_dimset_t *dims);
    int (*factor)(cvx_kktsolver_t *S, cvx_scaling_t *W, cvx_matrix_t *H, cvx_matrix_t *Df);
    int (*solve)(cvx_kktsolver_t *S, cvx_matrix_t *x, cvx_matrix_t *y, cvx_matgrp_t *z_g);
    cvx_size_t (*bytes)(int n, int m, const cvx_dimset_t *dims);
    cvx_size_t (*make)(cvx_kktsolver_t *S, cvx_problem_t *cp, int n, int m, const cvx_dimset_t *dims, void *ptr, cvx_size_t nbytes);
    void (*free)(cvx_kktsolver_t *S);
} cvx_kktfuncs_t;

#if 0
typedef struct cvx_kktsolver {
    cvx_kktfuncs_t *fnc;
    cvx_problem_t *cp;
} cvx_kktsolver_t;
#endif

typedef struct cvx_ldlsolver {
    cvx_kktfuncs_t *fnc;
    cvx_problem_t *cp;
    cvx_matrix_t K;
    cvx_matrix_t u;
    cvx_matrix_t g;
    cvx_memblk_t work;
    cvx_scaling_t *W;
    cvx_matrix_t *A;
    cvx_matrix_t *G;
    cvx_matrix_t *Df;
    const cvx_dimset_t *dims;
    size_t ldK;
    int *ipiv;
    size_t p;
    size_t n;
    size_t mnl;
} cvx_ldlsolver_t;

typedef struct cvx_qrsolver {
    cvx_kktfuncs_t *fnc;
    cvx_matrix_t QA;
    cvx_matrix_t tauA;
    cvx_matrix_t Gs;
    cvx_matrix_t tauG;
    cvx_matrix_t u;
    cvx_matrix_t vv;
    cvx_matrix_t w;
    cvx_memblk_t work;
    cvx_scaling_t *W;
    cvx_matrix_t *A;
    cvx_matrix_t *G;
    cvx_problem_t *cp;
    const cvx_dimset_t *dims;
    size_t p;
    size_t n;
} cvx_qrsolver_t;

typedef struct cvx_chainsolver {
    cvx_kktfuncs_t *fnc;
    cvx_problem_t *cp;
    cvx_kktsolver_t *next;
} cvx_chainsolver_t;

struct cvx_kktsolver {
    union {
        cvx_kktfuncs_t *fnc;
        cvx_ldlsolver_t ldl;
        cvx_chainsolver_t next;
    } u;
    int debug;
};

extern cvx_kktfuncs_t *cvx_ldlload(void *ptr);

__CVX_INLINE
void cvx_ldlsolver_init(cvx_kktsolver_t *kkt,
                        cvx_problem_t *cp,
                        int n, int m, const cvx_dimset_t *dims)
{
    cvx_kktfuncs_t *fnc = cvx_ldlload((void *)0);
    fnc->init(kkt, cp, n, m, dims);
}

__CVX_INLINE
int cvx_kktinit(cvx_kktsolver_t *S,
                cvx_problem_t *cp,
                int n, int m, const cvx_dimset_t *dims)
{
    return (*S->u.fnc->init)(S, cp, n, m, dims);
}

__CVX_INLINE
int cvx_kktfactor(cvx_kktsolver_t *S,
                  cvx_scaling_t *W, cvx_matrix_t *H, cvx_matrix_t *Df)
{
    return (*S->u.fnc->factor)(S, W, H, Df);
}

__CVX_INLINE
int cvx_kktsolve(cvx_kktsolver_t *S,
                 cvx_matrix_t *x, cvx_matrix_t *y, cvx_matgrp_t *z_g)
{
    return (*S->u.fnc->solve)(S, x, y, z_g);
}

typedef struct cvx_solution {
    int status;
    cvx_matrix_t *x;
    cvx_matrix_t *s;
    cvx_matrix_t *y;
    cvx_matrix_t *z;
    cvx_float_t primal_objective;
    cvx_float_t dual_objective;
    cvx_float_t gap;
    cvx_float_t relative_gap;
    cvx_float_t primal_infeasibility;
    cvx_float_t dual_infeasibility;
    cvx_float_t primal_slack;
    cvx_float_t dual_slack;
    cvx_float_t primal_residual_cert;
    cvx_float_t dual_residual_cert;
    int iterations;
} cvx_solution_t;

typedef struct cvx_solopts {
    cvx_float_t abstol;         ///< Absolute tolerance
    cvx_float_t reltol;         ///< Relative tolerance
    cvx_float_t feastol;        ///< Feasibility tolerance
    int max_iter;               ///< Maximum iterations
    int debug;                  ///< Debug
    int refinement;             ///< Refinement count
    int show_progress;          ///< Show progress of the iteration
    int kkt_solver_name;        ///< KKT solver function
} cvx_solopts_t;

typedef struct cvx_stats {
    cvx_float_t resx0;
    cvx_float_t resy0;
    cvx_float_t resz0;
    cvx_float_t resx;
    cvx_float_t resy;
    cvx_float_t resz;
    cvx_float_t hresx;
    cvx_float_t hresy;
    cvx_float_t hresz;
    cvx_float_t cx;
    cvx_float_t by;
    cvx_float_t hz;
    cvx_float_t rt;
    cvx_float_t dres;
    cvx_float_t pres;
    cvx_float_t dinfres;
    cvx_float_t pinfres;
    cvx_float_t gap;
    cvx_float_t relgap;
    cvx_float_t pcost;
    cvx_float_t dcost;
} cvx_stats_t;


typedef struct cvx_conelp_internal {
    cvx_float_t tau, kappa;
    cvx_float_t dkappa, dtau;
    cvx_float_t wkappa3;
    cvx_float_t nrms, nrmz;
    cvx_float_t dg, dgi;
    cvx_float_t ts, tz;                 // step sizes

    // internal result
    cvx_matrix_t x, y, s, z;

    cvx_matrix_t dx, dy, ds, dz;
    cvx_matrix_t x1, y1, z1;
    cvx_matrix_t rx, ry, rz;
    cvx_matrix_t hrx, hry, hrz;
    cvx_matrix_t sigs, sigz;
    cvx_matrix_t lmbda, lmbdasq;

    cvx_matrix_t wx, wy, ws, wz;
    cvx_matrix_t wx2, wy2, ws2, wz2;
    cvx_matrix_t ws3, wz3;
    cvx_matrix_t th;

    // matrix groups with correct indexing
    cvx_matgrp_t h_g;
    cvx_matgrp_t s_g, z_g;
    cvx_matgrp_t ds_g, dz_g;
    cvx_matgrp_t hrz_g;
    cvx_matgrp_t rz_g;
    cvx_matgrp_t ws_g, wz_g;
    cvx_matgrp_t ws2_g, wz2_g;
    cvx_matgrp_t ws3_g, wz3_g;
    cvx_matgrp_t th_g;
    cvx_matgrp_t z1_g;
    cvx_matgrp_t lmbda_g, lmbdasq_g;
    cvx_matgrp_t sigs_g, sigz_g;

    cvx_index_t index_full;             // indexing to full matrix group
    cvx_index_t index_packed;           // indexing to matrix group with packed 'S' space
    cvx_index_t index_diag;             // indexing to matrix group with diagonal 'S' space
    cvx_index_t index_sig;              // indexing to matrix group with only diagonal 'S" space

} cvx_conelp_internal_t;

typedef struct cvx_gpindex {
    cvx_size_t p;
    cvx_size_t *index;
    void *__bytes;
} cvx_gpindex_t;

typedef struct cvx_gp_params {
    cvx_matrix_t *F;
    cvx_matrix_t *g;
    cvx_matrix_t y;
    cvx_gpindex_t gpi;
} cvx_gp_params_t;

typedef struct cvx_gp_program {
    cvx_convex_program_t gp;
    cvx_gp_params_t gp_params;
} cvx_gp_program_t;



typedef struct cvx_cpl_internal {
    cvx_float_t tau, kappa;
    cvx_float_t dkappa, dtau;
    cvx_float_t wkappa3;
    cvx_float_t nrms, nrmz;
    cvx_float_t dg, dgi;
    cvx_float_t ts, tz;                 // step sizes
    cvx_float_t phi;
    cvx_float_t dphi;
    cvx_float_t phi0;
    cvx_float_t dphi0;
    cvx_float_t mu;
    cvx_float_t sigma;
    cvx_float_t sigma0;
    cvx_float_t eta;
    cvx_float_t eta0;
    cvx_float_t theta1;
    cvx_float_t theta2;
    cvx_float_t theta3;
    cvx_float_t step;
    cvx_float_t step0;
    cvx_float_t gap;
    cvx_float_t gap0;
    cvx_float_t dsdz;
    cvx_float_t dsdz0;

    cvx_float_t relgap;
    cvx_float_t resznl;
    cvx_float_t reszl;
    cvx_float_t resrznl;
    cvx_float_t resrzl;
    cvx_float_t relgap0;
    cvx_float_t resznl0;
    cvx_float_t reszl0;
    cvx_float_t resrznl0;
    cvx_float_t resrzl0;

    cvx_float_t resx0;
    cvx_float_t resy0;
    cvx_float_t resz0;
    cvx_float_t resx;
    cvx_float_t resy;
    cvx_float_t resz;
    cvx_float_t hresx;
    cvx_float_t hresy;
    cvx_float_t hresz;
    cvx_float_t cx;
    cvx_float_t by;
    cvx_float_t hz;
    cvx_float_t rt;
    cvx_float_t dres;
    cvx_float_t pres;
    cvx_float_t dinfres;
    cvx_float_t pinfres;
    cvx_float_t pcost;
    cvx_float_t dcost;
    cvx_float_t dres0;
    cvx_float_t pres0;

    // KKT solver for CPL and CP problems; chaning to proper KKT solver
    cvx_chainsolver_t cp_solver;

    // solution matrix
    cvx_matrix_t x, y, s, z;  // ok

    cvx_matrix_t c0;   // internal vector for non-linear target function
    cvx_matrix_t z_mnl;

    // f is point in domain of F (mln,1); Df gradient (mnl,n)
    // H is Hessian (n, n);
    cvx_matrix_t f, Df, H;
    cvx_matrix_t newf, newDf; //, H;

    cvx_matrix_t dx, dy, ds, dz;
    cvx_matrix_t dx0, dy0, ds0, dz0;
    cvx_matrix_t ds2, dz2, ds20, dz20;

    cvx_matrix_t x0, y0, s0, z0;
    cvx_matrix_t x1, y1, s1, z1;

    cvx_matrix_t rx, ry, rz;
    cvx_matrix_t rx0, ry0, rz0;
    cvx_matrix_t newx, newy, newz, news, newrx;
    cvx_matrix_t newrz0;
    cvx_matrix_t hrx, hry, hrz;
    cvx_matrix_t sigs, sigz;
    cvx_matrix_t lmbda, lmbdasq;
    cvx_matrix_t lmbda0;
    cvx_matrix_t lmbdasq0;

    cvx_matrix_t wx, wy, ws, wz;
    cvx_matrix_t wx2, wy2, ws2, wz2;
    cvx_matrix_t ws3, wz3;

    // matrix groups with correct indexing
    cvx_matgrp_t h_g;
    cvx_matgrp_t s_g, z_g;
    cvx_matgrp_t s0_g, z0_g;
    cvx_matgrp_t ds_g, dz_g;
    cvx_matgrp_t ds0_g, dz0_g;
    cvx_matgrp_t ds2_g, dz2_g;
    cvx_matgrp_t ds20_g, dz20_g;
    cvx_matgrp_t newz_g;
    cvx_matgrp_t newrz_g;
    cvx_matgrp_t news_g;
    cvx_matgrp_t rz_g;
    cvx_matgrp_t rzl_g, rznl_g;
    cvx_matgrp_t ws_g, wz_g;
    cvx_matgrp_t ws2_g, wz2_g;
    cvx_matgrp_t ws3_g, wz3_g;
    cvx_matgrp_t th_g;
    cvx_matgrp_t z1_g;
    cvx_matgrp_t lmbda_g, lmbdasq_g;
    cvx_matgrp_t sigs_g, sigz_g;

    cvx_index_t index_full;             // indexing to full matrix group
    cvx_index_t index_packed;           // indexing to matrix group with packed 'S' space
    cvx_index_t index_diag;             // indexing to matrix group with diagonal 'S' space
    cvx_index_t index_sig;              // indexing to matrix group with only diagonal 'S" space
    cvx_index_t index_cpt;              // indexing for G/h matrix with convex target function
    cvx_scaling_t W0;
    cvx_gp_program_t gp;                // Internal GP program structure.
} cvx_cpl_internal_t;

typedef struct cvx_problem {
    cvx_matrix_t *c;                    ///< Cost function coefficients
    cvx_matrix_t *G;                    ///< Inequality constraint coefficients
    cvx_matrix_t *h;                    ///< Inequality constraint limits
    cvx_matrix_t *A;                    ///< Equality constraint coefficients
    cvx_matrix_t *b;                    ///< Equality constraint limits
    // const cvx_dimset_t *dims;           ///< Problems dimensions
    cvx_kktsolver_t *solver;            ///< KKT solver

    // user defined starting points for CONELP solver
    cvx_matrix_t *primal_x;             ///< User defined starting point for primal
    cvx_matrix_t *primal_s;             ///< User defined starintg point for primal slacks
    cvx_matrix_t *dual_y;               ///< User defined starting point for dual
    cvx_matrix_t *dual_z;               ///< User defined starting point for dual slacks

    // Convex program for CP/CPL solver
    cvx_convex_program_t *F;

    // result
    cvx_matrix_t *x;
    cvx_matrix_t *y;
    cvx_matrix_t *s;
    cvx_matrix_t *z;

    // pointers to
    cvx_matrix_t *f;
    cvx_matrix_t *Df;
    cvx_matrix_t *H;

    cvx_index_t *index_g;
    int error;                          ///< Last error

    cvx_solution_t solution;            ///< Solution data

    // internal
    cvx_kktsolver_t __S;        // internal solver structure

    cvx_size_t cdim, cdim_diag;

    // statistics
    cvx_stats_t stats;
    // solver options
    cvx_solopts_t *solopts;

    cvx_scaling_t W;                    // scaling matrix group

    // solver internal variables
    union {
        cvx_conelp_internal_t conelp;
        cvx_cpl_internal_t cpl;
    } u;

    cvx_memblk_t work;                  // workspace
    cvx_size_t mlen;
    cvx_float_t *memory;

} cvx_problem_t;


extern cvx_float_t cvx_max_step(cvx_matgrp_t *x, cvx_matgrp_t *sigs, cvx_memblk_t *wrk);


extern cvx_size_t cvx_conelp_bytes(int n, int m, const cvx_dimset_t *dims);
extern cvx_size_t cvx_conelp_make(cvx_problem_t *prob,
                                  int n, int m, const cvx_dimset_t *dims,
                                  void *memory, cvx_size_t nbytes);

extern int cvx_conelp_isok(const cvx_matrix_t *c, const cvx_matrix_t *G, const cvx_matrix_t *h,
                           const cvx_matrix_t *A, const cvx_matrix_t *b, const cvx_dimset_t *dims);


extern cvx_size_t
cvx_conelp_setup(cvx_problem_t *prob,
                 cvx_matrix_t *c, cvx_matrix_t *G, cvx_matrix_t *h,
                 cvx_matrix_t *A, cvx_matrix_t *b, cvx_dimset_t *dims,
                 cvx_kktsolver_t *kktsolver);
extern void
cvx_conelp_set_start(cvx_problem_t *prob,
                     cvx_matrix_t *primal_x, cvx_matrix_t *primal_s,
                     cvx_matrix_t *dual_y, cvx_matrix_t *dual_z);
extern int
cvx_conelp_compute_start(cvx_problem_t *prob);

extern int
cvx_conelp_solve(cvx_problem_t *prob, cvx_solopts_t *opts);

extern int
cvx_cpl_compute_start(cvx_problem_t *cp);

extern int
cvx_cpl_solve(cvx_problem_t *cp, cvx_solopts_t *opts);

extern cvx_size_t
cvx_cpl_setup(cvx_problem_t *cp, cvx_convex_program_t *F,
              cvx_matrix_t *c, cvx_matrix_t *G, cvx_matrix_t *h, cvx_matrix_t *A,
              cvx_matrix_t *b,  const cvx_dimset_t *dims, cvx_kktsolver_t *kktsolver);

extern int cvx_cpl_isok(const cvx_matrix_t *c,
                        const cvx_convex_program_t *F,
                        const cvx_matrix_t *G,
                        const cvx_matrix_t *h,
                        const cvx_matrix_t *A,
                        const cvx_matrix_t *b,
                        const cvx_dimset_t *dims);

extern int cvx_cpl_setvars(cvx_problem_t *cp,
                           cvx_convex_program_t *F,
                           cvx_size_t n, cvx_size_t m,
                           cvx_matrix_t *c,
                           cvx_matrix_t *G,
                           cvx_matrix_t *h,
                           cvx_matrix_t *A,
                           cvx_matrix_t *b,
                           const cvx_dimset_t *dims,
                           cvx_kktsolver_t *kktsolver);

extern cvx_size_t cvx_cpl_allocate(cvx_problem_t *cp,
                                   int nl,
                                   cvx_size_t n,
                                   cvx_size_t m,
                                   cvx_size_t extra,
                                   const cvx_dimset_t *dims);
extern cvx_size_t cvx_cpl_make(cvx_problem_t *cp,
                               int n,
                               int m,
                               const cvx_dimset_t *dims,
                               int nl,
                               void *memory,
                               cvx_size_t nbytes);

extern cvx_size_t cvx_cpl_bytes(int n, int m, const cvx_dimset_t *dims, int nonlinear);


extern int
cvx_cp_isok(const cvx_convex_program_t *F,
            const cvx_matrix_t *G,
            const cvx_matrix_t *h,
            const cvx_matrix_t *A,
            const cvx_matrix_t *b,
            const cvx_dimset_t *dims);

extern int
cvx_cp_compute_start(cvx_problem_t *cp);

extern int
cvx_cp_solve(cvx_problem_t *cp, cvx_solopts_t *opts);

extern cvx_size_t
cvx_cp_setup(cvx_problem_t *cp, cvx_convex_program_t *F,
             cvx_matrix_t *G, cvx_matrix_t *h, cvx_matrix_t *A,
             cvx_matrix_t *b,  const cvx_dimset_t *dims, cvx_kktsolver_t *kktsolver);

extern int
cvx_cp_setvars(cvx_problem_t *cp, cvx_convex_program_t *F,
               cvx_size_t n, cvx_size_t m,
               cvx_matrix_t *G, cvx_matrix_t *h, cvx_matrix_t *A,
               cvx_matrix_t *b,  const cvx_dimset_t *dims, cvx_kktsolver_t *kktsolver);

extern cvx_size_t
cvx_cp_make(cvx_problem_t *cp, int n, int m,
            const cvx_dimset_t *dims, void *memory, cvx_size_t nbytes);

extern cvx_size_t
cvx_cp_bytes(int n, int m, const cvx_dimset_t *dims);


extern int
cvx_gp_setup(cvx_problem_t *cp, cvx_size_t *K,
             cvx_matrix_t *F, cvx_matrix_t *g,
             cvx_matrix_t *G, cvx_matrix_t *h, cvx_matrix_t *A,
             cvx_matrix_t *b,  cvx_dimset_t *dims, cvx_kktsolver_t *kktsolver);

extern int
cvx_gp_set_start(cvx_problem_t *cp,
                 cvx_matrix_t *x0,
                 cvx_matrix_t *s0,
                 cvx_matrix_t *y0,
                 cvx_matrix_t *z0);

extern int
cvx_gp_solve(cvx_problem_t *cp, cvx_solopts_t *opts);



#endif // __CVX_CONVEX_H

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:


