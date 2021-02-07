/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */
#ifndef __CVXC_CONVEX_H
#define __CVXC_CONVEX_H 1

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

#if ENABLE_FLOAT32
/* float32 */
#define SQRT   sqrtf
#define POW    powf
#define MAXF   fmaxf
#define MINF   fminf


#else
/* float64 */
#define SQRT   sqrt
#define POW    pow
#define MAXF   fmax
#define MINF   fmin

#endif

#define SQRT2            1.41421356237309504880
#define CVXC_ABSTOL      1e-7
#define CVXC_RELTOL      1e-6
#define CVXC_FEASTOL     1e-7
#define CVXC_MAXITER     100
#define CVXC_STEP        0.99

#define CVXC_STAT_OPTIMAL                0
#define CVXC_STAT_UNKNOWN                1
#define CVXC_STAT_PRIMAL_INFEASIBLE      2
#define CVXC_STAT_DUAL_INFEASIBLE        3
#define CVXC_STAT_SINGULAR               4
#define CVXC_STAT_MAX                    4

#define CVXC_ERR_NULLCOST        1
#define CVXC_ERR_DIMC            2
#define CVXC_ERR_DIMH            3
#define CVXC_ERR_DIMG            4
#define CVXC_ERR_DIMA            5
#define CVXC_ERR_DIMB            6
#define CVXC_ERR_RANK            7
#define CVXC_ERR_MEMORY          8
#define CVXC_ERR_NEG_INITIAL_S   9
#define CVXC_ERR_NEG_INITIAL_Z   10
#define CVXC_ERR_MAXITER         11
#define CVXC_ERR_SINGULAR        12

#include "cvxm.h"

static inline
cvxc_float_t __NaN()
{
#ifdef NAN
    return NAN;
#else
    return SQRT(-1.0);
#endif
}

static inline
int __isnull(const void *ptr) {
    return ptr == (void *)0;
}

static inline
cvxc_size_t __aligned128(cvxc_size_t n)
{
    return (n & 0xF) != 0 ? 16 - (n & 0xF) : 0;
}

static inline
cvxc_size_t __aligned64(cvxc_size_t n)
{
    return (n & 0x7) != 0 ? 8 - (n & 0x7) : 0;
}

static inline
cvxc_float_t cvxc_maxvec(int n, const cvxc_float_t *vec)
{
    cvxc_float_t r = vec[0];
    for (int k = 1; k < n; k++) {
        r = MAXF(vec[k], r);
    }
    return r;
}

static inline
cvxc_float_t cvxc_minvec(int n, const cvxc_float_t *vec)
{
    cvxc_float_t r = vec[0];
    for (int k = 1; k < n; k++) {
        r = MINF(vec[k], r);
    }
    return r;
}

#define __MAX2(a, b)  __maxvec(2, (cvxc_float_t []){(a), (b)})
#define __MIN2(a, b)  __minvec(2, (cvxc_float_t []){(a), (b)})

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
} cvxc_dim_enum;

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
} cvxc_mset_enum;

typedef enum  {
    CVXC_KKTSOL_LDL = 1,
    CVXC_KKTSOL_LDL2 = 2,
    CVXC_KKTSOL_QR = 3,
    CVXC_KKTSOL_CHOL = 4,
    CVXC_KKTSOL_CHOL2 = 5
} cvxc_kktsolver_name_t;

typedef enum  {
    CVXC_INDEX_NORMAL = 0,
    CVXC_INDEX_PACKED = 1,
    CVXC_INDEX_DIAG = 2,
    CVXC_INDEX_SIGS = 3
} cvxc_index_type;

// -------------------------------------------------------------------------------------

// very simple
typedef struct cvxc_memblk {
    cvxc_float_t *memory;
    cvxc_size_t mlen;
    void *__bytes;
} cvxc_memblk_t;

static inline
void cvxc_mblk_empty(cvxc_memblk_t *m)
{
    if (m) {
        m->memory = (cvxc_float_t *)0;
        m->mlen = 0;
    }
}
static inline
void cvxc_mblk_init(cvxc_memblk_t *m, cvxc_size_t n)
{
    m->memory = (cvxc_float_t *)0;
    m->mlen = 0;
    cvxc_float_t *space = (cvxc_float_t *)calloc(n, sizeof(cvxc_float_t));
    if (space) {
        m->memory = space;
        m->mlen = n;
        m->__bytes = space;
    }
}

static inline
cvxc_size_t cvxc_mblk_make(cvxc_memblk_t *m, cvxc_size_t mlen, void *ptr, cvxc_size_t nbytes)
{
    if (!ptr || nbytes < mlen)
        return 0;
    m->memory = (cvxc_float_t *)ptr;
    m->mlen = mlen/sizeof(cvxc_float_t);
    m->__bytes = (void *)0;
    return mlen;
}

static inline
void cvxc_mblk_release(cvxc_memblk_t *m)
{
    if (m->__bytes)
        free(m->__bytes);
    m->memory = (cvxc_float_t *)0;
    m->mlen = 0;
    m->__bytes = (void *)0;
}

static inline
cvxc_float_t *cvxc_mblk_offset(cvxc_memblk_t *m, cvxc_size_t off)
{
    if (off >= m->mlen) {
        abort();
        return (cvxc_float_t *)0;
    }
    return &m->memory[off];
}

static inline
void cvxc_mblk_subblk(cvxc_memblk_t *d, cvxc_memblk_t *m, cvxc_size_t off)
{
    if (off >= m->mlen) {
        d->memory = (cvxc_float_t *)0;
        d->mlen = 0;
    } else {
        d->memory = &m->memory[off];
        d->mlen = m->mlen - off;
    }
}

static inline
void cvxc_mblk_clear(cvxc_memblk_t *m)
{
    if (m && m->mlen > 0)
        memset(m->memory, 0, m->mlen*sizeof(cvxc_float_t));
}

//----------------------------------------------------------------------------------------

typedef struct cvxc_dimset {
    cvxc_size_t *qdims;         ///< Second Order Cone constraint dimensions
    cvxc_size_t *sdims;         ///< SDP constraint dimensions
    cvxc_size_t qlen;           ///< Number of SOCP dimensions
    cvxc_size_t slen;           ///< Number of SDP dimensions
    cvxc_size_t ldim;           ///< Linear constraints
    cvxc_size_t mnl;            ///< Non-linear constraints
    cvxc_size_t iscpt;          ///< Target function is non-linear convex (CP)
} cvxc_dimset_t;

extern cvxc_dimset_t *cvxc_dimset_new(cvxc_size_t nonlinear, cvxc_size_t linear, cvxc_size_t socp, cvxc_size_t sdp);
extern void cvxc_dimset_free(cvxc_dimset_t *dims);
extern cvxc_size_t cvxc_dimset_count(const cvxc_dimset_t *dims, cvxc_dim_enum name);
extern cvxc_dimset_t *cvxc_dimset_create(cvxc_dimset_t *dims, cvxc_size_t nonlinear, cvxc_size_t linear, cvxc_size_t socp, cvxc_size_t sdp);
extern cvxc_dimset_t *cvxc_dimset_alloc(cvxc_dimset_t *dims, cvxc_size_t linear, const cvxc_size_t *socp, const cvxc_size_t *sdp);
extern void cvxc_dimset_release(cvxc_dimset_t *dims);
extern cvxc_size_t cvxc_dimset_max(const cvxc_dimset_t *dims, cvxc_dim_enum name);
extern cvxc_size_t cvxc_dimset_sum(const cvxc_dimset_t *dims, cvxc_dim_enum name);
extern cvxc_size_t cvxc_dimset_sum_squared(const cvxc_dimset_t *dims, cvxc_dim_enum name);
extern cvxc_size_t cvxc_dimset_sum_packed(const cvxc_dimset_t *dims, cvxc_dim_enum name);

typedef struct cvxc_index {
    cvxc_size_t *index;          ///< Index data
    cvxc_size_t *indq;           ///< SOCP vector
    cvxc_size_t *inds;           ///< Indexes of SDP constraints
    cvxc_size_t *indl;           ///< Index to linear
    cvxc_size_t *indnl;          ///< Index to non-linear
    cvxc_size_t *indnlt;         ///< Index to non-linear target
    cvxc_size_t indlen;          ///< Total number of  index elements
    cvxc_size_t qlen;            ///< Number of SOCP dimensions
    cvxc_size_t slen;            ///< Number of SDP dimensions
    cvxc_index_type type;        ///< Index type
    void *__bytes;
} cvxc_index_t;

extern cvxc_size_t cvxc_index_bytes(const cvxc_dimset_t *dims, int kind);
extern cvxc_size_t cvxc_index_count(const cvxc_index_t *index, cvxc_dim_enum name);
extern cvxc_size_t cvxc_index_make(cvxc_index_t *ind, const cvxc_dimset_t *dims, int kind, void *buf, cvxc_size_t nbytes);
extern cvxc_index_t *cvxc_index_init(cvxc_index_t *index, const cvxc_dimset_t *dims, int packed);
extern void cvxc_index_release(cvxc_index_t *dims);
extern cvxc_size_t cvxc_index_elem(cvxc_matrix_t *x, const cvxc_matrix_t *y, const cvxc_index_t *ind, cvxc_dim_enum name, int k);
extern cvxc_size_t cvxc_index_length(const cvxc_index_t *ind, cvxc_dim_enum name);

extern void cvxc_index_create(cvxc_matrix_t *x, cvxc_index_t *index, const cvxc_dimset_t *dims, cvxc_index_type kind);
extern void cvxc_subindex(cvxc_index_t *ind, const cvxc_index_t *src, int parts);

typedef struct cvxc_scaling {
    cvxc_float_t *data;          // Scaling matrix data space, one big block
    cvxc_size_t *indexes;        // Data space for V, R, RTI sizes and offsets
    cvxc_float_t *dnlt;          // Pointer to DNLT vector
    cvxc_float_t *dnlti;         // Pointer to DNLTI vector
    int dnltsz;                 // Length of DNLT/DNLTI vectors (0 or 1)
    cvxc_float_t *dnl;           // Pointer to DNL vector
    cvxc_float_t *dnli;          // Pointer to DNLI vector
    int dnlsz;                  // Length of DNL/DNLI vectors; equal to dims->mnl
    cvxc_float_t *d;             // Pointer to D vector
    cvxc_float_t *di;            // Pointer to DI vector
    int dsz;                    // Length of D/DI vectors; equal to dims->ldim
    cvxc_float_t *beta;          // Pointer to BETA vector
    cvxc_size_t *indv;           // Offsets to V vector;
    int vcount;                 // # of V vectors; equal to dims->qlen
    cvxc_size_t *indr;           // Offsets to R matrices
    cvxc_size_t *indrti;         // Offsets to RTI matrices
    int rcount;                 // # of R/RTI matrices, equal to dims->slen
    unsigned char *__bytes;
    cvxc_size_t nbytes;
} cvxc_scaling_t;

extern cvxc_size_t cvxc_scaling_bytes(cvxc_size_t *isize, const cvxc_dimset_t *dims);
extern cvxc_size_t cvxc_scaling_make(cvxc_scaling_t *W, const cvxc_dimset_t *dims, void *mem, size_t nbytes);

extern cvxc_size_t cvxc_scaling_elem(cvxc_matrix_t *A, const cvxc_scaling_t *W, cvxc_mset_enum name, int ind);
extern cvxc_scaling_t *cvxc_scaling_init(cvxc_scaling_t *W, const cvxc_dimset_t *dims);
extern void cvxc_scaling_initial_value(cvxc_scaling_t *W);
extern void cvxc_scaling_release(cvxc_scaling_t *W);
extern void cvxc_scaling_printf(FILE *f, const char *format, cvxc_scaling_t *W, const char *s);
extern void cvxc_scaling_elem_printf(FILE *f, const char *format,
                                    const cvxc_scaling_t *W, cvxc_mset_enum name, int index, const char *s);

extern int cvxc_scaling_copy(cvxc_scaling_t *W, const cvxc_scaling_t *Ws);


/*
 * Matrix group stored in a continuous memory matrix.
 */
typedef struct cvxc_matgrp {
    cvxc_matrix_t *mat;
    cvxc_index_t *index;
} cvxc_matgrp_t;

#define __nilgrp (cvxc_matgrp_t *)0

static inline
void cvxc_mgrp_init(cvxc_matgrp_t *g, cvxc_matrix_t *m, cvxc_index_t *ind)
{
    g->mat = m; g->index = ind;
}

static inline
cvxc_size_t cvxc_mgrp_elem(cvxc_matrix_t *x, cvxc_matgrp_t *g, cvxc_dim_enum name, int k)
{
    return cvxc_index_elem(x, g->mat, g->index, name, k);
}

static inline
int cvxc_mgrp_count(cvxc_matgrp_t *g, cvxc_dim_enum name)
{
    return cvxc_index_count(g->index, name);
}

extern void cvxc_mgrp_copy(cvxc_matgrp_t *x_g, cvxc_matgrp_t *y_g, cvxc_dim_enum name);
extern void cvxc_mgrp_copy_lambda(cvxc_matgrp_t *ds_g, cvxc_matgrp_t *lmbda);
extern void cvxc_mgrp_scale_sz(cvxc_matgrp_t *ds_g, cvxc_float_t v, int flags);
extern void cvxc_mgrp_update_sz(cvxc_matgrp_t *ds_g, cvxc_float_t v, int flags);
extern void cvxc_mgrp_axpby_sz(cvxc_float_t beta, cvxc_matgrp_t *y_g, cvxc_float_t alpha, cvxc_matgrp_t *x_g, int flags);
extern void cvxc_mgrp_initial_value(cvxc_matgrp_t *x_g, int flags);
extern void cvxc_mgrp_printf(FILE *f, const char *format, cvxc_matgrp_t *g, const char *s);
extern void cvxc_mat_printf(FILE *f, const char *format, cvxc_matrix_t *g, const char *s);
extern void cvxc_mat_test_nan(const char *name, const cvxc_matrix_t *A);
extern void cvxc_mat_print_ifenv(const char *name, const cvxc_matrix_t *A, const char *s);

// -------------------------------------------------------------------------------------

// \brief hyperbolic ||x||_2
extern cvxc_float_t cvxc_jnrm2(cvxc_matrix_t *x);
// \brief hyperbolic x'*y
extern cvxc_float_t cvxc_jdot(cvxc_matrix_t *x, cvxc_matrix_t *y);
extern cvxc_float_t cvxc_snrm2(cvxc_matgrp_t *x);
extern cvxc_float_t cvxc_sdot(cvxc_matgrp_t *x, cvxc_matgrp_t *y);
extern cvxc_float_t cvxc_sdot_elem(cvxc_matgrp_t *x, cvxc_matgrp_t *y, cvxc_dim_enum name);
extern cvxc_float_t cvxc_snrm2_elem(cvxc_matgrp_t *x, cvxc_dim_enum name);

extern void cvxc_pack(cvxc_matrix_t *y, cvxc_matrix_t *x, const cvxc_index_t *index);
extern void cvxc_unpack(cvxc_matrix_t *y, cvxc_matrix_t *x, const cvxc_index_t *index);
extern void cvxc_pack2(cvxc_matrix_t *x, const cvxc_index_t *index, cvxc_memblk_t *work);

extern int cvxc_scale(cvxc_matgrp_t *x, cvxc_scaling_t *W, int flags, cvxc_memblk_t *work);
extern int cvxc_scale_part(cvxc_matgrp_t *x, cvxc_scaling_t *W, int flags, int parts, cvxc_memblk_t *work);
extern int cvxc_scale2(cvxc_matgrp_t *x_g, cvxc_matgrp_t *lmbda_g, int flags, cvxc_memblk_t *work);
extern int cvxc_sprod(cvxc_matgrp_t *x, cvxc_matgrp_t *y, int flags, cvxc_memblk_t *wrk);
extern int cvxc_ssqr(cvxc_matgrp_t *x, cvxc_matgrp_t *y);
extern int cvxc_sinv(cvxc_matgrp_t *x, cvxc_matgrp_t *y, cvxc_memblk_t *work);

extern int cvxc_sgemv(cvxc_float_t beta, cvxc_matrix_t *y,
                     cvxc_float_t alpha, const cvxc_matrix_t *A, cvxc_matgrp_t *x, int flags);
extern int cvxc_sgemv2(cvxc_float_t beta, cvxc_matrix_t *y, cvxc_float_t alpha,
                      const cvxc_matrix_t *A, const cvxc_matrix_t *B, cvxc_matgrp_t *x, int flags);
extern int cvxc_triusc(cvxc_matgrp_t *x);
extern int cvxc_trisc(cvxc_matgrp_t *x);
extern int cvxc_mksymm(cvxc_matgrp_t *x);

extern int
cvxc_compute_scaling(cvxc_scaling_t *W, cvxc_matgrp_t *s_g, cvxc_matgrp_t *z_g, cvxc_matgrp_t *l_g, cvxc_memblk_t *wrk);
extern int
cvxc_update_scaling(cvxc_scaling_t *W, cvxc_matgrp_t *s_g, cvxc_matgrp_t *z_g, cvxc_matgrp_t *l_g, cvxc_memblk_t *wrk);

// ------------------------------------------------------------------------------------------
struct cvxc_gpindex;

typedef struct cvxc_solopts {
    cvxc_float_t abstol;         ///< Absolute tolerance
    cvxc_float_t reltol;         ///< Relative tolerance
    cvxc_float_t feastol;        ///< Feasibility tolerance
    int max_iter;               ///< Maximum iterations
    int debug;                  ///< Debug
    int refinement;             ///< Refinement count
    int show_progress;          ///< Show progress of the iteration
    int kkt_solver_name;        ///< KKT solver function
} cvxc_solopts_t;

typedef struct cvxc_params
{
    cvxc_dimset_t *dims;
    cvxc_matrix_t *c;
    cvxc_matrix_t *G;
    cvxc_matrix_t *h;
    cvxc_matrix_t *A;
    cvxc_matrix_t *b;
    cvxc_matrix_t *F;
    struct cvxc_gpindex *K;
    struct cvxc_solopts *opts;
    char *module;
    char *args;
} cvxc_params_t;

struct cvxc_solution;

extern int cvxc_json_matrix_read(cvxc_matrix_t **A, cvxc_stream_t *ios);
extern int cvxc_json_matrix_write(cvxc_stream_t *ios, const cvxc_matrix_t *A);
extern int cvxc_json_write_token(cvxc_stream_t *ios, int tok, const void *val, size_t len);
extern int cvxc_json_write_simple_token(cvxc_stream_t *ios, int tok);
extern int cvxc_json_read_token(char *iob, size_t len, cvxc_stream_t *ios);
extern int cvxc_json_dimset_write(cvxc_stream_t *ios, const cvxc_dimset_t *dims);
extern int cvxc_json_dimset_read(cvxc_dimset_t **dims, cvxc_stream_t *ios);
extern int cvxc_json_gpindex_read(struct cvxc_gpindex **gpi, cvxc_stream_t *ios);
extern int cvxc_json_gpindex_write(cvxc_stream_t *ios, const struct cvxc_gpindex *gpi);
extern int cvxc_json_intarray_read(cvxc_size_t *array, cvxc_size_t len, cvxc_stream_t *ios);
extern int cvxc_json_intarray_write(cvxc_stream_t *ios, const cvxc_size_t *array, cvxc_size_t len);
extern int cvxc_json_write_solution(cvxc_stream_t *ios, const struct cvxc_solution *sol);
extern int cvxc_json_write_result(cvxc_stream_t *ios, const struct cvxc_solution *sol);
extern int cvxc_json_read_params(cvxc_params_t **pars, cvxc_stream_t *ios);
extern int cvxc_json_write_params(cvxc_stream_t *ios, cvxc_params_t *pars);
extern int cvxc_json_read_options(cvxc_solopts_t **opts, cvxc_stream_t *ios);
extern int cvxc_json_write_options(cvxc_stream_t *ios, const cvxc_solopts_t *opts, const char *kkt);

extern void cvxc_file_stream(cvxc_stream_t *ios, FILE *fp);
extern void cvxc_str_stream(cvxc_stream_t *ios, const char *s, int len);

#define ONERR(exp) \
    do { if ((exp) < 0) return -1; } while (0)


// ------------------------------------------------------------------------------------------
// solver things

typedef int (*cvxc_convex_f0)(cvxc_matrix_t *x0, void *user);
typedef int (*cvxc_convex_f1)(cvxc_matrix_t *f, cvxc_matrix_t *Df, const cvxc_matrix_t *x, void *user);
typedef int (*cvxc_convex_f2)(cvxc_matrix_t *f, cvxc_matrix_t *Df, cvxc_matrix_t *H,
                             const cvxc_matrix_t *x, const cvxc_matrix_t *z, void *user);

typedef struct cvxc_convex_program {
    cvxc_convex_f2 F;
    void *user;
}  cvxc_convex_program_t;

typedef int (*cvxc_cp_initfunc)(cvxc_convex_program_t *, void *);

static inline
void cvxc_convex_program_init(cvxc_convex_program_t *cp, cvxc_convex_f2 F, void *user)
{
    cp->F = F;
    cp->user = user;
}

static inline
int F0(cvxc_convex_program_t *cp, cvxc_matrix_t *x0)
{
    return cp->F(x0, __cvxnil, __cvxnil, __cvxnil, __cvxnil, cp->user);
}

static inline
int F1(cvxc_convex_program_t *cp, cvxc_matrix_t *f, cvxc_matrix_t *Df, const cvxc_matrix_t *x)
{
    int e = cp->F(f, Df, __cvxnil, x, __cvxnil, cp->user);
    if (e == 0 && cvxm_isepi(x) && f) {
        cvxm_set(f, 0, 0, cvxm_get(f, 0, 0) - x->t);
    }
    return e;
}

static inline
int F2(cvxc_convex_program_t *cp, cvxc_matrix_t *f, cvxc_matrix_t *Df, cvxc_matrix_t *H,
       const cvxc_matrix_t *x, const cvxc_matrix_t *z)
{
    int e = cp->F(f, Df, H, x, z, cp->user);
    if (e == 0 && cvxm_isepi(x) && f) {
        cvxm_set(f, 0, 0, cvxm_get(f, 0, 0) - x->t);
    }
    return e;
}


/* Forward declarations */
typedef struct cvxc_problem cvxc_problem_t;
typedef struct cvxc_kktsolver cvxc_kktsolver_t;

typedef struct cvxc_kktfuncs {
    cvxc_kktsolver_t *(*new)(cvxc_problem_t *cp, int n, int m, const cvxc_dimset_t *dims);
    int (*init)(cvxc_kktsolver_t *S, cvxc_problem_t *cp, int n, int m, const cvxc_dimset_t *dims);
    int (*factor)(cvxc_kktsolver_t *S, cvxc_scaling_t *W, cvxc_matrix_t *H, cvxc_matrix_t *Df);
    int (*solve)(cvxc_kktsolver_t *S, cvxc_matrix_t *x, cvxc_matrix_t *y, cvxc_matgrp_t *z_g);
    cvxc_size_t (*bytes)(int n, int m, const cvxc_dimset_t *dims);
    cvxc_size_t (*make)(cvxc_kktsolver_t *S, cvxc_problem_t *cp, int n, int m, const cvxc_dimset_t *dims, void *ptr, cvxc_size_t nbytes);
    void (*free)(cvxc_kktsolver_t *S);
} cvxc_kktfuncs_t;

#if 0
typedef struct cvxc_kktsolver {
    cvxc_kktfuncs_t *fnc;
    cvxc_problem_t *cp;
} cvxc_kktsolver_t;
#endif

typedef struct cvxc_ldlsolver {
    cvxc_kktfuncs_t *fnc;
    cvxc_problem_t *cp;
    cvxc_matrix_t K;
    cvxc_matrix_t u;
    cvxc_matrix_t g;
    cvxc_memblk_t work;
    cvxc_scaling_t *W;
    cvxc_matrix_t *A;
    cvxc_matrix_t *G;
    cvxc_matrix_t *Df;
    const cvxc_dimset_t *dims;
    size_t ldK;
    int *ipiv;
    size_t p;
    size_t n;
    size_t mnl;
} cvxc_ldlsolver_t;

typedef struct cvxc_qrsolver {
    cvxc_kktfuncs_t *fnc;
    cvxc_matrix_t QA;
    cvxc_matrix_t tauA;
    cvxc_matrix_t Gs;
    cvxc_matrix_t tauG;
    cvxc_matrix_t u;
    cvxc_matrix_t vv;
    cvxc_matrix_t w;
    cvxc_memblk_t work;
    cvxc_scaling_t *W;
    cvxc_matrix_t *A;
    cvxc_matrix_t *G;
    cvxc_problem_t *cp;
    const cvxc_dimset_t *dims;
    size_t p;
    size_t n;
} cvxc_qrsolver_t;

typedef struct cvxc_chainsolver {
    cvxc_kktfuncs_t *fnc;
    cvxc_problem_t *cp;
    cvxc_kktsolver_t *next;
} cvxc_chainsolver_t;

struct cvxc_kktsolver {
    union {
        cvxc_kktfuncs_t *fnc;
        cvxc_ldlsolver_t ldl;
        cvxc_chainsolver_t next;
    } u;
    int debug;
};

extern cvxc_kktfuncs_t *cvxc_ldlload(void *ptr);

__CVXC_INLINE
void cvxc_ldlsolver_init(cvxc_kktsolver_t *kkt,
                        cvxc_problem_t *cp,
                        int n, int m, const cvxc_dimset_t *dims)
{
    cvxc_kktfuncs_t *fnc = cvxc_ldlload((void *)0);
    fnc->init(kkt, cp, n, m, dims);
}

__CVXC_INLINE
int cvxc_kktinit(cvxc_kktsolver_t *S,
                cvxc_problem_t *cp,
                int n, int m, const cvxc_dimset_t *dims)
{
    return (*S->u.fnc->init)(S, cp, n, m, dims);
}

__CVXC_INLINE
int cvxc_kktfactor(cvxc_kktsolver_t *S,
                  cvxc_scaling_t *W, cvxc_matrix_t *H, cvxc_matrix_t *Df)
{
    return (*S->u.fnc->factor)(S, W, H, Df);
}

__CVXC_INLINE
int cvxc_kktsolve(cvxc_kktsolver_t *S,
                 cvxc_matrix_t *x, cvxc_matrix_t *y, cvxc_matgrp_t *z_g)
{
    return (*S->u.fnc->solve)(S, x, y, z_g);
}

typedef struct cvxc_solution {
    int status;
    cvxc_matrix_t *x;
    cvxc_matrix_t *s;
    cvxc_matrix_t *y;
    cvxc_matrix_t *z;
    cvxc_float_t primal_objective;
    cvxc_float_t dual_objective;
    cvxc_float_t gap;
    cvxc_float_t relative_gap;
    cvxc_float_t primal_infeasibility;
    cvxc_float_t dual_infeasibility;
    cvxc_float_t primal_slack;
    cvxc_float_t dual_slack;
    cvxc_float_t primal_residual_cert;
    cvxc_float_t dual_residual_cert;
    int iterations;
} cvxc_solution_t;

typedef struct cvxc_stats {
    cvxc_float_t resx0;
    cvxc_float_t resy0;
    cvxc_float_t resz0;
    cvxc_float_t resx;
    cvxc_float_t resy;
    cvxc_float_t resz;
    cvxc_float_t hresx;
    cvxc_float_t hresy;
    cvxc_float_t hresz;
    cvxc_float_t cx;
    cvxc_float_t by;
    cvxc_float_t hz;
    cvxc_float_t rt;
    cvxc_float_t dres;
    cvxc_float_t pres;
    cvxc_float_t dinfres;
    cvxc_float_t pinfres;
    cvxc_float_t gap;
    cvxc_float_t relgap;
    cvxc_float_t pcost;
    cvxc_float_t dcost;
} cvxc_stats_t;


typedef struct cvxc_conelp_internal {
    cvxc_float_t tau, kappa;
    cvxc_float_t dkappa, dtau;
    cvxc_float_t wkappa3;
    cvxc_float_t nrms, nrmz;
    cvxc_float_t dg, dgi;
    cvxc_float_t ts, tz;                 // step sizes

    // internal result
    cvxc_matrix_t x, y, s, z;

    cvxc_matrix_t dx, dy, ds, dz;
    cvxc_matrix_t x1, y1, z1;
    cvxc_matrix_t rx, ry, rz;
    cvxc_matrix_t hrx, hry, hrz;
    cvxc_matrix_t sigs, sigz;
    cvxc_matrix_t lmbda, lmbdasq;

    cvxc_matrix_t wx, wy, ws, wz;
    cvxc_matrix_t wx2, wy2, ws2, wz2;
    cvxc_matrix_t ws3, wz3;
    cvxc_matrix_t th;

    // matrix groups with correct indexing
    cvxc_matgrp_t h_g;
    cvxc_matgrp_t s_g, z_g;
    cvxc_matgrp_t ds_g, dz_g;
    cvxc_matgrp_t hrz_g;
    cvxc_matgrp_t rz_g;
    cvxc_matgrp_t ws_g, wz_g;
    cvxc_matgrp_t ws2_g, wz2_g;
    cvxc_matgrp_t ws3_g, wz3_g;
    cvxc_matgrp_t th_g;
    cvxc_matgrp_t z1_g;
    cvxc_matgrp_t lmbda_g, lmbdasq_g;
    cvxc_matgrp_t sigs_g, sigz_g;

    cvxc_index_t index_full;             // indexing to full matrix group
    cvxc_index_t index_packed;           // indexing to matrix group with packed 'S' space
    cvxc_index_t index_diag;             // indexing to matrix group with diagonal 'S' space
    cvxc_index_t index_sig;              // indexing to matrix group with only diagonal 'S" space

} cvxc_conelp_internal_t;

typedef struct cvxc_gpindex {
    cvxc_size_t p;
    cvxc_size_t *index;
    void *__bytes;
} cvxc_gpindex_t;

typedef struct cvxc_gp_params {
    cvxc_matrix_t *F;
    cvxc_matrix_t *g;
    cvxc_matrix_t y;
    cvxc_matrix_t u;
    cvxc_matrix_t Fs;
    cvxc_gpindex_t gpi;
} cvxc_gp_params_t;

typedef struct cvxc_gp_program {
    cvxc_convex_program_t gp;
    cvxc_gp_params_t gp_params;
} cvxc_gp_program_t;



typedef struct cvxc_cpl_internal {
    cvxc_float_t tau, kappa;
    cvxc_float_t dkappa, dtau;
    cvxc_float_t wkappa3;
    cvxc_float_t nrms, nrmz;
    cvxc_float_t dg, dgi;
    cvxc_float_t ts, tz;                 // step sizes
    cvxc_float_t phi;
    cvxc_float_t dphi;
    cvxc_float_t phi0;
    cvxc_float_t dphi0;
    cvxc_float_t mu;
    cvxc_float_t sigma;
    cvxc_float_t sigma0;
    cvxc_float_t eta;
    cvxc_float_t eta0;
    cvxc_float_t theta1;
    cvxc_float_t theta2;
    cvxc_float_t theta3;
    cvxc_float_t step;
    cvxc_float_t step0;
    cvxc_float_t gap;
    cvxc_float_t gap0;
    cvxc_float_t dsdz;
    cvxc_float_t dsdz0;

    cvxc_float_t relgap;
    cvxc_float_t resznl;
    cvxc_float_t reszl;
    cvxc_float_t resrznl;
    cvxc_float_t resrzl;
    cvxc_float_t relgap0;
    cvxc_float_t resznl0;
    cvxc_float_t reszl0;
    cvxc_float_t resrznl0;
    cvxc_float_t resrzl0;

    cvxc_float_t resx0;
    cvxc_float_t resy0;
    cvxc_float_t resz0;
    cvxc_float_t resx;
    cvxc_float_t resy;
    cvxc_float_t resz;
    cvxc_float_t hresx;
    cvxc_float_t hresy;
    cvxc_float_t hresz;
    cvxc_float_t cx;
    cvxc_float_t by;
    cvxc_float_t hz;
    cvxc_float_t rt;
    cvxc_float_t dres;
    cvxc_float_t pres;
    cvxc_float_t dinfres;
    cvxc_float_t pinfres;
    cvxc_float_t pcost;
    cvxc_float_t dcost;
    cvxc_float_t dres0;
    cvxc_float_t pres0;

    // KKT solver for CPL and CP problems; chaning to proper KKT solver
    cvxc_chainsolver_t cp_solver;

    // solution matrix
    cvxc_matrix_t x, y, s, z;  // ok

    cvxc_matrix_t c0;   // internal vector for non-linear target function
    cvxc_matrix_t z_mnl;

    // f is point in domain of F (mln,1); Df gradient (mnl,n)
    // H is Hessian (n, n);
    cvxc_matrix_t f, Df, H;
    cvxc_matrix_t newf, newDf; //, H;

    cvxc_matrix_t dx, dy, ds, dz;
    cvxc_matrix_t dx0, dy0, ds0, dz0;
    cvxc_matrix_t ds2, dz2, ds20, dz20;

    cvxc_matrix_t x0, y0, s0, z0;
    cvxc_matrix_t x1, y1, s1, z1;

    cvxc_matrix_t rx, ry, rz;
    cvxc_matrix_t rx0, ry0, rz0;
    cvxc_matrix_t newx, newy, newz, news, newrx;
    cvxc_matrix_t newrz0;
    cvxc_matrix_t hrx, hry, hrz;
    cvxc_matrix_t sigs, sigz;
    cvxc_matrix_t lmbda, lmbdasq;
    cvxc_matrix_t lmbda0;
    cvxc_matrix_t lmbdasq0;

    cvxc_matrix_t wx, wy, ws, wz;
    cvxc_matrix_t wx2, wy2, ws2, wz2;
    cvxc_matrix_t ws3, wz3;

    // matrix groups with correct indexing
    cvxc_matgrp_t h_g;
    cvxc_matgrp_t s_g, z_g;
    cvxc_matgrp_t s0_g, z0_g;
    cvxc_matgrp_t ds_g, dz_g;
    cvxc_matgrp_t ds0_g, dz0_g;
    cvxc_matgrp_t ds2_g, dz2_g;
    cvxc_matgrp_t ds20_g, dz20_g;
    cvxc_matgrp_t newz_g;
    cvxc_matgrp_t newrz_g;
    cvxc_matgrp_t news_g;
    cvxc_matgrp_t rz_g;
    cvxc_matgrp_t rzl_g, rznl_g;
    cvxc_matgrp_t ws_g, wz_g;
    cvxc_matgrp_t ws2_g, wz2_g;
    cvxc_matgrp_t ws3_g, wz3_g;
    cvxc_matgrp_t th_g;
    cvxc_matgrp_t z1_g;
    cvxc_matgrp_t lmbda_g, lmbdasq_g;
    cvxc_matgrp_t sigs_g, sigz_g;

    cvxc_index_t index_full;             // indexing to full matrix group
    cvxc_index_t index_packed;           // indexing to matrix group with packed 'S' space
    cvxc_index_t index_diag;             // indexing to matrix group with diagonal 'S' space
    cvxc_index_t index_sig;              // indexing to matrix group with only diagonal 'S" space
    cvxc_index_t index_cpt;              // indexing for G/h matrix with convex target function
    cvxc_scaling_t W0;
    cvxc_gp_program_t gp;                // Internal GP program structure.
} cvxc_cpl_internal_t;

typedef struct cvxc_problem {
    cvxc_matrix_t *c;                    ///< Cost function coefficients
    cvxc_matrix_t *G;                    ///< Inequality constraint coefficients
    cvxc_matrix_t *h;                    ///< Inequality constraint limits
    cvxc_matrix_t *A;                    ///< Equality constraint coefficients
    cvxc_matrix_t *b;                    ///< Equality constraint limits
    // const cvxc_dimset_t *dims;           ///< Problems dimensions
    cvxc_kktsolver_t *solver;            ///< KKT solver

    // user defined starting points for CONELP solver
    cvxc_matrix_t *primal_x;             ///< User defined starting point for primal
    cvxc_matrix_t *primal_s;             ///< User defined starintg point for primal slacks
    cvxc_matrix_t *dual_y;               ///< User defined starting point for dual
    cvxc_matrix_t *dual_z;               ///< User defined starting point for dual slacks

    // Convex program for CP/CPL solver
    cvxc_convex_program_t *F;

    // result
    cvxc_matrix_t *x;
    cvxc_matrix_t *y;
    cvxc_matrix_t *s;
    cvxc_matrix_t *z;

    // pointers to
    cvxc_matrix_t *f;
    cvxc_matrix_t *Df;
    cvxc_matrix_t *H;

    cvxc_index_t *index_g;
    int error;                          ///< Last error

    cvxc_solution_t solution;            ///< Solution data

    // internal
    cvxc_kktsolver_t __S;        // internal solver structure

    cvxc_size_t cdim, cdim_diag;

    // statistics
    cvxc_stats_t stats;
    // solver options
    cvxc_solopts_t *solopts;

    cvxc_scaling_t W;                    // scaling matrix group

    // solver internal variables
    union {
        cvxc_conelp_internal_t conelp;
        cvxc_cpl_internal_t cpl;
    } u;

    cvxc_memblk_t work;                  // workspace
    cvxc_size_t mlen;
    cvxc_float_t *memory;

} cvxc_problem_t;


extern cvxc_float_t cvxc_max_step(cvxc_matgrp_t *x, cvxc_matgrp_t *sigs, cvxc_memblk_t *wrk);


extern cvxc_size_t cvxc_conelp_bytes(int n, int m, const cvxc_dimset_t *dims);
extern cvxc_size_t cvxc_conelp_make(cvxc_problem_t *prob,
                                  int n, int m, const cvxc_dimset_t *dims,
                                  void *memory, cvxc_size_t nbytes);

extern int cvxc_conelp_isok(const cvxc_matrix_t *c, const cvxc_matrix_t *G, const cvxc_matrix_t *h,
                           const cvxc_matrix_t *A, const cvxc_matrix_t *b, const cvxc_dimset_t *dims);


extern cvxc_size_t
cvxc_conelp_setup(cvxc_problem_t *prob,
                 cvxc_matrix_t *c, cvxc_matrix_t *G, cvxc_matrix_t *h,
                 cvxc_matrix_t *A, cvxc_matrix_t *b, cvxc_dimset_t *dims,
                 cvxc_kktsolver_t *kktsolver);
extern void
cvxc_conelp_set_start(cvxc_problem_t *prob,
                     cvxc_matrix_t *primal_x, cvxc_matrix_t *primal_s,
                     cvxc_matrix_t *dual_y, cvxc_matrix_t *dual_z);
extern int
cvxc_conelp_compute_start(cvxc_problem_t *prob);

extern int
cvxc_conelp_solve(cvxc_problem_t *prob, cvxc_solopts_t *opts);

extern int
cvxc_cpl_compute_start(cvxc_problem_t *cp);

extern int
cvxc_cpl_solve(cvxc_problem_t *cp, cvxc_solopts_t *opts);

extern cvxc_size_t
cvxc_cpl_setup(cvxc_problem_t *cp, cvxc_convex_program_t *F,
              cvxc_matrix_t *c, cvxc_matrix_t *G, cvxc_matrix_t *h, cvxc_matrix_t *A,
              cvxc_matrix_t *b,  const cvxc_dimset_t *dims, cvxc_kktsolver_t *kktsolver);

extern int cvxc_cpl_isok(const cvxc_matrix_t *c,
                        const cvxc_convex_program_t *F,
                        const cvxc_matrix_t *G,
                        const cvxc_matrix_t *h,
                        const cvxc_matrix_t *A,
                        const cvxc_matrix_t *b,
                        const cvxc_dimset_t *dims);

extern int cvxc_cpl_setvars(cvxc_problem_t *cp,
                           cvxc_convex_program_t *F,
                           cvxc_size_t n, cvxc_size_t m,
                           cvxc_matrix_t *c,
                           cvxc_matrix_t *G,
                           cvxc_matrix_t *h,
                           cvxc_matrix_t *A,
                           cvxc_matrix_t *b,
                           const cvxc_dimset_t *dims,
                           cvxc_kktsolver_t *kktsolver);

extern cvxc_size_t cvxc_cpl_allocate(cvxc_problem_t *cp,
                                   int nl,
                                   cvxc_size_t n,
                                   cvxc_size_t m,
                                   cvxc_size_t extra,
                                   const cvxc_dimset_t *dims);
extern cvxc_size_t cvxc_cpl_make(cvxc_problem_t *cp,
                               int n,
                               int m,
                               const cvxc_dimset_t *dims,
                               int nl,
                               void *memory,
                               cvxc_size_t nbytes);

extern cvxc_size_t cvxc_cpl_bytes(int n, int m, const cvxc_dimset_t *dims, int nonlinear);


extern int
cvxc_cp_isok(const cvxc_convex_program_t *F,
            const cvxc_matrix_t *G,
            const cvxc_matrix_t *h,
            const cvxc_matrix_t *A,
            const cvxc_matrix_t *b,
            const cvxc_dimset_t *dims);

extern int
cvxc_cp_compute_start(cvxc_problem_t *cp);

extern int
cvxc_cp_solve(cvxc_problem_t *cp, cvxc_solopts_t *opts);

extern cvxc_size_t
cvxc_cp_setup(cvxc_problem_t *cp, cvxc_convex_program_t *F,
             cvxc_matrix_t *G, cvxc_matrix_t *h, cvxc_matrix_t *A,
             cvxc_matrix_t *b,  const cvxc_dimset_t *dims, cvxc_kktsolver_t *kktsolver);

extern int
cvxc_cp_setvars(cvxc_problem_t *cp, cvxc_convex_program_t *F,
               cvxc_size_t n, cvxc_size_t m,
               cvxc_matrix_t *G, cvxc_matrix_t *h, cvxc_matrix_t *A,
               cvxc_matrix_t *b,  const cvxc_dimset_t *dims, cvxc_kktsolver_t *kktsolver);

extern cvxc_size_t
cvxc_cp_make(cvxc_problem_t *cp, int n, int m,
            const cvxc_dimset_t *dims, void *memory, cvxc_size_t nbytes);

extern cvxc_size_t
cvxc_cp_bytes(int n, int m, const cvxc_dimset_t *dims);


extern int
cvxc_gp_setup(cvxc_problem_t *cp, cvxc_size_t *K,
             cvxc_matrix_t *F, cvxc_matrix_t *g,
             cvxc_matrix_t *G, cvxc_matrix_t *h, cvxc_matrix_t *A,
             cvxc_matrix_t *b, cvxc_kktsolver_t *kktsolver);

extern int
cvxc_gp_compute_start(cvxc_problem_t *cp);

extern int
cvxc_gp_solve(cvxc_problem_t *cp, cvxc_solopts_t *opts);



#endif // __CVXC_CONVEX_H
