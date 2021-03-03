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

#if ! defined(__CVXC_INLINE)
#  if __GNUC__
#    if ! __STDC_VERSION__
#      define __CVXC_INLINE extern inline
#    else
#      define __CVXC_INLINE inline
#      define __ARMAS_INLINE inline
#    endif
#  else
#    define __CVXC_INLINE inline
#  endif
#endif

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
    cvxc_float_t *memory;   ///< Array of floats
    cvxc_size_t mlen;       ///< Array length
    void *__bytes;          ///< Allocates bytes. If null then memory array provided externally
} cvxc_memblk_t;

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
    cvxc_matrix_t *g;
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

typedef int (*cvxc_convex_func)(cvxc_matrix_t *f, cvxc_matrix_t *Df, cvxc_matrix_t *H,
                             const cvxc_matrix_t *x, const cvxc_matrix_t *z, void *user);

typedef struct cvxc_convex_program {
    cvxc_convex_func F;
    void *user;
}  cvxc_convex_program_t;

typedef int (*cvxc_cp_initfunc)(cvxc_convex_program_t *, void *);

static inline
void cvxc_convex_program_init(cvxc_convex_program_t *cp, cvxc_convex_func F, void *user)
{
    cp->F = F;
    cp->user = user;
}


/* Forward declarations */
typedef struct cvxc_problem cvxc_problem_t;
typedef struct cvxc_kktsolver cvxc_kktsolver_t;

typedef struct cvxc_kktfuncs {
    int (*init)(cvxc_kktsolver_t *S, cvxc_problem_t *cp, int n, int m, const cvxc_dimset_t *dims);
    int (*factor)(cvxc_kktsolver_t *S, cvxc_scaling_t *W, cvxc_matrix_t *H, cvxc_matrix_t *Df);
    int (*solve)(cvxc_kktsolver_t *S, cvxc_matrix_t *x, cvxc_matrix_t *y, cvxc_matgrp_t *z_g);
    cvxc_size_t (*bytes)(int n, int m, const cvxc_dimset_t *dims);
    cvxc_size_t (*make)(cvxc_kktsolver_t *S, cvxc_problem_t *cp, int n, int m, const cvxc_dimset_t *dims, void *ptr, cvxc_size_t nbytes);
    void (*release)(cvxc_kktsolver_t *S);
} cvxc_kktfuncs_t;

typedef struct cvxc_kktsolver {
    cvxc_kktfuncs_t *vtable;
    void *private;
    cvxc_kktsolver_t *next;
} cvxc_kktsolver_t;


extern void cvxc_kktldl_load(cvxc_kktsolver_t *);

__CVXC_INLINE
int cvxc_kktinit(cvxc_kktsolver_t *S,
                cvxc_problem_t *cp,
                int n, int m, const cvxc_dimset_t *dims)
{
    return (*S->vtable->init)(S, cp, n, m, dims);
}

__CVXC_INLINE
int cvxc_kktfactor(cvxc_kktsolver_t *S,
                  cvxc_scaling_t *W, cvxc_matrix_t *H, cvxc_matrix_t *Df)
{
    return (*S->vtable->factor)(S, W, H, Df);
}

__CVXC_INLINE
int cvxc_kktsolve(cvxc_kktsolver_t *S,
                 cvxc_matrix_t *x, cvxc_matrix_t *y, cvxc_matgrp_t *z_g)
{
    return (*S->vtable->solve)(S, x, y, z_g);
}

__CVXC_INLINE
void cvxc_ldlsolver_init(cvxc_kktsolver_t *kkt,
                        cvxc_problem_t *cp,
                        int n, int m, const cvxc_dimset_t *dims)
{
    cvxc_kktldl_load(kkt);
    cvxc_kktinit(kkt, cp, n, m, dims);
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

typedef struct cvxc_gpindex {
    cvxc_size_t p;
    cvxc_size_t *index;
    void *__bytes;
} cvxc_gpindex_t;

extern int cvxc_gpi_setup(cvxc_gpindex_t *gpi, const cvxc_size_t *K, cvxc_size_t p);
extern int cvxc_gpi_init(cvxc_gpindex_t *gpi, const cvxc_size_t *K, cvxc_size_t p);
extern cvxc_size_t cvxc_gpi_make(cvxc_gpindex_t *gpi, cvxc_size_t p, void *memory, cvxc_size_t nbytes);
extern void cvxc_gpi_release(cvxc_gpindex_t *gpi);
extern void cvxc_gpi_free(cvxc_gpindex_t *gpi);
extern cvxc_size_t cvxc_gpi_elem(const cvxc_gpindex_t *gpi, cvxc_matrix_t *e, const cvxc_matrix_t *Fg, cvxc_size_t n);
extern cvxc_size_t cvxc_gpi_length(const cvxc_gpindex_t *gpi, cvxc_size_t n);

struct cvxc_conelp_internal;
struct cvxc_cpl_internal;

typedef struct cvxc_problem {
    cvxc_matrix_t *c;                    ///< Cost function coefficients
    cvxc_matrix_t *G;                    ///< Inequality constraint coefficients
    cvxc_matrix_t *h;                    ///< Inequality constraint limits
    cvxc_matrix_t *A;                    ///< Equality constraint coefficients
    cvxc_matrix_t *b;                    ///< Equality constraint limits
    cvxc_kktsolver_t *solver;            ///< KKT solver
    cvxc_convex_program_t *F;            ///< Convex program for CP/CPL solver

    cvxc_index_t *index_g;
    int error;                           ///< Last error

    // solver internal variables
    cvxc_kktsolver_t __S;
    cvxc_size_t cdim, cdim_diag;

    union {
        unsigned char *space;
        struct cvxc_conelp_internal *conelp;
        struct cvxc_cpl_internal *cpl;
    } u;

    cvxc_memblk_t *work;                  // workspace
    cvxc_size_t nbytes;
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
extern int
cvxc_conelp_compute_start_with(cvxc_problem_t *prob,
                               cvxc_matrix_t *primal_x, cvxc_matrix_t *primal_s,
                               cvxc_matrix_t *dual_y, cvxc_matrix_t *dual_z);
extern int
cvxc_conelp_compute_start(cvxc_problem_t *prob);

extern int
cvxc_conelp_solve(cvxc_solution_t *sol, cvxc_problem_t *prob, cvxc_solopts_t *opts);

extern int
cvxc_cpl_compute_start(cvxc_problem_t *cp);

extern int
cvxc_cpl_solve(cvxc_solution_t *sol, cvxc_problem_t *cp, cvxc_solopts_t *opts);

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
cvxc_cp_solve(cvxc_solution_t *sol, cvxc_problem_t *cp, cvxc_solopts_t *opts);

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
cvxc_gp_setup(cvxc_problem_t *cp, cvxc_gpindex_t *K,
             cvxc_matrix_t *F, cvxc_matrix_t *g,
             cvxc_matrix_t *G, cvxc_matrix_t *h, cvxc_matrix_t *A,
             cvxc_matrix_t *b, cvxc_kktsolver_t *kktsolver);

extern int
cvxc_gp_compute_start(cvxc_problem_t *cp);

extern int
cvxc_gp_solve(cvxc_solution_t *sol, cvxc_problem_t *cp, cvxc_solopts_t *opts);



#endif // __CVXC_CONVEX_H
