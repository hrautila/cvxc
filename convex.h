
#ifndef __CVX_CONVEX_H
#define __CVX_CONVEX_H 1

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <math.h>

#define SQRT sqrt
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
    L = 1,
    S = 2,
    Q = 4,
    NL = 8,
    CVXDIM_LINEAR = 1,
    CVXDIM_SDP = 2,
    CVXDIM_SOCP = 4,
    CVXDIM_NONLINEAR = 8
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
    CVXWS_D = 0,
    CVXWS_DI = 1,
    CVXWS_BETA = 2,
    CVXWS_V = 3,
    CVXWS_R = 4,
    CVXWS_RTI = 5,
    CVXWS_DNL = 6,
    CVXWS_DNLI = 7
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

typedef struct cvx_memblk {
    cvx_float_t *memory;
    cvx_size_t mlen;
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
    }
    
}

static inline
void __mblk_release(cvx_memblk_t *m)
{
    if (m->memory)
        free(m->memory);
    m->memory = (cvx_float_t *)0;
    m->mlen = 0;
}

static inline
cvx_float_t *__mblk_offset(cvx_memblk_t *m, cvx_size_t off)
{
    if (off >= m->mlen)
        return (cvx_float_t *)0;
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
    int *qdims;         ///< Second Order Cone constraint dimensions
    int *sdims;         ///< SDP constraint dimensions
    int qlen;           ///< Number of SOCP dimensions
    int slen;           ///< Number of SDP dimensions
    int ldim;           ///< Linear constraints
    int mnl;            ///< Non-linear constraints
} cvx_dimset_t;

extern cvx_dimset_t *cvx_dimset_alloc(cvx_dimset_t *dims, int linear, int *socp, int *sdp);
extern void cvx_dimset_release(cvx_dimset_t *dims);
extern cvx_size_t cvx_dimset_max(cvx_dimset_t *dims, cvx_dim_enum name);
extern cvx_size_t cvx_dimset_sum(cvx_dimset_t *dims, cvx_dim_enum name);
extern cvx_size_t cvx_dimset_sum_squared(cvx_dimset_t *dims, cvx_dim_enum name);
extern cvx_size_t cvx_dimset_sum_packed(cvx_dimset_t *dims, cvx_dim_enum name);

typedef struct cvx_index {
    cvx_size_t *index;         ///< Index data
    cvx_size_t *indq;          ///< SOCP vector
    cvx_size_t *inds;          ///< Indexes of SDP constraints
    cvx_size_t *indl;          ///< Index to linear
    cvx_size_t *indnl;         ///< Index to non-linear
    cvx_dimset_t *dims;
} cvx_index_t;

extern cvx_size_t cvx_index_count(cvx_index_t *index, cvx_dim_enum name);
extern cvx_index_t *cvx_index_init(cvx_index_t *index, cvx_dimset_t *dims, int packed);
extern void cvx_index_release(cvx_index_t *dims);
extern cvx_size_t cvx_index_elem(cvx_matrix_t *x, cvx_matrix_t *y, cvx_index_t *ind, cvx_dim_enum name, int k);
extern void cvx_index_create(cvx_matrix_t *x, cvx_index_t *index, cvx_dimset_t *dims, cvx_index_type kind);

typedef struct cvx_scaling {
    cvx_float_t *data;          // Scaling matrix data space, one big block
    cvx_size_t *indexes;        // Data space for V, R, RTI sizes and offsets
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
    cvx_dimset_t *dims;
    unsigned char *__bytes;
    cvx_size_t nbytes;
} cvx_scaling_t;

extern cvx_size_t cvx_scaling_elem(cvx_matrix_t *A, cvx_scaling_t *W, cvx_mset_enum name, int ind);
extern cvx_scaling_t *cvx_scaling_alloc(cvx_scaling_t *W, cvx_dimset_t *dims);
extern void cvx_scaling_release(cvx_scaling_t *W);


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

extern void cvx_mgrp_printf(FILE *f, const char *format, cvx_matgrp_t *g, const char *s);
extern void cvx_mat_printf(FILE *f, const char *format, cvx_matrix_t *g, const char *s);
extern void cvx_scaling_printf(FILE *f, const char *format, cvx_scaling_t *W, const char *s);

// -------------------------------------------------------------------------------------

// \brief hyperbolic ||x||_2 
extern cvx_float_t cvx_jnrm2(cvx_matrix_t *x);
// \brief hyperbolic x'*y
extern cvx_float_t cvx_jdot(cvx_matrix_t *x, cvx_matrix_t *y);

extern int cvx_pack(cvx_matrix_t *y, cvx_matrix_t *x, cvx_index_t *index);
extern int cvx_unpack(cvx_matrix_t *y, cvx_matrix_t *x, cvx_index_t *index);

extern int cvx_scale(cvx_matgrp_t *x, cvx_scaling_t *W, int flags, cvx_memblk_t *work);
extern int cvx_scale2(cvx_matgrp_t *x_g, cvx_matgrp_t *lmbda_g, int flags, cvx_memblk_t *work);
extern int cvx_sprod(cvx_matgrp_t *x, cvx_matgrp_t *y, int flags, cvx_memblk_t *wrk);
extern int cvx_ssqr(cvx_matgrp_t *x, cvx_matgrp_t *y);
extern int cvx_sinv(cvx_matgrp_t *x, cvx_matgrp_t *y, cvx_memblk_t *work);

extern int cvx_sgemv(cvx_float_t beta, cvx_matrix_t *y, cvx_float_t alpha, cvx_matrix_t *A, cvx_matgrp_t *x, 
                     int flags);
extern int cvx_triusc(cvx_matgrp_t *x);
extern int cvx_trisc(cvx_matgrp_t *x);
extern int cvx_mksymm(cvx_matgrp_t *x);
extern cvx_float_t cvx_snrm2(cvx_matgrp_t *x);
extern cvx_float_t cvx_sdot(cvx_matgrp_t *x, cvx_matgrp_t *y);

extern int
cvx_compute_scaling(cvx_scaling_t *W, cvx_matgrp_t *s_g, cvx_matgrp_t *z_g, cvx_matgrp_t *l_g, cvx_memblk_t *wrk);
extern int
cvx_update_scaling(cvx_scaling_t *W, cvx_matgrp_t *s_g, cvx_matgrp_t *z_g, cvx_matgrp_t *l_g, cvx_memblk_t *wrk);

// ------------------------------------------------------------------------------------------
// solver things

   
/* Forward declarations */
typedef struct cvx_conelp_problem cvx_conelp_problem_t;
typedef struct cvx_kktsolver_s cvx_kktsolver_t;

typedef struct cvx_kktfuncs_s {
    int (*factor)(cvx_kktsolver_t *S, cvx_scaling_t *W, cvx_matrix_t *H, cvx_matrix_t *Df);
    int (*solve)(cvx_kktsolver_t *S, cvx_matrix_t *x, cvx_matrix_t *y, cvx_matgrp_t *z_g);
} cvx_kktfuncs_t;
    
typedef struct cvx_ldlsolver_s {
    cvx_kktfuncs_t fnc;
    cvx_matrix_t K;
    cvx_matrix_t u;
    cvx_matrix_t g;
    cvx_memblk_t work;
    cvx_scaling_t *W;
    cvx_matrix_t *A;
    cvx_matrix_t *G;
    cvx_dimset_t *dims;
    size_t ldK;
    int *ipiv;
    size_t p;
    size_t n;
    size_t mnl;
} cvx_ldlsolver_t;


struct cvx_kktsolver_s {
    union {
        cvx_kktfuncs_t fnc;
        cvx_ldlsolver_t ldl;
    } u;
    cvx_conelp_problem_t *cp;
    int debug;
};

//extern void cvx_ldlsolver_init(cvx_kktsolver_t *, cvx_dimset_t *, cvx_matrix_t *G, cvx_matrix_t *A, int n);
extern void cvx_ldlsolver_init(cvx_kktsolver_t *kkt, cvx_conelp_problem_t *cp, cvx_dimset_t *dims, int mnl);

static inline
int cvx_kktfactor(cvx_kktsolver_t *S, cvx_scaling_t *W, cvx_matrix_t *H, cvx_matrix_t *Df)
{
    return (*S->u.fnc.factor)(S, W, H, Df);
}

static inline
int cvx_kktsolve(cvx_kktsolver_t *S, cvx_matrix_t *x, cvx_matrix_t *y, cvx_matgrp_t *z_g)
{
    return (*S->u.fnc.solve)(S, x, y, z_g);
}

typedef struct __solution_s {
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
        
typedef struct __solopts_s {
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
    cvx_float_t pcost;
    cvx_float_t dcost;
    cvx_float_t relgap;
} cvx_stats_t;

typedef struct cvx_conelp_problem {
    cvx_matrix_t *c;                    ///< Cost function coefficients
    cvx_matrix_t *G;                    ///< Inequality constraint coefficients
    cvx_matrix_t *h;                    ///< Inequality constraint limits
    cvx_matrix_t *A;                    ///< Equality constraint coefficients
    cvx_matrix_t *b;                    ///< Equality constraint limits
    cvx_dimset_t *dims;                 ///< Problems dimensions
    cvx_kktsolver_t *solver;            ///< KKT solver

    // user defined starting points
    cvx_matrix_t *primal_x;             ///< User defined starting point for primal
    cvx_matrix_t *primal_s;             ///< User defined starintg point for primal slacks
    cvx_matrix_t *dual_y;               ///< User defined starting point for dual
    cvx_matrix_t *dual_z;               ///< User defined starting point for dual slacks
    
    int error;                          ///< Last error
    
    cvx_solution_t solution;            ///< Solution data
    
    // internal 
    cvx_kktsolver_t __S;        // internal solver structure

    cvx_size_t cdim, cdim_diag;
    // solution matrix
    cvx_matrix_t x, y, s, z;

    cvx_float_t tau, kappa;
    
    cvx_matrix_t dx, dy, ds, dz;
    cvx_float_t dkappa, dtau;

    cvx_matrix_t x1, y1, z1;
    
    cvx_matrix_t rx, ry, rz, rs;
    cvx_matrix_t hrx, hry, hrz, hrs;
    cvx_matrix_t sigs, sigz;
    cvx_matrix_t lmbda, lmbdasq;
    
    cvx_matrix_t wx, wy, ws, wz;
    cvx_matrix_t wx2, wy2, ws2, wz2;
    cvx_matrix_t wx3, wy3, ws3, wz3;
    cvx_matrix_t th;
    cvx_float_t wkappa3;
    
    // matrix groups with correct indexing
    cvx_matgrp_t h_g;
    cvx_matgrp_t s_g, z_g;
    cvx_matgrp_t ds_g, dz_g;
    cvx_matgrp_t hrs_g, hrz_g;
    cvx_matgrp_t rs_g, rz_g;
    cvx_matgrp_t ws_g, wz_g;
    cvx_matgrp_t ws2_g, wz2_g;
    cvx_matgrp_t ws3_g, wz3_g;
    cvx_matgrp_t th_g;
    cvx_matgrp_t z1_g;
    cvx_matgrp_t lmbda_g, lmbdasq_g;
    cvx_matgrp_t sigs_g, sigz_g;
    
    // statistics
    cvx_stats_t stats;
    // solver options
    cvx_solopts_t *solopts;

    cvx_float_t nrms, nrmz;
    cvx_float_t dg, dgi;
    cvx_float_t ts, tz;                 // step sizes

    cvx_scaling_t W;                    // scaling matrix group

    cvx_index_t index_full;             // indexing to full matrix group
    cvx_index_t index_packed;           // indexing to matrix group with packed 'S' space
    cvx_index_t index_diag;             // indexing to matrix group with diagonal 'S' space
    cvx_index_t index_sig;              // indexing to matrix group with only diagonal 'S" space

    cvx_memblk_t work;                  // workspace
    
    cvx_size_t mlen;
    cvx_float_t *memory;
} cvx_conelp_problem_t;

extern cvx_conelp_problem_t *
cvx_conelp_setup(cvx_conelp_problem_t *prob,
                 cvx_matrix_t *c, cvx_matrix_t *G, cvx_matrix_t *h,
                 cvx_matrix_t *A, cvx_matrix_t *b, cvx_dimset_t *dims,
                 cvx_kktsolver_t *kktsolver);
extern void
cvx_conelp_set_start(cvx_conelp_problem_t *prob,
                     cvx_matrix_t *primal_x, cvx_matrix_t *primal_s,
                     cvx_matrix_t *dual_y, cvx_matrix_t *dual_z);
extern int
cvx_conelp_compute_start(cvx_conelp_problem_t *prob);

extern int
cvx_conelp_solve(cvx_conelp_problem_t *prob, cvx_solopts_t *opts);





#endif // __CVX_CONVEX_H

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:


