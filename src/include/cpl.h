
#ifndef __CVX_CPL_H
#define __CVX_CPL_H

#include "epi.h"
#include "convex.h"

#define MAX_RELAXED_ITERS 5
#ifndef __ZERO
#define __ZERO 0.0
#endif

typedef struct cvx_cpl_stats {
    cvx_float_t resx0;
    cvx_float_t resy0;
    cvx_float_t resz0;
    cvx_float_t resx;
    cvx_float_t resy;
    cvx_float_t resz;
    cvx_float_t resznl;
    cvx_float_t reszl;
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
} cvx_cpl_stats_t;

struct cvx_problem;

typedef int (*cvx_convex_program_t)(struct cvx_problem *cp,
                                    const cvx_matrix_t *x, const cvx_matrix_t *z);

#if 0
typedef struct cvx_convex_program_t cvx_convex_program_t;

struct cvx_convex_program_t {
    int (*F)(cvx_convex_program_t *cp, const cvx_matrix_t *f, const cvx_matrix_t *z);
    cvx_matrix_t *x0;
    cvx_matrix_t *f;
    cvx_matrix_t *Df;
    cvx_matrix_t *H;
};

static inline
int F0(cvx_convex_program_t *cp) {
    return cp->F(cp, (cvx_matrix_t *)0, (cvx_matrix_t *)0);
}

static inline
int F1(cvx_convex_program_t *cp, const cvx_matrix_t *x) {
    return cp->F(cp, x, (cvx_matrix_t *)0);
}

static inline
int F2(cvx_convex_program_t *cp, const cvx_matrix_t *x, const cvx_matrix_t *z) {
    return cp->F(cp, x, z);
}
#endif


extern cvx_size_t cvx_cpl_setup(cvx_problem_t *prob,
                                cvx_convex_program_t *F,
                                cvx_matrix_t *c,
                                cvx_matrix_t *G,
                                cvx_matrix_t *h,
                                cvx_matrix_t *A,
                                cvx_matrix_t *b,
                                cvx_dimset_t *dims,
                                cvx_kktsolver_t *kktsolver);


#endif // __CVX_CPL_H

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
