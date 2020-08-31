
#ifndef __CVXC_CPL_H
#define __CVXC_CPL_H

#include "epi.h"
#include "cvxc.h"

#define MAX_RELAXED_ITERS 5
#ifndef __ZERO
#define __ZERO 0.0
#endif

typedef struct cvxc_cpl_stats {
    cvxc_float_t resx0;
    cvxc_float_t resy0;
    cvxc_float_t resz0;
    cvxc_float_t resx;
    cvxc_float_t resy;
    cvxc_float_t resz;
    cvxc_float_t resznl;
    cvxc_float_t reszl;
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
    cvxc_float_t pcost;
    cvxc_float_t dcost;
    cvxc_float_t relgap;
} cvxc_cpl_stats_t;

struct cvxc_problem;

typedef int (*cvxc_convex_program_t)(struct cvxc_problem *cp,
                                    const cvxc_matrix_t *x, const cvxc_matrix_t *z);

#if 0
typedef struct cvxc_convex_program_t cvxc_convex_program_t;

struct cvxc_convex_program_t {
    int (*F)(cvxc_convex_program_t *cp, const cvxc_matrix_t *f, const cvxc_matrix_t *z);
    cvxc_matrix_t *x0;
    cvxc_matrix_t *f;
    cvxc_matrix_t *Df;
    cvxc_matrix_t *H;
};

static inline
int F0(cvxc_convex_program_t *cp) {
    return cp->F(cp, (cvxc_matrix_t *)0, (cvxc_matrix_t *)0);
}

static inline
int F1(cvxc_convex_program_t *cp, const cvxc_matrix_t *x) {
    return cp->F(cp, x, (cvxc_matrix_t *)0);
}

static inline
int F2(cvxc_convex_program_t *cp, const cvxc_matrix_t *x, const cvxc_matrix_t *z) {
    return cp->F(cp, x, z);
}
#endif


extern cvxc_size_t cvxc_cpl_setup(cvxc_problem_t *prob,
                                cvxc_convex_program_t *F,
                                cvxc_matrix_t *c,
                                cvxc_matrix_t *G,
                                cvxc_matrix_t *h,
                                cvxc_matrix_t *A,
                                cvxc_matrix_t *b,
                                cvxc_dimset_t *dims,
                                cvxc_kktsolver_t *kktsolver);


#endif // __CVXC_CPL_H

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
