
#ifndef __CVX_EPI_H
#define __CVX_EPI_H

#include "cvxm.h"

typedef struct cvx_epi_t {
    cvx_matrix_t *mat;
    cvx_int_t val;
} cvx_epi_t;

__CVX_INLINE
void cvx_epi_init(cvx_epi_t *e, cvx_matrix_t *m, cvx_int_t v)
{
    if (e) {
        e->mat = m;
        e->val = v;
    }
}

__CVX_INLINE
cvx_float_t cvx_epi_dot(const cvx_epi_t *a, cvx_epi_t *b)
{
    cvx_float_t r = cvxm_dot(a->mat, b->mat);
    r += a->val * b->val;
    return r;
}

__CVX_INLINE
void cvx_epi_scale(cvx_epi_t *e, cvx_float_t c)
{
    cvxm_scale(e->mat, c, 0);
    e->val *= c;
}

__CVX_INLINE
void cvx_epi_axpy(cvx_epi_t *y, cvx_float_t c, const cvx_epi_t *x)
{
    cvxm_axpy(y->mat, c, x->mat);
    y->val += c * x->val;
}

__CVX_INLINE
void cvx_epi_copy(cvx_epi_t *y, const cvx_epi_t *x)
{
    cvxm_copy(y->mat, x->mat, 0);
    y->val = x->val;
}



#endif // __CVX_EPI_H

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
