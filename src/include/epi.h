/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */

#ifndef __CVXC_EPI_H
#define __CVXC_EPI_H

#include "cvxm.h"

typedef struct cvxc_epi_t {
    cvxc_matrix_t *mat;
    cvxc_int_t val;
} cvxc_epi_t;

__CVXC_INLINE
void cvxc_epi_init(cvxc_epi_t *e, cvxc_matrix_t *m, cvxc_int_t v)
{
    if (e) {
        e->mat = m;
        e->val = v;
    }
}

__CVXC_INLINE
cvxc_float_t cvxc_epi_dot(const cvxc_epi_t *a, cvxc_epi_t *b)
{
    cvxc_float_t r = cvxm_dot(a->mat, b->mat);
    r += a->val * b->val;
    return r;
}

__CVXC_INLINE
void cvxc_epi_scale(cvxc_epi_t *e, cvxc_float_t c)
{
    cvxm_scale(e->mat, c, 0);
    e->val *= c;
}

__CVXC_INLINE
void cvxc_epi_axpy(cvxc_epi_t *y, cvxc_float_t c, const cvxc_epi_t *x)
{
    cvxm_axpy(y->mat, c, x->mat);
    y->val += c * x->val;
}

__CVXC_INLINE
void cvxc_epi_copy(cvxc_epi_t *y, const cvxc_epi_t *x)
{
    cvxm_copy(y->mat, x->mat, 0);
    y->val = x->val;
}

#endif // __CVXC_EPI_H
