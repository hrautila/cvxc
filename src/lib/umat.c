/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */
#if ! defined(__STDC_VERSION__)
#  if defined(__CVXC_INLINE)
#    undef __CVXC_INLINE
#  endif
#  define __CVXC_INLINE
#endif

#include <ctype.h>
#include <errno.h>

#include "cvxc.h"

#if __STDC_VERSION__
// inline function reference
void cvxc_umat_init(cvxc_umatrix_t *G, struct cvxc_umatrix_funcs *vt, const void *user);
void cvxc_umat_clear(cvxc_umatrix_t *G);
int cvxc_umat_mvmult(
    cvxc_float_t beta, cvxc_matrix_t *y, cvxc_float_t alpha, const cvxc_umatrix_t *G, cvxc_matrix_t *x, int flags);
#endif /* __STRC_VERSION__ */

static
int stdmat_mvmult(
    const void *U, cvxc_float_t beta, cvxc_matrix_t *y, cvxc_float_t alpha, cvxc_matrix_t *x, int flags)
{
    return cvxm_mvmult(beta, y, alpha, (const cvxc_matrix_t *)U, x, flags);
}

/*
 * User matrix operator table for standard matrices.
 */
static struct cvxc_umatrix_funcs std = {
    .mvmult = stdmat_mvmult
};

/**
 * @brief Make user defined matrix out of given standard matrix.
 */
void cvxc_umat_make(cvxc_umatrix_t *U, const cvxc_matrix_t *A)
{
    cvxc_umat_init(U, &std, A);
}
