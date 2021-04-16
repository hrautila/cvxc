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
void cvxc_ldlsolver_init(cvxc_kktsolver_t *S,
                        cvxc_problem_t *cp, int n, int m, const cvxc_dimset_t *dims);

int cvxc_kktinit(cvxc_kktsolver_t *S,
                cvxc_problem_t *cp, int n, int m, const cvxc_dimset_t *dims);

int cvxc_kktfactor(cvxc_kktsolver_t *S,
                  cvxc_scaling_t *W, cvxc_matrix_t *H, cvxc_matrix_t *Df);
int cvxc_kktsolve(cvxc_kktsolver_t *S, cvxc_matrix_t *x, cvxc_matrix_t *y, cvxc_matgrp_t *z_g);
void cvxc_kktrelease(cvxc_kktsolver_t *kkt);

#endif
