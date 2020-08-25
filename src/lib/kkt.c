
// Copyright: Harri Rautila, 2018 <harri.rautila@gmail.com>

#if ! defined(__STDC_VERSION__)
#  if defined(__CVX_INLINE)
#    undef __CVX_INLINE
#  endif
#  define __CVX_INLINE
#endif

#include <ctype.h>
#include <errno.h>

#include "convex.h"

#if __STDC_VERSION__
// inline function reference
void cvx_ldlsolver_init(cvx_kktsolver_t *S,
                        cvx_problem_t *cp, int n, int m, const cvx_dimset_t *dims);

int cvx_kktinit(cvx_kktsolver_t *S,
                cvx_problem_t *cp, int n, int m, const cvx_dimset_t *dims);

int cvx_kktfactor(cvx_kktsolver_t *S,
                  cvx_scaling_t *W, cvx_matrix_t *H, cvx_matrix_t *Df);
int cvx_kktsolve(cvx_kktsolver_t *S, cvx_matrix_t *x, cvx_matrix_t *y, cvx_matgrp_t *z_g);
#endif
