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
void cvxc_kktsolver_init(cvxc_kktsolver_t *kkt, cvxc_kktfuncs_t *vtable, void *private);
void cvxc_convex_program_init(cvxc_convex_program_t *cp, cvxc_convex_func F, void *user);

#endif
