
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include "cvxc.h"

// forward declarations
static
int chol2_init(cvxc_kktsolver_t *S,
             cvxc_conelp_problem_t *cp,
             int n,
             int m,
             const cvxc_dimset_t *dims);

static
int chol2_factor(cvxc_kktsolver_t *S,
               cvxc_scaling_t *W,
               cvxc_matrix_t *H,
               cvxc_matrix_t *Df);
static
int chol2_solve(cvxc_kktsolver_t *S,
              cvxc_matrix_t *x,
              cvxc_matrix_t *y,
              cvxc_matgrp_t *z_g);

static
cvxc_size_t chol2_bytes(int n, int m, const cvxc_dimset_t *dims);

static
cvxc_size_t chol2_make(cvxc_kktsolver_t *kkt,
                    cvxc_conelp_problem_t *cp,
                    int n,
                    int m,
                    const cvxc_dimset_t *dims,
                    void *mem,
                    cvxc_size_t nbytes);

static
cvxc_kktsolver_t *chol2_new(cvxc_conelp_problem_t *cp,
                         int n,
                         int m,
                         const cvxc_dimset_t *dims);

static
void chol2_free(cvxc_kktsolver_t *S);

// function table
static cvxc_kktfuncs_t chol2_functions = {
    .new    = chol2_new,
    .factor = chol2_factor,
    .solve  = chol2_solve,
    .init   = chol2_init,
    .bytes  = chol2_bytes,
    .make   = chol2_make,
    .free   = chol2_free
};


static
int chol2_factor(cvxc_kktsolver_t *S, cvxc_scaling_t *W, cvxc_matrix_t *H, cvxc_matrix_t *Df)
{
}

// Solve
//
//     [ H          A'   GG'*W^{-1} ]   [ ux   ]   [ bx        ]
//     [ A          0    0          ] * [ uy   [ = [ by        ]
//     [ W^{-T}*GG  0   -I          ]   [ W*uz ]   [ W^{-T}*bz ]
//
// and return ux, uy, W*uz.
//
// On entry, x, y, z contain bx, by, bz.  On exit, they contain
// the solution ux, uy, W*uz.

static
int chol2_solve(cvxc_kktsolver_t *S, cvxc_matrix_t *x, cvxc_matrix_t *y, cvxc_matgrp_t *z_g)
{
    return 0;
}


static
int chol2_init(cvxc_kktsolver_t *S, cvxc_conelp_problem_t *cp, int n, int m, const cvxc_dimset_t *dims)
{
    return 0;
}


static
cvxc_size_t chol2_bytes(int n, int m, const cvxc_dimset_t *dims)
{
    return 0;
}

static
cvxc_size_t chol2_make(cvxc_kktsolver_t *kkt,
                     cvxc_conelp_problem_t *cp,
                     int n,
                     int m,
                     const cvxc_dimset_t *dims,
                     void *mem,
                     cvxc_size_t nbytes)
{
    return 0;
}

static
cvxc_kktsolver_t *chol2_new(cvxc_conelp_problem_t *cp,
                          int n,
                          int m,
                          const cvxc_dimset_t *dims)
{
    return (cvxc_kktsolver_t *)0;
}

static
void chol2_free(cvxc_kktsolver_t *kkt)
{
    if (!kkt)
        return;
    
}

cvxc_kktfuncs_t *cvxc_chol2_load()
{
    return &chol2_functions;
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
