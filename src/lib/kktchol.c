
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include "convex.h"

// forward declarations
static
int chol_init(cvx_kktsolver_t *S,
             cvx_conelp_problem_t *cp,
             int n,
             int m,
             const cvx_dimset_t *dims);

static
int chol_factor(cvx_kktsolver_t *S,
               cvx_scaling_t *W,
               cvx_matrix_t *H,
               cvx_matrix_t *Df);
static
int chol_solve(cvx_kktsolver_t *S,
              cvx_matrix_t *x,
              cvx_matrix_t *y,
              cvx_matgrp_t *z_g);

static
cvx_size_t chol_bytes(int n, int m, const cvx_dimset_t *dims);

static
cvx_size_t chol_make(cvx_kktsolver_t *kkt,
                    cvx_conelp_problem_t *cp,
                    int n,
                    int m,
                    const cvx_dimset_t *dims,
                    void *mem,
                    cvx_size_t nbytes);

static
cvx_kktsolver_t *chol_new(cvx_conelp_problem_t *cp,
                         int n,
                         int m,
                         const cvx_dimset_t *dims);

static
void chol_free(cvx_kktsolver_t *S);

// function table
static cvx_kktfuncs_t chol_functions = {
    .new    = chol_new,
    .factor = chol_factor,
    .solve  = chol_solve,
    .init   = chol_init,
    .bytes  = chol_bytes,
    .make   = chol_make,
    .free   = chol_free
};


static
int chol_factor(cvx_kktsolver_t *S, cvx_scaling_t *W, cvx_matrix_t *H, cvx_matrix_t *Df)
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
int chol_solve(cvx_kktsolver_t *S, cvx_matrix_t *x, cvx_matrix_t *y, cvx_matgrp_t *z_g)
{
    return 0;
}


static
int chol_init(cvx_kktsolver_t *S, cvx_conelp_problem_t *cp, int n, int m, const cvx_dimset_t *dims)
{
    return 0;
}


static
cvx_size_t chol_bytes(int n, int m, const cvx_dimset_t *dims)
{
    return 0;
}

static
cvx_size_t chol_make(cvx_kktsolver_t *kkt,
                     cvx_conelp_problem_t *cp,
                     int n,
                     int m,
                     const cvx_dimset_t *dims,
                     void *mem,
                     cvx_size_t nbytes)
{
    return 0;
}

static
cvx_kktsolver_t *chol_new(cvx_conelp_problem_t *cp,
                          int n,
                          int m,
                          const cvx_dimset_t *dims)
{
    return (cvx_kktsolver_t *)0;
}

static
void chol_free(cvx_kktsolver_t *kkt)
{
    if (!kkt)
        return;
    
}

cvx_kktfuncs_t *cvx_chol_load()
{
    return &chol_functions;
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
