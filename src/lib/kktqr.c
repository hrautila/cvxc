
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include "convex.h"

// forward declarations
static
int qr_init(cvx_kktsolver_t *S,
            cvx_conelp_problem_t *cp,
            int n,
            int m,
            const cvx_dimset_t *dims);

static
int qr_factor(cvx_kktsolver_t *S,
              cvx_scaling_t *W,
              cvx_matrix_t *H,
              cvx_matrix_t *Df);
static
int qr_solve(cvx_kktsolver_t *S,
             cvx_matrix_t *x,
             cvx_matrix_t *y,
             cvx_matgrp_t *z_g);

static
cvx_size_t qr_bytes(int n, int m, const cvx_dimset_t *dims);

static
cvx_size_t qr_make(cvx_kktsolver_t *kkt,
                   cvx_conelp_problem_t *cp,
                   int n,
                   int m,
                   const cvx_dimset_t *dims,
                   void *mem,
                   cvx_size_t nbytes);

static
cvx_kktsolver_t *qr_new(cvx_conelp_problem_t *cp,
                        int n,
                        int m,
                        const cvx_dimset_t *dims);

static
void qr_free(cvx_kktsolver_t *S);

// function table
static cvx_kktfuncs_t qrfunctions = {
    .new    = qr_new,
    .factor = qr_factor,
    .solve  = qr_solve,
    .init   = qr_init,
    .bytes  = qr_bytes,
    .make   = qr_make,
    .free   = qr_free
};


// Solution of KKT equations with zero 1,1 block, by eliminating the
// equality constraints via a QR factorization, and solving the
// reduced KKT system by another QR factorization.
//    
// Computes the QR factorization
//    
//        A' = [Q1, Q2] * [R1; 0]
//    
// and returns a function that (1) computes the QR factorization 
//    
//        W^{-T} * G * Q2 = Q3 * R3
//   
// (with columns of W^{-T}*G in packed storage), and (2) returns a function for solving 
//    
//        [ 0    A'   G'    ]   [ ux ]   [ bx ]
//        [ A    0    0     ] * [ uy ] = [ by ].
//        [ G    0   -W'*W  ]   [ uz ]   [ bz ]
//    
// A is p x n and G is N x n where N = dims['l'] + sum(dims['q']) + 
// sum( k**2 for k in dims['s'] ).
//


static
int qr_factor(cvx_kktsolver_t *S,
              cvx_scaling_t *W,
              cvx_matrix_t *H,
              cvx_matrix_t *Df)
{
    cvx_qrsolver_t *qr = (cvx_qrsolver_t *)S;
    cvx_conelp_problem_t *cp = qr->cp;
    cvx_matgrp_t Gs_g;
    
    qr->W = W;

    // Gs = W^{-T}*G in packed storage
    cvxm_copy(&qr->Gs, cp->G, 0);
    cvx_mgrp_init(&Gs_g, &qr->Gs, &cp->index_full);
    cvx_scale(&Gs_g, W, CVX_INV|CVX_TRANS, &cp->work);
    cvx_pack2(&Gs_g);

    // Gs = [Gs1, Gs2] = Gs * [Q1; Q2]^T
    cvxm_lqmult(&Gs, &qr.QA, &qr.tauA, &qr.work, CVX_RIGHT|CVX_TRANS);
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
int qr_solve(cvx_kktsolver_t *S,
             cvx_matrix_t *x,
             cvx_matrix_t *y,
             cvx_matgrp_t *z_g)
{
    cvx_matrix_t vv1, Gs1, Gs2;
    cvx_qrsolver_t *qr = (cvx_qrsolver_t *)S;
    cvx_conelp_problem_t *cp = qr->cp;
    
#if 0    
    cvxm_view_map(&Gs1, &qr->Gs, 0, qr->p, cp->cdim_packed, qr->n - qr->p);
    cvxm_view_map(&Gs2, &qr->Gs, 0, qr->p, cp->cdim_packed, qr->n - qr->p);
    
    cvx_scale(z_g, &cp->W, CVX_INV|CVX_TRANS, &cp->work);
    cvx_pack(z_g);

    // vv := [ vv0, vv1 ] = [Q1*bx, R3^{-T}*Q2*bx]
    cvxm_copy(&qr->vv, x);
    cvxm_lqmult(&qr->vv, &qr->QA, &qr->tauA, 0);
    // 
    cvxm_view_map(&vv1, &qr->vv, qr->p, 0, qr->n-qr->p, 1);
    cvxm_mvsolve(&vv1, 1.0, &Gs2, CVX_LOWER);

    // x[:p] = R1^{-T} * by
    if (qr->p > 0) {
        cvxm_copy(x, y, 0);
        cvxm_mvsolve(x0, 1.0, &qr->QA, CVX_LOWER);
    }
    // w = w - Gs1 * x[:p]  = W^{-T}*bz - Gs1*y
    cvxm_mvmult(1.0, &qr->w, -1.0, &Gs1, x);
#endif
    return 0;
}


static
int qr_init(cvx_kktsolver_t *S,
            cvx_conelp_problem_t *cp,
            int n,
            int m,
            const cvx_dimset_t *dims)
{
    cvx_qrsolver_t *qr = (cvx_qrsolver_t *)S;
    qr->cp = cp;
    qr->n = n;
    qr->p = m;

    cvxm_init(&qr->QA, m, n);
    cvxm_init(&qr->tauA, n, 1);
    cvxm_init(&qr->tauG, cp->cdim, 1);
    cvxm_init(&qr->Gs, cp->cdim, n);
    cvxm_init(&qr->u, cp->cdim_packed, 1);
    cvxm_init(&qr->w, cp->cdim_packed, 1);
    cvxm_init(&qr->vv, m, 1);
    
    // QR(A^T) = LQ(A) = [Q1; Q2]^T
    cvxm_lqfactor(&qr->QA, &qr->tauA, &qr->work);
    return 0;
}


static
cvx_size_t qr_bytes(int n, int m, const cvx_dimset_t *dims)
{
    return 0;
}

static
cvx_size_t qr_make(cvx_kktsolver_t *kkt,
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
cvx_kktsolver_t *qr_new(cvx_conelp_problem_t *cp,
                        int n,
                        int m,
                        const cvx_dimset_t *dims)
{
    return (cvx_kktsolver_t *)0;
}

static
void qr_free(cvx_kktsolver_t *kkt)
{
    if (!kkt)
        return;
    
}

cvx_kktfuncs_t *cvx_qrload()
{
    return &qrfunctions;
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
