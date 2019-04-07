
// Copyright: Harri Rautila, 2018 <harri.rautila@gmail.com>

#include "convex.h"


static inline
void mksymm(cvx_matrix_t *x, int n)
{
    cvx_matrix_t xcol, xrow;
    for (int i = 0; i < n-1; i++) {
        // i'th column and i'th row
        cvxm_view_map(&xcol, x, i+1, i, n-1-i, 1);
        cvxm_view_map(&xrow, x, i, i+1, 1, n-1-i);
        cvxm_copy(&xrow, &xcol, CVX_ALL);
    }
}
/*
 * Convert lower triangular matrices in 'S' to symmetric. Fill in the
 * strictly upper triangular part of the symmetric matrix in x with
 * strictly lower triangular part
 */
int cvx_mksymm(cvx_matgrp_t *x_g)
{
    cvx_matrix_t xs;

    for (int k = 0; k < cvx_mgrp_count(x_g, CVXDIM_SDP); k++) {
        int n = cvx_mgrp_elem(&xs, x_g, CVXDIM_SDP, k);
        mksymm(&xs, n);
    }
    return 0;
}

/*
 *   Sets upper triangular part of the 'S' components of x equal to zero
 *   and scales the strictly lower triangular part by 2.0.
 */
int cvx_trisc(cvx_matgrp_t *x_g)
{
    cvx_matrix_t xs;
    int k;

    for (k = 0; k < cvx_mgrp_count(x_g, CVXDIM_SDP); k++) {
        cvx_mgrp_elem(&xs, x_g, CVXDIM_SDP, k);
        cvxm_make_trm(&xs, CVX_LOWER);
        cvxm_scale(&xs, (cvx_float_t)2.0, CVX_LOWER|CVX_UNIT);
        //cvx_mat_printf(stdout, "%13.6e", &xs, "trisc(x)");
    }
    return 0;
}


/*
 * Scales the strictly lower triangular part of the 's' components of x  by 0.5.
 */
int cvx_triusc(cvx_matgrp_t *x_g)
{
    cvx_matrix_t xs;
    int k;

    for (k = 0; k < cvx_mgrp_count(x_g, CVXDIM_SDP); k++) {
        cvx_mgrp_elem(&xs, x_g, CVXDIM_SDP, k);
        cvxm_scale(&xs, 0.5, CVX_LOWER|CVX_UNIT);
    }
    return 0;
}

/*
 *   Matrix-vector multiplication.
 *
 *   A is a matrix of size (m, n) where 
 *
 *       m = dims['l'] + sum(dims['q']) + sum( k**2 for k in dims['s'] ) 
 *
 *   representing a mapping from R^n to S.  
 *
 *   If no transpose defined, computes
 *
 *       y := alpha*A*x + beta * y   
 *
 *   x is a vector of length n.  y is a vector of length m.
 *
 *   If transpose flag CVX_TRANS defined, computes
 *
 *       y := alpha*A^T*x + beta * y 
 *
 *   x is a vector of length m.  y is a vector of length n.
 *
 *   The 's' components in S are stored in unpacked lower triangular storage.
 */
int cvx_sgemv(cvx_float_t beta,
              cvx_matrix_t *y,
              cvx_float_t alpha,
              const cvx_matrix_t *A,
              cvx_matgrp_t *x_g,
              int flags)
{
    if (flags & CVX_TRANS) {
        //cvx_mat_printf(stdout, "%e", x_g->mat, "sgemv x");
        cvx_trisc(x_g);
    }

    int err = cvxm_mvmult(beta, y, alpha, A, x_g->mat, flags);

    if (flags & CVX_TRANS) {
        cvx_triusc(x_g);
    }
    return err;
}

/*
 *   If no transpose defined, computes
 *
 *       y := alpha*[A; B]*x + beta * y
 *
 *   x is a vector of length n.  y is a vector of length m.
 *
 *   If transpose flag CVX_TRANS defined, computes
 *
 *       y := alpha*[A; B]^T*x + beta * y
 *
 *   x is a vector of length m.  y is a vector of length n.
 *
 *   The 's' components in S are stored in unpacked lower triangular storage.
 *
 *   A = Df[mnl, n]
 *   B =  G[M,   n]  M = dims[L] + sum(dims[Q_i}) + sum (dims[S_i]*dims[S_i])
 *
 * If CP problem computes:
 *        v  = beta * v  + alpha * Df * epi(u)    if NOTRANS
 *    epi(v) = beta * epi(v) + alpha * Df^T * u   if TRANS
 */
int cvx_sgemv2(cvx_float_t beta,
               cvx_matrix_t *y,
               cvx_float_t alpha,
               const cvx_matrix_t *A,
               const cvx_matrix_t *B,
               cvx_matgrp_t *x_g,
               int flags)
{
    cvx_matrix_t t0, t1;
    int err = 0;
    cvx_float_t v;
    cvx_size_t ar, ac, yr, yc, xr, xc;

    cvxm_size(&ar, &ac, A);
    cvxm_size(&yr, &yc, y);
    cvxm_size(&xr, &xc, x_g->mat);

    if (flags & CVX_TRANS) {
        cvx_trisc(x_g);
        cvx_size_t br, bc;
        cvxm_size(&br, &bc, B);
        cvxm_view_map(&t0, x_g->mat, 0,  0, ar, 1);
        cvxm_view_map(&t1, x_g->mat, ar, 0, xr-ar, 1);
        if (A) {
            err = cvxm_mvmult(beta, y, alpha, A, &t0, flags);
            if (cvxm_isepi(y)) {
                v = beta*cvxm_get_epi(y) - alpha*cvxm_get(x_g->mat, 0, 0);
                cvxm_set_epi(y, v);
            }
        }
        //cvx_mat_printf(stderr, "%e", y, "y.0");
        if (err == 0 && B)  {
            err = cvxm_mvmult(beta, y, alpha, B, &t1, flags);
            if (cvxm_isepi(y)) {
                cvxm_set_epi(y, beta*cvxm_get_epi(y));
            }
        }
        cvx_triusc(x_g);
        return err;
    }

    cvxm_view_map(&t0, y, 0,  0, ar, 1);
    cvxm_view_map(&t1, y, ar, 0, yr-ar, 1);
    err = cvxm_mvmult(beta, &t0, alpha, A, x_g->mat, 0);
    if (err != 0)
        return err;
    if (cvxm_isepi(x_g->mat)) {
        v = cvxm_get(y, 0, 0) - alpha * cvxm_get_epi(x_g->mat);
        cvxm_set(y, 0, 0, v);
    }
    err = cvxm_mvmult(beta, &t1, alpha, B, x_g->mat, 0);

    return err;
}

/*
 * The inverse product x := y o\ x, when 'S' components of y are diagonal.
 */
int cvx_sinv(cvx_matgrp_t *x_g,
             cvx_matgrp_t *y_g,
             cvx_memblk_t *work)
{
    cvx_matrix_t xk, yk, x1, y1;
    cvx_size_t m;
    cvx_float_t aa, ee, cc, dd;
    int k;

    // the nonlinear and 'l' blocks
    //   xk = yk o\ xk = yk / xk

    if (cvx_mgrp_count(x_g, CVXDIM_NLTARGET) > 0) {
        cvx_mgrp_elem(&xk, x_g, CVXDIM_NLTARGET, 0);
        cvx_mgrp_elem(&yk, y_g, CVXDIM_NLTARGET, 0);
        cvxm_solve_diag(&xk, 1.0, &yk, 0);
    }
    // Non-linear blocks
    if (cvx_mgrp_count(x_g, CVXDIM_NONLINEAR) > 0) {
        cvx_mgrp_elem(&xk, x_g, CVXDIM_NONLINEAR, 0);
        cvx_mgrp_elem(&yk, y_g, CVXDIM_NONLINEAR, 0);
        cvxm_solve_diag(&xk, 1.0, &yk, 0);
    }

    // 'L' blocks
    if (cvx_mgrp_count(x_g, CVXDIM_LINEAR) > 0) {
        cvx_mgrp_elem(&xk, x_g, CVXDIM_LINEAR, 0);
        cvx_mgrp_elem(&yk, y_g, CVXDIM_LINEAR, 0);
        cvxm_solve_diag(&xk, 1.0, &yk, 0);
    }

    // For the 'Q' (SOCP) blocks:
    //
    //                        [ y0   -y1'              ]
    //     yk o\ xk = 1/a^2 * [                        ] * xk
    //                        [ -y1  (a*I + y1*y1')/y0 ]
    //
    // where yk = (y0, y1) and a = y0^2 - y1'*y1.
    for (k = 0; k < cvx_mgrp_count(x_g, CVXDIM_SOCP); k++) {
        m = cvx_mgrp_elem(&xk, x_g, CVXDIM_SOCP, k);
        cvx_mgrp_elem(&yk, y_g, CVXDIM_SOCP, k);

        cvxm_view_map(&x1, &xk, 1, 0, m-1, 1);
        cvxm_view_map(&y1, &yk, 1, 0, m-1, 1);

        aa = cvxm_nrm2(&y1);
        ee = cvxm_get(&yk, 0, 0);
        cc = cvxm_get(&xk, 0, 0);
        dd = cvxm_dot(&x1, &y1);
        aa = (ee + aa)*(ee - aa);
        cvxm_set(&xk, 0, 0, cc*ee-dd);
        cvxm_scale(&x1, aa/ee, CVX_ALL);
        cvxm_axpy(&x1, dd/ee-cc, &y1);

        cvxm_scale(&xk, 1.0/aa, CVX_ALL);
    }

    // For the 'S' (SDP) blocks:
    //
    //     yk o\ xk =  xk ./ gamma
    //
    // where gamma_{i,j} = .5 * (yk_i + yk_j).
    for (k = 0; k < cvx_mgrp_count(x_g, S); k++) {
        cvx_matrix_t u, xc, yc;
        m = cvx_mgrp_elem(&xk, x_g, CVXDIM_SDP, k);
        cvx_mgrp_elem(&yk, y_g, CVXDIM_SDP, k);
        for (cvx_size_t i = 0; i < m; i++) {
            cvxm_map_data(&u, m-i, 1, __mblk_offset(work, 0));
            // column of x
            cvxm_view_map(&xc, &xk, i, i, m-i, 1);
            cvxm_view_map(&yc, &yk, i, 0, m-i, 1);
            // u := yk[i:] .+ yk[i]
            cvxm_copy(&u, &yc, CVX_ALL);
            cvxm_add(&u, cvxm_get(&yc, 0, 0), CVX_ALL);
            cvxm_scale(&u, 0.5, CVX_ALL);
            //cvx_mat_printf(stdout, "%e", &u, "sinv: u");
            //cvx_mat_printf(stdout, "%e", &xc, "sinv: xc");
            // x[:,i] := 2.0*u*x[:,i] ; multiply by 2 is equal to scaling u by 0.5
            cvxm_solve_diag(&xc, 1.0, &u, CVX_RIGHT);
            //cvx_mat_printf(stdout, "%e", &xc, "sinv: xc");
        }
    }
    return 0;
}


/*
 * The product x := (y o x).  If flag bit CVX_DIAG is set, the 's' part of y is
 * diagonal and only the diagonal is stored.
 */
int cvx_sprod(cvx_matgrp_t *x_g,
              cvx_matgrp_t *y_g,
              int flags,
              cvx_memblk_t *work)
{
    cvx_matrix_t xk, yk, x1, y1, A;
    cvx_size_t m;
    cvx_float_t dd, y0, x0;
    // cvx_index_t *index = x_g->index;
    int k;

    // For the nonlinear and 'l' blocks:
    //
    //     yk o xk = yk .* xk.
    // the nonlinear and 'l' blocks
    //   xk = yk o\ xk = yk / xk

    if (cvx_mgrp_count(x_g, CVXDIM_NLTARGET) > 0) {
        cvx_mgrp_elem(&xk, x_g, CVXDIM_NLTARGET, 0);
        cvx_mgrp_elem(&yk, y_g, CVXDIM_NLTARGET, 0);
        cvxm_mult_diag(&xk, 1.0, &yk, 0);
    }
    // Non-linear blocks
    if (cvx_mgrp_count(x_g, CVXDIM_NONLINEAR) > 0) {
        cvx_mgrp_elem(&xk, x_g, CVXDIM_NONLINEAR, 0);
        cvx_mgrp_elem(&yk, y_g, CVXDIM_NONLINEAR, 0);
        cvxm_mult_diag(&xk, 1.0, &yk, 0);
    }

    // 'L' blocks
    if (cvx_mgrp_count(x_g, CVXDIM_LINEAR) > 0) {
        cvx_mgrp_elem(&xk, x_g, CVXDIM_LINEAR, 0);
        cvx_mgrp_elem(&yk, y_g, CVXDIM_LINEAR, 0);
        cvxm_mult_diag(&xk, 1.0, &yk, 0);
    }

    // For 'q' blocks:
    //
    //               (y0   y1' )        (y0  y1' ) (x0)   (y0*x0+y1'*x1 )
    //     yk o xk = (         ) * xk = (        ) (  ) = (             )
    //               (y1   y0*I)        (y1  y0*I) (x1)   (y1*x0+y0*I*x1)
    //
    // where yk = (y0, y1).
    for (k = 0; k < cvx_mgrp_count(x_g, CVXDIM_SOCP); k++) {
        m = cvx_mgrp_elem(&xk, x_g, CVXDIM_SOCP, k);
        cvx_mgrp_elem(&yk, y_g, CVXDIM_SOCP, k);

        // ind points to start of xk, yk in x and y
        dd  = cvxm_dot(&yk, &xk);
        y0  = cvxm_get(&yk, 0, 0);
        x0  = cvxm_get(&xk, 0, 0);
        // x0 = yk^T*xk = y0*x0 + y1'*x1
        cvxm_set(&xk, 0, 0, dd);
        // remap  xk, yk = x1, y1
        cvxm_view_map(&x1, &xk, 1, 0, m-1, 1);
        cvxm_view_map(&y1, &yk, 1, 0, m-1, 1);
        // x1 = y0*x1
        cvxm_scale(&x1, y0, CVX_ALL);
        // x1 += y1*x0
        cvxm_axpy(&x1, x0, &y1);
    }

    // For the 's' blocks:
    //
    //    yk o xk = .5 * ( Yk * mat(xk) + mat(xk) * Yk ) ; symmetric rank update?
    //
    // where Yk = diag(yk) if CVX_DIAG is set in flags and Yk = mat(yk) otherwise
    for (k = 0; k < cvx_mgrp_count(x_g, CVXDIM_SDP); k++) {
        // .5 * (diag(yk) * mat(xk) + mat(xk) * diag(yk)) ==
        m = cvx_mgrp_elem(&xk, x_g, CVXDIM_SDP, k);
        cvx_mgrp_elem(&yk, y_g, CVXDIM_SDP, k);

        if (flags & CVX_DIAG) {
            cvx_matrix_t u, xc, yc;
            for (cvx_size_t i = 0; i < m; i++) {
                cvxm_map_data(&u, m-i, 1, __mblk_offset(work, 0));
                // column of x
                cvxm_view_map(&xc, &xk, i, i, m-i, 1);
                cvxm_view_map(&yc, &yk, i, 0, m-i, 1);
                // u := yk[i:] .+ yk[i]
                cvxm_copy(&u, &yc, CVX_ALL);
                cvxm_add(&u, cvxm_get(&yc, 0, 0), CVX_ALL);
                // x[:,i] := 0.5*u*x[:,i]
                cvxm_mult_diag(&xc, 0.5, &u, CVX_LEFT);
            }
        } else {
            cvxm_map_data(&A, m, m, __mblk_offset(work, 0));
            cvxm_copy(&A, &xk, CVX_LOWER);
            cvxm_mksymm(&A, m);
            cvxm_mksymm(&yk, m);
            //cvx_mat_printf(stdout, "%13.6e", &A, "symm(A)");
            //cvx_mat_printf(stdout, "%13.6e", &yk, "symm(yk)");
            // xk = 0.5*(yk*A^T + A*yk^T)  = 0.5*(A*yk^T + A*yk^T)
            //    =      yk*A ; mat(yk), A symmetric
            // would equal to x = tril(gemm(yk, A)), gemm implememtation maybe faster
            cvxm_update2_sym(0.0, &xk, 0.5, &A, &yk, 0);
        }
        // we don't care about the strictly upper tridiagonal part (could be zeroed)
        //cvxm_make_trm(&xk, CVX_LOWER);
    }
    //cvx_mat_printf(stdout, "%13.6e", x_g->mat, "sprod(x)");
    return 0;
}


/*
 * The product x := y o y.   The 's' components of y are diagonal and
 * only the diagonals of x and y are stored.
 */
int cvx_ssqr(cvx_matgrp_t *x_g,
             cvx_matgrp_t *y_g)
{
    cvx_matrix_t xk, yk, x1, y1;
    cvx_size_t ind, m;
    cvx_float_t aa;
    cvx_index_t *index = x_g->index;
    int k;

    ind = 0;

    // Non-linear target function
    if (index->indnlt) {
        cvx_mgrp_elem(&xk, x_g, CVXDIM_NLTARGET, 0);
        cvx_mgrp_elem(&yk, y_g, CVXDIM_NLTARGET, 0);
        cvxm_copy(&xk, &yk, CVX_ALL);
        cvxm_mult_diag(&xk, 1.0, &yk, 0);
        ind += 1;
    }

    // Non-linear blocks
    if (index->indnl) {
        cvx_mgrp_elem(&xk, x_g, CVXDIM_NONLINEAR, 0);
        cvx_mgrp_elem(&yk, y_g, CVXDIM_NONLINEAR, 0);
        cvxm_copy(&xk, &yk, CVX_ALL);
        cvxm_mult_diag(&xk, 1.0, &yk, 0);
        ind += cvx_index_length(index, CVXDIM_NONLINEAR);
    }

    // 'L' blocks
    if (index->indl) {
        cvx_mgrp_elem(&xk, x_g, CVXDIM_LINEAR, 0);
        cvx_mgrp_elem(&yk, y_g, CVXDIM_LINEAR, 0);
        cvxm_copy(&xk, &yk, CVX_ALL);
        cvxm_mult_diag(&xk, 1.0, &yk, 0);
        ind += cvx_index_length(index, CVXDIM_LINEAR);
    }

    //  'Q' part
    for (k = 0; k < cvx_mgrp_count(x_g, CVXDIM_SOCP); k++) {
        m = cvx_mgrp_elem(&yk, y_g, CVXDIM_SOCP, k);
        cvx_mgrp_elem(&xk, x_g, CVXDIM_SOCP, k);
        cvxm_copy(&xk, &yk, CVX_ALL);

        aa = cvxm_nrm2(&yk);

        cvxm_set(&xk, 0, 0, aa*aa);
        cvxm_view_map(&x1, &xk, 1, 0, m-1, 1);
        cvxm_view_map(&y1, &yk, 1, 0, m-1, 1);
        cvxm_scale(&x1, 2.0*cvxm_get(&yk, 0, 0), CVX_ALL);
        ind += m;
    }
    // 'S' part
    m = cvx_index_length(index, CVXDIM_SDP);
    cvxm_view_map(&xk, x_g->mat, ind, 0, m, 1);
    cvxm_view_map(&yk, y_g->mat, ind, 0, m, 1);
    cvxm_copy(&xk, &yk, CVX_ALL);
    cvxm_mult_diag(&xk, 1.0, &yk, 0);
    return 0;
}


// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
