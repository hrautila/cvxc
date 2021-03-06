/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */
#include "internal.h"

#ifndef SQRT2
#define SQRT2    1.41421356237309504880D
#endif
#define INVSQRT2 0.70710678118654752440D

/*
 * Element wise operators
 */
static
cvxc_float_t fsqrt(cvxc_float_t a)
{
    return SQRT(a);
}

static
cvxc_float_t finv(cvxc_float_t a)
{
    return 1.0/a;
}
static
int cvxc_scale_vector(cvxc_matgrp_t *x_g, cvxc_scaling_t *W, int flags, cvxc_memblk_t *work);

/*
 *  \brief Applies Nesterov-Todd scaling or its inverse.
 *
 *   Computes
 *
 *        x := W*x        (trans is false 'N', inverse = false 'N')
 *        x := W^T*x      (trans is true  'T', inverse = false 'N')
 *        x := W^{-1}*x   (trans is false 'N', inverse = true  'T')
 *        x := W^{-T}*x   (trans is true  'T', inverse = true  'T').
 *
 *   x is a dense m-by-n real matrix
 *
 *   W is a MatrixSet with entries:
 *
 *   - W['dnl']: positive vector
 *   - W['dnli']: componentwise inverse of W['dnl']
 *   - W['d']: positive vector
 *   - W['di']: componentwise inverse of W['d']
 *   - W['v']: lists of 2nd order cone vectors with unit hyperbolic norms
 *   - W['beta']: list of positive numbers
 *   - W['r']: list of square matrices
 *   - W['rti']: list of square matrices.  rti[k] is the inverse transpose
 *     of r[k].
 *
 *   The 'dnl' and 'dnli' entries are optional, and only present when the
 *   function is called from the nonlinear solver.
 *
 */
int cvxc_scale(cvxc_matgrp_t *x_g,
              cvxc_scaling_t *W,
              int flags,
              cvxc_memblk_t *work)
{
    cvxc_matrix_t xk, xk0, Ak, w, v, dnl, dl, r, beta;
    cvxc_size_t j, k, m, nrows, ncols, ind = 0;
    cvxc_matrix_t *x = x_g->mat;

    int inverse = (flags & CVXC_INV) != 0;

    cvxc_mblk_clear(work);

    ncols = 1;
    cvxm_size(&nrows, &ncols, x);

    if (ncols == 1)
        return cvxc_scale_vector(x_g, W, flags, work);

    // NOTE:
    // x is M-by-ncols matrix; special case if ncols == 1
    // - linear and non-linear case update each x column
    // - Q, S space as m-by-ncol matrix

    // Scaling for nonlinear component xk is xk := dnl .* xk;
    // inverse scaling is xk ./ dnl = dnli .* xk,
    // where dnl = W['dnl'], dnli = W['dnli'].
    if (W->dnltsz > 0 && cvxc_mgrp_count(x_g, CVXDIM_NLTARGET) > 0) {
        m = cvxc_scaling_elem(&dnl, W, (inverse ? CVXWS_DNLTI : CVXWS_DNLT), 0);
        cvxm_view_map(&xk, x, ind, 0, m, ncols);
        cvxm_mult_diag(&xk, 1.0, &dnl, CVXC_LEFT);
        ind += m;
    }
    if (W->dnlsz > 0 && cvxc_mgrp_count(x_g, CVXDIM_NONLINEAR) > 0) {
        m = cvxc_scaling_elem(&dnl, W, (inverse ? CVXWS_DNLI : CVXWS_DNL), 0);
        cvxm_view_map(&xk, x, ind, 0, m, ncols);
        cvxm_mult_diag(&xk, 1.0, &dnl, CVXC_LEFT);
        ind += m;
    }

    // Scaling for linear 'l' component xk is xk := d .* xk;
    // inverse scaling is xk ./ d = di .* xk,
    // where d = W['d'], di = W['di'].
    if (W->dsz > 0 && cvxc_mgrp_count(x_g, CVXDIM_LINEAR) > 0) {
        m = cvxc_scaling_elem(&dl, W, (inverse ? CVXWS_DI : CVXWS_D), 0);
        cvxm_view_map(&xk, x, ind, 0, m, ncols);
        cvxm_mult_diag(&xk, 1.0, &dl, CVXC_LEFT);
        ind += m;
    }

    // Scaling for 'q' component is
    //
    //    xk := beta * (2*v*v' - J) * xk
    //        = beta * (2*v*(xk'*v)' - J*xk)
    //
    // Inverse scaling for 'q' component is
    //
    //    xk := 1/beta * (2*J*v*v'*J - J) * xk
    //        = 1/beta * (-J) * (2*v*((-J*xk)'*v)' + xk).
    //
    // where beta = W['beta'][k], v = W['v'][k], J = [1, 0; 0, -I].
    //
    if (W->vcount > 0 && cvxc_mgrp_count(x_g, CVXDIM_SOCP) > 0) {
        cvxc_float_t bk;
        cvxm_map_data(&w, ncols, 1, cvxc_mblk_offset(work, 0));
        cvxc_scaling_elem(&beta, W, CVXWS_BETA, 0);
        for (k = 0; k < W->vcount; k++) {
            m = cvxc_scaling_elem(&v, W, CVXWS_V, k);
            cvxm_view_map(&xk, x, ind, 0, m, ncols);
            cvxm_view_map(&xk0, x, ind, 0, 1, ncols);
            if (inverse) {
                cvxm_scale(&xk0, -1.0, CVXC_ALL);
            }
            // w = xk'*v
            cvxm_mvmult(0.0, &w, 1.0, &xk, &v, CVXC_TRANS);
            // xk0= -xk0 ( -J*xk )
            cvxm_scale(&xk0, -1.0, CVXC_ALL);
            // xk = xk + 2.0*v*w
            cvxm_mvupdate(&xk, 2.0, &v, &w);
            if (inverse) {
                cvxm_scale(&xk0, -1.0, CVXC_ALL);
                bk = 1.0 / cvxm_get(&beta, k, 0);
            } else {
                bk = cvxm_get(&beta, k, 0);
            }
            // TODO: xk is m-x-ncols matrix; scale is vector operator
            cvxm_scale(&xk, bk, CVXC_ALL);
            // next
            ind += m;
        }
    }

    // Scaling for 's' component xk is
    //
    //     xk := vec( r' * mat(xk) * r )  if trans = 'N'
    //     xk := vec( r * mat(xk) * r' )  if trans = 'T'.
    //
    // r is kth element of W['r'].
    //
    // Inverse scaling for 's' component xk is
    //
    //     xk := vec( rti * mat(xk) * rti' )  if trans = 'N'
    //     xk := vec( rti' * mat(xk) * rti )  if trans = 'T'.
    //
    // rti is kth element of W['rti'].
    if (W->rcount > 0 && cvxc_mgrp_count(x_g, CVXDIM_SDP) > 0) {
        cvxc_matrix_t xd;
        int rflags = inverse ? flags & CVXC_TRANS : (flags ^ CVXC_TRANS);
        for (k = 0; k < W->rcount; k++) {
            m = cvxc_scaling_elem(&r, W, (inverse ? CVXWS_RTI : CVXWS_R), k);
            // workspace Ak = mat(xk)
            cvxm_map_data(&Ak, m, m, cvxc_mblk_offset(work, 0));
            // each column of x is matrix in lower triangular storage format
            for (j = 0; j < ncols; j++) {
                cvxm_map_data(&xk, m, m, cvxm_data(x, ind+j*nrows));
                // scale diagonal by 0.5
                cvxm_view_diag(&xd, &xk, 0);
                cvxm_scale(&xd, 0.5, CVXC_ALL);
                //
                cvxm_copy(&Ak, &r, CVXC_ALL);
                if ((rflags & CVXC_TRANS) == 0) {
                    cvxm_mult_trm(&Ak, 1.0, &xk, CVXC_RIGHT|CVXC_LOWER);
                    // x = 0.0*x + (r*a' + a*r') (lower triangular storage)
                    cvxm_update2_sym(0.0, &xk, 1.0, &r, &Ak, rflags|CVXC_LOWER);
                } else {
                    cvxm_mult_trm(&Ak, 1.0, &xk, CVXC_LEFT|CVXC_LOWER);
                    // x = 0.0*x + (r'*a + a'*r) (lower triangular storage)
                    cvxm_update2_sym(0.0, &xk, 1.0, &r, &Ak, rflags|CVXC_LOWER);
                }
            }
            ind += m*m;
        }
    }
    return 0;
}

// X is column vector
static
int cvxc_scale_vector(cvxc_matgrp_t *x_g,
                     cvxc_scaling_t *W,
                     int flags,
                     cvxc_memblk_t *work)
{
    cvxc_matrix_t xk, v, r, dnl, dl;
    cvxc_size_t m;

    int inverse = (flags & CVXC_INV) != 0;

    // Scaling for nonlinear component xk is xk := dnl .* xk; where dnl = W['dnl'],
    if (W->dnltsz > 0 && cvxc_mgrp_count(x_g, CVXDIM_NLTARGET) > 0) {
        cvxc_scaling_elem(&dnl, W, (inverse ? CVXWS_DNLTI : CVXWS_DNLT), 0);
        cvxc_mgrp_elem(&xk, x_g, CVXDIM_NLTARGET, 0);
        cvxm_mult_diag(&xk, 1.0, &dnl, CVXC_LEFT);
    }
    if (W->dnlsz > 0 && cvxc_mgrp_count(x_g, CVXDIM_NONLINEAR) > 0) {
        cvxc_scaling_elem(&dnl, W, (inverse ? CVXWS_DNLI : CVXWS_DNL), 0);
        cvxc_mgrp_elem(&xk, x_g, CVXDIM_NONLINEAR, 0);
        cvxm_mult_diag(&xk, 1.0, &dnl, CVXC_LEFT);
    }

    // Scaling for linear 'l' component xk is xk := d .* xk; where d = W['d']
    if (W->dsz > 0 && cvxc_mgrp_count(x_g, CVXDIM_LINEAR) > 0) {
        cvxc_scaling_elem(&dl, W, (inverse ? CVXWS_DI : CVXWS_D), 0);
        cvxc_mgrp_elem(&xk, x_g, CVXDIM_LINEAR, 0);
        cvxm_mult_diag(&xk, 1.0, &dl, CVXC_LEFT);
    }

    // Scaling for 'q' component is
    //
    //    xk := beta * (2*v*v' - J) * xk
    //        = beta * (2*v*(xk'*v)' - J*xk)
    //        = beta * (2*w*v + [-xk0; xk1])  where w = xk'*v
    //
    // Inverse scaling is
    //
    //    xk := 1/beta * (2*J*v*v'*J - J) * xk
    //        = 1/beta * (-J) * (2*v*((-J*xk)'*v)' + xk).
    //
    // where beta = beta(W, k), v = v(W, k), J = [1, 0; 0, -I].
    if (W->vcount > 0 && cvxc_mgrp_count(x_g, CVXDIM_SOCP) > 0) {
        cvxc_matrix_t beta;
        cvxc_float_t w, bk;
        cvxc_scaling_elem(&beta, W, CVXWS_BETA, 0);
        for (int k = 0; k < W->vcount; k++) {
            m = cvxc_scaling_elem(&v, W, CVXWS_V, k);
            cvxc_mgrp_elem(&xk, x_g, CVXDIM_SOCP, k);
            if (inverse) {
                // xk = -J*xk
                cvxm_set(&xk, 0, 0, -cvxm_get(&xk, 0, 0));
            }
            // w = 2*xk'*v or w = 2*(-J*xk)'*v
            w = 2.0 * cvxm_dot(&xk, &v);
            cvxm_set(&xk, 0, 0, -cvxm_get(&xk, 0, 0));
            cvxm_axpy(&xk, w, &v);
            bk = cvxm_get(&beta, k, 0);
            if (inverse) {
                // xk = -J*xk
                cvxm_set(&xk, 0, 0, -cvxm_get(&xk, 0, 0));
                bk = 1.0 / cvxm_get(&beta, k, 0);
            }
            cvxm_scale(&xk, bk, CVXC_ALL);
        }
    }

    // Scaling for 's' component xk is
    //
    //     xk := vec( r' * mat(xk) * r )  if trans = 'N'
    //     xk := vec( r * mat(xk) * r' )  if trans = 'T'.
    //
    // r is kth element of W['r'].
    if (W->rcount > 0 && cvxc_mgrp_count(x_g, CVXDIM_SDP) > 0) {
        cvxc_matrix_t Ak, xd;
        int rflags = inverse ? flags & CVXC_TRANS : (flags ^ CVXC_TRANS);
        for (int k = 0; k < W->rcount; k++) {
            m = cvxc_scaling_elem(&r, W, (inverse ? CVXWS_RTI : CVXWS_R), k); 
            // workspace Ak = mat(xk)
            cvxm_map_data(&Ak, m, m, cvxc_mblk_offset(work, 0));
            // x is matrix in lower triangular storage format
            cvxc_mgrp_elem(&xk, x_g, CVXDIM_SDP, k);
            // scale diagonal by 0.5
            cvxm_view_diag(&xd, &xk, 0);
            cvxm_scale(&xd, 0.5, CVXC_ALL);
            cvxm_copy(&Ak, &r, CVXC_ALL);
            if ((rflags & CVXC_TRANS) == 0) {
                cvxm_mult_trm(&Ak, 1.0, &xk, CVXC_RIGHT|CVXC_LOWER);
                // x = 0.0*x + (r*a' + a*r') (lower triangular storage)
                cvxm_update2_sym(0.0, &xk, 1.0, &r, &Ak, rflags|CVXC_LOWER);
            } else {
                cvxm_mult_trm(&Ak, 1.0, &xk, CVXC_LEFT|CVXC_LOWER);
                // x = 0.0*x + (r'*a + a'*r) (lower triangular storage)
                cvxm_update2_sym(0.0, &xk, 1.0, &r, &Ak, rflags|CVXC_LOWER);
            }
        }
    }
    return 0;
}

/*
 *   Returns the Nesterov-Todd scaling W at points s and z, and stores the 
 *   scaled variable in lmbda. 
 *
 *       W * z = W^{-T} * s = lmbda. 
 *
 *   W is a MatrixSet with entries:
 *
 *   - W[dnl]: positive vector
 *   - W[dnli]: componentwise inverse of W[dnl]
 *   - W[d]: positive vector
 *   - W[di]: componentwise inverse of W['d']
 *   - W[v]: lists of 2nd order cone vectors with unit hyperbolic norms
 *   - W[beta]: list of positive numbers
 *   - W[r]: list of square matrices 
 *   - W[rti]: list of square matrices.  rti[k] is the inverse transpose
 *     of r[k].
 *
 */
int cvxc_compute_scaling(cvxc_scaling_t *W,
                        cvxc_matgrp_t *s_g,
                        cvxc_matgrp_t *z_g,
                        cvxc_matgrp_t *lmbda_g,
                        cvxc_memblk_t *work)
{
    cvxc_matrix_t sk, zk, lk;
    cvxc_size_t k, m, ind = 0;

    // For the nonlinear target block
    // TODO: optimize to singleton matrix
    if (W->dnltsz > 0) {
        cvxc_matrix_t dnlt, dnlti;
        cvxc_mgrp_elem(&sk, s_g, CVXDIM_NLTARGET, 0);
        cvxc_mgrp_elem(&zk, z_g, CVXDIM_NLTARGET, 0);
        cvxc_scaling_elem(&dnlt, W, CVXWS_DNLT, 0);
        cvxc_scaling_elem(&dnlti, W, CVXWS_DNLTI, 0);
        // dnl = sk/zk
        cvxm_copy(&dnlt, &sk, CVXC_ALL);
        cvxm_solve_diag(&dnlt, 1.0, &zk, CVXC_LEFT);
        cvxm_apply(&dnlt, fsqrt, CVXC_ALL);
        // dnli = 1/dnl
        cvxm_copy(&dnlti, &dnlt, CVXC_ALL);
        cvxm_apply(&dnlti, finv, CVXC_ALL);
        // lmblda[:mnl]
        cvxc_mgrp_elem(&lk, lmbda_g, CVXDIM_NLTARGET, 0);
        cvxm_copy(&lk, &sk, CVXC_ALL);
        cvxm_mult_diag(&lk, 1.0, &zk, CVXC_LEFT);  // element wise multiply
        cvxm_apply(&lk, fsqrt, CVXC_ALL);
        ind = W->dnltsz;
    }

    // For the nonlinear block:
    //
    //     W[dnl] = sqrt( s[:mnl] ./ z[:mnl] )
    //     W[dnli] = sqrt( z[:mnl] ./ s[:mnl] )
    //     lambda[:mnl] = sqrt( s[:mnl] .* z[:mnl] )
    if (W->dnlsz > 0) {
        cvxc_matrix_t dnl, dnli;
        cvxc_mgrp_elem(&sk, s_g, CVXDIM_NONLINEAR, 0);
        cvxc_mgrp_elem(&zk, z_g, CVXDIM_NONLINEAR, 0);
        cvxc_scaling_elem(&dnl, W, CVXWS_DNL, 0);
        cvxc_scaling_elem(&dnli, W, CVXWS_DNLI, 0);
        // dnl = sk/zk
        cvxm_copy(&dnl, &sk, CVXC_ALL);
        cvxm_solve_diag(&dnl, 1.0, &zk, CVXC_LEFT);
        cvxm_apply(&dnl, fsqrt, CVXC_ALL);
        // dnli = 1/dnl
        cvxm_copy(&dnli, &dnl, CVXC_ALL);
        cvxm_apply(&dnli, finv, CVXC_ALL);
        // lmblda[:mnl]
        cvxc_mgrp_elem(&lk, lmbda_g, CVXDIM_NONLINEAR, 0);
        cvxm_copy(&lk, &sk, CVXC_ALL);
        cvxm_mult_diag(&lk, 1.0, &zk, CVXC_LEFT);  // element wise multiply
        cvxm_apply(&lk, fsqrt, CVXC_ALL);
        ind = W->dnlsz;
    }

    // For the 'l' block:
    //
    //     W[d] = sqrt( sk ./ zk )
    //     W[di] = sqrt( zk ./ sk )
    //     lambdak = sqrt( sk .* zk )
    //
    // where sk and zk are the first dims['l'] entries of s and z.
    // lambda_k is stored in the first dims['l'] positions of lmbda.
    if (W->dsz > 0) {
        cvxc_matrix_t d, di;
        cvxc_mgrp_elem(&sk, s_g, CVXDIM_LINEAR, 0);
        cvxc_mgrp_elem(&zk, z_g, CVXDIM_LINEAR, 0);
        cvxc_scaling_elem(&d, W, CVXWS_D, 0);
        cvxc_scaling_elem(&di, W, CVXWS_DI, 0);
        // d = sk/zk
        cvxm_copy(&d, &sk, CVXC_ALL);
        cvxm_solve_diag(&d, 1.0, &zk, CVXC_LEFT);
        cvxm_apply(&d, fsqrt, CVXC_ALL);
        // di = 1/d
        cvxm_copy(&di, &d, CVXC_ALL);
        cvxm_apply(&di, finv, CVXC_ALL);
        // lmblda[mnl:l]
        //cvxm_view_map(&lk, lmbda, ind, 0, W->dsz, 1);
        cvxc_mgrp_elem(&lk, lmbda_g, CVXDIM_LINEAR, 0);
        cvxm_copy(&lk, &sk, CVXC_ALL);
        cvxm_mult_diag(&lk, 1.0, &zk, CVXC_LEFT);  // element wise multiply
        cvxm_apply(&lk, fsqrt, CVXC_ALL);
        ind += W->dsz;
    }

    // For the 'q' blocks, compute lists 'v', 'beta'.
    //
    // The vector v[k] has unit hyperbolic norm:
    //
    //     (sqrt( v[k]' * J * v[k] ) = 1 with J = [1, 0; 0, -I]).
    //
    // beta[k] is a positive scalar.
    //
    // The hyperbolic Householder matrix H = 2*v[k]*v[k]' - J
    // defined by v[k] satisfies
    //
    //    (beta[k] * H) * zk  = (beta[k] * H) \ sk = lambda_k
    //
    // where sk = s[indq[k]:indq[k+1]], zk = z[indq[k]:indq[k+1]].
    //
    //  lambda_k is stored in lmbda[indq[k]:indq[k+1]].
    if (W->vcount > 0) {
        cvxc_matrix_t beta, v;
        cvxc_float_t a, b, c, dd, zz, ss;

        cvxc_scaling_elem(&beta, W, CVXWS_BETA, 0);

        for (k = 0; k < W->vcount; k++) {
            m = cvxc_scaling_elem(&v, W, CVXWS_V, k);

            cvxc_mgrp_elem(&sk, s_g, CVXDIM_SOCP, k);
            cvxc_mgrp_elem(&zk, z_g, CVXDIM_SOCP, k);

            a = cvxc_jnrm2(&sk);
            b = cvxc_jnrm2(&zk);
            cvxm_set(&beta, k, 0, SQRT(a/b));
            // c = sqrt(dot(sk/a, zk/b) + 1) / sqrt(2)
            c = cvxm_dot(&sk, &zk);
            c = SQRT(c/(2.0*a*b) + 0.5);
            // vk = 1/(2*c) * ((sk/a) + J*(zk/b)) [J*x == (x0; -x1)]
            cvxm_copy(&v, &zk, CVXC_ALL);
            // v = J*(zk/b)
            cvxm_scale(&v, -1.0/b, CVXC_ALL);
            cvxm_set(&v, 0, 0, -cvxm_get(&v, 0, 0));
            // v = (1/2c)*(v + (sk/a))
            cvxm_axpy(&v, 1.0/a, &sk);
            cvxm_scale(&v, 0.5/c, CVXC_ALL);
            // v = 1/sqrt(2*(vk_0 + 1)) * ( vk + e ),  e = [1; 0]
            cvxm_set(&v, 0, 0, 1.0+cvxm_get(&v, 0, 0));
            // v = (1/sqrt2)*(1/sqrt2(vk_0))
            cvxm_scale(&v, INVSQRT2/SQRT(cvxm_get(&v, 0, 0)), CVXC_ALL);
            // To get the scaled variable lambda_k
            //
            //     d =  sk0/a + zk0/b + 2*c
            //     lambda_k = [ c; (c + zk0/b)/d * sk1/a + (c + sk0/a)/d * zk1/b ]
            //              = [ c; (c + zk0/b)/(d*a) * sk1 + (c + sk0/a)/(d*b) * zk1 ]
            //              = [ c; z * sk1 + s * zk1 ]
            //     lambda_k *= sqrt(a * b)
            cvxc_mgrp_elem(&lk, lmbda_g, CVXDIM_SOCP, k);
            cvxm_copy(&lk, &sk, CVXC_ALL);
            // do it for full lambda and set lmbda[0] at the end
            dd = 2*c + cvxm_get(&sk, 0, 0)/a + cvxm_get(&zk, 0, 0)/b;
            zz = (c + cvxm_get(&zk, 0, 0)/b) / (dd*a);
            ss = (c + cvxm_get(&sk, 0, 0)/a) / (dd*b);
            cvxm_scale(&lk, zz, CVXC_ALL);
            cvxm_axpy(&lk, ss, &zk);
            // lmbda[0] = c
            cvxm_set(&lk, 0, 0, c);
            cvxm_scale(&lk, SQRT(a*b), CVXC_ALL);

            ind += m;
        }
    }

    // For the 's' blocks: compute two lists 'r' and 'rti'.
    //
    //    r[k]' * sk^{-1} * r[k] = diag(lambda_k)^{-1}
    //    r[k]' * zk * r[k] = diag(lambda_k)
    //
    // where sk and zk are the entries inds[k] : inds[k+1] of
    // s and z, reshaped into symmetric matrices.
    //
    // rti[k] is the inverse of r[k]', so
    //
    //   rti[k]' * sk * rti[k] = diag(lambda_k)^{-1}
    //   rti[k]' * zk^{-1} * rti[k] = diag(lambda_k).
    //
    // The vectors lambda_k are stored in
    //
    //   lmbda[ dims['l'] + sum(dims['q']) : -1 ]
    if (W->rcount > 0) {
        cvxc_matrix_t Ls, Lz, wrk, r, rti;
        cvxc_memblk_t Tmp;
        cvxc_size_t j, ind2 = ind;
        cvxc_float_t a;

        // here work is big enough to hold 4*SIZEMAX(S) matrices
        for (k = 0; k < W->rcount; k++) {
            m = cvxc_scaling_elem(&r, W, CVXWS_R, k);
            cvxc_scaling_elem(&rti, W, CVXWS_RTI, k);

            cvxc_mblk_clear(work);
            // space for intermediate matrices
            cvxm_map_data(&wrk, m, m, cvxc_mblk_offset(work, 0));
            cvxm_map_data(&Ls,  m, m, cvxc_mblk_offset(work, m*m));
            cvxm_map_data(&Lz,  m, m, cvxc_mblk_offset(work, 2*m*m));
            cvxc_mblk_subblk(&Tmp, work, 3*m*m);
            //
            cvxc_mgrp_elem(&sk, s_g, CVXDIM_SDP, k);
            cvxc_mgrp_elem(&zk, z_g, CVXDIM_SDP, k);
            cvxc_mgrp_elem(&lk, lmbda_g, CVXDIM_SDP, k);
            // factor sk = Ls*Ls'
            cvxm_copy(&Ls, &sk, CVXC_LOWER);
            cvxm_cholfactor(&Ls, CVXC_LOWER);
            // factor zk = Lz*Lz'
            cvxm_copy(&Lz, &zk, CVXC_LOWER);
            cvxm_cholfactor(&Lz, CVXC_LOWER);
            // SVD(Lz'*Ls) = U*diag(lmbda_k)*V'  (r holds U, V' is not formed)
            cvxm_make_trm(&Ls, CVXC_LOWER);
            cvxm_copy(&wrk, &Ls, 0);
            cvxm_mult_trm(&wrk, 1.0, &Lz, CVXC_LEFT|CVXC_LOWER|CVXC_TRANS);
            cvxm_svd(&lk, &r, (cvxc_matrix_t *)0, &wrk, CVXC_WANTU/*|ARMAS_FORWARD*/, &Tmp);
            cvxm_copy(&rti, &r, CVXC_ALL);
            // r = Lz^-T * U
            cvxm_solve_trm(&r, 1.0, &Lz, CVXC_LEFT|CVXC_TRANS|CVXC_LOWER);
            // rti = Lz* U
            cvxm_mult_trm(&rti, 1.0, &Lz, CVXC_LEFT|CVXC_LOWER);
            // r = r * diag(sqrt(lmbda_k))
            // rti = rti * diag(1/sqrt(lmbda_k))
            for (j = 0; j < m; j++) {
                cvxc_matrix_t rk;
                a = SQRT(cvxm_get(&lk, j, 0));
                cvxm_view_map(&rk, &r, 0, j, m, 1);
                cvxm_scale(&rk, a, CVXC_ALL);
                cvxm_view_map(&rk, &rti, 0, j, m, 1);
                cvxm_scale(&rk, 1.0/a, CVXC_ALL);
            }
            //
            ind2 += m*m;
            ind  += m;
        }
    }
    return 0;
}

/*
 *   Updates the Nesterov-Todd scaling matrix W and the scaled variable 
 *   lmbda so that on exit
 *
 *         W * zt = W^{-T} * st = lmbda.
 *
 *   On entry, the nonlinear, 'l' and 'q' components of the arguments s 
 *   and z contain W^{-T}*st and W*zt, i.e, the new iterates in the current 
 *   scaling.
 *
 *   The 's' components contain the factors Ls, Lz in a factorization of 
 *   the new iterates in the current scaling, W^{-T}*st = Ls*Ls',   
 *   W*zt = Lz*Lz'.
 *
 */
int cvxc_update_scaling(cvxc_scaling_t *W,
                       cvxc_matgrp_t *lmbda_g,
                       cvxc_matgrp_t *s_g,
                       cvxc_matgrp_t *z_g,
                       cvxc_memblk_t *work)
{
    cvxc_size_t m; //ind = 0;
    cvxc_matrix_t sk, zk, lk, wrk, tmp;
    cvxc_matrix_t dnl, dnli, dl, dli;

    //   Nonlinear and 'l' blocks
    //
    //    d :=  d .* sqrt( s ./ z )
    //    lmbda := lmbda .* sqrt(s) .* sqrt(z)
    //
#if 1
    if  (W->dnltsz > 0) {
        cvxc_scaling_elem(&dnl,  W, CVXWS_DNLT, 0);
        cvxc_scaling_elem(&dnli, W, CVXWS_DNLTI, 0);
        // s = sqrt(s) ; z = sqrt(z);
        cvxc_mgrp_elem(&sk, s_g, CVXDIM_NLTARGET, 0);
        cvxc_mgrp_elem(&zk, z_g, CVXDIM_NLTARGET, 0);
        cvxm_apply(&sk, fsqrt, CVXC_ALL);
        cvxm_apply(&zk, fsqrt, CVXC_ALL);
        // lmbda = sqrt(s) .* sqrt(z)
        cvxc_mgrp_elem(&lk, lmbda_g, CVXDIM_NLTARGET, 0);
        cvxm_copy(&lk, &sk, CVXC_ALL);
        cvxm_mult_diag(&lk, 1.0, &zk, CVXC_LEFT);
        // dnl = dnl .* s .* z
        cvxm_mult_diag(&dnl, 1.0, &sk, CVXC_LEFT);
        cvxm_solve_diag(&dnl, 1.0, &zk, CVXC_LEFT);
        // dnli = 1.0 ./ dnl
        cvxm_copy(&dnli, &dnl, CVXC_ALL);
        cvxm_apply(&dnli, finv, CVXC_ALL);
    }
    if  (W->dnlsz > 0) {
        cvxc_scaling_elem(&dnl,  W, CVXWS_DNL, 0);
        cvxc_scaling_elem(&dnli, W, CVXWS_DNLI, 0);
        // s = sqrt(s) ; z = sqrt(z);
        cvxc_mgrp_elem(&sk, s_g, CVXDIM_NONLINEAR, 0);
        cvxc_mgrp_elem(&zk, z_g, CVXDIM_NONLINEAR, 0);
        cvxm_apply(&sk, fsqrt, CVXC_ALL);
        cvxm_apply(&zk, fsqrt, CVXC_ALL);
        // lmbda = sqrt(s) .* sqrt(z)
        cvxc_mgrp_elem(&lk, lmbda_g, CVXDIM_NONLINEAR, 0);
        cvxm_copy(&lk, &sk, CVXC_ALL);
        cvxm_mult_diag(&lk, 1.0, &zk, CVXC_LEFT);
        // dnl = dnl .* s .* z
        cvxm_mult_diag(&dnl, 1.0, &sk, CVXC_LEFT);
        cvxm_solve_diag(&dnl, 1.0, &zk, CVXC_LEFT);
        // dnli = 1.0 ./ dnl
        cvxm_copy(&dnli, &dnl, CVXC_ALL);
        cvxm_apply(&dnli, finv, CVXC_ALL);
    }
    if (W->dsz > 0) {
        cvxc_scaling_elem(&dl,  W, CVXWS_D, 0);
        cvxc_scaling_elem(&dli, W, CVXWS_DI, 0);
        // s = sqrt(s) ; z = sqrt(z);
        cvxc_mgrp_elem(&sk, s_g, CVXDIM_LINEAR, 0);
        cvxc_mgrp_elem(&zk, z_g, CVXDIM_LINEAR, 0);
        cvxm_apply(&sk, fsqrt, CVXC_ALL);
        cvxm_apply(&zk, fsqrt, CVXC_ALL);
        // lmbda = sqrt(s) .* sqrt(z)
        cvxc_mgrp_elem(&lk, lmbda_g, CVXDIM_LINEAR, 0);
        cvxm_copy(&lk, &sk, CVXC_ALL);
        cvxm_mult_diag(&lk, 1.0, &zk, CVXC_LEFT);
        // dl = dl .* s ./ z
        cvxm_mult_diag(&dl, 1.0, &sk, CVXC_LEFT);
        cvxm_solve_diag(&dl, 1.0, &zk, CVXC_LEFT);
        // dli = 1.0 ./ dl
        cvxm_copy(&dli, &dl, CVXC_ALL);
        cvxm_apply(&dli, finv, CVXC_ALL);
    }
#endif
#if 0
    cvxc_matrix_t *s = s_g->mat, *z = z_g->mat, *lmbda = lmbda_g->mat;
    if (W->dsz + W->dnlsz > 0) {
        if (W->dnlsz > 0) {
            cvxc_scaling_elem(&dnl,  W, CVXWS_DNL, 0);
            cvxc_scaling_elem(&dnli, W, CVXWS_DNLI, 0);
        }
        cvxc_scaling_elem(&dl,  W, CVXWS_D, 0);
        cvxc_scaling_elem(&dli, W, CVXWS_DI, 0);
        // s = sqrt(s) ; z = sqrt(z);
        cvxm_view_map(&sk, s, 0, 0, W->dnlsz + W->dsz, 1);
        cvxm_view_map(&zk, z, 0, 0, W->dnlsz + W->dsz, 1);
        cvxm_apply(&sk, fsqrt, CVXC_ALL);
        cvxm_apply(&zk, fsqrt, CVXC_ALL);
        // lmbda = sqrt(s) .* sqrt(z)
        cvxm_view_map(&lk, lmbda, 0, 0, W->dnlsz + W->dsz, 1);
        cvxm_copy(&lk, &sk, CVXC_ALL);
        cvxm_mult_diag(&lk, 1.0, &zk, CVXC_LEFT);
        // update scaling matrices d, di, dnl, dnli
        if (W->dnlsz > 0) {
            cvxc_mgrp_elem(&sk, s_g, CVXDIM_NONLINEAR, 0);
            cvxc_mgrp_elem(&zk, z_g, CVXDIM_NONLINEAR, 0);
            // dnl = dnl .* s .* z
            cvxm_mult_diag(&dnl, 1.0, &sk, CVXC_LEFT);
            cvxm_solve_diag(&dnl, 1.0, &zk, CVXC_LEFT);
            // dnli = 1.0 ./ dnl
            cvxm_copy(&dnli, &dnl, CVXC_ALL);
            cvxm_apply(&dnli, finv, CVXC_ALL);
            //ind = W->dnlsz;
        }
        cvxc_mgrp_elem(&sk, s_g, CVXDIM_LINEAR, 0);
        cvxc_mgrp_elem(&zk, z_g, CVXDIM_LINEAR, 0);
        // dl = dl .* s ./ z
        cvxm_mult_diag(&dl, 1.0, &sk, CVXC_LEFT);
        cvxm_solve_diag(&dl, 1.0, &zk, CVXC_LEFT);
        // dli = 1.0 ./ dl
        cvxm_copy(&dli, &dl, CVXC_ALL);
        cvxm_apply(&dli, finv, CVXC_ALL);
        //ind += W->dsz;
    }
#endif
    // 'q' blocks.
    // Let st and zt be the new variables in the old scaling:
    //
    //     st = s_k,   zt = z_k
    //
    // and a = sqrt(st' * J * st),  b = sqrt(zt' * J * zt).
    //
    // 1. Compute the hyperbolic Householder transformation 2*q*q' - J
    //    that maps st/a to zt/b.
    //
    //        c = sqrt( (1 + st'*zt/(a*b)) / 2 )
    //        q = (st/a + J*zt/b) / (2*c).
    //
    //    The new scaling point is
    //
    //        wk := betak * sqrt(a/b) * (2*v[k]*v[k]' - J) * q
    //
    //    with betak = W['beta'][k].
    //
    // 3. The scaled variable:
    //
    //        lambda_k0 = sqrt(a*b) * c
    //        lambda_k1 = sqrt(a*b) * ( (2vk*vk' - J) * (-d*q + u/2) )_1
    //
    //    where
    //
    //        u = st/a - J*zt/b
    //        d = ( vk0 * (vk'*u) + u0/2 ) / (2*vk0 *(vk'*q) - q0 + 1).
    //
    // 4. Update scaling
    //
    //        v[k] := wk^1/2
    //              = 1 / sqrt(2*(wk0 + 1)) * (wk + e).
    //        beta[k] *=  sqrt(a/b)
    if (W->vcount > 0) {
        cvxc_matrix_t v, beta;
        cvxc_float_t a, b, c, d, vs, vz, vq, vu, wk0, v0;

        cvxc_scaling_elem(&beta, W, CVXWS_BETA, 0);
        for (int k = 0; k < W->vcount; k++) {
            m = cvxc_scaling_elem(&v, W, CVXWS_V, k);
            cvxc_mgrp_elem(&sk, s_g, CVXDIM_SOCP, k);
            cvxc_mgrp_elem(&zk, z_g, CVXDIM_SOCP, k);

            // a = sqrt( sk' * J * sk ); sk = sk/a
            a = cvxc_jnrm2(&sk);
            cvxm_scale(&sk, 1.0/a, CVXC_ALL);
            // b = sqrt( zk' * J * zk ); zk = zk/b
            b = cvxc_jnrm2(&zk);
            cvxm_scale(&zk, 1.0/b, CVXC_ALL);
            // c = dot(sk', zk)
            c = cvxm_dot(&sk, &zk);
            c = SQRT((1.0 + c)/2.0);
            // vs = v' * sk/a
            vs = cvxm_dot(&v, &sk);
            // vz = v' * J *zk/b
            vz = cvxc_jdot(&v, &zk);
            // vq = v' * q where q = (st/a + J * zt/b) / (2 * c)
            vq = (vs + vz) / (2.0*c);
            vu = vs - vz;

            cvxc_mgrp_elem(&lk, lmbda_g, CVXDIM_SOCP, k);

            // wk0 = 2 * vk0 * (vk' * q) - q0 
            wk0 = 2.0*cvxm_get(&v, 0, 0)*vq
                - (cvxm_get(&sk, 0, 0) + cvxm_get(&zk, 0, 0))/(2.0*c);
            // d = (v[0] * (vk' * u) - u0/2) / (wk0 + 1)
            d = (vu*cvxm_get(&v, 0, 0)
                 - (cvxm_get(&sk, 0, 0) - cvxm_get(&zk, 0, 0))/2.0) / (wk0 + 1.0);

            // lambda_k1 = 2 * v_k1 * vk' * (-d*q + u/2) - d*q1 + u1/2
            cvxm_copy(&lk, &v, CVXC_ALL);
            cvxm_scale(&lk, (2.0 * (-d*vq + 0.5*vu)), CVXC_ALL);
            cvxm_axpy(&lk, 0.5*(1.0 - d/c), &sk);
            cvxm_axpy(&lk, 0.5*(1.0 + d/c), &zk);
            // set lmbda[0] = c
            cvxm_set(&lk, 0, 0, c);
            // Scale so that sqrt(lambda_k' * J * lambda_k) = sqrt(aa*bb).
            cvxm_scale(&lk, SQRT(a*b), CVXC_ALL);

            // v := (2*v*v' - J) * q 
            //    = 2 * (v'*q) * v' - (J* st/a + zt/b) / (2*c)
            cvxm_scale(&v, 2.0*vq, CVXC_ALL);
            v0 = cvxm_get(&v, 0, 0);
            cvxm_axpy(&v, 0.5/c, &sk);
            cvxm_set(&v, 0, 0, v0 - (cvxm_get(&sk, 0, 0)/(2.0*c)));
            cvxm_axpy(&v, -0.5/c, &zk);

            // v := v^{1/2} = 1/sqrt(2 * (v0 + 1)) * (v + e)
            v0 = cvxm_get(&v, 0, 0) + 1.0;
            cvxm_set(&v, 0, 0, v0);
            // v := 1.0/SQRT(2*v) * v
            cvxm_scale(&v, INVSQRT2/SQRT(v0), CVXC_ALL);

            // beta[k] *= ( aa / bb )**1/2
            cvxm_set(&beta, k, 0, SQRT(a/b)*cvxm_get(&beta, k, 0));

            //ind += m;
       }
    }

    // 's' blocks
    //
    // Let st, zt be the updated variables in the old scaling:
    //
    //     st = Ls * Ls', zt = Lz * Lz'.
    //
    // where Ls and Lz are the 's' components of s, z.
    //
    // 1.  SVD Lz'*Ls = Uk * lambda_k^+ * Vk'.
    //
    // 2.  New scaling is
    //
    //         r[k] := r[k] * Ls * Vk * diag(lambda_k^+)^{-1/2}
    //         rti[k] := r[k] * Lz * Uk * diag(lambda_k^+)^{-1/2}.
    //

    if (W->rcount) {
        cvxc_matrix_t r, rti, Ls, Lz;
        cvxc_float_t a;
        cvxc_size_t k;
        cvxc_memblk_t Tmp;

        for (k = 0; k < W->rcount; k++) {
            cvxc_size_t j;
            m = cvxc_scaling_elem(&r, W, CVXWS_R, k);
            cvxc_scaling_elem(&rti, W, CVXWS_RTI, k);

            // Ls, Lz triangular
            cvxc_mgrp_elem(&Ls, s_g, CVXDIM_SDP, k);
            cvxc_mgrp_elem(&Lz, z_g, CVXDIM_SDP, k);
            cvxc_mgrp_elem(&lk, lmbda_g, CVXDIM_SDP, k);
            // workspaces
            cvxc_mblk_clear(work);
            cvxm_map_data(&wrk, m, m, cvxc_mblk_offset(work, 0));
            cvxc_mblk_subblk(&Tmp, work, m*m);
            // t = r*sk = r*Ls
            cvxm_mult(0.0, &wrk, 1.0, &r, &Ls, 0);
            cvxm_copy(&r, &wrk, CVXC_ALL);
            // rti = rti*Lz
            cvxm_mult(0.0, &wrk, 1.0, &rti, &Lz, 0);
            cvxm_copy(&rti, &wrk, CVXC_ALL);
            // SVD(Lz'*Ls) = U * lmbds^+ * V'; store U in sk and V' in zk. '
            cvxm_mult(0.0, &wrk, 1.0, &Lz, &Ls, CVXC_TRANS);
            // U = s; Vt = z
            cvxm_scale(&Ls, 0.0, 0);
            cvxm_scale(&Lz, 0.0, 0);
            cvxm_svd(&lk, &Ls, &Lz, &wrk, CVXC_WANTU|CVXC_WANTV, &Tmp);
            // r = r*V
            cvxm_mult(0.0, &wrk, 1.0, &r, &Lz, CVXC_TRANSB);
            cvxm_copy(&r, &wrk, CVXC_ALL);
            // rti = rti*U
            cvxm_mult(0.0, &wrk, 1.0, &rti, &Ls, 0);
            cvxm_copy(&rti, &wrk, CVXC_ALL);
            // scale r, rti with diag(1.0./sqrt(lmbda)); columnwise
            for (j = 0; j < m; j++) {
                a = 1.0 / SQRT(cvxm_get(&lk, j, 0));
                cvxm_view_map(&tmp, &r, 0, j, m, 1);
                cvxm_scale(&tmp, a, CVXC_ALL);
                cvxm_view_map(&tmp, &rti, 0, j, m, 1);
                cvxm_scale(&tmp, a, CVXC_ALL);
            }
        }
    }
    return 0;
}


//   Evaluates
//
//       x := H(lambda^{1/2}) * x   (CVXC_INV not set in flags)
//       x := H(lambda^{-1/2}) * x  (CVXC_INV is set in flags).
//
//   H is the Hessian of the logarithmic barrier.
int cvxc_scale2(cvxc_matgrp_t *x_g,
               cvxc_matgrp_t *lmbda_g,
               int flags,
               cvxc_memblk_t *work)
{
    cvxc_matrix_t xk, lk, x1, l1, lc, xcol;
    cvxc_index_t *index = x_g->index;
    cvxc_size_t m;
    cvxc_float_t lx, a, c, x0, l0;

    // For the nonlinear and 'l' blocks,
    //
    //     xk := xk ./ lk   (CVXC_INV is not set)
    //     xk := xk .* lk   (CVXC_INV is set)
    //
    if (cvxc_mgrp_count(x_g, CVXDIM_NLTARGET) > 0) {
        m = cvxc_mgrp_elem(&xk, x_g, CVXDIM_NLTARGET, 0);
        cvxc_mgrp_elem(&lk, lmbda_g, CVXDIM_NLTARGET, 0);
        if ((flags & CVXC_INV) != 0) {
            cvxm_mult_diag(&xk, 1.0, &lk, 0);
        } else {
            cvxm_solve_diag(&xk, 1.0, &lk, 0);
        }
    }
 #if 0
    if (index->indnlt) {
        m = cvxc_mgrp_elem(&xk, x_g, CVXDIM_NLTARGET, 0);
        cvxc_mgrp_elem(&lk, lmbda_g, CVXDIM_NLTARGET, 0);
        if ((flags & CVXC_INV) != 0) {
            cvxm_mult_diag(&xk, 1.0, &lk, 0);
        } else {
            cvxm_solve_diag(&xk, 1.0, &lk, 0);
        }
    }
#endif
    if (index->indnl) {
        m = cvxc_mgrp_elem(&xk, x_g, CVXDIM_NONLINEAR, 0);
        cvxc_mgrp_elem(&lk, lmbda_g, CVXDIM_NONLINEAR, 0);
        if ((flags & CVXC_INV) != 0) {
            cvxm_mult_diag(&xk, 1.0, &lk, 0);
        } else {
            cvxm_solve_diag(&xk, 1.0, &lk, 0);
        }
    }

    if (index->indl) {
        m = cvxc_mgrp_elem(&xk, x_g, CVXDIM_LINEAR, 0);
        cvxc_mgrp_elem(&lk, lmbda_g, CVXDIM_LINEAR, 0);
        if ((flags & CVXC_INV) != 0) {
            cvxm_mult_diag(&xk, 1.0, &lk, 0);
        } else {
            cvxm_solve_diag(&xk, 1.0, &lk, 0);
        }
    }

    // For 'q' blocks, if CVXC_INV is not set,
    //
    //     xk := 1/a * [ l'*J*xk;
    //         xk[1:] - (xk[0] + l'*J*xk) / (l[0] + 1) * l[1:] ].
    //
    // If CVXC_INV is set,
    //
    //     xk := a * [ l'*xk;
    //         xk[1:] + (xk[0] + l'*xk) / (l[0] + 1) * l[1:] ].
    //
    // a = sqrt(lambda_k' * J * lambda_k), l = lambda_k / a.
    if (index->indq) {
        for (int k = 0; k < cvxc_mgrp_count(x_g, CVXDIM_SOCP); k++) {
            m = cvxc_mgrp_elem(&xk, x_g, CVXDIM_SOCP, k);
            cvxc_mgrp_elem(&lk, lmbda_g, CVXDIM_SOCP, k);
            a = cvxc_jnrm2(&lk);
            if ((flags & CVXC_INV) == 0) {
                lx = cvxc_jdot(&lk, &xk);
                lx /= a;
            } else {
                lx = cvxm_dot(&lk, &xk);
                lx /= a;
            }
            x0 = cvxm_get(&xk, 0, 0);
            l0 = cvxm_get(&lk, 0, 0);
            cvxm_set(&xk, 0, 0, lx);
            c = (lx + x0) / (l0/a + 1.0) / a;
            if ((flags & CVXC_INV) == 0) {
                c = -c;
                a = 1.0 /a;
            }
            cvxm_view_map(&x1, &xk, 1, 0, m-1, 1);
            cvxm_view_map(&l1, &lk, 1, 0, m-1, 1);
            cvxm_axpy(&x1, c, &l1);
            cvxm_scale(&xk, a, CVXC_ALL);
        }
    }

    // For the 's' blocks, if CVXC_INV is not set
    //
    //     xk := vec( diag(l)^{-1/2} * mat(xk) * diag(k)^{-1/2}).
    //
    // If CVXC_INV is set,
    //
    //     xk := vec( diag(l)^{1/2} * mat(xk) * diag(k)^{1/2}).
    //
    // where l is k'th block of lambda.
    // 
    // We scale upper and lower triangular part of mat(xk) because the
    // inverse operation will be applied to nonsymmetric matrices.
    if (index->inds) {
        for (int k = 0; k < cvxc_mgrp_count(x_g, CVXDIM_SDP); k++) {
            m = cvxc_mgrp_elem(&xk, x_g, CVXDIM_SDP, k);
            cvxc_mgrp_elem(&lk, lmbda_g, CVXDIM_SDP, k);
            cvxm_map_data(&lc, m, 1, cvxc_mblk_offset(work, 0));
            for (cvxc_size_t i = 0; i < m; i++) {
                lx = SQRT(cvxm_get(&lk, i, 0));
                for (cvxc_size_t j = 0; j < m; j++) {   // changed:
                    c = lx * SQRT(cvxm_get(&lk, j, 0));
                    cvxm_set(&lc, j, 0, c);
                }
                cvxm_view_map(&xcol, &xk, 0, i, m, 1);
                if ((flags & CVXC_INV) == 0) {
                    cvxm_solve_diag(&xcol, 1.0, &lc, CVXC_RIGHT);
                } else {
                    cvxm_mult_diag(&xcol, 1.0, &lc, CVXC_RIGHT);
                }
            }
        }
    }
    return 0;
}
