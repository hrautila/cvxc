
// Copyright: Harri Rautila, 2018 <harri.rautila@gmail.com>

#include "convex.h"

void cvx_mgrp_copy(cvx_matgrp_t *x_g, cvx_matgrp_t *y_g, cvx_dim_enum name)
{
    cvx_matrix_t xs, ys;
    int k;

    switch (name) {
    case CVXDIM_NLTARGET:
        cvx_mgrp_elem(&xs, x_g, CVXDIM_NLTARGET, 0);
        cvx_mgrp_elem(&ys, y_g, CVXDIM_NLTARGET, 0);
        cvxm_copy(&xs, &ys, 0);
        break;
        break;
    case CVXDIM_NONLINEAR:
        cvx_mgrp_elem(&xs, x_g, CVXDIM_LINEAR, 0);
        cvx_mgrp_elem(&ys, y_g, CVXDIM_LINEAR, 0);
        cvxm_copy(&xs, &ys, 0);
        break;
    case CVXDIM_LINEAR:
        cvx_mgrp_elem(&xs, x_g, CVXDIM_LINEAR, 0);
        cvx_mgrp_elem(&ys, y_g, CVXDIM_LINEAR, 0);
        cvxm_copy(&xs, &ys, 0);
        break;
    case CVXDIM_SOCP:
        for (k = 0; k < cvx_mgrp_count(x_g, CVXDIM_SOCP); k++) {
            cvx_mgrp_elem(&xs, x_g, CVXDIM_SOCP, k);
            cvx_mgrp_elem(&ys, y_g, CVXDIM_SOCP, k);
            cvxm_copy(&xs, &ys, 0);
        }
        break;
    case CVXDIM_SDP:
        for (k = 0; k < cvx_mgrp_count(x_g, CVXDIM_SDP); k++) {
            cvx_mgrp_elem(&xs, x_g, CVXDIM_SDP, k);
            cvx_mgrp_elem(&ys, y_g, CVXDIM_SDP, k);
            cvxm_copy(&xs, &ys, 0);
        }
        break;
    case CVXDIM_CONVEX:
        cvx_mgrp_copy(x_g, y_g, CVXDIM_NLTARGET);
        cvx_mgrp_copy(x_g, y_g, CVXDIM_NONLINEAR);
        break;
    case CVXDIM_CONELP:
        cvx_mgrp_copy(x_g, y_g, CVXDIM_LINEAR);
        cvx_mgrp_copy(x_g, y_g, CVXDIM_SOCP);
        cvx_mgrp_copy(x_g, y_g, CVXDIM_SDP);
        break;
    case CVXDIM_CONVEXLP:
        cvx_mgrp_copy(x_g, y_g, CVXDIM_NONLINEAR);
        cvx_mgrp_copy(x_g, y_g, CVXDIM_LINEAR);
        cvx_mgrp_copy(x_g, y_g, CVXDIM_SOCP);
        cvx_mgrp_copy(x_g, y_g, CVXDIM_SDP);
        break;
    case CVXDIM_CONVEXPROG:
        cvx_mgrp_copy(x_g, y_g, CVXDIM_NLTARGET);
        cvx_mgrp_copy(x_g, y_g, CVXDIM_NONLINEAR);
        cvx_mgrp_copy(x_g, y_g, CVXDIM_LINEAR);
        cvx_mgrp_copy(x_g, y_g, CVXDIM_SOCP);
        cvx_mgrp_copy(x_g, y_g, CVXDIM_SDP);
        break;
    default:
        break;
    }
}

/*
 * Convert lower triangular matrices in 'S' to symmetric. Fill in the strictly upper triangular part
 * of the symmetric matrix in x with strictly lower triangular part
 */
int cvx_mgrp_mksymm(cvx_matgrp_t *x_g)
{
    cvx_matrix_t xs;

    for (int k = 0; k < cvx_mgrp_count(x_g, CVXDIM_SDP); k++) {
        int n = cvx_mgrp_elem(&xs, x_g, CVXDIM_SDP, k);
        cvxm_mksymm(&xs, n);
    }
    return 0;
}

// \brief Copy lambda vector with diagonal 'S' space to s/z with vectors with standard 'S' storage
void cvx_mgrp_copy_lambda(cvx_matgrp_t *ds_g, cvx_matgrp_t *lmbda_g)
{
    cvx_matrix_t dk, lk, d;

    const cvx_dimset_t *dims = ds_g->index->dims;
    cvx_size_t n =
        cvx_dimset_sum(dims, CVXDIM_NLTARGET) +
        cvx_dimset_sum(dims, CVXDIM_NONLINEAR) +
        cvx_dimset_sum(dims, CVXDIM_LINEAR) +
        cvx_dimset_sum(dims, CVXDIM_SOCP);

    cvxm_map_data(&dk, n, 1, cvxm_data(ds_g->mat, 0));
    cvxm_map_data(&lk, n, 1, cvxm_data(lmbda_g->mat, 0));
    cvxm_copy(&dk, &lk, CVX_ALL);

    for (int k = 0; k < cvx_mgrp_count(ds_g, CVXDIM_SDP); k++) {
        cvx_mgrp_elem(&dk, ds_g, CVXDIM_SDP, k);
        cvx_mgrp_elem(&lk, lmbda_g, CVXDIM_SDP, k);
        // clear it
        cvxm_mkconst(&dk, 0.0);
        // copy to diagonal
        cvxm_view_diag(&d, &dk, 0);
        cvxm_copy(&d, &lk, 0);
    }
}

// \brief Update matrix group by adding parameter value to its elements
void cvx_mgrp_update_sz(cvx_matgrp_t *ds_g, cvx_float_t val, int flags)
{
    cvx_matrix_t x, xd;
    if (flags == 0 || (flags & CVXDIM_NLTARGET) != 0) {
        // non-linear target function part; NLT => NLT +/ val
        cvx_mgrp_elem(&x, ds_g, CVXDIM_NLTARGET, 0);
        cvxm_add(&x, val, 0);
    }
    if (flags == 0 || (flags & CVXDIM_NONLINEAR) != 0) {
        // non-linear part; NL => NL +/ val
        cvx_mgrp_elem(&x, ds_g, CVXDIM_NONLINEAR, 0);
        cvxm_add(&x, val, 0);
    }
    if (flags == 0 || (flags & CVXDIM_LINEAR) != 0) {
        // linear part; L => L +/ val
        cvx_mgrp_elem(&x, ds_g, CVXDIM_LINEAR, 0);
        cvxm_add(&x, val, 0);
    }

    if (flags == 0 || (flags & CVXDIM_SOCP) != 0) {
        // SOCP part; (y0, y1) => (y0+val, y1)
        for (int k = 0; k < cvx_mgrp_count(ds_g, CVXDIM_SOCP); k++) {
            cvx_mgrp_elem(&x, ds_g, CVXDIM_SOCP, k);
            cvxm_set(&x, 0, 0, cvxm_get(&x, 0, 0) + val);
        }
    }

    if (flags == 0 || (flags & CVXDIM_SDP) != 0) {
        // SDP part; diag(S) => diag(S) +/ val
        for (int k = 0; k < cvx_mgrp_count(ds_g, CVXDIM_SDP); k++) {
            cvx_mgrp_elem(&x, ds_g, CVXDIM_SDP, k);
            cvxm_view_diag(&xd, &x, 0);
            cvxm_add(&xd, val, 0);
        }
    }
}

// \brief Update matrix group by scaling its elements by parameter value
void cvx_mgrp_scale_sz(cvx_matgrp_t *ds_g, cvx_float_t val, int flags)
{
    cvx_matrix_t x, xd;
   if (flags == 0 || (flags & CVXDIM_NLTARGET) != 0) {
        // non-linear part; NL => NL * val
        cvx_mgrp_elem(&x, ds_g, CVXDIM_NLTARGET, 0);
        cvxm_scale(&x, val, 0);
    }
   if (flags == 0 || (flags & CVXDIM_NONLINEAR) != 0) {
        // non-linear part; NL => NL * val
        cvx_mgrp_elem(&x, ds_g, CVXDIM_NONLINEAR, 0);
        cvxm_scale(&x, val, 0);
    }
    if (flags == 0 || (flags & CVXDIM_LINEAR) != 0) {
        // linear part; L => L * val
        cvx_mgrp_elem(&x, ds_g, CVXDIM_LINEAR, 0);
        cvxm_scale(&x, val, 0);
    }

    if (flags == 0 || (flags & CVXDIM_SOCP) != 0) {
        // SOCP part; (y0, y1) => (y0*val, y1*val)
        for (int k = 0; k < cvx_mgrp_count(ds_g, CVXDIM_SOCP); k++) {
            cvx_mgrp_elem(&x, ds_g, CVXDIM_SOCP, k);
            cvxm_scale(&x, val, 0);
        }
    }

    if (flags == 0 || (flags & CVXDIM_SDP) != 0) {
        // SDP part; diag(S) => diag(S) */ val
        for (int k = 0; k < cvx_mgrp_count(ds_g, CVXDIM_SDP); k++) {
            cvx_mgrp_elem(&x, ds_g, CVXDIM_SDP, k);
            cvxm_view_diag(&xd, &x, 0);
            cvxm_scale(&xd, val, 0);
        }
    }
}

// \brief Update matrix group: $f y_g = beta*y_g + alpha*x_g $f
void cvx_mgrp_axpby_sz(cvx_float_t beta, cvx_matgrp_t *y_g, cvx_float_t alpha, cvx_matgrp_t *x_g, int flags)
{
    cvx_matrix_t y, x;
    if (flags == 0 || flags == CVXDIM_CONVEXLP) {
        cvx_index_elem(&x, x_g->mat, x_g->index, CVXDIM_CONVEXLP, 0);
        cvx_index_elem(&y, y_g->mat, y_g->index, CVXDIM_CONVEXLP, 0);
        cvxm_axpby(beta, &y, alpha, &x);
        return;
    }
    if (flags == 0 || (flags & CVXDIM_NLTARGET) != 0) {
        // non-linear part; NL => NL * val
        cvx_mgrp_elem(&x, x_g, CVXDIM_NLTARGET, 0);
        cvx_mgrp_elem(&y, y_g, CVXDIM_NLTARGET, 0);
        cvxm_axpby(beta, &y, alpha, &x);
    }
    if (flags == 0 || (flags & CVXDIM_NONLINEAR) != 0) {
        // non-linear part; NL => NL * val
        cvx_mgrp_elem(&x, x_g, CVXDIM_NONLINEAR, 0);
        cvx_mgrp_elem(&y, y_g, CVXDIM_NONLINEAR, 0);
        cvxm_axpby(beta, &y, alpha, &x);
    }
    if (flags == 0 || (flags & CVXDIM_LINEAR) != 0) {
        // linear part; L => L * val
        cvx_mgrp_elem(&x, x_g, CVXDIM_LINEAR, 0);
        cvx_mgrp_elem(&y, y_g, CVXDIM_LINEAR, 0);
        cvxm_axpby(beta, &y, alpha, &x);
    }

    if (flags == 0 || (flags & CVXDIM_SOCP) != 0) {
        // SOCP part; (y0, y1) => (y0*val, y1*val)
        for (int k = 0; k < cvx_mgrp_count(x_g, CVXDIM_SOCP); k++) {
            cvx_mgrp_elem(&x, x_g, CVXDIM_SOCP, k);
            cvx_mgrp_elem(&y, y_g, CVXDIM_SOCP, k);
            cvxm_axpby(beta, &y, alpha, &x);
        }
    }

    if (flags == 0 || (flags & CVXDIM_SDP) != 0) {
        // SDP part; S := S * val
        cvx_matrix_t xv, yv;
        cvx_size_t m;
        for (int k = 0; k < cvx_mgrp_count(x_g, CVXDIM_SDP); k++) {
            m = cvx_mgrp_elem(&x, x_g, CVXDIM_SDP, k);
            cvx_mgrp_elem(&y, y_g, CVXDIM_SDP, k);
            cvxm_map_data(&xv, m*m, 1, cvxm_data(&x, 0));
            cvxm_map_data(&yv, m*m, 1, cvxm_data(&y, 0));
            cvxm_axpby(beta, &yv, alpha, &xv);
        }
    }
}

void cvx_mgrp_initial_value(cvx_matgrp_t *x_g, int flags)
{
    cvx_matrix_t x, xd;
    if (flags == 0 || (flags & CVXDIM_NLTARGET) != 0) {
        cvx_mgrp_elem(&x, x_g, CVXDIM_NLTARGET, 0);
        cvxm_set_all(&x, 1.0);
    }
    if (flags == 0 || (flags & CVXDIM_NONLINEAR) != 0) {
        // nonlinear part; L := (1, 1, ...)
        cvx_mgrp_elem(&x, x_g, CVXDIM_NONLINEAR, 0);
        cvxm_set_all(&x, 1.0);
    }
    if (flags == 0 || (flags & CVXDIM_LINEAR) != 0) {
        // linear part; L := (1, 1, ...)
        cvx_mgrp_elem(&x, x_g, CVXDIM_LINEAR, 0);
        cvxm_set_all(&x, 1.0);
    }
    if (flags == 0 || (flags & CVXDIM_SOCP) != 0) {
        // SOCP part; (y0, y1) := (1, 0, ...)
        for (int k = 0; k < cvx_mgrp_count(x_g, CVXDIM_SOCP); k++) {
            cvx_mgrp_elem(&x, x_g, CVXDIM_SOCP, k);
            cvxm_unit_vector(&x);
        }
    }

    if (flags == 0 || (flags & CVXDIM_SDP) != 0) {
        // SDP part; diag(S) := diag(1) 
        for (int k = 0; k < cvx_mgrp_count(x_g, CVXDIM_SDP); k++) {
            cvx_mgrp_elem(&x, x_g, CVXDIM_SDP, k);
            cvxm_view_diag(&xd, &x, 0);
            cvxm_set_all(&xd, 1.0);
        }
    }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
