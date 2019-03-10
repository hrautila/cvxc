
// Copyright: Harri Rautila, 2018 <harri.rautila@gmail.com>

#include <stdio.h>
#include "convex.h"


/*
 * Compute x^T * J * y, where J = [1, 0; 0, -I]
 *
 *   (x0 X1^T) (1  0) (y0)  = (x0 -X1^T) (y0)  = x0*y0 - X1^T*Y1
 *             (0 -I) (Y1)               (Y1)
 */
cvx_float_t cvx_jdot(cvx_matrix_t *x, cvx_matrix_t *y)
{
    cvx_float_t a, x0, y0;
    cvx_matrix_t x1, y1;
    cvx_size_t nr, nc;

    cvxm_size(&nr, &nc, x);

    cvxm_view_map(&x1, x, 1, 0, nr-1, 1);
    cvxm_view_map(&y1, y, 1, 0, nr-1, 1);

    a  = cvxm_dot(&x1, &y1);
    x0 = cvxm_get(x, 0, 0);
    y0 = cvxm_get(y, 0, 0);

    return x0*y0 - a;
}


/*
 * Compute sqrt(x^T * J * x), where J = [1, 0; 0, -I]
 *
 *   (x0 X1^T) (1  0) (x0)  = (x0 -X1^T) (x0)  = sqrt(x0*x0 - X1^T*X1) = 
 *             (0 -I) (X1)               (X1)
 *
 *   sqrt( x0*x0 - |X1|*|X1| ) = sqrt(x0 - |X1|)*sqrt(x0 + |X1|)
 */
cvx_float_t cvx_jnrm2(cvx_matrix_t *x)
{
    cvx_float_t a, x0;
    cvx_matrix_t x1;
    cvx_size_t nr, nc;

    cvxm_size(&nr, &nc, x);
    cvxm_view_map(&x1, x, 1, 0, nr-1, 1);
    a  = cvxm_nrm2(&x1);
    x0 = cvxm_get(x, 0, 0);
    return SQRT(x0 - a)*SQRT(x0 + a);
}


/*
 * Compute inner product of two vectors in S; S part is assumed to be in lower
 * storage format.
 */
cvx_float_t cvx_sdot(cvx_matgrp_t *x_g, cvx_matgrp_t *y_g)
{
    cvx_size_t ind;
    cvx_matrix_t x0, y0, xs, ys, *x = x_g->mat, *y = y_g->mat;
    cvx_float_t sdot = 0.0;
    int k;
    cvx_size_t j;

    ind = x_g->index->dims->mnl +
        cvx_dimset_sum(x_g->index->dims, CVXDIM_LINEAR) +
        cvx_dimset_sum(x_g->index->dims, CVXDIM_SOCP);

    cvxm_view_map(&x0, x, 0, 0, ind, 1);
    cvxm_view_map(&y0, y, 0, 0, ind, 1);
    sdot += cvxm_dot(&x0, &y0);

    for (k = 0; k < cvx_mgrp_count(x_g, CVXDIM_SDP); k++) {
        cvx_size_t m = cvx_mgrp_elem(&xs, x_g, CVXDIM_SDP, k);
        cvx_mgrp_elem(&ys, y_g, CVXDIM_SDP, k);

        cvxm_view_diag(&x0, &xs, 0);
        cvxm_view_diag(&y0, &ys, 0);
        sdot += cvxm_dot(&x0, &y0);

        for (j = 1; j < m; j++) {
            // j'th subdiagonal
            cvxm_view_diag(&x0, &xs, -j);
            cvxm_view_diag(&y0, &ys, -j);
            sdot += 2.0 * cvxm_dot(&x0, &y0);
        }
    }
    return sdot;
}

/*
 * Compute the norm of vector in S
 */
cvx_float_t cvx_snrm2(cvx_matgrp_t *x_g)
{
    return sqrt(cvx_sdot(x_g, x_g));
}

/*
 * Compute inner product of two vectors in subset of S; SDP part is assumed to be in lower
 * storage format.
 */
cvx_float_t cvx_sdot_elem(cvx_matgrp_t *x_g, cvx_matgrp_t *y_g, cvx_dim_enum name)
{
    cvx_matrix_t x0, y0, xs, ys;
    cvx_float_t sdot = 0.0;
    int k;
    cvx_size_t j;

    switch (name) {
    case CVXDIM_NLTARGET:
    case CVXDIM_CONVEX:
        break;
    case CVXDIM_NONLINEAR:
        cvx_mgrp_elem(&xs, x_g, CVXDIM_LINEAR, 0);
        cvx_mgrp_elem(&ys, y_g, CVXDIM_LINEAR, 0);
        sdot += cvxm_dot(&xs, &ys);
        break;
    case CVXDIM_LINEAR:
        cvx_mgrp_elem(&xs, x_g, CVXDIM_LINEAR, 0);
        cvx_mgrp_elem(&ys, y_g, CVXDIM_LINEAR, 0);
        sdot += cvxm_dot(&xs, &ys);
        break;
    case CVXDIM_SOCP:
        for (k = 0; k < cvx_mgrp_count(x_g, CVXDIM_SOCP); k++) {
            cvx_mgrp_elem(&xs, x_g, CVXDIM_SOCP, k);
            cvx_mgrp_elem(&ys, y_g, CVXDIM_SOCP, k);
            sdot += cvxm_dot(&xs, &ys);
        }
        break;
    case CVXDIM_SDP:
        for (k = 0; k < cvx_mgrp_count(x_g, CVXDIM_SDP); k++) {
            cvx_size_t m = cvx_mgrp_elem(&xs, x_g, CVXDIM_SDP, k);
            cvx_mgrp_elem(&ys, y_g, CVXDIM_SDP, k);

            cvxm_view_diag(&x0, &xs, 0);
            cvxm_view_diag(&y0, &ys, 0);
            sdot += cvxm_dot(&x0, &y0);

            for (j = 1; j < m; j++) {
                // j'th subdiagonal
                cvxm_view_diag(&x0, &xs, -j);
                cvxm_view_diag(&y0, &ys, -j);
                sdot += 2.0 * cvxm_dot(&x0, &y0);
            }
        }
        break;
    case CVXDIM_CONELP:
        sdot += cvx_sdot_elem(x_g, y_g, CVXDIM_LINEAR);
        sdot += cvx_sdot_elem(x_g, y_g, CVXDIM_SOCP);
        sdot += cvx_sdot_elem(x_g, y_g, CVXDIM_SDP);
        break;
    case CVXDIM_CONVEXPROG:
        sdot += cvx_sdot_elem(x_g, y_g, CVXDIM_NLTARGET);
    case CVXDIM_CONVEXLP:
        sdot += cvx_sdot_elem(x_g, y_g, CVXDIM_NONLINEAR);
        sdot += cvx_sdot_elem(x_g, y_g, CVXDIM_LINEAR);
        sdot += cvx_sdot_elem(x_g, y_g, CVXDIM_SOCP);
        sdot += cvx_sdot_elem(x_g, y_g, CVXDIM_SDP);
        break;
    }
    return sdot;
}

cvx_float_t cvx_snrm2_elem(cvx_matgrp_t *x_g, cvx_dim_enum name)
{
    return sqrt(cvx_sdot_elem(x_g, x_g, name));
}


// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
