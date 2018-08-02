
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

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
 *   \brief y = pack(x), copy x to y using packed storage.
 *
 *   The vector x is an element of S, with the 's' components stored in 
 *   unpacked storage.  On return, x is copied to y with the 's' components
 *   stored in packed storage and the off-diagonal entries scaled by 
 *   sqrt(2).
 *
 *   \param y
 *      On exit, copy of x with S component in packed storage; in vector format
 *   \param x
 *      Source data; in vector format
 *   \param dims
 *      Problems dimensions 
 *   \param offsetx 
 *      Offset to X vector
 *
 *   Standard storage to packed storage
 *
 *    a00 
 *    a10 a11     -->  a00 a10 a20 a11 a21 a22 
 *    a20 a21 a22
 *
 *    N*N elems        N*(N+1)/2 elems
 */
int cvx_pack(cvx_matrix_t *y, cvx_matrix_t *x, cvx_index_t *index)
{
    cvx_size_t nlq, ip, n, np;
    int j, k;
    cvx_matrix_t x1, y1, xk;
    
    nlq = cvx_dimset_sum(index->dims, CVXDIM_NONLINEAR) +
        cvx_dimset_sum(index->dims, CVXDIM_LINEAR) +
        cvx_dimset_sum(index->dims, CVXDIM_SOCP);

    if (nlq > 0) {
        cvxm_view_map(&x1, x, 0, 0, nlq, 1);
        cvxm_view_map(&y1, y, 0, 0, nlq, 1);
        cvxm_copy(&y1, &x1, CVX_ALL);
    }

    ip = nlq;
    
    for (j = 0; j < cvx_index_count(index, CVXDIM_SDP); j++) {
        n = cvx_index_elem(&xk, x, index, CVXDIM_SDP, j);
        for (k = 0; k < n; k++) {
            // copy column; entries on diagonal and below
            cvxm_view_map(&x1, &xk, k*(n+1), 0, n-k, 1);
            cvxm_view_map(&y1, y, ip, 0, n-k, 1);

            cvxm_copy(&y1, &x1, CVX_ALL);
            // divide diagonal entry by sqrt(2)
            cvxm_set(&y1, 0, 0, cvxm_get(&y1, 0, 0)/SQRT2);
            ip += n - k;
        }
    }
    np = ip - nlq;
    // scale all of 's' with sqrt(2) ('unscaling' the diagonal entries)
    cvxm_view_map(&y1, y, nlq, 0, np, 1);
    cvxm_scale(&y1, SQRT2, CVX_ALL);
    return 0;
}

/*
 * \brief y = unpack(x); inverse of cvx_pack operation.
 *
 *  The vector x is an element of S, with the 's' components stored
 *  in packed storage and off-diagonal entries scaled by sqrt(2).
 *  On return, x is copied to y with the 's' components stored in 
 *  unpacked storage and off-diagonal entries scaled by 1/sqrt(2).
 */
int cvx_unpack(cvx_matrix_t *y, cvx_matrix_t *x, cvx_index_t *index)
{
    cvx_size_t nlq, ip, n;
    int j, k;
    cvx_matrix_t x1, yk, y1;
    
    nlq = cvx_dimset_sum(index->dims, CVXDIM_NONLINEAR) +
        cvx_dimset_sum(index->dims, CVXDIM_LINEAR) +
        cvx_dimset_sum(index->dims, CVXDIM_SOCP);

    if (nlq > 0) {
        cvxm_view_map(&x1, x, 0, 0, nlq, 1);
        cvxm_view_map(&y1, y, 0, 0, nlq, 1);
        cvxm_copy(&y1, &x1, CVX_ALL);
    }

    ip = nlq;
    //iu = nlq;

    for (j = 0; j < cvx_index_count(index, CVXDIM_SDP); j++) {
        n = cvx_index_elem(&yk, y, index, CVXDIM_SDP, j);
        for (k = 0; k < n; k++) {
            // create submatrix view
            cvxm_view_map(&x1, x, ip, 0, n-k, 1);
            cvxm_view_map(&y1, &yk, k, k, n-k, 1);

            cvxm_copy(&y1, &x1, CVX_ALL);
            //cvxm_set(&y1, 0, 0, cvxm_get(&y1, 0, 0)/SQRT2);
            ip += n - k;

            // scale column off-diagonal entries
            cvxm_view_map(&y1, &yk, k+1, k, n-k-1, 1);
            cvxm_scale(&y1, 1.0/SQRT2, CVX_ALL);
        }
        //iu += n*n;
    }
    return 0;
}

/*
 * \brief x = pack(x)
 *
 * In-place version of pack(), which also accepts matrix arguments x.
 * The columns of x are elements of S, with the 's' components stored
 * in unpacked storage.  On return, the 's' components are stored in
 * packed storage and the off-diagonal entries are scaled by sqrt(2).
 */
int cvx_pack2(cvx_matrix_t *x, cvx_index_t *index, cvx_memblk_t *work)
{
    cvx_size_t iu, ip, n, nr, nc;
    int i, j, k;
    cvx_matrix_t xk, row;
    
    iu = ip = cvx_dimset_sum(index->dims, CVXDIM_NONLINEAR) +
        cvx_dimset_sum(index->dims, CVXDIM_LINEAR) +
        cvx_dimset_sum(index->dims, CVXDIM_SOCP);

    cvxm_size(&nr, &nc, x);
    cvxm_map_data(&row, nc, 1, __mblk_offset(work, 0));

    for (j = 0; j < cvx_index_count(index, CVXDIM_SDP); j++) {
        // get row count for j'th element
        n = cvx_index_elem(__cvxnil, x, index, CVXDIM_SDP, j);

        for (k = 0; k < n; k++) {
            // k'th diagonal element on x
            cvxm_view_map(&xk, x, iu + k*(n+1), 0, 1, nc);
            cvxm_view_map(&row, x, ip, 0, 1, nc);
            cvxm_copy(&row, &xk, CVX_ALL);
            // off diagonal entries of k'th column
            for (i = 0; i < n-k; i++) {
                cvxm_view_map(&xk, x, iu + k*(n+1)+i, 0, 1, nc);
                cvxm_view_map(&row, x, ip+i, 0, 1, nc);
                cvxm_copy(&row, &xk, CVX_ALL);
                cvxm_scale(&row, SQRT2, CVX_ALL);
            }
            ip += n - k;
        }
        iu += n*n;
    }

    return 0;
}


// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
