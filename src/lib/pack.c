/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */
#include <stdio.h>
#include "internal.h"

/**
 *   @brief y = pack(x), copy x to y using packed storage.
 *
 *   The vector x is an element of S, with the 's' components stored in 
 *   unpacked storage.  On return, x is copied to y with the 's' components
 *   stored in packed storage and the off-diagonal entries scaled by 
 *   sqrt(2).
 *
 *   @param y
 *      On exit, copy of x with S component in packed storage; in vector format
 *   @param x
 *      Source data; in vector format
 *   @param dims
 *      Problems dimensions 
 *   @param offsetx 
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
void cvxc_pack(cvxc_matrix_t *y, cvxc_matrix_t *x, const cvxc_index_t *index)
{
    cvxc_size_t nlq, ip, n, np;
    cvxc_size_t j, k;
    cvxc_matrix_t x1, y1, xk;

    nlq = cvxc_index_length(index, CVXDIM_NONLINEAR) +
        cvxc_index_length(index, CVXDIM_LINEAR) +
        cvxc_index_length(index, CVXDIM_SOCP);

    if (nlq > 0) {
        cvxm_view_map(&x1, x, 0, 0, nlq, 1);
        cvxm_view_map(&y1, y, 0, 0, nlq, 1);
        cvxm_copy(&y1, &x1, CVXC_ALL);
    }
    ip = nlq;

    for (j = 0; j < cvxc_index_count(index, CVXDIM_SDP); j++) {
        n = cvxc_index_elem(&xk, x, index, CVXDIM_SDP, j);
        for (k = 0; k < n; k++) {
            // copy column; entries on diagonal and below
            cvxm_view_map(&x1, &xk, k*(n+1), 0, n-k, 1);
            cvxm_view_map(&y1, y, ip, 0, n-k, 1);

            cvxm_copy(&y1, &x1, CVXC_ALL);
            // divide diagonal entry by sqrt(2)
            cvxm_set(&y1, 0, 0, cvxm_get(&y1, 0, 0)/SQRT2);
            ip += n - k;
        }
    }
    np = ip - nlq;
    // scale all of 's' with sqrt(2) ('unscaling' the diagonal entries)
    cvxm_view_map(&y1, y, nlq, 0, np, 1);
    cvxm_scale(&y1, SQRT2, CVXC_ALL);
}

/**
 * @brief y = unpack(x); inverse of cvxc_pack operation.
 *
 *  The vector x is an element of S, with the 's' components stored
 *  in packed storage and off-diagonal entries scaled by sqrt(2).
 *  On return, x is copied to y with the 's' components stored in 
 *  unpacked storage and off-diagonal entries scaled by 1/sqrt(2).
 */
void cvxc_unpack(cvxc_matrix_t *y, cvxc_matrix_t *x, const cvxc_index_t *index)
{
    cvxc_size_t nlq, ip, n;
    cvxc_size_t j, k;
    cvxc_matrix_t x1, yk, y1;

    nlq = cvxc_index_length(index, CVXDIM_NONLINEAR) +
        cvxc_index_length(index, CVXDIM_LINEAR) +
        cvxc_index_length(index, CVXDIM_SOCP);

    if (nlq > 0) {
        cvxm_view_map(&x1, x, 0, 0, nlq, 1);
        cvxm_view_map(&y1, y, 0, 0, nlq, 1);
        cvxm_copy(&y1, &x1, CVXC_ALL);
    }

    ip = nlq;
    //iu = nlq;

    for (j = 0; j < cvxc_index_count(index, CVXDIM_SDP); j++) {
        n = cvxc_index_elem(&yk, y, index, CVXDIM_SDP, j);
        for (k = 0; k < n; k++) {
            // create submatrix view
            cvxm_view_map(&x1, x, ip, 0, n-k, 1);
            cvxm_view_map(&y1, &yk, k, k, n-k, 1);

            cvxm_copy(&y1, &x1, CVXC_ALL);
            //cvxm_set(&y1, 0, 0, cvxm_get(&y1, 0, 0)/SQRT2);
            ip += n - k;

            // scale column off-diagonal entries
            cvxm_view_map(&y1, &yk, k+1, k, n-k-1, 1);
            cvxm_scale(&y1, 1.0/SQRT2, CVXC_ALL);
        }
        //iu += n*n;
    }
}

/*
 * \brief x = pack(x)
 *
 * In-place version of pack(), which also accepts matrix arguments x.
 * The columns of x are elements of S, with the 's' components stored
 * in unpacked storage.  On return, the 's' components are stored in
 * packed storage and the off-diagonal entries are scaled by sqrt(2).
 */
void cvxc_pack2(cvxc_matrix_t *x, const cvxc_index_t *index, cvxc_memblk_t *work)
{
    cvxc_size_t iu, ip, n, nr, nc;
    cvxc_size_t i, j, k;
    cvxc_matrix_t xk, row;

    iu = ip =
        cvxc_index_length(index, CVXDIM_NONLINEAR) +
        cvxc_index_length(index, CVXDIM_LINEAR) +
        cvxc_index_length(index, CVXDIM_SOCP);

    cvxm_size(&nr, &nc, x);
    cvxm_map_data(&row, nc, 1, cvxc_mblk_offset(work, 0));

    for (j = 0; j < cvxc_index_count(index, CVXDIM_SDP); j++) {
        // get row count for j'th element
        n = cvxc_index_elem(__cvxnil, x, index, CVXDIM_SDP, j);

        for (k = 0; k < n; k++) {
            // k'th diagonal element on x
            cvxm_view_map(&xk, x, iu + k*(n+1), 0, 1, nc);
            cvxm_view_map(&row, x, ip, 0, 1, nc);
            cvxm_copy(&row, &xk, CVXC_ALL);
            // off diagonal entries of k'th column
            for (i = 0; i < n-k; i++) {
                cvxm_view_map(&xk, x, iu + k*(n+1)+i, 0, 1, nc);
                cvxm_view_map(&row, x, ip+i, 0, 1, nc);
                cvxm_copy(&row, &xk, CVXC_ALL);
                cvxm_scale(&row, SQRT2, CVXC_ALL);
            }
            ip += n - k;
        }
        iu += n*n;
    }
}
