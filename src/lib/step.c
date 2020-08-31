
// Copyright: Harri Rautila, 2018 <harri.rautila@gmail.com>

#include "cvxc.h"

/**
 * @brief Compute \$ min {t | x + t*e >= 0} \$
 *
 * @param[in,out] x_g
 *     Argument vector group. If sigma_g is non-null then on exit SDP components 
 *     of x_g hold the eigenvectors of the component.
 * @param[out]    sigma_g
 *     If non-null then eigenvalues of SDP components is return
 * @param         wrk
 *    Workspace
 *
 * Computes min {t | x + t*e >= 0}, where e is defined as follows
 *
 *  - For the nonlinear and 'L' blocks: e is the vector of ones.
 *  - For the 'Q' blocks: e is the first unit vector.
 *  - For the 'S' blocks: e is the identity matrix.
 *
 */
cvxc_float_t cvxc_max_step(cvxc_matgrp_t *x_g,
                         cvxc_matgrp_t *sigma_g,
                         cvxc_memblk_t *wrk)
{
    cvxc_float_t v, tmax = -1e18;
    cvxc_matrix_t u, u1;
    cvxc_index_t *index = x_g->index;
    cvxc_size_t m;
    // static int svd = 0;

    if (index->indnlt) {
        m = cvxc_mgrp_elem(&u, x_g, CVXDIM_NLTARGET, 0);
        for (cvxc_size_t k = 0; k < m; k++) {
            if (tmax < -cvxm_get(&u, k, 0)) {
                tmax = -cvxm_get(&u, k, 0);
            }
        }
    }
    if (index->indnl) {
        m = cvxc_mgrp_elem(&u, x_g, CVXDIM_NONLINEAR, 0);
        for (cvxc_size_t k = 0; k < m; k++) {
            if (tmax < -cvxm_get(&u, k, 0)) {
                tmax = -cvxm_get(&u, k, 0);
            }
        }
    }
    if (index->indl) {
        m = cvxc_mgrp_elem(&u, x_g, CVXDIM_LINEAR, 0);
        for (cvxc_size_t k = 0; k < m; k++) {
            if (tmax < -cvxm_get(&u, k, 0)) {
                tmax = -cvxm_get(&u, k, 0);
            }
        }
    }
    if (index->indq) {
        for (int k = 0; k < cvxc_mgrp_count(x_g, CVXDIM_SOCP); k++) {
            m = cvxc_mgrp_elem(&u, x_g, CVXDIM_SOCP, k);
            cvxm_view_map(&u1, &u, 1, 0, m-1, 1);
            v = cvxm_nrm2(&u1);
            v -= cvxm_get(&u, 0, 0);
            if (tmax < v) {
                tmax = v;
            }
        }
    }
    if (index->inds) {
        cvxc_matrix_t Q, w, lk;
        cvxc_memblk_t mem;
        for (int k = 0; k < cvxc_mgrp_count(x_g, CVXDIM_SDP); k++) {
            m = cvxc_mgrp_elem(&u, x_g,  CVXDIM_SDP, k);
#if 0
            if (svd > 3 && svd < 7) {
                cvxc_mat_printf(stderr, "%.12e", &u, "svd.u");
            }
#endif
            if (sigma_g) {
                cvxc_mgrp_elem(&lk, sigma_g, CVXDIM_SDP, k);
                cvxm_evd_sym(&lk, &u, CVXC_WANTV|CVXC_LOWER|ARMAS_FORWARD, wrk);
                v = - cvxm_get(&lk, 0, 0);
#if 0
                if (svd > 3 && svd < 7) {
                    cvxc_mat_printf(stderr, "%.12e", &lk, "svd.lk");
                }
#endif
            } else {
                cvxm_map_data(&Q, m, m, __mblk_offset(wrk, 0));
                cvxm_map_data(&w, m, 1, __mblk_offset(wrk, m*m));
                __mblk_subblk(&mem, wrk, m*m +m);
                cvxm_copy(&Q, &u, 0);
                cvxm_evd_sym(&w, &Q, CVXC_LOWER, &mem);

                v = - cvxm_get(&w, 0, 0);
#if 0
                if (svd > 3 && svd < 7) {
                    cvxc_mat_printf(stderr, "%.12e", &w, "svd.w");
                }
#endif
            }
            if (tmax < v)
                tmax = v;
        }
        //svd++;
    }
    return tmax;
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
