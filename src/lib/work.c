
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include "cvxc.h"
#include "cvxm.h"

// \brief Allocate work space for scaling computation, when 
void cvxc_work_alloc(cvxc_workspace_t *work, cvxc_dims_t *dims)
{
    cvxc_size_t maxs = cvxc_dims_max(dims, S);
    // temp buffer for scaling computation in Q
    cvxc_work_stemp(work, 3*maxs*maxs);
    // work space for SVD computation
    cvxc_matrix_t Tmp;
    armas_wbuf_t wb;
    wb->bytes = 0;
    cvxm_map_data(&Tmp, maxs, maxs, (cvxc_float_t *)0);
    if (armas_svd_w(&Tmp, CVXC_WANTU, &wb, NULLCONF) < 0) {
        // error;
    }
    cvxm_init(&work->work, wb->bytes, 1);
}

void cvxc_work_free(cvxc_workspace_t *work)
{
    if (!work)
        return;
    if (work->stmp)
        free(work->stmp);
    cvxm_free(&work->work);
}
