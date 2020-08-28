
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include "cvxc.h"
#include "cvxm.h"

// \brief Allocate work space for scaling computation, when 
void cvx_work_alloc(cvx_workspace_t *work, cvx_dims_t *dims)
{
    cvx_size_t maxs = cvx_dims_max(dims, S);
    // temp buffer for scaling computation in Q
    cvx_work_stemp(work, 3*maxs*maxs);
    // work space for SVD computation
    cvx_matrix_t Tmp;
    cvxm_map_data(&Tmp, maxs, maxs, (cvx_float_t *)0);
    cvx_size_t svdsz = armas_d_svd_work(&Tmp, CVX_WANTU, NULLCONF);
    cvxm_init(&work->work, svdsz, 1);
}

void cvx_work_free(cvx_workspace_t *work)
{
    if (!work)
        return;
    if (work->stmp)
        free(work->stmp);
    cvxm_free(&work->work);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
