
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
    cvxm_map_data(&Tmp, maxs, maxs, (cvxc_float_t *)0);
    cvxc_size_t svdsz = armas_d_svd_work(&Tmp, CVXC_WANTU, NULLCONF);
    cvxm_init(&work->work, svdsz, 1);
}

void cvxc_work_free(cvxc_workspace_t *work)
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
