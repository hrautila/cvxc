
// Copyright: Harri Rautila, 2018 <harri.rautila@gmail.com>

#include "convex.h"
#include "cvxm.h"


cvx_dimset_t *cvx_dimset_alloc(cvx_dimset_t *dims, int linear, const int *socp, const int *sdp)
{
    int n;
    if (!dims)
        return dims;
    for (n = 0; socp && socp[n] > 0; n++);
    dims->qlen = n;
    
    for (n = 0; sdp && sdp[n] > 0; n++);
    dims->slen = n;

    // allocate memory to hold both arrays
    dims->qdims = (int *)calloc(dims->slen + dims->qlen, sizeof(int));
    if (! dims->qdims) {
        dims->qlen = dims->slen = 0;
        return (cvx_dimset_t *)0;
    }
    // SDP dimension after SOCP dimension
    dims->sdims = dims->slen > 0 ? &dims->qdims[dims->qlen] : (int *)0;
    dims->mnl = 0;
    dims->ldim = linear;
    // copy SOCP dimensions
    for (n = 0; n < dims->qlen; dims->qdims[n] = socp[n], n++);
    // copy SDP dimensions
    for (n = 0; n < dims->slen; dims->sdims[n] = sdp[n], n++);
    return dims;
}

void cvx_dimset_release(cvx_dimset_t *dims)
{
    if (!dims)
        return;

    if (dims->qdims)
        free(dims->qdims);
    dims->qlen = dims->slen = 0;
    dims->qdims = (int *)0;
    dims->sdims = (int *)0;
    dims->ldim = 0;
}

cvx_size_t cvx_dimset_max(const cvx_dimset_t *dims, cvx_dim_enum name)
{
    if (!dims)
        return 0;
    
    cvx_size_t m = 0;

    switch (name) {
    case CVXDIM_NONLINEAR:
        return dims->mnl;
    case CVXDIM_LINEAR:
        return dims->ldim;
    case CVXDIM_SOCP:
        for (int k = 0; k < dims->qlen; k++) {
            if (dims->qdims[k] > m)
                m = dims->qdims[k];
        }
        break;
    case CVXDIM_SDP:
        for (int k = 0; k < dims->slen; k++) {
            if (dims->sdims[k] > m)
                m = dims->sdims[k];
        }
        break;
    default:
        break;
    }
    return m;
}

cvx_size_t cvx_dimset_sum(const cvx_dimset_t *dims, cvx_dim_enum name)
{
    if (!dims)
        return 0;
    
    cvx_size_t sum = 0;

    switch (name) {
    case CVXDIM_NONLINEAR:
        return dims->mnl;
    case CVXDIM_LINEAR:
        return dims->ldim;
    case CVXDIM_SOCP:
        for (int k = 0; k < dims->qlen; k++)
            sum += dims->qdims[k];
        break;
    case CVXDIM_SDP:
        for (int k = 0; k < dims->slen; k++)
            sum += dims->sdims[k];
        break;
    case CVXDIM_CONELP:
        sum = dims->ldim;
        for (int k = 0; k < dims->qlen; k++)
            sum += dims->qdims[k];
        for (int k = 0; k < dims->slen; k++)
            sum += dims->sdims[k];
        break;
    default:
        break;
    }
    return sum;
}

cvx_size_t cvx_dimset_sum_squared(const cvx_dimset_t *dims, cvx_dim_enum name)
{
    if (!dims)
        return 0;
    
    cvx_size_t sum = 0;

    switch (name) {
    case CVXDIM_NONLINEAR:
        return dims->mnl;
    case CVXDIM_LINEAR:
        return dims->ldim;
    case CVXDIM_SOCP:
        for (int k = 0; k < dims->qlen; k++)
            sum += dims->qdims[k];
        break;
    case CVXDIM_SDP:
        for (int k = 0; k < dims->slen; k++)
            sum += dims->sdims[k] * dims->sdims[k];
        break;
    case CVXDIM_CONELP:
        sum = dims->ldim;
        for (int k = 0; k < dims->qlen; k++)
            sum += dims->qdims[k];
        for (int k = 0; k < dims->slen; k++)
            sum += dims->sdims[k] * dims->sdims[k];
        break;
    default:
        break;
    }

    return sum;
}

cvx_size_t cvx_dimset_sum_packed(const cvx_dimset_t *dims, cvx_dim_enum name)
{
    if (!dims)
        return 0;
    
    cvx_size_t sum = 0;

    switch (name) {
    case CVXDIM_NONLINEAR:
    case CVXDIM_LINEAR:
    case CVXDIM_SOCP:
        return 0;
    case CVXDIM_SDP:
        for (int k = 0; k < dims->slen; k++)
            sum += dims->sdims[k] * (dims->sdims[k] + 1)/2;
        break;
    case CVXDIM_CONELP:
        sum = dims->ldim;
        for (int k = 0; k < dims->qlen; k++)
            sum += dims->qdims[k];
        for (int k = 0; k < dims->slen; k++)
            sum += dims->sdims[k] * (dims->sdims[k] + 1)/2;
        break;
    default:
        break;
    }
    return sum;
}


// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
