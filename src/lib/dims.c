
// Copyright: Harri Rautila, 2018 <harri.rautila@gmail.com>

#include "cvxc.h"
#include "cvxm.h"


cvx_dimset_t *cvx_dimset_alloc(cvx_dimset_t *dims, cvx_size_t linear, const cvx_size_t *socp, const cvx_size_t *sdp)
{
    cvx_size_t n;
    if (!dims)
        return dims;
    for (n = 0; socp && socp[n] > 0; n++);
    dims->qlen = n;

    for (n = 0; sdp && sdp[n] > 0; n++);
    dims->slen = n;

    // allocate memory to hold both arrays
    dims->qdims = (cvx_size_t *)calloc(dims->slen + dims->qlen, sizeof(cvx_size_t));
    if (! dims->qdims) {
        dims->qlen = dims->slen = 0;
        return (cvx_dimset_t *)0;
    }
    // SDP dimension after SOCP dimension
    dims->sdims = dims->slen > 0 ? &dims->qdims[dims->qlen] : (cvx_size_t *)0;
    dims->mnl = dims->iscpt = 0;
    dims->ldim = linear;
    // copy SOCP dimensions
    for (n = 0; n < dims->qlen; dims->qdims[n] = socp[n], n++);
    // copy SDP dimensions
    for (n = 0; n < dims->slen; dims->sdims[n] = sdp[n], n++);
    return dims;
}

cvx_dimset_t *cvx_dimset_create(cvx_dimset_t *dims, cvx_size_t nonlinear, cvx_size_t linear, cvx_size_t socp, cvx_size_t sdp)
{
    dims->iscpt = 0;
    dims->mnl = nonlinear;
    dims->ldim = linear;
    dims->qlen = socp;
    dims->slen = sdp;

    dims->qdims = (cvx_size_t *)0;
    dims->sdims = (cvx_size_t *)0;

    if (socp + sdp > 0) {
        // allocate memory to hold both arrays
        dims->qdims = (cvx_size_t *)calloc(dims->slen + dims->qlen, sizeof(cvx_size_t));
        if (! dims->qdims) {
            dims->qlen = dims->slen = 0;
            return (cvx_dimset_t *)0;
        }
        // SDP dimension after SOCP dimension
        dims->sdims = dims->slen > 0 ? &dims->qdims[dims->qlen] : (cvx_size_t *)0;
    }
    return dims;
}

cvx_dimset_t *cvx_dimset_new(cvx_size_t nonlinear, cvx_size_t linear, cvx_size_t socp, cvx_size_t sdp)
{
    cvx_dimset_t *dims = (cvx_dimset_t *)malloc(sizeof(cvx_dimset_t));
    if (!dims)
        return (cvx_dimset_t *)0;

    if (! cvx_dimset_create(dims, nonlinear, linear, socp, sdp)) {
        free(dims);
        return (cvx_dimset_t *)0;
    }
    return dims;
}

void cvx_dimset_release(cvx_dimset_t *dims)
{
    if (!dims)
        return;

    if (dims->qdims)
        free(dims->qdims);
    dims->qlen = dims->slen = dims->iscpt = 0;
    dims->qdims = (cvx_size_t *)0;
    dims->sdims = (cvx_size_t *)0;
    dims->ldim = 0;
}

void cvx_dimset_free(cvx_dimset_t *dims)
{
    if (!dims)
        return;
    cvx_dimset_release(dims);
    free(dims);
}

cvx_size_t cvx_dimset_count(const cvx_dimset_t *dims, cvx_dim_enum name)
{
    if (!dims)
        return 0;

    cvx_size_t m = 1;

    switch (name) {
    case CVXDIM_SOCP:
        return dims->qlen;
    case CVXDIM_SDP:
        return dims->slen;
    default:
        break;
    }
    return m;
}

cvx_size_t cvx_dimset_max(const cvx_dimset_t *dims, cvx_dim_enum name)
{
    if (!dims)
        return 0;

    cvx_size_t m = 0;

    switch (name) {
    case CVXDIM_NLTARGET:
        return dims->iscpt ? 1 : 0;
    case CVXDIM_NONLINEAR:
        return dims->mnl;
    case CVXDIM_LINEAR:
        return dims->ldim;
    case CVXDIM_SOCP:
        for (cvx_size_t k = 0; k < dims->qlen; k++) {
            if (dims->qdims[k] > m)
                m = dims->qdims[k];
        }
        break;
    case CVXDIM_SDP:
        for (cvx_size_t k = 0; k < dims->slen; k++) {
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
    case CVXDIM_NLTARGET:
        return dims->iscpt ? 1 : 0;
    case CVXDIM_NONLINEAR:
        return dims->mnl;
    case CVXDIM_LINEAR:
        return dims->ldim;
    case CVXDIM_SOCP:
        for (cvx_size_t k = 0; k < dims->qlen; k++)
            sum += dims->qdims[k];
        break;
    case CVXDIM_SDP:
        for (cvx_size_t k = 0; k < dims->slen; k++)
            sum += dims->sdims[k];
        break;
    case CVXDIM_CONVEXPROG:
        sum += dims->iscpt ? 1 : 0;
    case CVXDIM_CONVEXLP:
        sum += dims->mnl;
    case CVXDIM_CONELP:
        sum += dims->ldim;
        for (cvx_size_t k = 0; k < dims->qlen; k++)
            sum += dims->qdims[k];
        for (cvx_size_t k = 0; k < dims->slen; k++)
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
    case CVXDIM_NLTARGET:
        return dims->iscpt ? 1 : 0;
    case CVXDIM_NONLINEAR:
        return dims->mnl;
    case CVXDIM_LINEAR:
        return dims->ldim;
    case CVXDIM_SOCP:
        for (cvx_size_t k = 0; k < dims->qlen; k++)
            sum += dims->qdims[k];
        break;
    case CVXDIM_SDP:
        for (cvx_size_t k = 0; k < dims->slen; k++)
            sum += dims->sdims[k] * dims->sdims[k];
        break;
    case CVXDIM_CONVEXPROG:
        sum += dims->iscpt ? 1 : 0;
    case CVXDIM_CONVEXLP:
        sum += dims->mnl;
    case CVXDIM_CONELP:
        sum += dims->ldim;
        for (cvx_size_t k = 0; k < dims->qlen; k++)
            sum += dims->qdims[k];
        for (cvx_size_t k = 0; k < dims->slen; k++)
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
    case CVXDIM_NLTARGET:
    case CVXDIM_NONLINEAR:
    case CVXDIM_LINEAR:
    case CVXDIM_SOCP:
        return cvx_dimset_sum(dims, name);
    case CVXDIM_SDP:
        for (cvx_size_t k = 0; k < dims->slen; k++)
            sum += dims->sdims[k] * (dims->sdims[k] + 1)/2;
        break;
    case CVXDIM_CONVEXPROG:
        sum += dims->iscpt ? 1 : 0;
    case CVXDIM_CONVEXLP:
        sum += dims->mnl;
    case CVXDIM_CONELP:
        sum += dims->ldim;
        for (cvx_size_t k = 0; k < dims->qlen; k++)
            sum += dims->qdims[k];
        for (cvx_size_t k = 0; k < dims->slen; k++)
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
