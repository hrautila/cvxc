
// Copyright: Harri Rautila, 2018 <harri.rautila@gmail.com>

#include "convex.h"
#include "cvxm.h"


/**
 * @brief Allocate and initialize new indexing over spesified dimension set.
 *
 * @param[in] dims
 *    Dimension set to create indexing on.
 * @param[in] kind
 *    Indexing type for SDP-space sets. If dimension set does not include any 
 *    SDP sets then this parameter is not referenced.
 *
 * Parameter 'kind' spesifies indexing type for SDP ('S') space. 
 *  0  : standard storage indexing (m*m elements)
 *  1  : packed storage indexing (m*(m+1)/2 elements)
 *  2  : diagonal storage indexing (m elements)
 *  3  : only 'S' space diagonal storage indexing
 *
 * @return 
 *    New index set or null;
 */
cvx_index_t *cvx_index_new(cvx_dimset_t *dims, int kind)
{
    cvx_index_t *ind = (cvx_index_t *)malloc(sizeof(cvx_index_t));
    if (ind)
        return cvx_index_init(ind, dims, kind);
    return ind;           
}

cvx_size_t cvx_index_bytes(const cvx_dimset_t *dims, int kind)
{
    cvx_size_t n = 0;
    if (!dims)
        return n;

    if (kind != 3) {
        n = dims->slen +
            dims->qlen +
            (dims->ldim > 0 ? 1 : 0) +
            (dims->mnl > 0 ? 1 : 0) +
            1;
    } else {
        n = dims->slen + 1;
    }
    return n*sizeof(cvx_size_t);
}

/**
 * @brief Initialize new indexing over spesified dimension set.
 *
 * @param[in,out] ind
 *    On entry uninitialized indexing. On exit initialized indexing with
 *    proper memory allocation.
 * @param[in] dims
 *    Dimension set to create indexing on.
 * @param[in] kind
 *    Indexing type for SDP-space sets. 
 *
 * @return
 *    Pointer to initialized indexing or null if no space allocation done.
 *
 * Parameter 'kind' spesifies indexing type for SDP ('S') space. 
 *  0  : standard storage indexing (m*m elements)
 *  1  : packed storage indexing (m*(m+1)/2 elements)
 *  2  : diagonal storage indexing (m elements)
 *  3  : only 'S' space diagonal storage indexing
 *
 */
cvx_index_t *cvx_index_init(cvx_index_t *ind, const cvx_dimset_t *dims, int kind)
{
    int n;
    if (! ind)
        return ind;
    if (! dims)
        return (cvx_index_t *)0;

    // allocate for offsets
    if (kind != 3) {
        n = dims->slen +
            dims->qlen +
            (dims->ldim > 0 ? 1 : 0) +
            (dims->mnl > 0 ? 1 : 0) +
            1;
    } else {
        n = dims->slen + 1;
    }
    
    // allocate memory to hold both arrays
    ind->index = (cvx_size_t *)calloc(n, sizeof(cvx_size_t));
    if (! ind->index)
        return (cvx_index_t *)0;
    
    int k = 0;
    int off = 0;
    ind->indnl = ind->indl = ind->indq = ind->inds = (cvx_size_t *)0;
    ind->dims = dims;

    if (kind != 3) {
        if (dims->mnl > 0) {
            ind->indnl = &ind->index[k];
            ind->indnl[0] = off;
            k++;
            off += dims->mnl;
        }
        if (dims->ldim > 0) {
            ind->indl = &ind->index[k];
            ind->index[k] = off;
            k++;
            off += dims->ldim;
        }
        for (int j = 0; j < dims->qlen; j++) {
            if (j == 0)
                ind->indq = &ind->index[k];
            ind->index[k] = off;
            k++;
            off += dims->qdims[j];
        }
    }

    for (int j = 0; j < dims->slen; j++) {
        if (j == 0)
            ind->inds = &ind->index[k];
        ind->index[k] = off;
        k++;
        off += kind == 0 ?
            dims->sdims[j] * dims->sdims[j] :           // normal storage for S
            ( kind == 1 ?   
              dims->sdims[j] * (dims->sdims[j] + 1)/2 : // packed storage for S
              dims->sdims[j]);                          // diagonal storage for S (kind == 2|3)
    }
    // offset past the last entry;
    ind->index[k] = off;
    return ind;
}

void cvx_index_release(cvx_index_t *ind)
{
    if (! ind)
        return;
    if (ind->index)
        free(ind->index);
    ind->index = ind->indnl = ind->indl = ind->indq = ind->inds = (cvx_size_t *)0;
}

cvx_size_t cvx_index_count(const cvx_index_t *ind,
                           cvx_dim_enum name)
{
    switch (name) {
    case CVXDIM_SOCP:
        return ind->dims->qlen;
    case CVXDIM_SDP:
        return ind->dims->slen;
    case CVXDIM_LINEAR:
        return ind->indl ? 1 : 0;
    case CVXDIM_NONLINEAR:
    default:
        return ind->indnl ? 1 : 0;
    }
    return 0;
}

/*
 * @brief Make 'x' to present k'th element of y's dimension set.
 *
 * If x is null then returns dimensionality of k'th element in the spesified 
 * dimension set.
 */
cvx_size_t cvx_index_elem(cvx_matrix_t *x,
                          const cvx_matrix_t *y,
                          const cvx_index_t *ind,
                          cvx_dim_enum name,
                          int k)
{
    cvx_size_t n, m = 0;

    if (!y || !ind)
        return 0;

    if (x)
        cvxm_map_data(x, 0, 0, (cvx_float_t *)0);

    switch (name) {
    case CVXDIM_LINEAR:
        if (ind->indl) {
            m = ind->indl[1] - ind->indl[0];
            if (x)
                cvxm_map_data(x, m, 1, cvxm_data(y, ind->indl[0]));
        }
        break; 
    case CVXDIM_SOCP:
        if (ind->indq) {
            m = ind->indq[k+1] - ind->indq[k];
            if (x)
                cvxm_map_data(x, m, 1, cvxm_data(y, ind->indq[k]));
        }
        break; 
    case CVXDIM_SDP:
        if (ind->inds) {
            m = ind->dims->sdims[k];
            // TODO: what if indexing is not standard storage??
            if (x) {
                n = ind->inds[k+1] - ind->inds[k];
                n = n == m ? 1 : m;  // standard vs. diagonal storage
                cvxm_map_data(x, m, n, cvxm_data(y, ind->inds[k]));
            }
        }
        break;
    case CVXDIM_NONLINEAR:
        if (ind->indnl) {
            m = ind->indnl[1] - ind->indnl[0];
            if (x)
                cvxm_map_data(x, m, 1, cvxm_data(y, ind->indnl[0]));
        }
        break;
    }
    return m;
}


void cvx_index_create(cvx_matrix_t *x, cvx_index_t *index, const cvx_dimset_t *dims, cvx_index_type kind)
{
    cvx_size_t cdim = cvx_dimset_sum(dims, CVXDIM_NONLINEAR) +
        cvx_dimset_sum(dims, CVXDIM_LINEAR) +
        cvx_dimset_sum(dims, CVXDIM_SOCP);
    switch (kind) {
    case CVX_INDEX_PACKED:
        cdim += cvx_dimset_sum_packed(dims, CVXDIM_SDP);
        break;
    case CVX_INDEX_DIAG:
        cdim += cvx_dimset_sum(dims, CVXDIM_SDP);
        break;
    case CVX_INDEX_SIGS:
        cdim = cvx_dimset_sum(dims, CVXDIM_SDP);
        break;
    case CVX_INDEX_NORMAL:
    default:
        cdim += cvx_dimset_sum_squared(dims, CVXDIM_SDP);
        break;
    }
    cvxm_init(x, cdim, 1);
    cvx_index_init(index, dims, kind);
}


// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
