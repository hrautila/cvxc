
// Copyright: Harri Rautila, 2018 <harri.rautila@gmail.com>

#include "convex.h"
#include "cvxm.h"

/*
   Index layout:

       0       : 0
       1       : length of non-linear
       2       : ind[1] + length-of-linear

 */

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
cvx_index_t *cvx_index_new(const cvx_dimset_t *dims, int kind)
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
            (dims->iscpt > 0 ? 1 : 0) +
            1;
    } else {
        n = dims->slen + 1;
    }
    n *= sizeof(cvx_size_t);
    // align to 64bits
    n += (n & 0x7) != 0 ? 8 - (n & 0x7) : 0;
    return n;
}

/**
 * @brief Overlay index structure to provided memory block.
 */
cvx_size_t cvx_index_make(cvx_index_t *ind,
                          const cvx_dimset_t *dims,
                          int kind,
                          void *buf,
                          cvx_size_t nbytes)
{
    if (!ind || !dims)
        return 0;

    cvx_size_t n = cvx_index_bytes(dims, kind);
    if (nbytes < n)
        return 0;

    cvx_size_t k = 0;
    cvx_size_t off = 0;
    ind->indnlt = ind->indnl = ind->indl = ind->indq = ind->inds = (cvx_size_t *)0;
    ind->dims = dims;
    ind->index = (cvx_size_t *)buf;

    if (kind != CVX_INDEX_SIGS) {
        if (dims->iscpt > 0) {
            ind->indnlt = &ind->index[k];
            k++;
            off += 1;
        }
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
        for (cvx_size_t j = 0; j < dims->qlen; j++) {
            if (j == 0)
                ind->indq = &ind->index[k];
            ind->index[k] = off;
            k++;
            off += dims->qdims[j];
        }
    }

    for (cvx_size_t j = 0; j < dims->slen; j++) {
        if (j == 0)
            ind->inds = &ind->index[k];
        ind->index[k] = off;
        k++;
        off += kind == CVX_INDEX_NORMAL ?
            dims->sdims[j] * dims->sdims[j] :           // normal storage for S
            ( kind == CVX_INDEX_PACKED ?
              dims->sdims[j] * (dims->sdims[j] + 1)/2 : // packed storage for S
              dims->sdims[j]);                          // diagonal storage for S (kind == 2|3)
    }

    // offset past the last entry;
    ind->index[k] = off;
    ind->indlen = k;
    return n;
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
    if (! ind)
        return ind;
    if (! dims)
        return (cvx_index_t *)0;

    cvx_size_t nb = cvx_index_bytes(dims, kind);
    void *mem = calloc(nb, 1);
    if (!mem)
        return (cvx_index_t *)0;

    cvx_index_make(ind, dims, kind, mem, nb);
    ind->__bytes = mem;
    return ind;

}

void cvx_index_release(cvx_index_t *ind)
{
    if (! ind)
        return;
    if (ind->__bytes) {
        free(ind->__bytes);
        ind->__bytes = (void *)0;
    }
    ind->index = ind->indnlt = ind->indnl = ind->indl = ind->indq = ind->inds = (cvx_size_t *)0;
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
        return ind->indnl ? 1 : 0;
    case CVXDIM_NLTARGET:
        return ind->indnlt ? 1 : 0;
    default:
        break;
    }
    return 0;
}

/**
 * @brief      Create subindex with given parts from source index.
 *
 * @param      ind
 *     Result index
 * @param      src
 *     Source index
 * @param      parts
 *     Flags from cvx_dim_enum set to indicate which part are included.
 *
 */
void cvx_subindex(cvx_index_t *ind, const cvx_index_t *src, int parts)
{
    ind->__bytes = (void *)0;
    ind->index = ind->indnlt = ind->indnl = ind->indl = ind->indq = ind->inds = (cvx_size_t *)0;
    if ((parts & CVXDIM_NLTARGET) != 0) {
        ind->indnlt = src->indnlt;
    }
    if ((parts & CVXDIM_NONLINEAR) != 0) {
        ind->indnl = src->indnl;
    }
    if ((parts & CVXDIM_LINEAR) != 0) {
        ind->indl = src->indl;
    }
    if ((parts & CVXDIM_SOCP) != 0) {
        ind->indq = src->indq;
    }
    if ((parts & CVXDIM_SDP) != 0) {
        ind->inds = src->inds;
    }
    ind->index = src->index;
    ind->indlen = src->indlen;
}

/**
 * @brief Make 'x' to present k'th element of y's dimension set.
 *
 * @param[out] x
 *     If not null; on entry uninitialized matrix, on exit k'th element 
 *     of requested set.
 * @param[in]  y
 *     Source vector
 * @param[in] ind
 *     Index set associated with y
 * @param[in] name
 *     Index set name, CVXDIM_NONLINEAR (NL), CVXDIM_LINEAR (L),  CVXDIM_SOCP (Q)
 *     CVXDIM_SDP (S) or CVXDIM_CONELP.
 * @param[in] k
 *     Requested element. Non-zero value meaningfull only for SOCP and SDP sets.
 *
 * @return Length of the requested element;
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
    cvx_size_t n = 0, m = 0;

    if (!y || !ind)
        return 0;

    if (x)
        cvxm_map_data(x, 0, 0, (cvx_float_t *)0);

    switch (name) {
    case CVXDIM_CONVEX:
        if (ind->indnlt && ind->indnl) {
            m = ind->indnl[1] - ind->indnlt[0];
            n = ind->indnlt[0];
        } else if (ind->indnlt) {
            m = ind->indnlt[1] - ind->indnlt[0];
            n = ind->indnlt[0];
        } else if (ind->indnl) {
            m = ind->indnl[1] - ind->indnl[0];
            n = ind->indnl[0];
        }
        if (x)
            cvxm_map_data(x, m, 1, cvxm_data(y, n));
        break;
    case CVXDIM_NLTARGET:
        if (ind->indnlt) {
            m = ind->indnlt[1] - ind->indnlt[0];
            if (x)
                cvxm_map_data(x, m, 1, cvxm_data(y, ind->indnlt[0]));
        }
        break;
    case CVXDIM_NONLINEAR:
        if (ind->indnl) {
            m = ind->indnl[1] - ind->indnl[0];
            if (x)
                cvxm_map_data(x, m, 1, cvxm_data(y, ind->indnl[0]));
        }
        break;
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
    case CVXDIM_CONELP:
        // map linear,socp and sdp parts onto x
        if (ind->indl)
            n = ind->indl[0];  // have linear part
        else if (ind->indq)
            n = ind->indq[0];  // have socp part
        else if (ind->inds)
            n = ind->inds[0];  // last resort; must have sdp part
        m = ind->index[ind->indlen] - n;
        if (x)
            cvxm_map_data(x, m, 1, cvxm_data(y, n));
        break;
    case CVXDIM_CONVEXLP:
        // map all but non-linear target function
        if (ind->indnl)
            n = ind->indnl[0];  // have linear part
        else if (ind->indl)
            n = ind->indl[0];  // have linear part
        else if (ind->indq)
            n = ind->indq[0];  // have socp part
        else if (ind->inds)
            n = ind->inds[0];  // last resort; must have sdp part
        m = ind->index[ind->indlen] - n;
        if (x)
            cvxm_map_data(x, m, 1, cvxm_data(y, (ind->indnlt ? 1 : 0)));
        break;
    case CVXDIM_CONVEXPROG:
        // map all parts
        m = ind->index[ind->indlen];
        if (x)
            cvxm_map_data(x, m, 1, cvxm_data(y, 0));
        break;
    default:
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
