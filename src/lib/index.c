/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */
#include "cvxc.h"

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
cvxc_index_t *cvxc_index_new(const cvxc_dimset_t *dims, int kind)
{
    cvxc_index_t *ind = (cvxc_index_t *)malloc(sizeof(cvxc_index_t));
    if (ind)
        return cvxc_index_init(ind, dims, kind);
    return ind;
}

cvxc_size_t cvxc_index_bytes(const cvxc_dimset_t *dims, int kind)
{
    cvxc_size_t n = 0;
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
    n *= sizeof(cvxc_size_t);
    // align to 64bits
    n += (n & 0x7) != 0 ? 8 - (n & 0x7) : 0;
    return n;
}

/**
 * @brief Overlay index structure to provided memory block.
 */
cvxc_size_t cvxc_index_make(cvxc_index_t *ind,
                          const cvxc_dimset_t *dims,
                          int kind,
                          void *buf,
                          cvxc_size_t nbytes)
{
    if (!ind || !dims)
        return 0;

    cvxc_size_t n = cvxc_index_bytes(dims, kind);
    if (nbytes < n)
        return 0;

    cvxc_size_t k = 0;
    cvxc_size_t off = 0;
    ind->indnlt = ind->indnl = ind->indl = ind->indq = ind->inds = (cvxc_size_t *)0;
    //ind->dims = dims;
    ind->type = kind;
    ind->index = (cvxc_size_t *)buf;

    if (kind != CVXC_INDEX_SIGS) {
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
        for (cvxc_size_t j = 0; j < dims->qlen; j++) {
            if (j == 0)
                ind->indq = &ind->index[k];
            ind->index[k] = off;
            k++;
            off += dims->qdims[j];
        }
        ind->qlen = dims->qlen;
    }

    for (cvxc_size_t j = 0; j < dims->slen; j++) {
        if (j == 0)
            ind->inds = &ind->index[k];
        ind->index[k] = off;
        k++;
        off += kind == CVXC_INDEX_NORMAL ?
            dims->sdims[j] * dims->sdims[j] :           // normal storage for S
            ( kind == CVXC_INDEX_PACKED ?
              dims->sdims[j] * (dims->sdims[j] + 1)/2 : // packed storage for S
              dims->sdims[j]);                          // diagonal storage for S (kind == 2|3)
    }
    ind->slen = dims->slen;

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
cvxc_index_t *cvxc_index_init(cvxc_index_t *ind, const cvxc_dimset_t *dims, int kind)
{
    if (! ind)
        return ind;
    if (! dims)
        return (cvxc_index_t *)0;

    cvxc_size_t nb = cvxc_index_bytes(dims, kind);
    void *mem = calloc(nb, 1);
    if (!mem)
        return (cvxc_index_t *)0;

    cvxc_index_make(ind, dims, kind, mem, nb);
    ind->__bytes = mem;
    return ind;

}

void cvxc_index_release(cvxc_index_t *ind)
{
    if (! ind)
        return;
    if (ind->__bytes) {
        free(ind->__bytes);
        ind->__bytes = (void *)0;
    }
    ind->index = ind->indnlt = ind->indnl = ind->indl = ind->indq = ind->inds = (cvxc_size_t *)0;
}

cvxc_size_t cvxc_index_count(const cvxc_index_t *ind,
                           cvxc_dim_enum name)
{
    switch (name) {
    case CVXDIM_SOCP:
        return ind->qlen;
    case CVXDIM_SDP:
        return ind->slen;
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

cvxc_size_t cvxc_index_length(const cvxc_index_t *ind, cvxc_dim_enum name)
{
    switch (name) {
    case CVXDIM_SOCP:
        return ind->indq ? ind->indq[ind->qlen] - ind->indq[0] : 0;
    case CVXDIM_SDP:
        return ind->inds ? ind->inds[ind->slen] - ind->inds[0] : 0;
    case CVXDIM_LINEAR:
        return ind->indl ? ind->indl[1] -ind->indl[0] : 0;
    case CVXDIM_NONLINEAR:
        return ind->indnl ? ind->indnl[1] - ind->indnl[0] : 0;
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
 *     Flags from cvxc_dim_enum set to indicate which part are included.
 *
 */
void cvxc_subindex(cvxc_index_t *ind, const cvxc_index_t *src, int parts)
{
    ind->__bytes = (void *)0;
    ind->index = ind->indnlt = ind->indnl = ind->indl = ind->indq = ind->inds = (cvxc_size_t *)0;
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
cvxc_size_t cvxc_index_elem(cvxc_matrix_t *x,
                          const cvxc_matrix_t *y,
                          const cvxc_index_t *ind,
                          cvxc_dim_enum name,
                          int k)
{
    cvxc_size_t n = 0, m = 0;

    if (!y || !ind)
        return 0;

    if (x)
        cvxm_map_data(x, 0, 0, (cvxc_float_t *)0);

    switch (name) {
    case CVXDIM_CONVEX:
        // map non-linear space + non-linear target
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
            //m = ind->dims->sdims[k];
            n = ind->inds[k+1] - ind->inds[k];
            m = ind->type == CVXC_INDEX_NORMAL ? (cvxc_size_t)floor(sqrt(n)) : n;
            // TODO: what if indexing is not standard storage??
            if (x) {
                // standard vs. diagonal storage
                n = ind->type == CVXC_INDEX_NORMAL ? m : 1;
                cvxm_map_data(x, m, n, cvxm_data(y, ind->inds[k]));
            }
        }
        break;
    case CVXDIM_CONELP:
        // map linear,socp and sdp parts onto x; find start index
        if (ind->indl)
            n = ind->indl[0];  // have linear part
        else if (ind->indq)
            n = ind->indq[0];  // have socp part
        else if (ind->inds)
            n = ind->inds[0];  // have sdp part
        else
            n = ind->index[ind->indlen];
        m = ind->index[ind->indlen] - n;
        if (x)
            cvxm_map_data(x, m, 1, cvxm_data(y, n));
        break;
    case CVXDIM_CONVEXLP:
        // map all but non-linear target function; find start index
        if (ind->indnl)
            n = ind->indnl[0];  // have linear part
        else if (ind->indl)
            n = ind->indl[0];  // have linear part
        else if (ind->indq)
            n = ind->indq[0];  // have socp part
        else if (ind->inds)
            n = ind->inds[0];  // have sdp part
        else
            n = ind->index[ind->indlen];
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


void cvxc_index_create(cvxc_matrix_t *x, cvxc_index_t *index, const cvxc_dimset_t *dims, cvxc_index_type kind)
{
    cvxc_size_t cdim = cvxc_dimset_sum(dims, CVXDIM_NONLINEAR) +
        cvxc_dimset_sum(dims, CVXDIM_LINEAR) +
        cvxc_dimset_sum(dims, CVXDIM_SOCP);
    switch (kind) {
    case CVXC_INDEX_PACKED:
        cdim += cvxc_dimset_sum_packed(dims, CVXDIM_SDP);
        break;
    case CVXC_INDEX_DIAG:
        cdim += cvxc_dimset_sum(dims, CVXDIM_SDP);
        break;
    case CVXC_INDEX_SIGS:
        cdim = cvxc_dimset_sum(dims, CVXDIM_SDP);
        break;
    case CVXC_INDEX_NORMAL:
    default:
        cdim += cvxc_dimset_sum_squared(dims, CVXDIM_SDP);
        break;
    }
    cvxm_init(x, cdim, 1);
    cvxc_index_init(index, dims, kind);
}
