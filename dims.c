
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include "convex.h"
#include "cvxm.h"


cvx_dimset_t *cvx_dimset_alloc(cvx_dimset_t *dims, int linear, int *socp, int *sdp)
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

cvx_size_t cvx_dimset_max(cvx_dimset_t *dims, cvx_dim_enum name)
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
    }
    return m;
}

cvx_size_t cvx_dimset_sum(cvx_dimset_t *dims, cvx_dim_enum name)
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
    }
    return sum;
}

cvx_size_t cvx_dimset_sum_squared(cvx_dimset_t *dims, cvx_dim_enum name)
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
            sum += dims->sdims[k] * dims->sdims[k];
        break;
    }
    return sum;
}

cvx_size_t cvx_dimset_sum_packed(cvx_dimset_t *dims, cvx_dim_enum name)
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
    }
    return sum;
}

/*
 * \brief Compute indexing set from spesified dimension set.
 *
 * Parameter 'kind' spesifies indexing type for SDP ('S') space. 
 *  0  : standard storage indexing (m*m elements)
 *  1  : packed storage indexing (m*(m+1)/2 elements)
 *  2  : diagonal storage indexing (m elements)
 *  3  : only 'S' space diagonal storage indexing
 */
cvx_index_t *cvx_index_init(cvx_index_t *ind, cvx_dimset_t *dims, int kind)
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

cvx_size_t cvx_index_count(cvx_index_t *ind,
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
 * \brief Make 'x' to present k'th element of y's dimension set.
 *
 * If x is null then returns dimensionality of k'th element in the spesified 
 * dimension set.
 */
cvx_size_t cvx_index_elem(cvx_matrix_t *x,
                          cvx_matrix_t *y,
                          cvx_index_t *ind,
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


void cvx_index_create(cvx_matrix_t *x, cvx_index_t *index, cvx_dimset_t *dims, cvx_index_type kind)
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

/*
 *  Scaling matrix W, one big allocation block that is divided to index space
 *  and matrix data space.
 *
 *  Size of index space is: 
 *       (N_SOCP + 2*N_SDP)*sizeof(size_t) aligned to 64bit
 *  Size of matrix data space: 
 *       (2*(NONLINEAR + LINEAR) + sum(SOCP_i) + N_SOCP + 2*sum(SDP_i^2)) * sizeof(cvx_float_t)
 *
 */

cvx_scaling_t *cvx_scaling_alloc(cvx_scaling_t *W, cvx_dimset_t *dims)
{
    // SOCP space; for V matrices
    cvx_size_t vspace = cvx_dimset_sum(dims, CVXDIM_SOCP);

    // SDP space; for R and RTI matrices
    cvx_size_t rspace = cvx_dimset_sum_squared(dims, CVXDIM_SDP);

    // calculate index space length in bytes and align to 64bits;
    cvx_size_t itotal = (dims->qlen + 2*dims->slen)*sizeof(cvx_size_t);
    itotal += (itotal & 0x7) == 0 ? 8 - (itotal & 0x7) : 0;
    
    // matrix data space length
    cvx_size_t ntotal =
        2 * dims->mnl +         // Non-linear DNL and DNLI vectors
        2 * dims->ldim +        // Linear D and DI vectors
        dims->qlen + vspace +   // SOCP BETAs and V vectors
        2 * rspace;             // SDP R and RTI matrices
    ntotal *= sizeof(cvx_float_t);

    unsigned char *buf = calloc(itotal + ntotal, sizeof(char));
    if (!buf)
        return (cvx_scaling_t *)0;
    W->__bytes = buf;
    W->nbytes = itotal + ntotal;
#if 0    
    W->data = (cvx_float_t *)calloc(ntotal, sizeof(cvx_float_t));
    if (!W->data)
        return (cvx_scaling_t *)0;
#endif    
    W->indexes = (cvx_size_t *)W->__bytes;
    W->data    = (cvx_float_t *)&W->__bytes[itotal];
    if (dims->qlen > 0 || dims->slen > 0) {
#if 0
        // allocate space to save indexes to matrix data space
        W->indexes = (cvx_size_t *)calloc(dims->qlen + 2*dims->slen, sizeof(cvx_size_t));
        if (!W->indexes) {
            free(W->data);
            W->data = (cvx_float_t *)0;
            return (cvx_scaling_t *)0;
        }
#endif
        W->vcount = dims->qlen;
        W->rcount = dims->slen;
        if (dims->qlen > 0) {
            W->indv = W->indexes;
        }
        if (dims->slen > 0) {
            W->indr = &W->indexes[dims->qlen];
            W->indrti = &W->indexes[dims->qlen + dims->slen];
        }
    } else {
        W->indexes = W->indv = W->indr = W->indrti = (cvx_size_t *)0;
        W->vcount = W->rcount = 0;
    }
            
    W->dims = dims;
    
    // setup indexes
    cvx_size_t offset = 0;
    W->dnlsz = dims->mnl;
    if (dims->mnl > 0) {
        // Direct pointers to DNL/DNLI data space
        W->dnl = W->data;
        offset += dims->mnl;
        W->dnli = &W->data[offset];
        offset += dims->mnl;
    }
    W->dsz = dims->ldim;
    if (dims->ldim > 0) {
        // Direct pointers to D/DI data space
        W->d = &W->data[offset];
        offset += dims->ldim;
        W->di = &W->data[offset];
        offset += dims->ldim;
    }
    if (dims->qlen > 0) {
        // Direct pointer to BETA vector data space
        W->beta = &W->data[offset];
        offset += dims->qlen;
        // Sizes of and offsets to V vector data space
        for (int k = 0; k < dims->qlen; k++) {
            W->indv[k] = offset;
            offset += dims->qdims[k];
        }
    }
    if (dims->slen > 0) {
        // Save sizes of and offsets to R/RTI matrices data space
        for (int k = 0; k < dims->slen; k++) {
            W->indr[k] = offset;
            offset += dims->sdims[k] * dims->sdims[k];
            W->indrti[k] = offset;
            offset += dims->sdims[k] * dims->sdims[k];
        }
    }
    return W;
}

// \brief Get scaling element NAME at index ind; returns size of element
cvx_size_t cvx_scaling_elem(cvx_matrix_t *A, cvx_scaling_t *W, cvx_mset_enum name, int ind)
{
    if (!W || !A)
        return 0;

    cvxm_map_data(A, 0, 0, (cvx_float_t *)0);
    
    switch (name) {
    case CVXWS_DNL:
        if (W->dnlsz > 0) {
            cvxm_map_data(A, W->dnlsz, 1, W->dnl);
            return W->dnlsz;
        }
    case CVXWS_DNLI:
        if (W->dnlsz > 0) {
            cvxm_map_data(A, W->dnlsz, 1, W->dnli);
            return W->dnlsz;
        }
    case CVXWS_D:
        if (W->dsz > 0) {
            cvxm_map_data(A, W->dsz, 1, W->d);
            return W->dsz;
        }
    case CVXWS_DI:
        if (W->dsz > 0) {
            cvxm_map_data(A, W->dsz, 1, W->di);
            return W->dsz;
        }
    case CVXWS_BETA:
        if (W->vcount > 0) {
            cvxm_map_data(A, W->vcount, 1, W->beta);
            return W->vcount;
        }
    case CVXWS_V:
        if (ind < W->vcount) {
            cvxm_map_data(A, W->dims->qdims[ind], 1, &W->data[W->indv[ind]]);
            return W->indv[ind];
        }
    case CVXWS_R:
        if (ind < W->rcount) {
            cvxm_map_data(A, W->dims->sdims[ind], W->dims->sdims[ind], &W->data[W->indr[ind]]);
            return W->dims->sdims[ind];
        }
    case CVXWS_RTI:
        if (ind < W->rcount) {
            cvxm_map_data(A, W->dims->sdims[ind], W->dims->sdims[ind], &W->data[W->indrti[ind]]);
            return W->dims->sdims[ind];
        }
    }
    return 0;
}


void cvx_scaling_release(cvx_scaling_t *W)
{
    if (!W)
        return;
    if (W->__bytes) {
        free(W->__bytes);
        W->__bytes = (unsigned char *)0;
    }
        
#if 0
    if (W->data)
        free(W->data);
    if (W->indexes)
        free(W->indexes);
#endif
    W->data = W->dnl = W->dnli = W->d = W->di = W->beta = (cvx_float_t *)0;
    W->dnlsz = W->dsz = W->vcount = W->rcount = 0;
    W->indexes = W->indv = W->indr = W->indrti = (cvx_size_t *)0;
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
