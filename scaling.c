
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include "convex.h"
#include "cvxm.h"

/**
 * @brief Compute memory size of scaling matrix W for spesified dimension set.
 *
 * @param[out] isize
 *     If non-null pointer then set to number of bytes needed for scaling matrix indexes.
 * @param[in] dims
 *     Dimension set
 */
cvx_size_t cvx_scaling_bytes(cvx_size_t *isize, const cvx_dimset_t *dims)
{
    // SOCP space; for V matrices
    cvx_size_t vspace = cvx_dimset_sum(dims, CVXDIM_SOCP);

    // SDP space; for R and RTI matrices
    cvx_size_t rspace = cvx_dimset_sum_squared(dims, CVXDIM_SDP);

    // calculate index space length in bytes and align to 64bits;
    // (2 is for linear and non-linear spaces.)
    cvx_size_t itotal = (2 + dims->qlen + 2*dims->slen)*sizeof(cvx_size_t);
    itotal += (itotal & 0x7) != 0 ? 8 - (itotal & 0x7) : 0;
    
    // matrix data space length
    cvx_size_t ntotal =
        2 * dims->mnl +         // Non-linear DNL and DNLI vectors
        2 * dims->ldim +        // Linear D and DI vectors
        dims->qlen + vspace +   // SOCP BETAs and V vectors
        2 * rspace;             // SDP R and RTI matrices
    ntotal *= sizeof(cvx_float_t);

    // if non-null pointer return 
    if (isize)
        *isize = itotal;
    return ntotal + itotal;
}

/**
 * @brief Overlay scaling matrix W on provided memory block of size nbytes.
 *
 * @param[in,out] W
 *   On entry uninitialized scaling matrix. On exit initialized scaling.
 * @param[in] dims
 *   Dimension set.
 * @param[in] mem
 *   Pointer to memory.
 * @param[in] nbytes
 *   Size of memory block.
 *
 * @return
 *   Number of bytes used from memory.
 */
cvx_size_t cvx_scaling_make(cvx_scaling_t *W, const cvx_dimset_t *dims, void *mem, size_t nbytes)
{
    cvx_size_t ntotal, itotal;
    ntotal = cvx_scaling_bytes(&itotal, dims);

    if (nbytes < ntotal)
        return 0;

    // set __bytes to null to indicate that memory block is not owned by scaling matrix
    W->__bytes = (unsigned char *)0;
    unsigned char *buf = (unsigned char *)mem;
    W->nbytes = ntotal;
    W->indexes = (cvx_size_t *)buf;
    W->data    = (cvx_float_t *)&buf[itotal];
    
    // initialized pointers and indexes.    
    if (dims->qlen > 0 || dims->slen > 0) {
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
    W->dnl = W->dnli = (cvx_float_t *)0;
    if (dims->mnl > 0) {
        // Direct pointers to DNL/DNLI data space
        W->dnl = W->data;
        offset += dims->mnl;
        W->dnli = &W->data[offset];
        offset += dims->mnl;
    }
    W->dsz = dims->ldim;
    W->d = W->di = (cvx_float_t *)0;
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
    return ntotal;
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

/**
 * @brief Initialize scaling matrix for spesified dimension set.
 *
 * @return
 *    Pointer to scaling matrix or null if memory allocation failed.
 */
cvx_scaling_t *cvx_scaling_init(cvx_scaling_t *W, const cvx_dimset_t *dims)
{
    cvx_size_t ntotal, itotal;
    ntotal = cvx_scaling_bytes(&itotal, dims);
    
    unsigned char *buf = calloc(ntotal, sizeof(char));
    if (!buf)
        return (cvx_scaling_t *)0;
    cvx_scaling_make(W, dims, buf, ntotal);
    // set __bytes to indicate that memory is owned by this scaling matrix
    W->__bytes = buf;
    return W;
}

/**
 * @brief Create a new scaling matrix for dimension set.
 *
 * @return 
 *    Pointer to 
 */
cvx_scaling_t *cvx_scaling_new(const cvx_dimset_t *dims)
{
    cvx_scaling_t *W = (cvx_scaling_t *)malloc(sizeof(cvx_scaling_t));
    if (!W)
        return W;

    if (!cvx_scaling_init(W, dims)) {
        free(W);
        W = (cvx_scaling_t *)0;
    }
    return W;
}

/**
 * @brief Release resources reserved for scaling matrix and clear indexing.
 */
void cvx_scaling_release(cvx_scaling_t *W)
{
    if (!W)
        return;
    if (W->__bytes) {
        free(W->__bytes);
        W->__bytes = (unsigned char *)0;
    }
    W->nbytes = 0;
    W->data = W->dnl = W->dnli = W->d = W->di = W->beta = (cvx_float_t *)0;
    W->dnlsz = W->dsz = W->vcount = W->rcount = 0;
    W->indexes = W->indv = W->indr = W->indrti = (cvx_size_t *)0;
}

void cvx_scaling_free(cvx_scaling_t *W)
{
    cvx_scaling_release(W);
    if (W)
        free(W);
}

// \brief Get scaling element NAME at index ind; returns size of element
cvx_size_t cvx_scaling_elem(cvx_matrix_t *A, const cvx_scaling_t *W, cvx_mset_enum name, int ind)
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



// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
