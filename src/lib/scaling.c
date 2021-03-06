/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */
#include "cvxc.h"

/**
 * @brief Compute memory size of scaling matrix W for spesified dimension set.
 *
 * @param[out] isize
 *     If non-null pointer then set to number of bytes needed for scaling matrix indexes.
 * @param[in] dims
 *     Dimension set
 */
cvxc_size_t cvxc_scaling_bytes(cvxc_size_t *isize, const cvxc_dimset_t *dims)
{
    // SOCP space; for V matrices
    cvxc_size_t vspace = cvxc_dimset_sum(dims, CVXDIM_SOCP);

    // SDP space; for R and RTI matrices
    cvxc_size_t rspace = cvxc_dimset_sum_squared(dims, CVXDIM_SDP);

    // calculate index space length in bytes and align to 64bits;
    // (3 is for linear, non-linear and non-linear target spaces.)
    // SDP saved as tuple (offset, size)
    cvxc_size_t itotal = (3 + dims->qlen + 1 + 4*dims->slen)*sizeof(cvxc_size_t);
    itotal += (itotal & 0x7) != 0 ? 8 - (itotal & 0x7) : 0;

    // matrix data space length
    cvxc_size_t ntotal =
        2 * dims->iscpt +       // Non-linear target function
        2 * dims->mnl +         // Non-linear DNL and DNLI vectors
        2 * dims->ldim +        // Linear D and DI vectors
        dims->qlen + vspace +   // SOCP BETAs and V vectors
        2 * rspace;             // SDP R and RTI matrices
    ntotal *= sizeof(cvxc_float_t);

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
cvxc_size_t cvxc_scaling_make(cvxc_scaling_t *W, const cvxc_dimset_t *dims, void *mem, size_t nbytes)
{
    cvxc_size_t ntotal, itotal;
    ntotal = cvxc_scaling_bytes(&itotal, dims);

    if (nbytes < ntotal) {
        return 0;
    }

    // set __bytes to null to indicate that memory block is not owned by scaling matrix
    W->__bytes = (unsigned char *)0;
    unsigned char *buf = (unsigned char *)mem;
    W->nbytes = ntotal;
    W->indexes = (cvxc_size_t *)buf;
    W->data    = (cvxc_float_t *)&buf[itotal];

    // initialized pointers and indexes.
    if (dims->qlen > 0 || dims->slen > 0) {
        W->vcount = dims->qlen;
        W->rcount = dims->slen;
        if (dims->qlen > 0) {
            W->indv = W->indexes;
        }
        if (dims->slen > 0) {
            W->indr = &W->indexes[dims->qlen == 0 ? 0 : dims->qlen+1];
            W->indrti = &W->indexes[(dims->qlen == 0 ? 0 : dims->qlen+1) + 2*dims->slen];
        }
    } else {
        W->indexes = W->indv = W->indr = W->indrti = (cvxc_size_t *)0;
        W->vcount = W->rcount = 0;
    }

    // setup indexes
    cvxc_size_t offset = 0;
    W->dnltsz = dims->iscpt ? 1 : 0;
    W->dnlt = W->dnlti = (cvxc_float_t *)0;
    W->dnlsz = dims->mnl;
    W->dnl = W->dnli = (cvxc_float_t *)0;
    W->dsz = dims->ldim;
    W->d = W->di = (cvxc_float_t *)0;

    // Arrange NLTARGET, NONLINEAR and LINEAR spaces to continous
    // memory blocks for direct and inverse cases.
    // Memory layout:
    //   |DNLT|DNL|D|DNLTI|DNLI|DI|BETA|V[0]| ..|V[q]|R[0]|..|R[n]|RTI[0]|..|RTI[n]|
    //
    if (dims->iscpt > 0) {
        // Direct pointers to DNLT data space
        W->dnlt  = &W->data[offset];
        offset  += W->dnltsz;
    }
    if (dims->mnl > 0) {
        // Direct pointers to DNL data space
        W->dnl  = &W->data[offset];
        offset += dims->mnl;
    }
    if (dims->ldim > 0) {
        // Direct pointers to D data space
        W->d = &W->data[offset];
        offset += dims->ldim;
    }
    // Inverse pointers
    if (dims->iscpt > 0) {
        // Direct pointers to DNLTI data space
        W->dnlti = &W->data[offset];
        offset  += W->dnltsz;
    }
    if (dims->mnl > 0) {
        // Direct pointers to DNLI data space
        W->dnli = &W->data[offset];
        offset += dims->mnl;
    }
    if (dims->ldim > 0) {
        // Direct pointers to DI data space
        W->di = &W->data[offset];
        offset += dims->ldim;
    }

    if (dims->qlen > 0) {
        // Direct pointer to BETA vector data space
        W->beta = &W->data[offset];
        offset += dims->qlen;
        // Sizes of and offsets to V vector data space
        for (cvxc_size_t k = 0; k < dims->qlen; k++) {
            W->indv[k] = offset;
            offset += dims->qdims[k];
        }
        // offset past last entry to allow length computation
        W->indv[dims->qlen] = offset;
    }
    if (dims->slen > 0) {
        // Save sizes of and offsets to R/RTI matrices data space
        for (cvxc_size_t k = 0; k < 2*dims->slen; k += 2) {
            W->indr[k] = offset;
            W->indr[k+1] = dims->sdims[k];
            offset += dims->sdims[k] * dims->sdims[k];
        }
        for (cvxc_size_t k = 0; k < 2*dims->slen; k += 2) {
            W->indrti[k] = offset;
            W->indrti[k+1] = dims->sdims[k];
            offset += dims->sdims[k] * dims->sdims[k];
        }
    }
    return ntotal;
}

/*
 *  Scaling matrix W, one big allocation block that is divided to index space
 *  and matrix data space.
 */

/**
 * @brief Initialize scaling matrix for spesified dimension set.
 *
 * @return
 *    Pointer to scaling matrix or null if memory allocation failed.
 */
cvxc_scaling_t *cvxc_scaling_init(cvxc_scaling_t *W, const cvxc_dimset_t *dims)
{
    cvxc_size_t ntotal, itotal;
    ntotal = cvxc_scaling_bytes(&itotal, dims);

    unsigned char *buf = calloc(ntotal, sizeof(char));
    if (!buf)
        return (cvxc_scaling_t *)0;
    cvxc_scaling_make(W, dims, buf, ntotal);
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
cvxc_scaling_t *cvxc_scaling_new(const cvxc_dimset_t *dims)
{
    cvxc_scaling_t *W = (cvxc_scaling_t *)malloc(sizeof(cvxc_scaling_t));
    if (!W)
        return W;

    if (!cvxc_scaling_init(W, dims)) {
        free(W);
        W = (cvxc_scaling_t *)0;
    }
    return W;
}

/**
 * @brief Release resources reserved for scaling matrix and clear indexing.
 */
void cvxc_scaling_release(cvxc_scaling_t *W)
{
    if (!W)
        return;
    if (W->__bytes) {
        free(W->__bytes);
        W->__bytes = (unsigned char *)0;
    }
    W->nbytes = 0;
    W->data = W->dnl = W->dnli = W->d = W->di = W->beta = (cvxc_float_t *)0;
    W->dnlsz = W->dsz = W->vcount = W->rcount = 0;
    W->indexes = W->indv = W->indr = W->indrti = (cvxc_size_t *)0;
}

void cvxc_scaling_free(cvxc_scaling_t *W)
{
    cvxc_scaling_release(W);
    if (W)
        free(W);
}

int cvxc_scaling_copy(cvxc_scaling_t *W, const cvxc_scaling_t *Ws)
{
    if (W->nbytes != Ws->nbytes)
        return -1;
    if (W->__bytes && Ws->__bytes) {
        memcpy(W->__bytes, Ws->__bytes, W->nbytes);
    } else if (Ws->indexes && W->indexes) {
        // memory block not owned by structure; layout is [indexes; data]
        // copy nbytes from the start of indexes.
        memcpy(W->indexes, Ws->indexes, W->nbytes);
    }
    return 0;
}

// \brief Get scaling element NAME at index ind; returns size of element
cvxc_size_t cvxc_scaling_elem(cvxc_matrix_t *A, const cvxc_scaling_t *W, cvxc_mset_enum name, int ind)
{
    if (!W || !A)
        return 0;

    cvxm_map_data(A, 0, 0, (cvxc_float_t *)0);

    switch (name) {
    case CVXWS_DNLT:
        if (W->dnltsz > 0) {
            cvxm_map_data(A, W->dnltsz, 1, W->dnlt);
            return W->dnlsz;
        }
    case CVXWS_DNLTI:
        if (W->dnltsz > 0) {
            cvxm_map_data(A, W->dnltsz, 1, W->dnlti);
            return W->dnlsz;
        }
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
            cvxc_size_t n = W->indv[ind+1] - W->indv[ind];
            cvxm_map_data(A, n, 1, &W->data[W->indv[ind]]);
            return n;
        }
    case CVXWS_R:
        if (ind < W->rcount) {
            cvxc_size_t n = W->indr[2*ind+1];
            cvxm_map_data(A, n, n, &W->data[W->indr[2*ind]]);
            return n;
        }
    case CVXWS_RTI:
        if (ind < W->rcount) {
            cvxc_size_t n = W->indrti[2*ind+1];
            cvxm_map_data(A, n, n, &W->data[W->indrti[2*ind]]);
            return n;
        }
    }
    return 0;
}

/*
 * @brief Set scaling matrix to initial value
 *
 * d, di to ones, beta_j = 1, v_j = e_1, r_k, rti_k = diag(1)
 *
 */
void cvxc_scaling_initial_value(cvxc_scaling_t *W)
{
    cvxc_matrix_t x;
    // D = ones
    cvxc_scaling_elem(&x, W, CVXWS_D, 0);
    cvxm_mkconst(&x, 1.0);
    // DI = ones
    cvxc_scaling_elem(&x, W, CVXWS_DI, 0);
    cvxm_mkconst(&x, 1.0);
    if (W->vcount > 0) {
        // BETA = ones
        cvxc_scaling_elem(&x, W, CVXWS_BETA, 0);
        cvxm_mkconst(&x, 1.0);
        // V is unit vector
        for (int k = 0; k < W->vcount; k++) {
            cvxc_scaling_elem(&x, W, CVXWS_V, k);
            cvxm_mkident(&x);
        }
    }
    // R & RTI are identity
    for (int k = 0; k < W->rcount; k++) {
        cvxc_scaling_elem(&x, W, CVXWS_R, k);
        cvxm_mkident(&x);
        cvxc_scaling_elem(&x, W, CVXWS_RTI, k);
        cvxm_mkident(&x);
    }
}
