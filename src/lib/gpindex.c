/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */
#include "cvxc.h"

cvxc_size_t cvxc_gpi_make(cvxc_gpindex_t *gpi, cvxc_size_t p, void *memory, cvxc_size_t nbytes)
{
    if (nbytes < (p+1)*sizeof(cvxc_size_t))
        return 0;

    gpi->p = p;
    gpi->index = (cvxc_size_t *)memory;
    gpi->index[0] = 0;
    for (cvxc_size_t k = 0; k < p; k++) {
        gpi->index[k+1] = 0;
    }
    return (p+1)*sizeof(cvxc_size_t);
}

/**
 * @brief Get the n'th element in GP matrix F/g.
 *
 * @return Row index of element
 */
cvxc_size_t cvxc_gpi_elem(const cvxc_gpindex_t *gpi, cvxc_matrix_t *e, const cvxc_matrix_t *Fg, cvxc_size_t n)
{
    cvxc_size_t mS, nS;
    if (n < gpi->p) {
        cvxm_size(&mS, &nS, Fg);
        cvxm_view_map(e, Fg, gpi->index[n], 0, gpi->index[n+1]-gpi->index[n], nS);
        return gpi->index[n];
    }
    return 0;
}

/**
 * @brief Get length of n'th GP index element.
 */
cvxc_size_t cvxc_gpi_length(const cvxc_gpindex_t *gpi, cvxc_size_t n)
{
    return n < gpi->p ? gpi->index[n+1] - gpi->index[n] : 0;
}

int cvxc_gpi_setup(cvxc_gpindex_t *gpi, const cvxc_size_t *K, cvxc_size_t p)
{
    gpi->index[0] = 0;
    for (cvxc_size_t k = 0; k < p; k++) {
        gpi->index[k+1] = gpi->index[k] + K[k];
    }
    gpi->p = p;
    return 0;
}

int cvxc_gpi_init(cvxc_gpindex_t *gpi, const cvxc_size_t *K, cvxc_size_t p)
{
    gpi->__bytes = calloc(p + 1, sizeof(cvxc_size_t));
    if (!gpi->__bytes)
        return -1;

    gpi->index = (cvxc_size_t *)gpi->__bytes;
    return cvxc_gpi_setup(gpi, K, p);
}

void cvxc_gpi_release(cvxc_gpindex_t *gpi)
{
    if (!gpi)
        return;
    if (gpi->__bytes)
        free(gpi->__bytes);
    gpi->index = 0;
    gpi->__bytes = 0;
}

void cvxc_gpi_free(cvxc_gpindex_t *gpi)
{
    if (gpi) {
        cvxc_gpi_release(gpi);
        free(gpi);
    }
}
