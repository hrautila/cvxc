
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include <stdio.h>
#include "convex.h"


void test_dims()
{
    cvxc_dimset_t dims;
    cvxc_index_t index;
    
    cvxc_dimset_alloc(&dims, 4, (int[]){6, 0}, (int[]){4, 0});
    cvxc_index_init(&index, &dims, CVXC_INDEX_NORMAL);

    printf("sum(L) = %ld\n", cvxc_dimset_sum(&dims, L));
    printf("sum(Q) = %ld\n", cvxc_dimset_sum(&dims, Q));
    printf("sum(S) = %ld\n", cvxc_dimset_sum_squared(&dims, S));

    printf("offset(L) = %ld\n", cvxc_index_count(&index, CVXDIM_LINEAR) ? index.indl[0] : 0);
    for (int k = 0; k < cvxc_index_count(&index, CVXDIM_SOCP); k++) {
        printf("offset(Q) = %ld\n", index.indq[k]);
    }
    for (int k = 0; k < cvxc_index_count(&index, CVXDIM_SDP); k++) {
        printf("offset(S) = %ld\n", index.inds[k]);
    }

    cvxc_index_release(&index);
    cvxc_dimset_release(&dims);
}

void test_index()
{
    cvxc_dimset_t dims;
    cvxc_index_t index;
    cvxc_matrix_t x;
    
    cvxc_dimset_alloc(&dims, 4, (int[]){6, 0}, (int[]){4, 0});
    cvxc_index_create(&x, &index, &dims, CVXC_INDEX_NORMAL);

    printf("sum(L) = %ld\n", cvxc_dimset_sum(&dims, L));
    printf("sum(Q) = %ld\n", cvxc_dimset_sum(&dims, Q));
    printf("sum(S) = %ld\n", cvxc_dimset_sum_squared(&dims, S));

    printf("offset(L) = %ld\n", cvxc_index_count(&index, CVXDIM_LINEAR) ? index.indl[0] : 0);
    for (int k = 0; k < cvxc_index_count(&index, CVXDIM_SOCP); k++) {
        printf("offset(Q) = %ld\n", index.indq[k]);
    }
    for (int k = 0; k < cvxc_index_count(&index, CVXDIM_SDP); k++) {
        printf("offset(S) = %ld\n", index.inds[k]);
    }
    printf("size(x) = (%d, %d)\n", x.rows, x.cols);
    
    cvxc_index_release(&index);
    cvxc_dimset_release(&dims);
}

int main(int argc, char **argv)
{
    test_dims();
    test_index();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
