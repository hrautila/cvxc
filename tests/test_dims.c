
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include <stdio.h>
#include "convex.h"


void test_dims()
{
    cvx_dimset_t dims;
    cvx_index_t index;
    
    cvx_dimset_alloc(&dims, 4, (int[]){6, 0}, (int[]){4, 0});
    cvx_index_init(&index, &dims, CVX_INDEX_NORMAL);

    printf("sum(L) = %ld\n", cvx_dimset_sum(&dims, L));
    printf("sum(Q) = %ld\n", cvx_dimset_sum(&dims, Q));
    printf("sum(S) = %ld\n", cvx_dimset_sum_squared(&dims, S));

    printf("offset(L) = %ld\n", cvx_index_count(&index, CVXDIM_LINEAR) ? index.indl[0] : 0);
    for (int k = 0; k < cvx_index_count(&index, CVXDIM_SOCP); k++) {
        printf("offset(Q) = %ld\n", index.indq[k]);
    }
    for (int k = 0; k < cvx_index_count(&index, CVXDIM_SDP); k++) {
        printf("offset(S) = %ld\n", index.inds[k]);
    }

    cvx_index_release(&index);
    cvx_dimset_release(&dims);
}

void test_index()
{
    cvx_dimset_t dims;
    cvx_index_t index;
    cvx_matrix_t x;
    
    cvx_dimset_alloc(&dims, 4, (int[]){6, 0}, (int[]){4, 0});
    cvx_index_create(&x, &index, &dims, CVX_INDEX_NORMAL);

    printf("sum(L) = %ld\n", cvx_dimset_sum(&dims, L));
    printf("sum(Q) = %ld\n", cvx_dimset_sum(&dims, Q));
    printf("sum(S) = %ld\n", cvx_dimset_sum_squared(&dims, S));

    printf("offset(L) = %ld\n", cvx_index_count(&index, CVXDIM_LINEAR) ? index.indl[0] : 0);
    for (int k = 0; k < cvx_index_count(&index, CVXDIM_SOCP); k++) {
        printf("offset(Q) = %ld\n", index.indq[k]);
    }
    for (int k = 0; k < cvx_index_count(&index, CVXDIM_SDP); k++) {
        printf("offset(S) = %ld\n", index.inds[k]);
    }
    printf("size(x) = (%d, %d)\n", x.rows, x.cols);
    
    cvx_index_release(&index);
    cvx_dimset_release(&dims);
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
