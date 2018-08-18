
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include <stdio.h>
#include "convex.h"


void cvx_mgrp_printf(FILE *f, const char *format, cvx_matgrp_t *g, const char *s)
{
    cvx_matrix_t xk;
    const char *frm = format ? format : "%9.2e";
    
    if (cvx_mgrp_elem(__cvxnil, g, CVXDIM_LINEAR, 0) > 0) {
        cvx_mgrp_elem(&xk, g, CVXDIM_LINEAR, 0);
        if (s)
            fprintf(f, "L(%s)\n", s);
        else 
            fprintf(f, "L()\n");
        cvxm_printf(f, frm, &xk);
    }

    if (cvx_mgrp_count(g, CVXDIM_SOCP) > 0) {
        cvx_mgrp_elem(&xk, g, CVXDIM_SOCP, 0);
        if (s)
            fprintf(f, "Q(%s)\n", s);
        else
            fprintf(f, "Q()\n");
        cvxm_printf(f, frm, &xk);
    }

    if (cvx_mgrp_count(g, CVXDIM_SDP) > 0) {
            cvx_mgrp_elem(&xk, g, CVXDIM_SDP, 0);
        if (s)
            fprintf(f, "S(%s)\n", s);
        else
            fprintf(f, "S()\n");
        cvxm_printf(f, frm, &xk);
    }
}

void cvx_mat_printf(FILE *f, const char *format, cvx_matrix_t *g, const char *s)
{
    const char *frm = format ? format : "%9.2e";
    if (s)
        fprintf(f, "%s\n", s);
    cvxm_printf(f, frm, g);
}

void cvx_scaling_printf(FILE *f, const char *format, cvx_scaling_t *W, const char *s)
{
    cvx_matrix_t m;
    char buf[64];
    fprintf(f, "** start scaling %s\n", s);

    cvx_scaling_elem(&m, W, CVXWS_D, 0);
    cvx_mat_printf(f, format, &m, "D");
    cvx_scaling_elem(&m, W, CVXWS_DI, 0);
    cvx_mat_printf(f, format, &m, "DI");

    if (W->vcount > 0) {
        cvx_scaling_elem(&m, W, CVXWS_BETA, 0);   
        cvx_mat_printf(f, format, &m, "beta");
        for (int i = 0; i < W->vcount; i++) {
            cvx_scaling_elem(&m, W, CVXWS_V, i);   
            sprintf(buf, "V[%d]", i);
            cvx_mat_printf(f, format, &m, buf);
        }
    }

    if (W->rcount > 0) {
        for (int i = 0; i < W->rcount; i++) {
        }
    }
    fprintf(f, "** end scaling %s\n", s);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
