
// Copyright: Harri Rautila, 2018 <harri.rautila@gmail.com>

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctype.h>
#include <unistd.h>
#include "cvxc.h"

void cvxc_mgrp_printf(FILE *f, const char *format, cvxc_matgrp_t *g, const char *s)
{
    cvxc_matrix_t xk;
    const char *frm = format ? format : "%9.2e";

    if (cvxc_mgrp_elem(__cvxnil, g, CVXDIM_LINEAR, 0) > 0) {
        cvxc_mgrp_elem(&xk, g, CVXDIM_LINEAR, 0);
        if (s)
            fprintf(f, "L(%s)\n", s);
        else
            fprintf(f, "L()\n");
        cvxm_printf(f, frm, &xk);
    }

    if (cvxc_mgrp_count(g, CVXDIM_SOCP) > 0) {
        cvxc_mgrp_elem(&xk, g, CVXDIM_SOCP, 0);
        if (s)
            fprintf(f, "Q(%s)\n", s);
        else
            fprintf(f, "Q()\n");
        cvxm_printf(f, frm, &xk);
    }

    if (cvxc_mgrp_count(g, CVXDIM_SDP) > 0) {
            cvxc_mgrp_elem(&xk, g, CVXDIM_SDP, 0);
        if (s)
            fprintf(f, "S(%s)\n", s);
        else
            fprintf(f, "S()\n");
        cvxm_printf(f, frm, &xk);
    }
}

void cvxc_mat_printf(FILE *f, const char *format, cvxc_matrix_t *g, const char *s)
{
    cvxc_size_t rows, cols, i, j;
    const char *frm = format ? format : "%9.2e";
    if (s)
        fprintf(f, "%s\n", s);
    cvxm_size(&rows, &cols, g);
    for (i = 0; i < rows; i++) {
        fprintf(f, "[");
        for (j = 0; j < cols; j++) {
            if (j > 0)
                fprintf(f, ",");
            fprintf(f, frm, cvxm_get(g, i, j));
        }
        fprintf(f, "]\n");
    }
    if (cvxm_isepi(g)) {
        fprintf(f, "/");
        fprintf(f, format, g->t);
        fprintf(f, "/\n");
    }
    //cvxm_printf(f, frm, g);
}

void cvxc_scaling_printf(FILE *f, const char *format, cvxc_scaling_t *W, const char *s)
{
    cvxc_matrix_t m;
    char buf[64];
    fprintf(f, "** start scaling %s\n", s);

    if (W->dnlsz > 0) {
        cvxc_scaling_elem(&m, W, CVXWS_DNL, 0);
        cvxc_mat_printf(f, format, &m, "DNL");
        cvxc_scaling_elem(&m, W, CVXWS_DNLI, 0);
        cvxc_mat_printf(f, format, &m, "DNLI");
    }
    if (W->dsz > 0) {
        cvxc_scaling_elem(&m, W, CVXWS_D, 0);
        cvxc_mat_printf(f, format, &m, "D");
        cvxc_scaling_elem(&m, W, CVXWS_DI, 0);
        cvxc_mat_printf(f, format, &m, "DI");
    }
    if (W->vcount > 0) {
        cvxc_scaling_elem(&m, W, CVXWS_BETA, 0);
        cvxc_mat_printf(f, format, &m, "beta");
        for (int i = 0; i < W->vcount; i++) {
            cvxc_scaling_elem(&m, W, CVXWS_V, i);
            sprintf(buf, "V[%d]", i);
            cvxc_mat_printf(f, format, &m, buf);
        }
    }

    if (W->rcount > 0) {
        for (int i = 0; i < W->rcount; i++) {
            cvxc_scaling_elem(&m, W, CVXWS_R, i);
            sprintf(buf, "R[%d]", i);
            cvxc_mat_printf(f, format, &m, buf);

            cvxc_scaling_elem(&m, W, CVXWS_RTI, i);
            sprintf(buf, "RTI[%d]", i);
            cvxc_mat_printf(f, format, &m, buf);
        }
    }
    fprintf(f, "** end scaling %s\n", s);
}

void cvxc_scaling_elem_printf(FILE *f,
                             const char *format,
                             const cvxc_scaling_t *W,
                             cvxc_mset_enum name,
                             int index,
                             const char *s)
{
    cvxc_matrix_t e;
    cvxc_scaling_elem(&e, W, name, index);
    cvxc_mat_printf(f, format, &e, s);
}

static void
cvxm_print_cdata(FILE *f, const cvxc_matrix_t *A)
{
    cvxc_size_t n, i, j, nr, nc;
    char *rowmajor = getenv("ROWMAJOR");

    cvxm_size(&nr, &nc, A);
    n = 0;
    if (rowmajor && tolower(*rowmajor) == 'y') {
        fprintf(f, "/* row major order: [%ld %ld] */\n", nr, nc);
        fprintf(f, "{\n");
        for (i = 0; i < nr; i++) {
            fprintf(f, "/* row:%3ld */ ", i);
            for (j = 0; j < nc; j++) {
                if (n > 0)
                    fprintf(f, ",");
                fprintf(f, "%11.4e", cvxm_get(A, i, j));
                n++;
            }
            fprintf(f, "\n");
        }
    }
    else {
        fprintf(f, "/* column major order: [%ld %ld] */\n", nr, nc);
        fprintf(f, "{\n");
        for (j = 0; j < nc; j++) {
            fprintf(f, "/* col:%3ld */ ", j);
            for (i = 0; i < nr; i++) {
                if (n > 0)
                    fprintf(f, ",");
                fprintf(f, "%11.4e", cvxm_get(A, i, j));
                n++;
            }
            fprintf(f, "\n");
        }
    }
    fprintf(f, "}\n");
}

void cvxc_mat_print_ifenv(const char *name, const cvxc_matrix_t *A, const char *s)
{
    struct stat statbuf;
    char *p = getenv(name);
    char *path;
    if (!p)
        return;

    path = *p == '!' ? &p[1] : p;

    // if file exists don't write unless starts with '!'
    if (*p != '!' && stat(path, &statbuf) == 0)
        return;

    FILE *f = fopen(path, "w");
    if (!f)
        f = stderr;

    if (s)
        fprintf(f, "/* %s */\n", s);
    cvxm_print_cdata(f, A);
    if (f != stderr)
        fclose(f);
}

void cvxc_mat_test_nan(const char *name, const cvxc_matrix_t *A)
{
    cvxc_size_t i, j, nr, nc;
    cvxm_size(&nr, &nc, A);
    for (j = 0; j < nc; j++) {
        for (i = 0; i < nr; i++) {
            if (isnan(cvxm_get(A, i, j))) {
                cvxc_mat_print_ifenv(name, A, "nan-test");
                assert(0);
            }
        }
    }

}
