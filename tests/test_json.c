
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include <stdio.h>
#include "cvxm.h"
#include "convex.h"
#include "json.h"

typedef struct iobuf {
    char *data;
    int len;
    char *cp;
} iobuf_t;

int iob_getc(void *ptr)
{
    iobuf_t *b = (iobuf_t *)ptr;
    if (!b->cp)
        b->cp = b->data;
    if (! *b->cp)
        return -1;
    return *b->cp++;
}

void iob_ungetc(void *ptr, int c)
{
    iobuf_t *b = (iobuf_t *)ptr;
    if (b->cp)
        b->cp--;
}

int iob_putc(void *ptr, int c)
{
    iobuf_t *b = (iobuf_t *)ptr;
    if (!b->cp)
        b->cp = b->data;
    if ((b->cp - b->data) < b->len-1) {
        *b->cp++ = c;
        *b->cp = '\0';
        return c;
    }
    return -1;
}

static armas_iostream_vtable_t iob_vt = (armas_iostream_vtable_t){
    .get_char = iob_getc,
    .unget_char = iob_ungetc,
    .put_char = iob_putc
};

void iob_init(iobuf_t *iob, char *data, int len)
{
    iob->data = data;
    iob->len = len;
    iob->cp = iob->data;
}

void iob_stream_init(armas_iostream_t *ios, iobuf_t *iob)
{
    ios->uptr = iob;
    ios->vt = &iob_vt;
}


int file_getc(void *ptr)
{
    FILE *fp = (FILE *)ptr;
    return fgetc(fp);
}

void file_ungetc(void *ptr, int c)
{
    FILE *fp = (FILE *)ptr;
    ungetc(c, fp);
}

int file_putc(void *ptr, int c)
{
    FILE *fp = (FILE *)ptr;
    return fputc(c, fp);
}

static armas_iostream_vtable_t file_vt = (armas_iostream_vtable_t){
    .get_char = file_getc,
    .unget_char = file_ungetc,
    .put_char = file_putc
};

int file_stream_open(armas_iostream_t *ios, const char *path)
{
    FILE *fp;
    if (!(fp = fopen(path, "rw"))) {
        return -1;
    }
    ios->uptr = fp;
    ios->vt = &file_vt;
    return 0;
}

void file_stream_close(armas_iostream_t *ios)
{
    FILE *fp = (FILE *)ios->uptr;
    if (fp)
        fclose(fp);
}

void print_dims(const cvxc_dimset_t *d)
{
    printf("NL: %d\n", d->mnl);
    printf(" L: %d\n", d->ldim);
    printf(" Q: [");
    for (int i = 0; i < d->qlen; i++) {
        if (i > 0)
            printf(", ");
        printf("%d", d->qdims[i]);
    }
    printf("]\n");
    printf(" S: [");
    for (int i = 0; i < d->slen; i++) {
        if (i > 0)
            printf(", ");
        printf("%d", d->sdims[i]);
    }
    printf("]\n");
}

int test_dims(int n)
{
    cvxc_dimset_t dims, *dims2 = (cvxc_dimset_t *)0;
    iobuf_t iob;
    armas_iostream_t ios;
    char buf[512];

    iob_init(&iob, buf, sizeof(buf));
    iob_stream_init(&ios, &iob);
    switch (n) {
    case 0:
        cvxc_dimset_alloc(&dims, 5, (int *)0, (int*)0);
        break;
    case 1:
        cvxc_dimset_alloc(&dims, 5, (int[]){4, 6, 3, 0}, (int*)0);
        break;
    case 2:
        cvxc_dimset_alloc(&dims, 5, (int*)0, (int[]){6, 6, 0});
        break;
    default:
        cvxc_dimset_alloc(&dims, 5, (int[]){4, 6, 3, 0}, (int[]){6, 6, 0});
        break;
    }

    print_dims(&dims);
    cvxc_json_dimset_write(&ios, &dims);

    printf("JSON: %s\n", buf);
    iob.cp = iob.data;
    cvxc_json_dimset_read(&dims2, &ios);
    print_dims(dims2);

    return 0;
}

int test_opts(int n)
{
    iobuf_t iob;
    armas_iostream_t ios;
    char buf[512];
    cvxc_solopts_t *opts2, opts = (cvxc_solopts_t){
        .abstol = 1e-8, .reltol = 1e-9, .feastol = 1e-7, .max_iter = 100
    };

    iob_init(&iob, buf, sizeof(buf));
    iob_stream_init(&ios, &iob);

    opts2 = (cvxc_solopts_t *)0;
    
    cvxc_json_write_options(&ios, &opts, (char *)0);
    printf("JSON: %s\n", buf);
    iob.cp = iob.data;
    cvxc_json_read_options(&opts2, &ios);
    if (opts2) {
        printf("abstol: %e, reltol: %e, feastol: %e, maxiter: %d, kkt:%d\n",
               opts2->abstol, opts2->reltol, opts2->feastol,
               opts2->max_iter, opts2->kkt_solver_name);        
        free(opts2);
    }
    else
        printf("JSON parse error!\n");
    return 1;
}

int test_parms()
{
    iobuf_t iob;
    armas_iostream_t ios;
    char buf[512];

    cvxc_matrix_t c, G, h, A, b;
    cvxc_dimset_t dims;
    cvxc_params_t *pars;

    cvxc_float_t gdata[] = {
        2.0, 1.0, -1.0, 0.0,
        1.0, 2.0, 0.0, -1.0
    };
    cvxc_float_t cdata[] = {-4.0, -5.0};
    cvxc_float_t hdata[] = {3.0, 3.0, 0.0, 0.0};
    cvxc_solopts_t opts = (cvxc_solopts_t){
        .abstol = 1e-6,
        .reltol = 1e-7,
        .feastol = 1e-8,
        .max_iter = 100,
        .debug = 0,
        .refinement = 0,
        .kkt_solver_name = 0,
        .show_progress = 1
    };
    
    iob_init(&iob, buf, sizeof(buf));
    iob_stream_init(&ios, &iob);

    // 
    cvxm_map_data(&c, 2, 1, cdata);
    // inequality constraints; G*x <= h
    cvxm_map_data(&G, 4, 2, gdata);
    cvxm_map_data(&h, 4, 1, hdata);
    // equality constraints: Ax = b  (empty matrices if missing)
    cvxm_map_data(&A, 0, 2, (cvxc_float_t *)0);
    cvxm_map_data(&b, 0, 1, (cvxc_float_t *)0);

    cvxc_dimset_alloc(&dims, 4, (int *)0, (int *)0);

    cvxc_json_write_params(&ios, &opts, &dims, &c, &G, &h, &A, &b);
    printf("JSON: %s\n", buf);
    iob.cp = iob.data;
    memset(&pars, 0, sizeof(pars));
    cvxc_json_read_params(&pars, &ios);
    if (pars) {
        print_dims(pars->dims);
        printf("c\n"); cvxm_printf(stdout, "%9.2e", pars->c);
        printf("G\n"); cvxm_printf(stdout, "%9.2e", pars->G);
        printf("h\n"); cvxm_printf(stdout, "%9.2e", pars->h);
    }
    return 1;
}

int main(int argc, char **argv)
{
    int n = 0;
    if (argc > 1)
        n = atoi(argv[1]);
    //test_dims(n);
    //test_opts(n);
    test_parms();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
