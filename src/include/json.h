

#ifndef __CONVEX_JSON_H
#define __CONVEX_JSON_H 1

#include "cvxc.h"

#include <armas/armas.h>

typedef armas_iostream_t cvxc_stream_t;

enum json_tokens {
    CVXC_JSON_STRING = ARMAS_JSON_STRING,
    CVXC_JSON_INT = ARMAS_JSON_INT,
    CVXC_JSON_NUMBER = ARMAS_JSON_NUMBER,
    CVXC_JSON_NULL = ARMAS_JSON_NULL
};

typedef struct cvxc_params
{
    cvxc_dimset_t *dims;
    cvxc_matrix_t *c;
    cvxc_matrix_t *G;
    cvxc_matrix_t *h;
    cvxc_matrix_t *A;
    cvxc_matrix_t *b;
    cvxc_solopts_t *opts;
} cvxc_params_t;


extern int cvxc_json_matrix_read(cvxc_matrix_t **A, cvxc_stream_t *ios);
extern int cvxc_json_matrix_write(cvxc_stream_t *ios, const cvxc_matrix_t *A);
extern int cvxc_json_write_token(cvxc_stream_t *ios, int tok, const void *val, size_t len);
extern int cvxc_json_write_simple_token(cvxc_stream_t *ios, int tok);
extern int cvxc_json_read_token(char *iob, size_t len, cvxc_stream_t *ios);
extern int cvxc_json_dimset_write(cvxc_stream_t *ios, const cvxc_dimset_t *dims);
extern int cvxc_json_dimset_read(cvxc_dimset_t **dims, cvxc_stream_t *ios);
extern int cvxc_json_write_solution(cvxc_stream_t *ios, const cvxc_solution_t *sol);
extern int cvxc_json_write_result(cvxc_stream_t *ios, const cvxc_solution_t *sol);
extern int cvxc_json_read_params(cvxc_params_t **pars, cvxc_stream_t *ios);
int cvxc_json_write_params(cvxc_stream_t *ios, const cvxc_solopts_t *opts, const cvxc_dimset_t *dims,
                          const cvxc_matrix_t *c, const cvxc_matrix_t *G, const cvxc_matrix_t *h,
                          const cvxc_matrix_t *A, const cvxc_matrix_t *b);
extern int cvxc_json_read_options(cvxc_solopts_t **opts, cvxc_stream_t *ios);
extern int cvxc_json_write_options(cvxc_stream_t *ios, const cvxc_solopts_t *opts, const char *kkt);

#define ONERR(exp) \
    do { if ((exp) < 0) return -1; } while (0)



#endif
