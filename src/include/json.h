

#ifndef __CONVEX_JSON_H
#define __CONVEX_JSON_H 1

#include "convex.h"

#include <armas/armas.h>

typedef armas_iostream_t cvx_stream_t;

enum json_tokens {
    CVX_JSON_STRING = ARMAS_JSON_STRING,
    CVX_JSON_INT = ARMAS_JSON_INT,
    CVX_JSON_NUMBER = ARMAS_JSON_NUMBER,
    CVX_JSON_NULL = ARMAS_JSON_NULL
};

typedef struct cvx_params
{
    cvx_dimset_t *dims;
    cvx_matrix_t *c;
    cvx_matrix_t *G;
    cvx_matrix_t *h;
    cvx_matrix_t *A;
    cvx_matrix_t *b;
    cvx_solopts_t *opts;
} cvx_params_t;


extern int cvx_json_matrix_read(cvx_matrix_t **A, cvx_stream_t *ios);
extern int cvx_json_matrix_write(cvx_stream_t *ios, const cvx_matrix_t *A);
extern int cvx_json_write_token(cvx_stream_t *ios, int tok, const void *val, size_t len);
extern int cvx_json_write_simple_token(cvx_stream_t *ios, int tok);
extern int cvx_json_read_token(char *iob, size_t len, cvx_stream_t *ios);
extern int cvx_json_dimset_write(cvx_stream_t *ios, const cvx_dimset_t *dims);
extern int cvx_json_dimset_read(cvx_dimset_t **dims, cvx_stream_t *ios);
extern int cvx_json_write_solution(cvx_stream_t *ios, const cvx_solution_t *sol);
extern int cvx_json_write_result(cvx_stream_t *ios, const cvx_solution_t *sol);
extern int cvx_json_read_params(cvx_params_t **pars, cvx_stream_t *ios);
int cvx_json_write_params(cvx_stream_t *ios, const cvx_solopts_t *opts, const cvx_dimset_t *dims,
                          const cvx_matrix_t *c, const cvx_matrix_t *G, const cvx_matrix_t *h,
                          const cvx_matrix_t *A, const cvx_matrix_t *b);
extern int cvx_json_read_options(cvx_solopts_t **opts, cvx_stream_t *ios);
extern int cvx_json_write_options(cvx_stream_t *ios, const cvx_solopts_t *opts, const char *kkt);

#define ONERR(exp) \
    do { if ((exp) < 0) return -1; } while (0)



#endif
