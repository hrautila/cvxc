
// Copyright: Harri Rautila, 2018 <harri.rautila@gmail.com>

#include "json.h"
#if 0
#include "convex.h"

#include <armas/armas.h>

typedef armas_iostream_t cvx_stream_t;

enum json_tokens {
    CVX_JSON_STRING = ARMAS_JSON_STRING,
    CVX_JSON_INT = ARMAS_JSON_INT,
    CVX_JSON_NUMBER = ARMAS_JSON_NUMBER,
    CVX_JSON_NULL = ARMAS_JSON_NULL
};
#endif

int cvx_solver_number(const char *name)
{
    if (!name)
        return 0;
    if (strncmp(name, "ldl2", 4) == 0) {
        return CVX_KKTSOL_LDL2;
    } else if (strncmp(name, "chol2", 5) == 0) {
        return CVX_KKTSOL_CHOL2;
    } else if (strncmp(name, "chol", 4) == 0) {
        return CVX_KKTSOL_CHOL;
    } else if (strncmp(name, "qr", 2) == 0) {
        return CVX_KKTSOL_QR;
    } 
    return CVX_KKTSOL_LDL;
}

int cvx_json_write_token(cvx_stream_t *ios, int tok, const void *val, size_t len)
{
    return armas_json_write_token(tok, val, len, ios);
}

int cvx_json_write_simple_token(cvx_stream_t *ios, int tok)
{
    return armas_json_write_token(tok, (void *)0, 0, ios);
}

int cvx_json_read_token(char *iob, size_t len, cvx_stream_t *ios)
{
    return armas_json_read_token(iob, len, ios);
}

int cvx_json_matrix_read(cvx_matrix_t **A, cvx_stream_t *ios)
{
    return armas_d_json_read(A, ios);
}

int cvx_json_matrix_write(cvx_stream_t *ios, const cvx_matrix_t *A)
{
    return armas_d_json_write(ios, A, 0);
}

enum cvx_json_states {
    JSON_STATE_KEY = 100,
    JSON_STATE_SEP,
    JSON_STATE_NLSEP,
    JSON_STATE_NLVAL,
    JSON_STATE_LSEP,
    JSON_STATE_LVAL,
    JSON_STATE_QSEP,
    JSON_STATE_QVAL,
    JSON_STATE_SSEP,
    JSON_STATE_SVAL,
    JSON_STATE_QDIMS_SEP,
    JSON_STATE_QDIMS_VAL,
    JSON_STATE_SDIMS_SEP,
    JSON_STATE_SDIMS_VAL,
    JSON_STATE_HAVE_SIZES,
    JSON_STATE_HAVE_ALL,
    //
    JSON_STATE_OPTS_SEP,
    JSON_STATE_ATOL_SEP,
    JSON_STATE_ATOL_VAL,
    JSON_STATE_RTOL_SEP,
    JSON_STATE_RTOL_VAL,
    JSON_STATE_FTOL_SEP,
    JSON_STATE_FTOL_VAL,
    JSON_STATE_MAXI_SEP,
    JSON_STATE_MAXI_VAL,
    //
    JSON_STATE_KKT_SEP,
    JSON_STATE_KKT_VAL,
    JSON_STATE_DIMS_SEP,
    JSON_STATE_DIMS_VAL,
    JSON_STATE_MATC_SEP,
    JSON_STATE_MATC_VAL,
    JSON_STATE_MATG_SEP,
    JSON_STATE_MATG_VAL,
    JSON_STATE_MATH_SEP,
    JSON_STATE_MATH_VAL,
    JSON_STATE_MATA_SEP,
    JSON_STATE_MATA_VAL,
    JSON_STATE_MATB_SEP,
    JSON_STATE_MATB_VAL
};

enum cvx_json_flags {
    HAVE_SDIMS = 0x1,
    HAVE_QDIMS = 0x2
};

/**
 * @brief Write dimension set into stream.
 *
 * JSON: 
 *    {"nl:int, "l":int , "qn":int, "sn":int, "q":[ints], "s":[ints]}
 * 
 *    nl : non-linear count
 *     l : linear count
 *    qn : SOCP dimension count
 *     q : SOCP dimension sizes
 *    sn : SDP dimension count
 *     s : SDP dimension sizes
 */
int cvx_json_dimset_write(cvx_stream_t *ios, const cvx_dimset_t *dims)
{
    int lval;
    //double lval;
    ONERR(cvx_json_write_simple_token(ios, '{'));
    ONERR(cvx_json_write_token(ios, CVX_JSON_STRING, "nl", 2));
    ONERR(cvx_json_write_simple_token(ios, ':'));
    lval = dims->mnl;
    ONERR(cvx_json_write_token(ios, CVX_JSON_INT, &lval, sizeof(lval)));
    ONERR(cvx_json_write_simple_token(ios, ','));

    ONERR(cvx_json_write_token(ios, CVX_JSON_STRING, "l", 1));
    ONERR(cvx_json_write_simple_token(ios, ':'));
    lval = dims->ldim;
    ONERR(cvx_json_write_token(ios, CVX_JSON_INT, &lval, sizeof(lval)));
    ONERR(cvx_json_write_simple_token(ios, ','));

    ONERR(cvx_json_write_token(ios, CVX_JSON_STRING, "nq", 2));
    ONERR(cvx_json_write_simple_token(ios, ':'));
    lval = dims->qlen;
    ONERR(cvx_json_write_token(ios, CVX_JSON_INT, &lval, sizeof(lval)));
    ONERR(cvx_json_write_simple_token(ios, ','));

    ONERR(cvx_json_write_token(ios, CVX_JSON_STRING, "ns", 2));
    ONERR(cvx_json_write_simple_token(ios, ':'));
    lval = dims->slen;
    ONERR(cvx_json_write_token(ios, CVX_JSON_INT, &lval, sizeof(lval)));

    ONERR(cvx_json_write_simple_token(ios, ','));
    ONERR(cvx_json_write_token(ios, CVX_JSON_STRING, "q", 1));
    ONERR(cvx_json_write_simple_token(ios, ':'));
    ONERR(cvx_json_write_simple_token(ios, '['));
    for (int i = 0; i < dims->qlen; i++) {
        if (i > 0) {
            ONERR(cvx_json_write_simple_token(ios, ','));
        }
        lval = dims->qdims[i];
        ONERR(cvx_json_write_token(ios, CVX_JSON_INT, &lval, sizeof(lval)));
    }
    ONERR(cvx_json_write_simple_token(ios, ']'));

    ONERR(cvx_json_write_simple_token(ios, ','));
    ONERR(cvx_json_write_token(ios, CVX_JSON_STRING, "s", 1));
    ONERR(cvx_json_write_simple_token(ios, ':'));
    ONERR(cvx_json_write_simple_token(ios, '['));
    for (int i = 0; i < dims->slen; i++) {
        if (i > 0) {
            ONERR(cvx_json_write_simple_token(ios, ','));
        }
        lval = dims->sdims[i];
        ONERR(cvx_json_write_token(ios, CVX_JSON_INT, &lval, sizeof(lval)));
    }
    ONERR(cvx_json_write_simple_token(ios, ']'));
    // close json object
    ONERR(cvx_json_write_simple_token(ios, '}'));
    return 0;
}


/**
 * @brief Deserialize dimensiot set structure.
 *
 * JSON: 
 *   {"nl":int, "l":int, "q":int, "s":int, "qdims":[int,..], "sdims":[int,...]}
 */
    
int cvx_json_dimset_read(cvx_dimset_t **dims, cvx_stream_t *ios)
{
    char iob[64];
    int tok, ntok;
    int state = JSON_STATE_KEY;
    int nldim, ldim, qdim, sdim;;
    cvx_dimset_t *dd;

    nldim = ldim = qdim = sdim = -1;
    
    dd = *dims;
    tok = cvx_json_read_token(iob, sizeof(iob), ios);
    if (tok != '{') {
        if (tok == CVX_JSON_NULL) {
            if (dd)
                *dims = (cvx_dimset_t *)0;
            return 0;
        }
        return -1;
    }
    // expect: rows|cols|flags|nnz ':' val
    for (ntok = 0; state != JSON_STATE_HAVE_SIZES ; ntok++) {
        tok = cvx_json_read_token(iob, sizeof(iob), ios);
        if (ntok == 0 && tok == '}') {
            // null dimset if on first token
            if (dd)
                *dims = (cvx_dimset_t *)0;
            return 0;
        }
        switch (state) {
        case JSON_STATE_KEY:
            if (iob[0] == 'l' && iob[1] == '\0') {
                state = JSON_STATE_LSEP;
            } else if (strncmp(iob, "nq", 2) == 0) {
                state = JSON_STATE_QSEP;                
            } else if (strncmp(iob, "ns", 2) == 0) {
                state = JSON_STATE_SSEP;                
            } else if (strncmp(iob, "nl", 2) == 0) {
                state = JSON_STATE_NLSEP;                
            } else {
                return -1;
            }
            break;
        case JSON_STATE_SEP:
            if (tok != ',') return -1;
            state = JSON_STATE_KEY;
            break;
        case JSON_STATE_NLSEP:
            if (tok != ':') return -1;
            state = JSON_STATE_NLVAL;
            break;
        case JSON_STATE_LSEP:
            if (tok != ':') return -1;
            state = JSON_STATE_LVAL;
            break;
        case JSON_STATE_QSEP:
            if (tok != ':') return -1;
            state = JSON_STATE_QVAL;
            break;
        case JSON_STATE_SSEP:
            if (tok != ':') return -1;
            state = JSON_STATE_SVAL;
            break;

        case JSON_STATE_NLVAL:
            if (tok != CVX_JSON_NUMBER)
                return -1;
            nldim = atoi(iob);
            state = JSON_STATE_SEP;
            break;
        case JSON_STATE_LVAL:
            if (tok != CVX_JSON_NUMBER)
                return -1;
            ldim = atoi(iob);
            state = JSON_STATE_SEP;
            break;
        case JSON_STATE_QVAL:
            if (tok != CVX_JSON_NUMBER)
                return -1;
            qdim = atoi(iob);
            state = JSON_STATE_SEP;
            break;
        case JSON_STATE_SVAL:
            if (tok != CVX_JSON_NUMBER)
                return -1;
            sdim = atoi(iob);
            state = JSON_STATE_SEP;
            break;
        }
        // if all sizes then break out from this loop;
        if (nldim != -1 && ldim != -1 && qdim != -1 && sdim != -1) {
            state = JSON_STATE_HAVE_SIZES;
        }
    }

    int have_new = 0;
    if (dd)
        cvx_dimset_create(dd, nldim, ldim, qdim, sdim);
    else {
        have_new = 1;
        dd = cvx_dimset_new(nldim, ldim, qdim, sdim);
        if (!dd)
            return -1;
    }

    // read dimension arrays; TODO: error handling, release reservations on error?
    int read_bits = 0;
    state = JSON_STATE_SEP;
    
    for (; state != JSON_STATE_HAVE_ALL ;) {
        tok = cvx_json_read_token(iob, sizeof(iob), ios);
        switch (state) {
        case JSON_STATE_KEY:
            if (strncmp(iob, "q", 1) == 0) {
                state = JSON_STATE_QDIMS_SEP;
            } else if (strncmp(iob, "s", 1) == 0) {
                state = JSON_STATE_SDIMS_SEP;
            } else {
                goto error_exit;
            }
            break;
            
        case JSON_STATE_SEP:
            if (tok != ',') goto error_exit;
            state = JSON_STATE_KEY;
            break;

        case JSON_STATE_QDIMS_SEP:
            if (tok != ':') goto error_exit;
            state = JSON_STATE_QDIMS_VAL;
            break;
        case JSON_STATE_SDIMS_SEP:
            if (tok != ':') goto error_exit;
            state = JSON_STATE_SDIMS_VAL;
            break;
            
        case JSON_STATE_QDIMS_VAL:
            if (tok != '[')
                goto error_exit;
            for (int i = 0; i < dd->qlen; i++) {
                if (i > 0) {
                    if (cvx_json_read_token(iob, sizeof(iob), ios) != ',')
                        goto error_exit;
                }
                if (cvx_json_read_token(iob, sizeof(iob), ios) != CVX_JSON_NUMBER)
                    goto error_exit;
                dd->qdims[i] = atoi(iob);
            }
            if (cvx_json_read_token(iob, sizeof(iob), ios) != ']')
                goto error_exit;
            read_bits |= HAVE_QDIMS;
            state = JSON_STATE_SEP;
            break;
            
        case JSON_STATE_SDIMS_VAL:
            if (tok != '[')
                goto error_exit;

            for (int i = 0; i < dd->slen; i++) {
                if (i > 0) {
                    if (cvx_json_read_token(iob, sizeof(iob), ios) != ',')
                        goto error_exit;
                }
                if (cvx_json_read_token(iob, sizeof(iob), ios) != CVX_JSON_NUMBER)
                    goto error_exit;
                dd->sdims[i] = atoi(iob);
            }
            if (cvx_json_read_token(iob, sizeof(iob), ios) != ']')
                goto error_exit;
            read_bits |= HAVE_SDIMS;
            state = JSON_STATE_SEP;
            break;
        }
        if ((read_bits & (HAVE_QDIMS|HAVE_SDIMS)) == (HAVE_QDIMS|HAVE_SDIMS)) {
            state = JSON_STATE_HAVE_ALL;
        }               
    }
    // next need to be object closing brace
    if (cvx_json_read_token(iob, sizeof(iob), ios) != '}') {
        goto error_exit;
    }
    if (have_new)
        *dims = dd;
    return 0;
    
 error_exit:
    if (have_new)
        cvx_dimset_free(dd);
    else 
        cvx_dimset_release(dd);
    return -1;
}


/*
 *  JSON: problem parameters
 *   { "opts": { ...opts... },
 *     "dims": { ...dimset... }, 
 *        "c": { ...matrix... },
 *        "G": { ...matrix... },
 *        "h": { ...matrix... },
 *        "A": { ...matrix... },
 *        "b": { ...matrix... } }
 *
 *   {...matrix...}  = {"rows":m, "cols":n, "data":[...numbers...]} | null
 *   [...numbers...] = [] | [array of numbers] (length = m*n)
 *
 *   JSON: result 
 *   { "status": INT,
 *     "result": { ...result... },
 *          "x": { ...matrix... },
 *          "s": { ...matrix... },
 *          "y": { ...matrix... },
 *          "z": { ...matrix... } }
 */

enum cvx_param_bits {
    HAVE_C = 0x1,
    HAVE_G = 0x2,
    HAVE_H = 0x4,
    HAVE_A = 0x8,
    HAVE_B = 0x10,
    HAVE_DIMS = 0x20,
    HAVE_OPTS = 0x40,
    HAVE_ALL = 0x7F
};

int cvx_json_write_options(cvx_stream_t *ios, const cvx_solopts_t *opts, const char *kkt)
{
    ONERR(cvx_json_write_simple_token(ios, '{'));
    // "abstol": NUM
    ONERR(cvx_json_write_token(ios, CVX_JSON_STRING, "abstol", 6));
    ONERR(cvx_json_write_simple_token(ios, ':'));
    ONERR(cvx_json_write_token(ios, CVX_JSON_NUMBER, &opts->abstol, sizeof(opts->abstol)));
    ONERR(cvx_json_write_simple_token(ios, ','));

    // "reltol": NUM
    ONERR(cvx_json_write_token(ios, CVX_JSON_STRING, "reltol", 6));
    ONERR(cvx_json_write_simple_token(ios, ':'));
    ONERR(cvx_json_write_token(ios, CVX_JSON_NUMBER, &opts->reltol, sizeof(opts->reltol)));
    ONERR(cvx_json_write_simple_token(ios, ','));

    // "feastol": NUM
    ONERR(cvx_json_write_token(ios, CVX_JSON_STRING, "feastol", 7));
    ONERR(cvx_json_write_simple_token(ios, ':'));
    ONERR(cvx_json_write_token(ios, CVX_JSON_NUMBER, &opts->feastol, sizeof(opts->feastol)));
    ONERR(cvx_json_write_simple_token(ios, ','));

    // "max_iter": INT
    ONERR(cvx_json_write_token(ios, CVX_JSON_STRING, "maxiter", 7));
    ONERR(cvx_json_write_simple_token(ios, ':'));
    ONERR(cvx_json_write_token(ios, CVX_JSON_INT, &opts->max_iter, sizeof(opts->max_iter)));
    ONERR(cvx_json_write_simple_token(ios, ','));

    // "kkt": INT
    ONERR(cvx_json_write_token(ios, CVX_JSON_STRING, "kkt", 3));
    ONERR(cvx_json_write_simple_token(ios, ':'));
    int kktnum = cvx_solver_number(kkt);
    ONERR(cvx_json_write_token(ios, CVX_JSON_INT, &kktnum, sizeof(kktnum)));

    ONERR(cvx_json_write_simple_token(ios, '}'));
    return 0;
}


/**
 *
 */
int cvx_json_read_options(cvx_solopts_t **opts, cvx_stream_t *ios)
{
    int tok, ntok, state;
    char iob[64];
    cvx_solopts_t *lopts;
    
    if (*opts)
        lopts = *opts;
    else {
        lopts = (cvx_solopts_t *)malloc(sizeof(*lopts));
        if (!lopts)
            return -1;
    }
    memset(lopts, 0, sizeof(*lopts));
    
    tok = cvx_json_read_token(iob, sizeof(iob), ios);
    // first token ...
    switch (tok) {
    case CVX_JSON_NULL:
        *opts = (cvx_solopts_t *)0;
        return 0;
    case '{':
        break;
    default:
        return -1;
    }

    state = JSON_STATE_KEY;
    for (ntok = 0; state != JSON_STATE_HAVE_ALL; ntok++) {
        tok = cvx_json_read_token(iob, sizeof(iob), ios);
        switch (state) {
        case JSON_STATE_KEY:
            if (strncmp(iob, "abstol", 6) == 0) {
                state = JSON_STATE_ATOL_SEP;
            } else if (strncmp(iob, "reltol", 6) == 0) {
                state = JSON_STATE_RTOL_SEP;
            } else if (strncmp(iob, "feastol", 7) == 0) {
                state = JSON_STATE_FTOL_SEP;
            } else if (strncmp(iob, "maxiter", 7) == 0) {
                state = JSON_STATE_MAXI_SEP;
            } else if (strncmp(iob, "kkt", 3) == 0) {
                state = JSON_STATE_KKT_SEP;
            } else {
                return -1;
            }
            break;
        case JSON_STATE_SEP:
            if (tok == '}') {
                state = JSON_STATE_HAVE_ALL;
                break;
            }
            if (tok != ',') goto error;
            state = JSON_STATE_KEY;
            break;

        case JSON_STATE_ATOL_SEP:
            if (tok != ':') goto error;
            state = JSON_STATE_ATOL_VAL;
            break;
        case JSON_STATE_ATOL_VAL:
            if (tok != CVX_JSON_NUMBER) goto error;
            lopts->abstol = strtod(iob, (char **)0);
            state = JSON_STATE_SEP;
            break;

        case JSON_STATE_RTOL_SEP:
            if (tok != ':') goto error;
            state = JSON_STATE_RTOL_VAL;
            break;
        case JSON_STATE_RTOL_VAL:
            if (tok != CVX_JSON_NUMBER) goto error;
            lopts->reltol = strtod(iob, (char **)0);
            state = JSON_STATE_SEP;
            break;

        case JSON_STATE_FTOL_SEP:
            if (tok != ':') goto error;
            state = JSON_STATE_FTOL_VAL;
            break;
        case JSON_STATE_FTOL_VAL:
            if (tok != CVX_JSON_NUMBER) goto error;
            lopts->feastol = strtod(iob, (char **)0);
            state = JSON_STATE_SEP;
            break;

        case JSON_STATE_MAXI_SEP:
            if (tok != ':') goto error;
            state = JSON_STATE_MAXI_VAL;
            break;
        case JSON_STATE_MAXI_VAL:
            if (tok != CVX_JSON_NUMBER) goto error;
            lopts->max_iter = strtol(iob, (char **)0, 0);
            state = JSON_STATE_SEP;
            break;

        case JSON_STATE_KKT_SEP:
            if (tok != ':') goto error;
            state = JSON_STATE_KKT_VAL;
            break;
        case JSON_STATE_KKT_VAL:
            if (tok != CVX_JSON_NUMBER) goto error;
            lopts->kkt_solver_name = cvx_solver_number(iob);
            state = JSON_STATE_SEP;
            break;
        }
    }
    if (!*opts)
        *opts = lopts;
    return 0;
 error:
    if (!*opts)
        free(lopts);
    return -1;
}

int cvx_json_write_params(cvx_stream_t *ios,
                          const cvx_solopts_t *opts,
                          const cvx_dimset_t *dims,
                          const cvx_matrix_t *c,
                          const cvx_matrix_t *G,
                          const cvx_matrix_t *h,
                          const cvx_matrix_t *A,
                          const cvx_matrix_t *b)
{
    ONERR(cvx_json_write_simple_token(ios, '{'));
    // "status": INT
    ONERR(cvx_json_write_token(ios, CVX_JSON_STRING, "dims", 4));
    ONERR(cvx_json_write_simple_token(ios, ':'));
    ONERR(cvx_json_dimset_write(ios, dims));
    ONERR(cvx_json_write_simple_token(ios, ','));
    
    ONERR(cvx_json_write_token(ios, CVX_JSON_STRING, "c", 1));
    ONERR(cvx_json_write_simple_token(ios, ':'));
    ONERR(cvx_json_matrix_write(ios, c));
    ONERR(cvx_json_write_simple_token(ios, ','));

    ONERR(cvx_json_write_token(ios, CVX_JSON_STRING, "G", 1));
    ONERR(cvx_json_write_simple_token(ios, ':'));
    ONERR(cvx_json_matrix_write(ios, G));
    ONERR(cvx_json_write_simple_token(ios, ','));

    ONERR(cvx_json_write_token(ios, CVX_JSON_STRING, "h", 1));
    ONERR(cvx_json_write_simple_token(ios, ':'));
    ONERR(cvx_json_matrix_write(ios, h));
    ONERR(cvx_json_write_simple_token(ios, ','));

    ONERR(cvx_json_write_token(ios, CVX_JSON_STRING, "A", 1));
    ONERR(cvx_json_write_simple_token(ios, ':'));
    ONERR(cvx_json_matrix_write(ios, A));
    ONERR(cvx_json_write_simple_token(ios, ','));

    ONERR(cvx_json_write_token(ios, CVX_JSON_STRING, "b", 1));
    ONERR(cvx_json_write_simple_token(ios, ':'));
    ONERR(cvx_json_matrix_write(ios, b));
    ONERR(cvx_json_write_simple_token(ios, ','));

    ONERR(cvx_json_write_token(ios, CVX_JSON_STRING, "opts", 4));
    ONERR(cvx_json_write_simple_token(ios, ':'));
    ONERR(cvx_json_write_options(ios, opts, (char *)0));

    ONERR(cvx_json_write_simple_token(ios, '}'));
    return 0;
}

int cvx_json_read_params(cvx_params_t **pars, cvx_stream_t *ios)
{
    int tok, ntok, state, err, bits = 0;
    cvx_params_t *pptr;
    char iob[64];
    if (*pars)
        pptr = *pars;
    
    tok = cvx_json_read_token(iob, sizeof(iob), ios);
    if (tok != '{') {
        return -1;
    }
    if (!*pars) {
        pptr = (cvx_params_t *)calloc(1, sizeof(cvx_params_t));
        if (!pptr)
            return -1;
    }
    
    state = JSON_STATE_KEY;
    for (ntok = 0; state != JSON_STATE_HAVE_ALL; ntok++) {
        tok = cvx_json_read_token(iob, sizeof(iob), ios);
        switch (state) {
        case JSON_STATE_KEY:
            if (iob[0] == 'c' && iob[1] == '\0') {
                state = JSON_STATE_MATC_SEP;
            } else if (iob[0] == 'G' && iob[1] == '\0') {
                state = JSON_STATE_MATG_SEP;
            } else if (iob[0] == 'h' && iob[1] == '\0') {
                state = JSON_STATE_MATH_SEP;
            } else if (iob[0] == 'A' && iob[1] == '\0') {
                state = JSON_STATE_MATA_SEP;
            } else if (iob[0] == 'b' && iob[1] == '\0') {
                state = JSON_STATE_MATB_SEP;
            } else if (strncmp(iob, "opts", 4) == 0) {
                state = JSON_STATE_OPTS_SEP;
            } else if (strncmp(iob, "dims", 4) == 0) {
                state = JSON_STATE_DIMS_SEP;
            }
            break;

        case JSON_STATE_SEP:
            if (tok == '}') {
                state = JSON_STATE_HAVE_ALL;
                break;
            }
            if (tok != ',') goto error_exit;
            state = JSON_STATE_KEY;
            break;
            
        case JSON_STATE_MATC_SEP:
            if (tok != ':') goto error_exit;
            if ((err = cvx_json_matrix_read(&pptr->c, ios)) < 0)
                goto error_exit;
            state = JSON_STATE_SEP;
            bits |= HAVE_C;
            break;

        case JSON_STATE_MATG_SEP:
            if (tok != ':') goto error_exit;
            if ((err = cvx_json_matrix_read(&pptr->G, ios)) < 0)
                goto error_exit;
            state = JSON_STATE_SEP;
            bits |= HAVE_G;
            break;

        case JSON_STATE_MATH_SEP:
            if (tok != ':') goto error_exit;
            if ((err = cvx_json_matrix_read(&pptr->h, ios)) < 0)
                goto error_exit;
            state = JSON_STATE_SEP;
            bits |= HAVE_H;
            break;

        case JSON_STATE_MATA_SEP:
            if (tok != ':') goto error_exit;
            if ((err = cvx_json_matrix_read(&pptr->A, ios)) < 0)
                goto error_exit;
            state = JSON_STATE_SEP;
            bits |= HAVE_A;
            break;

        case JSON_STATE_MATB_SEP:
            if (tok != ':') goto error_exit;
            if ((err = cvx_json_matrix_read(&pptr->b, ios)) < 0)
                goto error_exit;
            state = JSON_STATE_SEP;
            bits |= HAVE_B;
            break;

        case JSON_STATE_DIMS_SEP:
            if (tok != ':') goto error_exit;
            if ((err = cvx_json_dimset_read(&pptr->dims, ios)) < 0)
                goto error_exit;
            state = JSON_STATE_SEP;
            bits |= HAVE_DIMS;
            break;

        case JSON_STATE_OPTS_SEP:
            if (tok != ':') goto error_exit;
            if (cvx_json_read_options(&pptr->opts, ios) < 0)
                goto error_exit;
            state = JSON_STATE_SEP;
            bits |= HAVE_OPTS;
            break;
        default:
            fprintf(stderr, "cvx_json_read_params: unknown state %d\n", state);
            goto error_exit;
        }
        if ((bits & HAVE_ALL) == HAVE_ALL) {
            state = JSON_STATE_HAVE_ALL;
        }
    }
    // must get ending '}'
    if (tok != '}') {
        tok = cvx_json_read_token(iob, sizeof(iob), ios);
    }
    if (tok == '}') {
        if (!*pars)
            *pars = pptr;
        return 0;
    }

 error_exit:
    // release all
    if (pptr->c)
        cvxm_free(pptr->c);
    if (pptr->G)
        cvxm_free(pptr->G);
    if (pptr->h)
        cvxm_free(pptr->h);
    if (pptr->A)
        cvxm_free(pptr->A);
    if (pptr->b)
        cvxm_free(pptr->b);
    if (pptr->dims)
        cvx_dimset_free(pptr->dims);
    if (pptr->opts) {
        free(pptr->opts);
    }
    if (! *pars)
        free(pptr);   
    return -1;
}

int cvx_json_write_result(cvx_stream_t *ios, const cvx_solution_t *sol)
{
    // write result fields from solution
    ONERR(cvx_json_write_simple_token(ios, '{'));

    ONERR(cvx_json_write_simple_token(ios, '}'));
    return 0;
}

int cvx_json_write_solution(cvx_stream_t *ios, const cvx_solution_t *sol)
{
    int lval;


    ONERR(cvx_json_write_simple_token(ios, '{'));
    // "status": INT
    ONERR(cvx_json_write_token(ios, CVX_JSON_STRING, "status", 6));
    ONERR(cvx_json_write_simple_token(ios, ':'));
    lval = sol->status;
    ONERR(cvx_json_write_token(ios, CVX_JSON_INT, &lval, sizeof(lval)));

    ONERR(cvx_json_write_simple_token(ios, ','));
    // "result": {result-data}
    ONERR(cvx_json_write_token(ios, CVX_JSON_STRING, "result", 6));
    ONERR(cvx_json_write_simple_token(ios, ':'));
    ONERR(cvx_json_write_result(ios, sol));

    ONERR(cvx_json_write_simple_token(ios, ','));

    int have_result = (sol->status == CVX_STAT_OPTIMAL || sol->status == CVX_STAT_UNKNOWN);

    // "x": {matrix}
    ONERR(cvx_json_write_token(ios, CVX_JSON_STRING, "x", 1));
    ONERR(cvx_json_write_simple_token(ios, ':'));
    if (have_result)
        ONERR(cvx_json_matrix_write(ios, sol->x));
    else
        ONERR(cvx_json_write_simple_token(ios, CVX_JSON_NULL));
        
    ONERR(cvx_json_write_simple_token(ios, ','));
        
    // "s": {matrix}
    ONERR(cvx_json_write_token(ios, CVX_JSON_STRING, "s", 1));
    ONERR(cvx_json_write_simple_token(ios, ':'));
    if (have_result)
        ONERR(cvx_json_matrix_write(ios, sol->s));
    else
        ONERR(cvx_json_write_simple_token(ios, CVX_JSON_NULL));
        
    ONERR(cvx_json_write_simple_token(ios, ','));

    // "y": {matrix}
    ONERR(cvx_json_write_token(ios, CVX_JSON_STRING, "y", 1));
    ONERR(cvx_json_write_simple_token(ios, ':'));
    if (have_result)
        ONERR(cvx_json_matrix_write(ios, sol->y));
    else
        ONERR(cvx_json_write_simple_token(ios, CVX_JSON_NULL));
    
    ONERR(cvx_json_write_simple_token(ios, ','));

    // "z": {matrix}
    ONERR(cvx_json_write_token(ios, CVX_JSON_STRING, "z", 1));
    ONERR(cvx_json_write_simple_token(ios, ':'));
    if (have_result)
        ONERR(cvx_json_matrix_write(ios, sol->z));
    else
        ONERR(cvx_json_write_simple_token(ios, CVX_JSON_NULL));
    
    ONERR(cvx_json_write_simple_token(ios, '}'));

    return 0;
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
