/*
 * Copyright by libcvxc authors. See AUTHORS file in this archive.
 *
 * This file is part of libcvxc library. It is free software,
 * distributed under the terms of GNU Lesser General Public License Version 3, or
 * any later version. See the COPYING file included in this archive.
 */
#define _ISOC99_SOURCE
#include <math.h>
#include "cvxc.h"

#define IOBLEN 128

enum json_tokens {
    CVXC_JSON_STRING = ARMAS_JSON_STRING,
    CVXC_JSON_INT = ARMAS_JSON_INT,
    CVXC_JSON_NUMBER = ARMAS_JSON_NUMBER,
    CVXC_JSON_NULL = ARMAS_JSON_NULL,
    CVXC_JSON_EOF = ARMAS_JSON_EOF,
    CVXC_JSON_TRUE = ARMAS_JSON_TRUE,
    CVXC_JSON_FALSE = ARMAS_JSON_FALSE
};

void cvxc_file_stream(cvxc_stream_t *ios, FILE *fp)
{
    armas_ios_file(ios, fp);
}

void cvxc_str_stream(cvxc_stream_t *ios, const char *s, int len)
{
    armas_ios_string(ios, s, len);
}

int cvxc_solver_number(const char *name)
{
    if (!name)
        return 0;
    if (strncmp(name, "ldl2", 4) == 0) {
        return CVXC_KKTSOL_LDL2;
    } else if (strncmp(name, "chol2", 5) == 0) {
        return CVXC_KKTSOL_CHOL2;
    } else if (strncmp(name, "chol", 4) == 0) {
        return CVXC_KKTSOL_CHOL;
    } else if (strncmp(name, "qr", 2) == 0) {
        return CVXC_KKTSOL_QR;
    }
    return CVXC_KKTSOL_LDL;
}

int cvxc_json_write_token(cvxc_stream_t *ios, int tok, const void *val, size_t len)
{
    return armas_json_write_token(tok, val, len, ios);
}

int cvxc_json_write_simple_token(cvxc_stream_t *ios, int tok)
{
    return armas_json_write_token(tok, (void *)0, 0, ios);
}

int cvxc_json_read_token(char *iob, size_t len, cvxc_stream_t *ios)
{
    return armas_json_read_token(iob, len, ios);
}

void cvxc_json_unget_simple_token(cvxc_stream_t *ios, int tok)
{
    armas_ios_ungetchar(ios, tok);
}

int cvxc_json_matrix_read(cvxc_matrix_t **A, cvxc_stream_t *ios)
{
    if (! *A) {
        cvxc_matrix_t *m = malloc(sizeof(cvxc_matrix_t));
        if (!m)
            return -1;
        m->t = 0.0;
        m->bits = 0;
        *A = m;
    }
    (*A)->t = 0.0;
    (*A)->bits = 0;
    armas_dense_t *amat = &(*A)->data;
    return armas_json_read(&amat, ios);
}

int cvxc_json_matrix_write(cvxc_stream_t *ios, const cvxc_matrix_t *A)
{
    return armas_json_write(ios, &A->data, 0);
}

enum cvxc_json_states {
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
    JSON_STATE_MATB_VAL,
    JSON_STATE_MOD_SEP,
    JSON_STATE_ARGS_SEP
};

enum cvxc_json_flags {
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
int cvxc_json_dimset_write(cvxc_stream_t *ios, const cvxc_dimset_t *dims)
{
    int lval;
    if (!dims) {
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_NULL, 0, 0));
        return 0;
    }
    //double lval;
    ONERR(cvxc_json_write_simple_token(ios, '{'));
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "nl", 2));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    lval = dims->mnl;
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_INT, &lval, sizeof(lval)));
    ONERR(cvxc_json_write_simple_token(ios, ','));

    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "l", 1));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    lval = dims->ldim;
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_INT, &lval, sizeof(lval)));
    ONERR(cvxc_json_write_simple_token(ios, ','));

    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "nq", 2));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    lval = dims->qlen;
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_INT, &lval, sizeof(lval)));
    ONERR(cvxc_json_write_simple_token(ios, ','));

    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "ns", 2));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    lval = dims->slen;
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_INT, &lval, sizeof(lval)));

    ONERR(cvxc_json_write_simple_token(ios, ','));
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "q", 1));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    ONERR(cvxc_json_write_simple_token(ios, '['));
    for (int i = 0; i < dims->qlen; i++) {
        if (i > 0) {
            ONERR(cvxc_json_write_simple_token(ios, ','));
        }
        lval = dims->qdims[i];
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_INT, &lval, sizeof(lval)));
    }
    ONERR(cvxc_json_write_simple_token(ios, ']'));

    ONERR(cvxc_json_write_simple_token(ios, ','));
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "s", 1));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    ONERR(cvxc_json_write_simple_token(ios, '['));
    for (int i = 0; i < dims->slen; i++) {
        if (i > 0) {
            ONERR(cvxc_json_write_simple_token(ios, ','));
        }
        lval = dims->sdims[i];
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_INT, &lval, sizeof(lval)));
    }
    ONERR(cvxc_json_write_simple_token(ios, ']'));
    // close json object
    ONERR(cvxc_json_write_simple_token(ios, '}'));
    return 0;
}


/**
 * @brief Deserialize dimensiot set structure.
 *
 * JSON: 
 *   {"nl":int, "l":int, "q":int, "s":int, "qdims":[int,..], "sdims":[int,...]}
 */
int cvxc_json_dimset_read(cvxc_dimset_t **dims, cvxc_stream_t *ios)
{
    char iob[IOBLEN];
    int tok, ntok;
    int state = JSON_STATE_KEY;
    int nldim, ldim, qdim, sdim;;
    cvxc_dimset_t *dd;

    nldim = ldim = qdim = sdim = -1;

    dd = *dims;
    tok = cvxc_json_read_token(iob, sizeof(iob), ios);
    if (tok != '{') {
        if (tok == CVXC_JSON_NULL) {
            if (dd)
                *dims = (cvxc_dimset_t *)0;
            return 0;
        }
        return -1;
    }
    // assume sizeof attributes l, nl, nq, ns first
    for (ntok = 0; state != JSON_STATE_HAVE_SIZES ; ntok++) {
        tok = cvxc_json_read_token(iob, sizeof(iob), ios);
        if (ntok == 0 && tok == '}') {
            // null dimset if on first token
            if (dd)
                *dims = (cvxc_dimset_t *)0;
            return 0;
        }
        switch (state) {
        case JSON_STATE_KEY:
            if (tok == CVXC_JSON_STRING) {
                if (iob[0] == 'l' && iob[1] == '\0') {
                    state = JSON_STATE_LSEP;
                } else if (strncmp(iob, "nq", 2) == 0) {
                    state = JSON_STATE_QSEP;
                } else if (strncmp(iob, "ns", 2) == 0) {
                    state = JSON_STATE_SSEP;
                } else if (strncmp(iob, "nl", 2) == 0) {
                    state = JSON_STATE_NLSEP;
                } else {
                    if (*iob != 's' && *iob != 'q')
                        return -1;
                    state = JSON_STATE_HAVE_SIZES;
                }
            } else if (tok == '}') {
                state = JSON_STATE_HAVE_SIZES;
                break;
            } else {
                return -1;
            }
            break;
        case JSON_STATE_SEP:
            if (tok != ',' && tok != '}') return -1;
            state = tok == ',' ? JSON_STATE_KEY : JSON_STATE_HAVE_SIZES;
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
            if (tok != CVXC_JSON_NUMBER)
                return -1;
            nldim = atoi(iob);
            state = JSON_STATE_SEP;
            break;
        case JSON_STATE_LVAL:
            if (tok != CVXC_JSON_NUMBER)
                return -1;
            ldim = atoi(iob);
            state = JSON_STATE_SEP;
            break;
        case JSON_STATE_QVAL:
            if (tok != CVXC_JSON_NUMBER)
                return -1;
            qdim = atoi(iob);
            state = JSON_STATE_SEP;
            break;
        case JSON_STATE_SVAL:
            if (tok != CVXC_JSON_NUMBER)
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

    if (nldim == -1)
        nldim = 0;
    if (ldim == -1)
        ldim = 0;
    if (qdim == -1)
        qdim = 0;
    if (sdim == -1)
        sdim = 0;

    int have_new = 0;
    if (dd)
        cvxc_dimset_create(dd, nldim, ldim, qdim, sdim);
    else {
        have_new = 1;
        dd = cvxc_dimset_new(nldim, ldim, qdim, sdim);
        if (!dd)
            return -1;
    }

    // if last token was object ending curly brace, return here
    if (tok == '}') {
        if (have_new)
            *dims = dd;
        return 0;
    }

    // read dimension arrays; TODO: error handling, release reservations on error?
    int read_bits = 0;
    state = JSON_STATE_SEP;
    if (tok == CVXC_JSON_STRING) {
        // last token was either 's' or 'q' key;
        state = *iob == 's' ? JSON_STATE_SDIMS_SEP : JSON_STATE_QDIMS_SEP;
    }

    for (; state != JSON_STATE_HAVE_ALL;) {
        tok = cvxc_json_read_token(iob, sizeof(iob), ios);
        switch (state) {
        case JSON_STATE_KEY:
            if (tok == CVXC_JSON_STRING) {
                if (strncmp(iob, "q", 1) == 0) {
                    state = JSON_STATE_QDIMS_SEP;
                } else if (strncmp(iob, "s", 1) == 0) {
                    state = JSON_STATE_SDIMS_SEP;
                } else {
                    goto error_exit;
                }
            } else {
                goto error_exit;
            }
            break;

        case JSON_STATE_SEP:
            if (tok != ',' && tok != '}') goto error_exit;
            state = tok == ',' ? JSON_STATE_KEY : JSON_STATE_HAVE_ALL;
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
                    if (cvxc_json_read_token(iob, sizeof(iob), ios) != ',')
                        goto error_exit;
                }
                if (cvxc_json_read_token(iob, sizeof(iob), ios) != CVXC_JSON_NUMBER)
                    goto error_exit;
                dd->qdims[i] = atoi(iob);
            }
            if (cvxc_json_read_token(iob, sizeof(iob), ios) != ']')
                goto error_exit;
            read_bits |= HAVE_QDIMS;
            state = JSON_STATE_SEP;
            break;

        case JSON_STATE_SDIMS_VAL:
            if (tok != '[')
                goto error_exit;

            for (int i = 0; i < dd->slen; i++) {
                if (i > 0) {
                    if (cvxc_json_read_token(iob, sizeof(iob), ios) != ',')
                        goto error_exit;
                }
                if (cvxc_json_read_token(iob, sizeof(iob), ios) != CVXC_JSON_NUMBER)
                    goto error_exit;
                dd->sdims[i] = atoi(iob);
            }
            if (cvxc_json_read_token(iob, sizeof(iob), ios) != ']')
                goto error_exit;
            read_bits |= HAVE_SDIMS;
            state = JSON_STATE_SEP;
            break;
        }
        if ((read_bits & (HAVE_QDIMS|HAVE_SDIMS)) == (HAVE_QDIMS|HAVE_SDIMS)) {
            state = JSON_STATE_HAVE_ALL;
        }
    }
    if (tok != '}') {
        // have not seen the ending curly brace; read it from the stream;
        if (cvxc_json_read_token(iob, sizeof(iob), ios) != '}') {
            goto error_exit;
        }
    }
    if (have_new)
        *dims = dd;
    return 0;

 error_exit:
    if (have_new)
        cvxc_dimset_free(dd);
    else
        cvxc_dimset_release(dd);
    return -1;
}


/*
 *  JSON: problem parameters
 *   {
 *     "opts": { ...opts... },
 *     "dims": { ...dimset... },
 *        "c": { ...matrix... },
 *        "G": { ...matrix... },
 *        "h": { ...matrix... },
 *        "A": { ...matrix... },
 *        "b": { ...matrix... },
 *   "module": STRING,
 *     "args": STRING
 *   }
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

enum cvxc_param_bits {
    HAVE_C = 0x1,
    HAVE_G = 0x2,
    HAVE_H = 0x4,
    HAVE_A = 0x8,
    HAVE_B = 0x10,
    HAVE_DIMS = 0x20,
    HAVE_OPTS = 0x40,
    HAVE_ALL = 0x7F
};

int cvxc_json_write_options(cvxc_stream_t *ios, const cvxc_solopts_t *opts, const char *kkt)
{
    if (!opts) {
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_NULL, 0, 0));
        return 0;
    }
    ONERR(cvxc_json_write_simple_token(ios, '{'));
    // "abstol": NUM
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "abstol", 6));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_NUMBER, &opts->abstol, sizeof(opts->abstol)));
    ONERR(cvxc_json_write_simple_token(ios, ','));

    // "reltol": NUM
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "reltol", 6));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_NUMBER, &opts->reltol, sizeof(opts->reltol)));
    ONERR(cvxc_json_write_simple_token(ios, ','));

    // "feastol": NUM
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "feastol", 7));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_NUMBER, &opts->feastol, sizeof(opts->feastol)));
    ONERR(cvxc_json_write_simple_token(ios, ','));

    // "max_iter": INT
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "maxiter", 7));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_INT, &opts->max_iter, sizeof(opts->max_iter)));
    ONERR(cvxc_json_write_simple_token(ios, ','));

    // "kkt": INT
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "kkt", 3));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    int kktnum = cvxc_solver_number(kkt);
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_INT, &kktnum, sizeof(kktnum)));

    ONERR(cvxc_json_write_simple_token(ios, '}'));
    return 0;
}


/**
 * @brief Deserialize convex program options
 */
int cvxc_json_read_options(cvxc_solopts_t **opts, cvxc_stream_t *ios)
{
    int tok, ntok, state;
    char iob[IOBLEN];
    cvxc_solopts_t *lopts;

    if (*opts)
        lopts = *opts;
    else {
        lopts = (cvxc_solopts_t *)malloc(sizeof(*lopts));
        if (!lopts)
            return -1;
    }
    memset(lopts, 0, sizeof(*lopts));

    tok = cvxc_json_read_token(iob, sizeof(iob), ios);
    // first token ...
    switch (tok) {
    case CVXC_JSON_NULL:
        *opts = (cvxc_solopts_t *)0;
        return 0;
    case '{':
        break;
    default:
        return -1;
    }

    state = JSON_STATE_KEY;
    for (ntok = 0; state != JSON_STATE_HAVE_ALL; ntok++) {
        tok = cvxc_json_read_token(iob, sizeof(iob), ios);
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
            if (tok != CVXC_JSON_NUMBER) goto error;
            lopts->abstol = strtod(iob, (char **)0);
            state = JSON_STATE_SEP;
            break;

        case JSON_STATE_RTOL_SEP:
            if (tok != ':') goto error;
            state = JSON_STATE_RTOL_VAL;
            break;
        case JSON_STATE_RTOL_VAL:
            if (tok != CVXC_JSON_NUMBER) goto error;
            lopts->reltol = strtod(iob, (char **)0);
            state = JSON_STATE_SEP;
            break;

        case JSON_STATE_FTOL_SEP:
            if (tok != ':') goto error;
            state = JSON_STATE_FTOL_VAL;
            break;
        case JSON_STATE_FTOL_VAL:
            if (tok != CVXC_JSON_NUMBER) goto error;
            lopts->feastol = strtod(iob, (char **)0);
            state = JSON_STATE_SEP;
            break;

        case JSON_STATE_MAXI_SEP:
            if (tok != ':') goto error;
            state = JSON_STATE_MAXI_VAL;
            break;
        case JSON_STATE_MAXI_VAL:
            if (tok != CVXC_JSON_NUMBER) goto error;
            lopts->max_iter = strtol(iob, (char **)0, 0);
            state = JSON_STATE_SEP;
            break;

        case JSON_STATE_KKT_SEP:
            if (tok != ':') goto error;
            state = JSON_STATE_KKT_VAL;
            break;
        case JSON_STATE_KKT_VAL:
            if (tok != CVXC_JSON_NUMBER) goto error;
            lopts->kkt_solver_name = cvxc_solver_number(iob);
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

int cvxc_json_write_param_items(
    cvxc_stream_t *ios,
    const cvxc_solopts_t *opts,
    const cvxc_dimset_t *dims,
    const cvxc_matrix_t *c,
    const cvxc_matrix_t *G,
    const cvxc_matrix_t *h,
    const cvxc_matrix_t *A,
    const cvxc_matrix_t *b,
    const char *module, const char *args)
{
    ONERR(cvxc_json_write_simple_token(ios, '{'));

    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "dims", 4));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    ONERR(cvxc_json_dimset_write(ios, dims));
    ONERR(cvxc_json_write_simple_token(ios, ','));

    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "c", 1));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    ONERR(cvxc_json_matrix_write(ios, c));
    ONERR(cvxc_json_write_simple_token(ios, ','));

    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "G", 1));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    ONERR(cvxc_json_matrix_write(ios, G));
    ONERR(cvxc_json_write_simple_token(ios, ','));

    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "h", 1));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    ONERR(cvxc_json_matrix_write(ios, h));
    ONERR(cvxc_json_write_simple_token(ios, ','));

    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "A", 1));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    ONERR(cvxc_json_matrix_write(ios, A));
    ONERR(cvxc_json_write_simple_token(ios, ','));

    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "b", 1));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    ONERR(cvxc_json_matrix_write(ios, b));
    ONERR(cvxc_json_write_simple_token(ios, ','));

    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "opts", 4));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    ONERR(cvxc_json_write_options(ios, opts, (char *)0));

    if (module) {
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "module", 6));
        ONERR(cvxc_json_write_simple_token(ios, ':'));
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, module, strlen(args)));
    }
    if (args) {
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "args", 4));
        ONERR(cvxc_json_write_simple_token(ios, ':'));
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, args, strlen(args)));
    }

    ONERR(cvxc_json_write_simple_token(ios, '}'));
    return 0;
}

int cvxc_json_write_params(cvxc_stream_t *ios, cvxc_params_t *pars)
{
    if (!pars) {
        ONERR(cvxc_json_write_simple_token(ios, '{'));
        ONERR(cvxc_json_write_simple_token(ios, '}'));
        return 0;
    }
    return cvxc_json_write_param_items(
        ios, pars->opts, pars->dims, pars->c, pars->G, pars->h,
        pars->A, pars->b, pars->module, pars->args);
}

int cvxc_json_read_params(cvxc_params_t **pars, cvxc_stream_t *ios)
{
    int tok, ntok, state, err, bits = 0;
    cvxc_params_t *pptr = (cvxc_params_t *)0;
    char iob[IOBLEN];
    if (*pars)
        pptr = *pars;

    tok = cvxc_json_read_token(iob, sizeof(iob), ios);
    if (tok != '{') {
        return -1;
    }
    if (!*pars) {
        pptr = (cvxc_params_t *)calloc(1, sizeof(cvxc_params_t));
        if (!pptr)
            return -1;
    } else {
        memset(*pars, 0, sizeof(cvxc_params_t));
    }

    state = JSON_STATE_KEY;
    for (ntok = 0; state != JSON_STATE_HAVE_ALL; ntok++) {
        tok = cvxc_json_read_token(iob, sizeof(iob), ios);
        switch (state) {
        case JSON_STATE_KEY:
            if (tok == CVXC_JSON_STRING) {
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
                } else if (strncmp(iob, "module", 6) == 0) {
                    state = JSON_STATE_MOD_SEP;
                } else if (strncmp(iob, "args", 4) == 0) {
                    state = JSON_STATE_ARGS_SEP;
                }
            } else if (tok == '}' && ntok == 0) {
                state = JSON_STATE_HAVE_ALL;
                break;
            } else {
                goto error_exit;
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
            if ((err = cvxc_json_matrix_read(&pptr->c, ios)) < 0)
                goto error_exit;
            state = JSON_STATE_SEP;
            bits |= HAVE_C;
            break;

        case JSON_STATE_MATG_SEP:
            if (tok != ':') goto error_exit;
            if ((err = cvxc_json_matrix_read(&pptr->G, ios)) < 0)
                goto error_exit;
            state = JSON_STATE_SEP;
            bits |= HAVE_G;
            break;

        case JSON_STATE_MATH_SEP:
            if (tok != ':') goto error_exit;
            if ((err = cvxc_json_matrix_read(&pptr->h, ios)) < 0)
                goto error_exit;
            state = JSON_STATE_SEP;
            bits |= HAVE_H;
            break;

        case JSON_STATE_MATA_SEP:
            if (tok != ':') goto error_exit;
            if ((err = cvxc_json_matrix_read(&pptr->A, ios)) < 0)
                goto error_exit;
            state = JSON_STATE_SEP;
            bits |= HAVE_A;
            break;

        case JSON_STATE_MATB_SEP:
            if (tok != ':') goto error_exit;
            if ((err = cvxc_json_matrix_read(&pptr->b, ios)) < 0)
                goto error_exit;
            state = JSON_STATE_SEP;
            bits |= HAVE_B;
            break;

        case JSON_STATE_DIMS_SEP:
            if (tok != ':') goto error_exit;
            if  (cvxc_json_dimset_read(&pptr->dims, ios) < 0)
                goto error_exit;
            state = JSON_STATE_SEP;
            bits |= HAVE_DIMS;
            break;

        case JSON_STATE_OPTS_SEP:
            if (tok != ':') goto error_exit;
            if (cvxc_json_read_options(&pptr->opts, ios) < 0)
                goto error_exit;
            state = JSON_STATE_SEP;
            bits |= HAVE_OPTS;
            break;

        case JSON_STATE_MOD_SEP:
            if (tok != ':') goto error_exit;
            if ((err = cvxc_json_read_token(iob, sizeof(iob), ios)) != CVXC_JSON_STRING)
                goto error_exit;
            pptr->module = malloc(strlen(iob) + 1);
            memcpy(pptr->module, iob, strlen(iob) + 1);
            state = JSON_STATE_SEP;
            break;

        case JSON_STATE_ARGS_SEP:
            if (tok != ':') goto error_exit;
            if ((err = cvxc_json_read_token(iob, sizeof(iob), ios)) != CVXC_JSON_STRING)
                goto error_exit;
            pptr->args = malloc(strlen(iob) + 1);
            memcpy(pptr->args, iob, strlen(iob) + 1);
            state = JSON_STATE_SEP;
            break;

        default:
            fprintf(stderr, "cvxc_json_read_params: unknown state %d\n", state);
            goto error_exit;
        }
        /* if ((bits & HAVE_ALL) == HAVE_ALL) { */
        /*     state = JSON_STATE_HAVE_ALL; */
        /* } */
    }
    // must get ending '}'
    if (tok != '}') {
        tok = cvxc_json_read_token(iob, sizeof(iob), ios);
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
        cvxc_dimset_free(pptr->dims);
    if (pptr->opts)
        free(pptr->opts);
    if (pptr->module)
        free(pptr->module);
    if (pptr->args)
        free(pptr->args);

    if (!*pars)
        free(pptr);
    return -1;
}

int cvxc_json_write_result(cvxc_stream_t *ios, const cvxc_solution_t *sol)
{
    // write result fields from solution
    ONERR(cvxc_json_write_simple_token(ios, '{'));

    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "primal_objective", 16));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_NUMBER,
                                &sol->primal_objective, sizeof(sol->primal_objective)));
    ONERR(cvxc_json_write_simple_token(ios, ','));

    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "dual_objective", 14));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_NUMBER,
                                &sol->dual_objective, sizeof(sol->dual_objective)));
    ONERR(cvxc_json_write_simple_token(ios, ','));

    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "gap", 3));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_NUMBER, &sol->gap, sizeof(sol->gap)));

    ONERR(cvxc_json_write_simple_token(ios, ','));
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "relative_gap", 12));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_NUMBER,
                                &sol->relative_gap, sizeof(sol->relative_gap)));

    ONERR(cvxc_json_write_simple_token(ios, ','));
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "primal_infeasibility", 19));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_NUMBER,
                                &sol->primal_infeasibility, sizeof(sol->primal_infeasibility)));

    ONERR(cvxc_json_write_simple_token(ios, ','));
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "dual_infeasibility", 17));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_NUMBER,
                                &sol->dual_infeasibility, sizeof(sol->dual_infeasibility)));

    ONERR(cvxc_json_write_simple_token(ios, ','));
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "primal_slack", 12));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_NUMBER,
                                &sol->primal_slack, sizeof(sol->primal_slack)));

    ONERR(cvxc_json_write_simple_token(ios, ','));
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "dual_slack", 10));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_NUMBER,
                                &sol->dual_slack, sizeof(sol->dual_slack)));

    if (isfinite(sol->primal_residual_cert)) {
        ONERR(cvxc_json_write_simple_token(ios, ','));
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "primal_residual_cert", 20));
        ONERR(cvxc_json_write_simple_token(ios, ':'));
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_NUMBER,
                                &sol->primal_residual_cert, sizeof(sol->primal_residual_cert)));
    }

    if (isfinite(sol->dual_residual_cert)) {
        ONERR(cvxc_json_write_simple_token(ios, ','));
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "dual_residual_cert", 18));
        ONERR(cvxc_json_write_simple_token(ios, ':'));
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_NUMBER,
                                &sol->dual_residual_cert, sizeof(sol->dual_residual_cert)));
    }

    ONERR(cvxc_json_write_simple_token(ios, '}'));
    return 0;
}

int cvxc_json_write_solution(cvxc_stream_t *ios, const cvxc_solution_t *sol)
{
    int lval;

    if (!sol) {
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_NULL, 0, 0));
        return 0;
    }

    ONERR(cvxc_json_write_simple_token(ios, '{'));
    // "status": INT
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "status", 6));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    lval = sol->status;
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_INT, &lval, sizeof(lval)));

    ONERR(cvxc_json_write_simple_token(ios, ','));
    // "result": {result-data}
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "result", 6));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    ONERR(cvxc_json_write_result(ios, sol));

    ONERR(cvxc_json_write_simple_token(ios, ','));

    int have_result = (sol->status == CVXC_STAT_OPTIMAL || sol->status == CVXC_STAT_UNKNOWN);

    // "x": {matrix}
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "x", 1));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    if (have_result)
        ONERR(cvxc_json_matrix_write(ios, sol->x));
    else
        ONERR(cvxc_json_write_simple_token(ios, CVXC_JSON_NULL));

    ONERR(cvxc_json_write_simple_token(ios, ','));

    // "s": {matrix}
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "s", 1));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    if (have_result)
        ONERR(cvxc_json_matrix_write(ios, sol->s));
    else
        ONERR(cvxc_json_write_simple_token(ios, CVXC_JSON_NULL));

    ONERR(cvxc_json_write_simple_token(ios, ','));

    // "y": {matrix}
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "y", 1));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    if (have_result)
        ONERR(cvxc_json_matrix_write(ios, sol->y));
    else
        ONERR(cvxc_json_write_simple_token(ios, CVXC_JSON_NULL));

    ONERR(cvxc_json_write_simple_token(ios, ','));

    // "z": {matrix}
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "z", 1));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    if (have_result)
        ONERR(cvxc_json_matrix_write(ios, sol->z));
    else
        ONERR(cvxc_json_write_simple_token(ios, CVXC_JSON_NULL));

    ONERR(cvxc_json_write_simple_token(ios, '}'));

    return 0;
}
