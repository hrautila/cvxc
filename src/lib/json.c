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
    armas_dense_t *amat;
    cvxc_matrix_t *m = 0;
    if (!*A) {
        m = malloc(sizeof(cvxc_matrix_t));
        if (!m)
            return -1;
        amat = m;
    } else {
        amat = *A;
    }
    if (armas_json_read(&amat, ios) < 0) {
        if (!*A)
            free(m);
        return -1;
    }
    if (!*A) {
        // If matrix is null then armas_json_read returns null pointer in matrix
        // argument, release allocated matrix and return null.
        if (!amat) {
            free(m);
        } else {
            *A = m;
        }
    }
    return 0;
}

int cvxc_json_matrix_write(cvxc_stream_t *ios, const cvxc_matrix_t *A)
{
    return armas_json_write(ios, A, 0);
}

enum cvxc_json_states {
    JSON_STATE_KEY = 100,
    JSON_STATE_SEP,
    JSON_STATE_VAL,
    JSON_STATE_KEY_SEP,
    JSON_STATE_KEY_VAL,
};

extern int cvxc_json_intarray_read(cvxc_size_t *array, cvxc_size_t len, cvxc_stream_t *ios);

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
    int keyid = -1;
    int ready = 0;
    for (ntok = 0; !ready ; ntok++) {
        tok = cvxc_json_read_token(iob, sizeof(iob), ios);
        switch (state) {
        case JSON_STATE_KEY:
            if (tok == CVXC_JSON_STRING) {
                state = JSON_STATE_KEY_SEP;
                if (iob[0] == 'l' && iob[1] == '\0') {
                    keyid = 0;
                } else if (strncmp(iob, "nq", 2) == 0) {
                    keyid = 1;
                } else if (strncmp(iob, "ns", 2) == 0) {
                    keyid = 2;
                } else if (strncmp(iob, "nl", 2) == 0) {
                    keyid = 3;
                } else {
                    if (*iob != 's' && *iob != 'q')
                        return -1;
                    ready = 1;
                }
            } else if (tok == '}') {
                ready = 1;
                break;
            } else {
                return -1;
            }
            break;

        case JSON_STATE_SEP:
            if (tok != ',' && tok != '}') return -1;
            state = JSON_STATE_KEY;
            if (tok == '}' || (nldim != -1 && ldim != -1 && qdim != -1 && sdim != -1)) {
                // if seen all sizes or end-of-object brace break out from this loop;
                ready = 1;
            }
            break;

        case JSON_STATE_KEY_SEP:
            if (tok != ':') return -1;
            tok = cvxc_json_read_token(iob, sizeof(iob), ios);
            if (tok != CVXC_JSON_NUMBER) {
                fprintf(stderr, "json(dimset): unexpected token: %d\n", tok);
                return -1;
            }
            state = JSON_STATE_SEP;
            switch (keyid) {
            case 0:
                ldim = atoi(iob);
                break;
            case 1:
                qdim = atoi(iob);
                break;
            case 2:
                sdim = atoi(iob);
                break;
            case 3:
                nldim = atoi(iob);
                break;
            }
            break;
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

    // here last token was ',' or STRING matching 's' or 'q' arrays.
    state = JSON_STATE_KEY;
    keyid = -1;
    if (tok == CVXC_JSON_STRING) {
        // last token was either 's' or 'q' key, we expect ':' token
        keyid = *iob == 'q' ? 0 : (*iob == 's' ? 1 : -1);
        state = JSON_STATE_KEY_SEP;
    }

    for (ready = 0; !ready;) {
        tok = cvxc_json_read_token(iob, sizeof(iob), ios);
        switch (state) {
        case JSON_STATE_KEY:
            if (tok == CVXC_JSON_STRING) {
                if (strncmp(iob, "q", 1) == 0) {
                    keyid = 0;
                } else if (strncmp(iob, "s", 1) == 0) {
                    keyid = 1;
                } else {
                    goto error_exit;
                }
                state = JSON_STATE_KEY_SEP;
            } else {
                goto error_exit;
            }
            break;

        case JSON_STATE_SEP:
            if (tok != ',' && tok != '}') goto error_exit;
            state = JSON_STATE_KEY;
            ready = tok == '}';
            break;

        case JSON_STATE_KEY_SEP:
            if (tok != ':') goto error_exit;
            state = JSON_STATE_SEP;
            switch (keyid) {
            case 0:
                if (cvxc_json_intarray_read(dd->qdims, dd->qlen, ios) < 0) {
                    goto error_exit;
                }
                break;
            case 1:
                if (cvxc_json_intarray_read(dd->sdims, dd->slen, ios) < 0) {
                    goto error_exit;
                }
                break;
            default:
                fprintf(stderr, "json(dimset) unexpected key value: %s\n", iob );
                goto error_exit;
            }
            break;
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

int cvxc_json_intarray_read(cvxc_size_t *array, cvxc_size_t len, cvxc_stream_t *ios)
{
    char iob[IOBLEN];
    int tok, ntok, cnt = 0;
    cvxc_size_t value;
    int state = JSON_STATE_VAL;

    tok = cvxc_json_read_token(iob, sizeof(iob), ios);
    if (tok != '[') {
        fprintf(stderr, "json(intarray): unexpected array start: %d\n", tok);
        return -1;
    }

    int ready = 0;
    for (ntok = 0; !ready; ntok++) {
        tok = cvxc_json_read_token(iob, sizeof(iob), ios);
        switch (state) {
        case JSON_STATE_VAL:
            if (ntok == 0 && tok == ']')
                return 0;
            if (tok == CVXC_JSON_NUMBER) {
                value = strtol(iob, (char **)0, 10);
                if (cnt < len && array)
                    array[cnt++] = value;
                state = JSON_STATE_SEP;
            } else {
                fprintf(stderr, "json(intarray): unexpected token: %d\n", tok);
                return -1;
            }
            break;

        case JSON_STATE_SEP:
            state = JSON_STATE_VAL;
            if (tok == ']') {
                ready = 1;
            } else  if (tok != ',') {
                fprintf(stderr, "json(intarray): unexpected token: %d\n", tok);
                return -1;
            }
            break;

        default:
            return -1;
        }
    }
    return cnt;
}

int cvxc_json_intarray_write(cvxc_stream_t *ios, const cvxc_size_t *array, cvxc_size_t len)
{
    ONERR(cvxc_json_write_simple_token(ios, '['));
    for (int k = 0; k < len; k++) {
        if (k > 0)
            ONERR(cvxc_json_write_simple_token(ios, ','));
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_NUMBER, &array[k], sizeof(*array)));
    }
    ONERR(cvxc_json_write_simple_token(ios, ']'));
    return 0;
}

/**
 * @brief Deserialize GP problem index set.
 *
 * JSON:
 *   {"len":int, "index":[int,..]}
 */
int cvxc_json_gpindex_read(cvxc_gpindex_t **gpi, cvxc_stream_t *ios)
{
    char iob[IOBLEN];
    int tok, ntok, ready, keyid;
    int state = JSON_STATE_KEY;
    cvxc_gpindex_t *lgpi;

    tok = cvxc_json_read_token(iob, sizeof(iob), ios);
    // first token ...
    switch (tok) {
    case CVXC_JSON_NULL:
        if (*gpi)
            *gpi = (cvxc_gpindex_t *)0;
        return 0;
    case '{':
        break;
    default:
        return -1;
    }

    if (*gpi)
       lgpi = *gpi;
    else {
        lgpi = (cvxc_gpindex_t *)malloc(sizeof(*lgpi));
        if (!lgpi)
            return -1;
    }
    memset(lgpi, 0, sizeof(*lgpi));

    state = JSON_STATE_KEY;
    keyid = -1;
    ready = 0;
    for (ntok = 0; !ready; ntok++) {
        tok = cvxc_json_read_token(iob, sizeof(iob), ios);
        switch (state) {
        case JSON_STATE_KEY:
            if (tok == CVXC_JSON_STRING) {
                if (strncmp(iob, "len", 3) == 0) {
                    keyid = 0;
                } else if (strncmp(iob, "sizes", 5) == 0) {
                    keyid = 1;
                } else {
                    fprintf(stderr, "json(gpindex): unexpected key: %s\n", iob);
                    goto error_exit;
                }
                state = JSON_STATE_KEY_SEP;
            } else if (ntok == 0 && tok == '}') {
                // empty structure
                if (*gpi)
                    *gpi = lgpi;
                return 0;
            } else {
                fprintf(stderr, "json(gpindex): unexpected token: %d\n", tok);
                goto error_exit;
            }
            break;

        case JSON_STATE_SEP:
            if (tok == '}') {
                ready = 1;
                break;
            }
            if (tok != ',') {
                fprintf(stderr, "json(gpindex): unexpected token: %d\n", tok);
                goto error_exit;
            }
            state = JSON_STATE_KEY;
            break;

        case JSON_STATE_KEY_SEP:
            if (tok != ':') {
                fprintf(stderr, "json(gpindex): unexpected token: %d\n", tok);
                goto error_exit;
            }
            switch (keyid) {
            case 0:
                tok = cvxc_json_read_token(iob, sizeof(iob), ios);
                if (tok != CVXC_JSON_NUMBER) goto error_exit;
                lgpi->p = strtol(iob, (char **)0, 10);
                break;
            case 1:
                if (!(lgpi->__bytes = malloc(lgpi->p * sizeof(*lgpi->index))))
                    goto error_exit;
                lgpi->index = (cvxc_size_t *)lgpi->__bytes;
                if (cvxc_json_intarray_read(lgpi->index, lgpi->p, ios) < 0) {
                    free(lgpi->__bytes);
                    goto error_exit;
                }
                break;
            default:
                goto error_exit;
            }
            state = JSON_STATE_SEP;
            break;

        default:
            goto error_exit;
        }
    }
    if (!*gpi)
        *gpi = lgpi;
    return 0;
error_exit:
    if (!*gpi)
        free(lgpi);
    return -1;
}

int cvxc_json_gpindex_write(cvxc_stream_t *ios, const cvxc_gpindex_t *gpi)
{
    if (!gpi) {
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_NULL, 0, 0));
        return 0;
    }
    ONERR(cvxc_json_write_simple_token(ios, '{'));
    // "abstol": NUM
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "len", 3));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    ONERR(cvxc_json_write_token(ios, CVXC_JSON_NUMBER, &gpi->p, sizeof(gpi->p)));
    ONERR(cvxc_json_write_simple_token(ios, ','));

    ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "sizes", 5));
    ONERR(cvxc_json_write_simple_token(ios, ':'));
    ONERR(cvxc_json_intarray_write(ios, gpi->index, gpi->p));

    ONERR(cvxc_json_write_simple_token(ios, '}'));
    return 0;
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
 *        "F": { ...matrix... },
 *        "K": { "length": int, "data": [] },
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
    cvxc_solopts_t *lopts = 0;

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

    if (*opts)
        lopts = *opts;
    else {
        lopts = (cvxc_solopts_t *)malloc(sizeof(*lopts));
        if (!lopts)
            return -1;
    }
    memset(lopts, 0, sizeof(*lopts));

    int keyid = -1;
    int ready = 0;
    state = JSON_STATE_KEY;
    for (ntok = 0; !ready; ntok++) {
        tok = cvxc_json_read_token(iob, sizeof(iob), ios);
        switch (state) {
        case JSON_STATE_KEY:
            if (strncmp(iob, "abstol", 6) == 0) {
                keyid = 0;
            } else if (strncmp(iob, "reltol", 6) == 0) {
                keyid = 1;
            } else if (strncmp(iob, "feastol", 7) == 0) {
                keyid = 2;
            } else if (strncmp(iob, "maxiter", 7) == 0) {
                keyid = 3;
            } else if (strncmp(iob, "kkt", 3) == 0) {
                keyid = 4;
            } else {
                goto error;
            }
            state = JSON_STATE_KEY_SEP;
            break;

        case JSON_STATE_SEP:
            if (tok == '}') {
                ready = 1;
                break;
            }
            if (tok != ',') {
                fprintf(stderr, "json(options): unexpected token: %d\n", tok);
                goto error;
            }
            state = JSON_STATE_KEY;
            break;

        case JSON_STATE_KEY_SEP:
            if (tok == ':') {
                fprintf(stderr, "json(options): unexpected token: %d\n", tok);
                goto error;
            }
            tok = cvxc_json_read_token(iob, sizeof(iob), ios);
            if (!(tok == CVXC_JSON_NUMBER || (keyid == 4 && tok == CVXC_JSON_STRING))) {
                fprintf(stderr, "json(options): unexpected token: %d\n", tok);
            }
            switch (keyid) {
            case 0:
                lopts->abstol = strtod(iob, (char **)0);
                break;
            case 1:
                lopts->reltol = strtod(iob, (char **)0);
                break;
            case 2:
                lopts->feastol = strtod(iob, (char **)0);
                break;
            case 3:
                lopts->max_iter = strtol(iob, (char **)0, 0);
                break;
            case 4:
                lopts->kkt_solver_name = cvxc_solver_number(iob);
                break;
            }
            state = JSON_STATE_KEY;
            break;

        default:
            goto error;
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

int cvxc_json_write_params(cvxc_stream_t *ios, cvxc_params_t *pars)
{
    int nelem = 0;
    if (!pars) {
        ONERR(cvxc_json_write_simple_token(ios, '{'));
        ONERR(cvxc_json_write_simple_token(ios, '}'));
        return 0;
    }

    ONERR(cvxc_json_write_simple_token(ios, '{'));

    if (pars->dims) {
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "dims", 4));
        ONERR(cvxc_json_write_simple_token(ios, ':'));
        ONERR(cvxc_json_dimset_write(ios, pars->dims));
        nelem++;
    }
    if (pars->c) {
        if (nelem)
            ONERR(cvxc_json_write_simple_token(ios, ','));
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "c", 1));
        ONERR(cvxc_json_write_simple_token(ios, ':'));
        ONERR(cvxc_json_matrix_write(ios, pars->c));
        ONERR(cvxc_json_write_simple_token(ios, ','));
        nelem++;
    }

    if (pars->G) {
        if (nelem)
            ONERR(cvxc_json_write_simple_token(ios, ','));
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "G", 1));
        ONERR(cvxc_json_write_simple_token(ios, ':'));
        ONERR(cvxc_json_matrix_write(ios, pars->G));
        nelem++;
    }

    if (pars->h) {
        if (nelem)
            ONERR(cvxc_json_write_simple_token(ios, ','));
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "h", 1));
        ONERR(cvxc_json_write_simple_token(ios, ':'));
        ONERR(cvxc_json_matrix_write(ios, pars->h));
        nelem++;
    }

    if (pars->A) {
        if (nelem)
            ONERR(cvxc_json_write_simple_token(ios, ','));
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "A", 1));
        ONERR(cvxc_json_write_simple_token(ios, ':'));
        ONERR(cvxc_json_matrix_write(ios, pars->A));
        nelem++;
    }

    if (pars->b) {
        if (nelem)
            ONERR(cvxc_json_write_simple_token(ios, ','));
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "b", 1));
        ONERR(cvxc_json_write_simple_token(ios, ':'));
        ONERR(cvxc_json_matrix_write(ios, pars->b));
        nelem++;
    }

    if (pars->F) {
        if (nelem)
            ONERR(cvxc_json_write_simple_token(ios, ','));
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "F", 1));
        ONERR(cvxc_json_write_simple_token(ios, ':'));
        ONERR(cvxc_json_matrix_write(ios, pars->F));
        nelem++;
    }

    if (pars->g) {
        if (nelem)
            ONERR(cvxc_json_write_simple_token(ios, ','));
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "g", 1));
        ONERR(cvxc_json_write_simple_token(ios, ':'));
        ONERR(cvxc_json_matrix_write(ios, pars->g));
        nelem++;
    }

    if (pars->K) {
        if (nelem)
            ONERR(cvxc_json_write_simple_token(ios, ','));
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "K", 1));
        ONERR(cvxc_json_write_simple_token(ios, ':'));
        ONERR(cvxc_json_gpindex_write(ios, pars->K));
        nelem++;
    }

    if (pars->opts) {
        if (nelem)
            ONERR(cvxc_json_write_simple_token(ios, ','));
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "opts", 4));
        ONERR(cvxc_json_write_simple_token(ios, ':'));
        ONERR(cvxc_json_write_options(ios, pars->opts, (char *)0));
        nelem++;
    }

    if (pars->module) {
        if (nelem)
            ONERR(cvxc_json_write_simple_token(ios, ','));
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "module", 6));
        ONERR(cvxc_json_write_simple_token(ios, ':'));
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, pars->module, strlen(pars->module)));
        nelem++;
    }
    if (pars->args) {
        if (nelem)
            ONERR(cvxc_json_write_simple_token(ios, ','));
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, "args", 4));
        ONERR(cvxc_json_write_simple_token(ios, ':'));
        ONERR(cvxc_json_write_token(ios, CVXC_JSON_STRING, pars->args, strlen(pars->args)));
        nelem++;
    }

    ONERR(cvxc_json_write_simple_token(ios, '}'));
    return 0;
}

int cvxc_json_read_params(cvxc_params_t **pars, cvxc_stream_t *ios)
{
    int tok, ntok, state, err;
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
    int keyid = -1;
    int ready = 0;
    for (ntok = 0; !ready; ntok++) {
        tok = cvxc_json_read_token(iob, sizeof(iob), ios);
        switch (state) {
        case JSON_STATE_KEY:
            if (tok == CVXC_JSON_STRING) {
                if (iob[0] == 'c' && iob[1] == '\0') {
                    keyid = 0;
                } else if (iob[0] == 'G' && iob[1] == '\0') {
                    keyid = 1;
                } else if (iob[0] == 'h' && iob[1] == '\0') {
                    keyid = 2;
                } else if (iob[0] == 'A' && iob[1] == '\0') {
                    keyid = 3;
                } else if (iob[0] == 'b' && iob[1] == '\0') {
                    keyid = 4;
                } else if (iob[0] == 'F' && iob[1] == '\0') {
                    keyid = 5;
                } else if (iob[0] == 'K' && iob[1] == '\0') {
                    keyid = 6;
                } else if (strncmp(iob, "opts", 4) == 0) {
                    keyid = 7;
                } else if (strncmp(iob, "dims", 4) == 0) {
                    keyid = 8;
                } else if (strncmp(iob, "module", 6) == 0) {
                    keyid = 9;
                } else if (strncmp(iob, "args", 4) == 0) {
                    keyid = 10;
                } else if (iob[0] == 'g' && iob[1] == '\0') {
                    keyid = 11;
                } else {
                    keyid = -1;
                }
                state = JSON_STATE_KEY_SEP;
            } else if (tok == '}' && ntok == 0) {
                ready = 1;
                break;
            } else {
                goto error_exit;
            }
            break;

        case JSON_STATE_SEP:
            if (tok == '}') {
                ready = 1;
                break;
            }
            if (tok != ',') goto error_exit;
            state = JSON_STATE_KEY;
            break;

        case JSON_STATE_KEY_SEP:
            if (tok != ':') goto error_exit;
            switch (keyid) {
            case 0:
                if ((err = cvxc_json_matrix_read(&pptr->c, ios)) < 0)
                    goto error_exit;
                break;
            case 1:
                if ((err = cvxc_json_matrix_read(&pptr->G, ios)) < 0)
                    goto error_exit;
                break;
            case 2:
                if ((err = cvxc_json_matrix_read(&pptr->h, ios)) < 0)
                    goto error_exit;
                break;
            case 3:
                if ((err = cvxc_json_matrix_read(&pptr->A, ios)) < 0)
                    goto error_exit;
                break;
            case 4:
                if ((err = cvxc_json_matrix_read(&pptr->b, ios)) < 0)
                    goto error_exit;
                break;
            case 5:
                if ((err = cvxc_json_matrix_read(&pptr->F, ios)) < 0)
                    goto error_exit;
                break;
            case 6:
                if ((err = cvxc_json_gpindex_read(&pptr->K, ios)) < 0)
                    goto error_exit;
                break;
            case 7:
                if (cvxc_json_read_options(&pptr->opts, ios) < 0)
                    goto error_exit;
                break;
            case 8:
                if  (cvxc_json_dimset_read(&pptr->dims, ios) < 0)
                    goto error_exit;
                break;
            case 9:
                if ((err = cvxc_json_read_token(iob, sizeof(iob), ios)) != CVXC_JSON_STRING)
                    goto error_exit;
                pptr->module = malloc(strlen(iob) + 1);
                memcpy(pptr->module, iob, strlen(iob) + 1);
                break;
            case 10:
                if ((err = cvxc_json_read_token(iob, sizeof(iob), ios)) != CVXC_JSON_STRING)
                    goto error_exit;
                pptr->args = malloc(strlen(iob) + 1);
                memcpy(pptr->args, iob, strlen(iob) + 1);
                break;
            case 11:
                if ((err = cvxc_json_matrix_read(&pptr->g, ios)) < 0)
                    goto error_exit;
                break;
            default:
                fprintf(stderr, "unexpected parameter name: %s", iob);
                break;
            }
            state = JSON_STATE_SEP;
            break;

        default:
            fprintf(stderr, "cvxc_json_read_params: unknown state %d\n", state);
            goto error_exit;
        }
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
    if (pptr->F)
        cvxm_free(pptr->F);
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

void cvxc_params_release(cvxc_params_t *pars)
{
    if (pars->c)
        cvxm_free(pars->c);
    if (pars->G)
        cvxm_free(pars->G);
    if (pars->h)
        cvxm_free(pars->h);
    if (pars->A)
        cvxm_free(pars->A);
    if (pars->b)
        cvxm_free(pars->b);
    if (pars->F)
        cvxm_free(pars->F);
    if (pars->dims)
        cvxc_dimset_free(pars->dims);
    if (pars->opts)
        free(pars->opts);
    if (pars->module)
        free(pars->module);
    if (pars->args)
        free(pars->args);
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
