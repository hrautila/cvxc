
// Copyright: Harri Rautila, 2018 <harri.rautila@gmail.com>

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include "cvxm.h"

/*
 * \brief Read matrix from string buffer
 *
 * Contents of string must be in format: '{' nrows ncol '[' elements']' '}'
 * Elements in column major order;
 *
 * \retval 0 OK
 * \retval -1 Parse error
 */
int cvxm_read_sbuffer(cvx_matrix_t *m, const char *s)
{
    const char *sp = s;
    char *endp;
    cvx_size_t nrows, ncols, nelems;
    
    while (*sp && *sp != '{')
        sp++;
    if (*sp == '\0')
        return -1;
    sp++;

    nrows = strtoul(sp, &endp, 0);
    for (sp = endp; *sp && *sp != ',' && !isspace(*sp); sp++)
    if (*sp == '\0')
        return -1;

    ncols = strtoul(sp, &endp, 0);
    for (sp = endp; *sp && *sp != ',' && !isspace(*sp); sp++)
    if (*sp == '\0')
        return -1;

    nelems = nrows*ncols;
    cvxm_init(m, nrows, ncols);
    
    if (nelems == 0)
        return 0;
    
    while (*sp && isspace(*sp) && *sp != '[')
        sp++;
    
    if (*sp == '\0')
        return -1;
    sp++;

    // parse data elements;
    cvx_float_t eval, *ep;
    int k = 0;
    ep = cvxm_data(m, 0);
    
    while (nelems > 0) {
        //printf("sp: %s\n", sp);
        errno = 0;
        eval = strtod(sp, &endp);
        if (eval == 0.0 && errno != 0) {
            // no conversion
            return -1;
        }
        ep[k] = eval;
        k++;
        nelems--;
        for (sp = endp; *sp && *sp != ',' && !isspace(*sp); sp++)
        if (*sp == ']')
            break;

    }
    return 0;
}

#define T_ERROR 256
#define T_INT   257
#define T_FLOAT 258
#define T_STRING 259

static
int get_tok(char *buf, size_t blen, FILE *fp)
{
    int c;
    size_t i;
    int exp_seen = 0, dot_seen = 0;
    do {
        c = fgetc(fp);
    } while (isspace(c) && c != EOF);

    switch (c) {
    case EOF:
    case '{':
    case '}':
    case '[':
    case ']':
    case ',':
    case ':':
        *buf = '\0';
        return c;

    case '\"':          // starts a string
        i = 0;
        c = fgetc(fp);
        while (c != '\"' && i < blen-1) {
            if (iscntrl(c)) {
                // unexpected control character
                return T_ERROR;
            }
            else if (c == '\\') {
                c = fgetc(fp);
                switch (c) {
                case 'u':
                    // unicode escape
                    break;
                case 'x':
                    // hex 
                    break;
                case '0':
                    // octal
                    break;
                default:
                    buf[i] = c;
                    break;
                }
            }
            else {
                buf[i] = c;
            }                
            c = fgetc(fp);
            i++;
        }
        buf[i] = '\0';
        // unfinish string
        return c == '\"' ? T_STRING : T_ERROR;

    default:
        // check if starts a number
        if (! isdigit(c) && c != '-' && c != '.') {
            return T_ERROR;
        }
        // number; int or float
        for (i = 0; i < blen-1; i++) {
            switch (c) {
            case '-':
                if (i != 0 && tolower(buf[i-1]) != 'e')
                    // unexpected '-'
                    return T_ERROR;
                break;
            case 'e':
            case 'E':
                if (exp_seen)
                    // unexpected 'e|E'
                    return T_ERROR;
                exp_seen = 1;
                break;
            case '.':
                if (dot_seen || exp_seen)
                    // unexpected '.'
                    return T_ERROR;
                dot_seen = 0;
                break;
            default:
                if (!isdigit(c)) {
                    goto done;
                }
                break;
            }
            buf[i] = c;
            c = fgetc(fp);
        } 
    done:
        buf[i] = '\0';
        return dot_seen || exp_seen ? T_FLOAT : T_INT;   
    }
    // should not come here at all
    return T_ERROR;
}

int cvxm_read_file(cvx_matrix_t *m, FILE *fp)
{
    int tok;
    double val;
    char buf[64], *endp;
    
    cvx_size_t nrows, ncols;
    
    tok = get_tok(buf, sizeof(buf), fp);
    if (tok != '{')
        return -1;
    
    tok = get_tok(buf, sizeof(buf), fp);
    if (tok != T_INT)
        return -1;
    nrows = strtol(buf, &endp, 0);

    tok = get_tok(buf, sizeof(buf), fp);
    if (tok == ',')
        tok = get_tok(buf, sizeof(buf), fp);
    if (tok != T_INT)
        return -1;
    ncols = strtol(buf, &endp, 0);

    tok = get_tok(buf, sizeof(buf), fp);
    if (tok == ',')
        tok = get_tok(buf, sizeof(buf), fp);
    if (tok != '[')
        return -1;

    cvxm_init(m, nrows, ncols);

    for (size_t i = 0; i < nrows*ncols; i++) {
        tok = get_tok(buf, sizeof(buf), fp);
        if (tok == ',')
            tok = get_tok(buf, sizeof(buf), fp);

        if (tok != T_FLOAT && tok != T_INT)
            goto release;
        val = strtod(buf, &endp);
        m->data.elems[i] = val;
    }
    tok = get_tok(buf, sizeof(buf), fp);
    if (tok != ']')
        goto release;
    tok = get_tok(buf, sizeof(buf), fp);
    if (tok != '}')
        goto release;

    return 0;

 release:
    cvxm_release(m);
    return -1;
}

int cvxm_write_file(FILE *fp, const cvx_matrix_t *m)
{
    int n = 0;
    n += fprintf(fp, "{%d, %d, [", m->data.rows, m->data.cols);
    for (int j = 0; j < m->data.cols; j++) {
        for (int i = 0; i < m->data.rows; i++) {
            if ((j == 0 && i > 0) || j > 0) {
                fputc(',', fp);             
                n++;
            } 
            n += fprintf(fp, "%.9e", cvxm_get((cvx_matrix_t *)m, i, j));
        }
    }
    n += fprintf(fp, "]}\n");
    return n;
}

int cvxm_json_write_file(FILE *fp, const cvx_matrix_t *m)
{
    int n = 0;
    n += fprintf(fp, "{\"rows\":%d, \"cols\": %d, \"data\":[", m->data.rows, m->data.cols);
    for (int j = 0; j < m->data.cols; j++) {
        for (int i = 0; i < m->data.rows; i++) {
            if ((j == 0 && i > 0) || j > 0) {
                fputc(',', fp);             
                n++;
            } 
            n += fprintf(fp, "%.9e", cvxm_get((cvx_matrix_t *)m, i, j));
        }
    }
    n += fprintf(fp, "]}\n");
    return n;
}


// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
