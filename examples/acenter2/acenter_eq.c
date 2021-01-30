
#include "cvxc.h"

/*
 * The analytic centering with equality constraints example in section 9.1 of CVXBook.
 */
typedef struct acenter {
    int rows;
    int cols;
} acenter_t;

static  acenter_t Acenter_eq = (acenter_t){ .rows = 0, .cols = 0};

static
cvxc_float_t inv(cvxc_float_t v)
{
    return 1.0/v;
}

static
cvxc_float_t inv_square(cvxc_float_t v)
{
    return (1.0/v)*(1.0/v);
}

static
int log_accumulate(cvxc_float_t e, void *p)
{
    cvxc_float_t *accum = (cvxc_float_t *)p;
    if (e <= 0.0)
        return -1;
    *accum += (cvxc_float_t)log((double)e);
    return 0;
}

static int cvxm_logsum(cvxc_float_t *result, const cvxc_matrix_t *x)
{
    cvxc_float_t sum = 0.0;
    if (armas_iterate(&x->data, log_accumulate, &sum, 0) < 0)
        return -1;
    *result = sum;
    return 0;
}

static
int acenter_eq_F(cvxc_matrix_t *f,
                 cvxc_matrix_t *Df,
                 cvxc_matrix_t *H,
                 const cvxc_matrix_t *x,
                 const cvxc_matrix_t *z,
                 void *user)
{
    cvxc_matrix_t d;
    cvxc_float_t logsum;
    if (!x) {
        cvxm_set_all(f, 1.0);
        return 0;
    }

    logsum = 0.0;
    if (cvxm_logsum(&logsum, x) < 0)
        return -1;
    if (f)
        cvxm_set(f, 0, 0, -logsum);

    // Df = 1 / x
    cvxm_copy(Df, x, 0);
    cvxm_apply(Df, inv, 0);
    cvxm_scale(Df, -1.0, 0);

    if (!z)
        return 0;

    // Hessian: diag(z[0] * x^-2)
    cvxm_view_diag(&d, H, 0);
    cvxm_copy(&d, x, 0);
    cvxm_apply(&d, inv_square, 0);
    cvxm_scale(&d, cvxm_get(z, 0, 0), 0);

    return 0;
 }

/*
 * Load convex program
 */
int cp_load(cvxc_convex_program_t *cp, void *data)
{
    cp->F = acenter_eq_F;
    cp->user = &Acenter_eq;
    // This is just an example how parameters can passed.
    // The size of the problem is implicitely defined by the argument matrices.
    if (data)
        sscanf((const char *)data, "m=%d,n=%d", &Acenter_eq.rows, &Acenter_eq.cols);
    return 0;
}

void cp_release(cvxc_convex_program_t *cp)
{
}
