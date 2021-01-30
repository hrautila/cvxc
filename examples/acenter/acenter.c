
#include "cvxc.h"

/*
 * The analytic centering with cone constraints example in section 9.1 of CVXBook.
 */
typedef struct acenter {
    int rows;
    int cols;
} acenter_t;

static  acenter_t Acenter = (acenter_t){ .rows = 3, .cols = 1};

static
int acenter_F(cvxc_matrix_t *f,
              cvxc_matrix_t *Df,
              cvxc_matrix_t *H,
              const cvxc_matrix_t *x,
              const cvxc_matrix_t *z,
              void *user)
{
    cvxc_size_t r, c;
    acenter_t *p = (acenter_t *)user;
    if (!x && !z) {
        cvxm_size(&r, &c, f);
        if ((int)r != p->rows || (int)c != p->cols)
            return -1;
        cvxm_set_all(f, 0.0);
        return 0;
    }
    if (cvxm_amax(x) >= 1.0)
        return -1;

    // compute u := 1 - x**2 ; f = u
    cvxc_float_t z0, xv, u, usum = 0.0;
    if (z) {
        z0 = cvxm_get(z, 0, 0);
        cvxm_set_all(H, 0.0);
    }

    //  f = [1, 1]; Df = [1, 3]; H = [3, 3]
    //
    //  f = - sum( log(1 - x_i^2) )
    // Df = 2*x*u = 2*x/(1 - x**2)
    for (int i = 0; i < p->rows; i++) {
        xv = cvxm_get(x, i, 0);
        u = 1.0 - xv*xv;
        cvxm_set(Df, 0, i, 2.0*xv/u);
        if (z) {
            // H[i,i] =  2*z0* (1 + u**2)/u**2 == 2*z0*(1 + 1/u**2)
            cvxm_set(H, i, i, 2*z0*((1.0 + u*u)/(u*u)));
        }
        usum += log(u);
    }
    if (f)
        cvxm_set(f, 0, 0, -usum);

    return 0;
 }

/*
 * Load convex program
 */
int cp_load(cvxc_convex_program_t *cp, void *data)
{
    cp->F = acenter_F;
    cp->user = &Acenter;
    return 0;
}

void cp_release(cvxc_convex_program_t *cp)
{
}
