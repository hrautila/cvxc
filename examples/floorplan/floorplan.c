
#include "cvxc.h"

typedef struct floorplan {
    cvxc_matrix_t Amin;
} floorplan_t;

static  cvxc_float_t amin_data[] = {100.0, 100.0, 100.0, 100.0, 100.0};

static
int floorplan_F(cvxc_matrix_t *f,
                cvxc_matrix_t *Df,
                cvxc_matrix_t *H,
                const cvxc_matrix_t *x,
                const cvxc_matrix_t *z,
                void *user)
{
    floorplan_t *p = (floorplan_t *)user;

    cvxc_matrix_t x1, x2, Df1, Df2, H1;
    cvxc_float_t x1val, x2val, aval, xval, zval;

    if (!x && !z) {
        /* F0;
         * assume len(f) is 5 last elements to 1.0
         */
        for (int i = 0; i < 5; i++)
            cvxm_set(f, 21-i, 0, 1.0);
        return 0;
    }

    /* F1 or F2 */
    cvxm_init(&x1, 0, 0);
    cvxm_init(&x2, 0, 0);
    cvxm_init(&Df1, 0, 0);
    cvxm_init(&Df2, 0, 0);

    cvxm_view_map(&x1, x, 12, 0, 5, 1);
    cvxm_view_map(&x2, x, 17, 0, 5, 1);
    cvxm_view_map(&Df1, Df, 0, 12, 5, 5);
    cvxm_view_map(&Df2, Df, 0, 17, 5, 5);

    for (int i = 0; i < 5; i++) {
        x1val = cvxm_get(&x1, i, 0);
        aval = cvxm_get(&p->Amin, i, 0);
        x2val = cvxm_get(&x2, i, 0);
        if (f)
            cvxm_set(f, i, 0, -x1val+aval/x2val);
        cvxm_set(&Df1, i, i, -1.0);
        cvxm_set(&Df2, i, i, -aval/(x2val*x2val));
    }

    if (!z)
        return 0;

    /* F2 */
    cvxm_init(&H1, 0, 0);
    cvxm_view_map(&H1, H, 17, 17, 5, 5);
    for (int i = 0; i < 5; i++) {
        aval = cvxm_get(&p->Amin, i, 0);
        xval = cvxm_get(&x2, i, 0);
        zval = cvxm_get(z, i, 0);
        // H_ii = 2*z_i * Amin_i / x_i^3
        cvxm_set(&H1, i, i, (2.0/xval)*(zval/xval)*(aval/xval));
    }
    return 0;
 }

/*
 * Load convex program
 */
int cp_load(cvxc_convex_program_t *cp, void *data)
{
    floorplan_t *plan = malloc(sizeof(*plan));
    if (!plan) {
        return -1;
    }
    cvxm_map_data(&plan->Amin, 5, 1, amin_data);
    cp->F = floorplan_F;
    cp->user = plan;
    return 0;
}

void cp_release(cvxc_convex_program_t *cp)
{
    if (cp->user)
        free(cp->user);
}
