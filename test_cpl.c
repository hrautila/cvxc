
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include <unistd.h>
#include "convex.h"


/*
 *     minimize    W+H
 *     subject to  Amin1 / h1 <= w1
 *                 Amin2 / h2 <= w2
 *                 Amin3 / h3 <= w3
 *                 Amin4 / h4 <= w4
 *                 Amin5 / h5 <= w5
 *                 x1 >= 0
 *                 x2 >= 0
 *                 x4 >= 0
 *                 x1 + w1 + rho <= x3
 *                 x2 + w2 + rho <= x3
 *                 x3 + w3 + rho <= x5
 *                 x4 + w4 + rho <= x5
 *                 x5 + w5 <= W
 *                 y2 >= 0
 *                 y3 >= 0
 *                 y5 >= 0
 *                 y2 + h2 + rho <= y1
 *                 y1 + h1 + rho <= y4
 *                 y3 + h3 + rho <= y4
 *                 y4 + h4 <= H
 *                 y5 + h5 <= H
 *                 h1/gamma <= w1 <= gamma*h1
 *                 h2/gamma <= w2 <= gamma*h2
 *                 h3/gamma <= w3 <= gamma*h3
 *                 h4/gamma <= w4 <= gamma*h4
 *                 h5/gamma <= w5 <= gamma*h5
 *
 * 22 Variables W, H, x (5), y (5), w (5), h (5).
 *
 * W, H:  scalars; bounding box width and height
 * x, y:  5-vectors; coordinates of bottom left corners of blocks
 * w, h:  5-vectors; widths and heigths of the 5 blocks
 *
 * Objective: 
 *      min W+H
 *
 * Five nonlinear constraints:
 *
 *     -w1 + Amin1 / h1 <= 0
 *     -w2 + Amin2 / h2 <= 0
 *     -w3 + Amin3 / h3 <= 0
 *     -w4 + Amin4 / h4 <= 0
 *     -w5 + Amin5 / h5 <= 0
 */

typedef struct floorplan {
    cvx_matrix_t Amin;
} floorplan_t;


int floorplan_init(cvx_convex_program_t *cp, void *data)
{
    return 0;
}

void floorplan_release(cvx_convex_program_t *cp)
{
}

int floorplan_F(cvx_matrix_t *f,
                cvx_matrix_t *Df,
                cvx_matrix_t *H,
                const cvx_matrix_t *x,
                const cvx_matrix_t *z,
                void *user)
{
    floorplan_t *p = (floorplan_t *)user;
    
    cvx_matrix_t x1, x2, Df1, Df2, H1;
    cvx_float_t x1val, x2val, aval, xval, zval;
    
    if (!x && !z) {
        // F0; assume len(f) is 
        // 5 last elements to 1.0
        for (int i = 0; i < 5; i++)
            cvxm_set(f, 21-i, 0, 1.0);
        return 0;
    }
    // F1
    
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
        
    cvxm_view_map(&H1, H, 17, 17, 5, 5);
    for (int i = 0; i < 5; i++) {
        aval = cvxm_get(&p->Amin, i, 0);
        xval = cvxm_get(&x2, i, 0);
        zval = cvxm_get(z, i, 0);
        // H_ii = 2*z_i * Amin_i / x_i^3
        cvxm_set(&H1, i, i, (2.0/xval)*(zval/xval)*(aval/xval));
    }
    // cvx_mat_printf(stdout, "%6.3f", H, "F2:H");
    return 0;
 }

void floorplan_constraints(cvx_matrix_t *c,
                           cvx_matrix_t *G,
                           cvx_matrix_t *h,
                           cvx_float_t rho,
                           cvx_float_t gamma)
{
    // [W, H, x(5), y(5), w(5), h(5)]
    cvxm_set(c, 0, 0, 1.0);
    cvxm_set(c, 1, 0, 1.0);

    // -x1 <= 0; -x2 <= 0; -x4 <= 0
    cvxm_set(G, 0, 2, -1.0);
    cvxm_set(G, 1, 3, -1.0);
    cvxm_set(G, 2, 5, -1.0);

    // x1 - x3 + w1 <= -rho
    cvxm_set(G, 3, 2,  1.0);
    cvxm_set(G, 3, 4, -1.0);
    cvxm_set(G, 3, 12, 1.0);
    cvxm_set(h, 3, 0, -rho);

    // x2 - x3 + w2 <= -rho
    cvxm_set(G, 4, 3,  1.0);
    cvxm_set(G, 4, 4, -1.0);
    cvxm_set(G, 4, 13, 1.0);
    cvxm_set(h, 4, 0, -rho);

    // x3 - x5 + w3 <= -rho
    cvxm_set(G, 5, 4,  1.0);
    cvxm_set(G, 5, 6, -1.0);
    cvxm_set(G, 5, 14, 1.0);
    cvxm_set(h, 5, 0, -rho);

    // x4 - x5 + w4 <= -rho
    cvxm_set(G, 6, 5,  1.0);
    cvxm_set(G, 6, 6, -1.0);
    cvxm_set(G, 6, 15, 1.0);
    cvxm_set(h, 6, 0, -rho);

    // -W + x5 + w5 <= 0  
    cvxm_set(G, 7, 0, -1.0);
    cvxm_set(G, 7, 6,  1.0);
    cvxm_set(G, 7, 16, 1.0);

    // -y2 <= 0  
    cvxm_set(G, 8, 8, -1.0);

    // -y3 <= 0  
    cvxm_set(G, 9, 9, -1.0);

    // -y5 <= 0  
    cvxm_set(G, 10, 11, -1.0);

    // -y1 + y2 + h2 <= -rho  
    cvxm_set(G, 11, 7, -1.0);
    cvxm_set(G, 11, 8,  1.0);
    cvxm_set(G, 11, 18, 1.0);
    cvxm_set(h, 11, 0, -rho);

    // y1 - y4 + h1 <= -rho  
    cvxm_set(G, 12, 7,  1.0);
    cvxm_set(G, 12, 10,-1.0);
    cvxm_set(G, 12, 17, 1.0);
    cvxm_set(h, 12, 0, -rho);

    // y3 - y4 + h3 <= -rho  
    cvxm_set(G, 13, 9,  1.0);
    cvxm_set(G, 13, 10,-1.0);
    cvxm_set(G, 13, 19, 1.0);
    cvxm_set(h, 13, 0, -rho);

    // -H + y4 + h4 <= 0  
    cvxm_set(G, 14, 1, -1.0);
    cvxm_set(G, 14, 10, 1.0);
    cvxm_set(G, 14, 20, 1.0);

    // -H + y5 + h5 <= 0  
    cvxm_set(G, 15, 1, -1.0);
    cvxm_set(G, 15, 11, 1.0);
    cvxm_set(G, 15, 21, 1.0);

    // -w1 + h1/gamma <= 0  
    cvxm_set(G, 16, 12,-1.0);
    cvxm_set(G, 16, 17, 1.0/gamma);

    // w1 - gamma * h1 <= 0  
    cvxm_set(G, 17, 12, 1.0);
    cvxm_set(G, 17, 17, -gamma);

    // -w2 + h2/gamma <= 0  
    cvxm_set(G, 18, 13,-1.0);
    cvxm_set(G, 18, 18, 1.0/gamma);

    //  w2 - gamma * h2 <= 0  
    cvxm_set(G, 19, 13, 1.0);
    cvxm_set(G, 19, 18, -gamma);

    // -w3 + h3/gamma <= 0  
    cvxm_set(G, 20, 14,-1.0);
    cvxm_set(G, 20, 18, 1.0/gamma);

    //  w3 - gamma * h3 <= 0  
    cvxm_set(G, 21, 14, 1.0);
    cvxm_set(G, 21, 19, -gamma);

    // -w4  + h4/gamma <= 0  
    cvxm_set(G, 22, 15,-1.0);
    cvxm_set(G, 22, 19, 1.0/gamma);

    //  w4 - gamma * h4 <= 0  
    cvxm_set(G, 23, 15, 1.0);
    cvxm_set(G, 23, 20, -gamma);

    // -w5 + h5/gamma <= 0  
    cvxm_set(G, 24, 16,-1.0);
    cvxm_set(G, 24, 21, 1.0/gamma);

    //  w5 - gamma * h5 <= 0.0  
    cvxm_set(G, 25, 16, 1.0);
    cvxm_set(G, 25, 21, -gamma);
}


int main(int argc, char **argv)
{
    cvx_matrix_t *c, *G, *h, A, b;
    cvx_problem_t cpl;
    cvx_dimset_t dims;
    cvx_float_t rho, gamma;
    int opt;
    
    cvx_float_t amin_data[] = {100.0, 100.0, 100.0, 100.0, 100.0};
    cvx_solopts_t opts = (cvx_solopts_t){
        .abstol = 0.0,
        .reltol = 0.0,
        .feastol = 0.0,
        .max_iter = 1,
        .debug = 0,
        .refinement = 0,
        .kkt_solver_name = 0,
        .show_progress = 1
    };
    floorplan_t Fplan;
    cvx_convex_program_t F;
    
    
    while ((opt = getopt(argc, argv, "N:")) != -1) {
        switch (opt) {
        case 'N':
            opts.max_iter = atoi(optarg);
            break;
        default:
            break;
        }
    }
    

    cvx_dimset_create(&dims, 5, 26, 0, 0);

    cvxm_map_data(&Fplan.Amin, 5, 1, amin_data);
    cvx_convex_program_init(&F, floorplan_F, &Fplan);
    
    rho = 1.0; gamma = 5.0;
    c = cvxm_new(22, 1);
    h = cvxm_new(26, 1);
    G = cvxm_new(26, 22);;
    cvxm_map_data(&A, 0, 22, (cvx_float_t *)0);
    cvxm_map_data(&b, 0, 1,  (cvx_float_t *)0);

    floorplan_constraints(c, G, h, rho, gamma);
    //cvx_mat_printf(stderr, "%9.2e", G, "G");

    if (opts.max_iter == 0)
        return 0;

    cvx_cpl_setup(&cpl, &F, c, G, h, &A, &b, &dims, (cvx_kktsolver_t *)0);
    cvx_cpl_compute_start(&cpl);
    cvx_cpl_solve(&cpl, &opts);
    return 0;
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
