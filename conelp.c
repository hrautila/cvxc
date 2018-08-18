
// Copyright: Harri Rautila, 2016 <harri.rautila@gmail.com>

#include "convex.h"

#define EFRM "%8.5f"
//#define EFRM "%9.2e"
#define STEP 0.99


static 
int cvx_res(cvx_conelp_problem_t *cp,
            cvx_matrix_t *ux,
            cvx_matrix_t *uy,
            cvx_matgrp_t *uz_g,
            double utau,
            cvx_matgrp_t *us_g,
            double ukappa,
            cvx_matrix_t *vx,
            cvx_matrix_t *vy,
            cvx_matgrp_t *vz_g,
            double *vtau,
            cvx_matgrp_t *vs_g,
            double *vkappa,
            cvx_scaling_t *W,
            double dg,
            cvx_matgrp_t *lmbda_g)
{
    cvx_matrix_t *us = us_g->mat, *uz = uz_g->mat;
    cvx_matrix_t *vs = vs_g->mat, *vz = vz_g->mat;
    cvx_matrix_t *lmbda = lmbda_g->mat;
    cvx_float_t tauplus, lscale, vkplus;
    cvx_matgrp_t ux_g = (cvx_matgrp_t){ .mat = ux, .index = (cvx_index_t *)0 };
    int err = 0;
    
    //printf("res:start utau=%e, ukappa=%e, vtau=%e, vkappa=%e, dg=%e\n", utau, ukappa, *vtau, *vkappa, dg);
    // vx = vx - A^T*uy - G^T*W^-1*uz - c*utau/dg
    cvxm_mult(1.0, vx, -1.0, cp->A, uy, CVX_TRANSA);
    cvxm_copy(&cp->wz3, uz, CVX_ALL);
    cvx_scale(&cp->wz3_g, &cp->W, CVX_INV, &cp->work);
    //cvx_mat_printf(stdout, "%e", &cp->wz3, "res_entry:wz3");
    cvx_sgemv(1.0, vx, -1.0, cp->G, &cp->wz3_g, CVX_TRANS);
    //cvx_mat_printf(stdout, "%e", vx, "res_entry.2:vx");
    cvxm_axpy(vx, -utau/dg, cp->c);

    // vy = vy + A*ux - b*utau/dg
    cvxm_mult(1.0, vy, 1.0, cp->A, ux, 0);
    cvxm_axpy(vy, -utau/dg, cp->b);

    // vz = vz + G*ux - h*utau/dg + W^T*us
    cvx_sgemv(1.0, vz, 1.0, cp->G, &ux_g, 0);
    cvxm_axpy(vz, -utau/dg, cp->h);
    cvxm_copy(&cp->ws3, us, 0);
    cvx_scale(&cp->ws3_g, &cp->W, CVX_TRANS, &cp->work);
    cvxm_axpy(vz, 1.0, &cp->ws3);
    //cvx_mat_printf(stdout, "%e", vz, "res_entry:update(vz)");

    // vtau : vtau + c'*ux + b'uy + h'W^-1*uz + dg*ukappa
    tauplus  = dg*ukappa;
    tauplus += cvxm_dot(cp->c, ux);
    tauplus += cvxm_dot(cp->b, uy);
    tauplus += cvx_sdot(&cp->h_g, &cp->wz3_g);
    *vtau   += tauplus;
    
    // vs = vs + lmbda o (uz + us)
    cvxm_copy(&cp->ws3, us, 0);
    cvxm_axpy(&cp->ws3, 1.0, uz);
    cvx_sprod(&cp->ws3_g, lmbda_g, CVX_DIAG, &cp->work);
    cvxm_axpy(vs, 1.0, &cp->ws3);

    // vkappa += vkappa + lmbdag * (utau + ukappa)
    lscale = cvxm_get(lmbda, cp->cdim_diag, 0);
    vkplus = lscale * (utau + ukappa);
    *vkappa += vkplus;
        
    //printf("res:end vtau=%e, vkappa=%e, tauplus=%e, lscale=%e, vkplus=%e\n", *vtau, *vkappa, tauplus, lscale, vkplus);
    
    return err;
}

// f6_no_ir(x, y, z, tau, s, kappa) solves
//
//    [ 0         ]   [  0   A'  G'  c ] [ ux        ]    [ bx   ]
//    [ 0         ]   [ -A   0   0   b ] [ uy        ]    [ by   ]
//    [ W'*us     ] - [ -G   0   0   h ] [ W^{-1}*uz ] = -[ bz   ]
//    [ dg*ukappa ]   [ -c' -b' -h'  0 ] [ utau/dg   ]    [ btau ]
//
//    lmbda o (uz + us) = -bs
//    lmbdag * (utau + ukappa) = -bkappa.
//
// On entry, x, y, z, tau, s, kappa contain bx, by, bz, btau, 
// bkappa.  On exit, they contain ux, uy, uz, utau, ukappa.
static 
int f6_no_ir(cvx_conelp_problem_t *cp,
             cvx_matrix_t *x,
             cvx_matrix_t *y,
             cvx_matgrp_t *z_g,
             cvx_float_t *tau,
             cvx_matgrp_t *s_g,
             cvx_float_t *kappa)
{
    cvx_matrix_t *s = s_g->mat, *z = z_g->mat;
    //cvx_matrix_t r;

    // Solve 
    // 
    // [  0   A'  G'    0   ] [ ux        ]   
    // [ -A   0   0     b   ] [ uy        ]  
    // [ -G   0   W'*W  h   ] [ W^{-1}*uz ] 
    // [ -c' -b' -h'    k/t ] [ utau/dg   ]
    // 
    //   [ bx                    ]
    //   [ by                    ]
    // = [ bz - W'*(lmbda o\ bs) ]
    //   [ btau - bkappa/tau     ]
    //
    // us = -lmbda o\ bs - uz
    // ukappa = -bkappa/lmbdag - utau.
    
    // First solve 
    // 
    // [ 0  A' G'   ] [ ux        ]   [  bx                    ]
    // [ A  0  0    ] [ uy        ] = [ -by                    ]
    // [ G  0 -W'*W ] [ W^{-1}*uz ]   [ -bz + W'*(lmbda o\ bs) ]
    
    //printf("f6_no_ir:start tau=%e, kappa=%e\n", *tau, *kappa);
    cvxm_scale(y, -1.0, 0);
    // s = -lmbda o\ s = - lmbda o\ bs
    //cvx_mat_printf(stdout, "%e", &cp->lmbda, "f6_no_ir:lmbda");
    cvx_sinv(s_g, &cp->lmbda_g, &cp->work);
    cvxm_scale(s, -1.0, 0);
    //cvx_mat_printf(stdout, "%e", s, "f6_no_ir:s");
    
    // z = -(z + W'*s) = -bz + W'(lmbda o\ bs)
    cvxm_copy(&cp->ws3, s, CVX_ALL);
    //cvx_mat_printf(stdout, "%e", &cp->ws3, "f6_no_ir:ws3");
    //cvx_scaling_elem(&r, &cp->W, CVXWS_R, 0);
    //cvx_mat_printf(stdout, "%13.6e", &r, "f6_no_ir:r");
    cvx_scale(&cp->ws3_g, &cp->W, CVX_TRANS, &cp->work);
    //cvx_mat_printf(stdout, "%e", &cp->ws3, "f6_no_ir:scaled(ws3)");
    cvxm_axpy(z, 1.0, &cp->ws3);
    cvxm_scale(z, -1.0, 0);
    //cvx_mat_printf(stdout, "%e", z, "f6_no_ir:z");

    //cvx_mat_printf(stdout, "%e", x, "f6_no_ir:x");
    cvx_kktsolve(cp->solver, x, y, z_g);
    //cvx_mat_printf(stdout, "%e", x, "f6_no_ir:solved(x)");
    //cvx_mat_printf(stdout, "%e", z, "f6_no_ir:solved(z)");
    // Combine with solution of 
    // 
    // [ 0   A'  G'    ] [ x1         ]          [ c ]
    // [-A   0   0     ] [ y1         ] = -dgi * [ b ]
    // [-G   0   W'*W  ] [ W^{-1}*dzl ]          [ h ]
    // 
    // to satisfy
    // 
    // -c'*x - b'*y - h'*W^{-1}*z + dg*tau = btau - bkappa/tau. '
    cvx_float_t lkappa = - (*kappa) / cvxm_get(&cp->lmbda, cp->cdim_diag, 0);
    cvx_float_t ltau  = *tau + lkappa/cp->dgi;
    //printf("f6_no_ir:tau.0=%e, kappa=%e\n", ltau, lkappa);
    ltau += cvxm_dot(cp->c, x);
    //printf("f6_no_ir:tau.1=%e, kappa=%e\n", ltau, lkappa);
    ltau += cvxm_dot(cp->b, y);
    //printf("f6_no_ir:tau.2=%e, kappa=%e\n", ltau, lkappa);
    ltau += cvx_sdot(&cp->th_g, z_g);  // TODO: check this
    //printf("f6_no_ir:tau.3=%e, kappa=%e\n", ltau, lkappa);
    ltau = cp->dgi * ltau  / (1.0 + cvx_sdot(&cp->z1_g, &cp->z1_g));
    //printf("f6_no_ir:tau=%e, kappa=%e\n", ltau, lkappa);

    *tau = ltau;
    cvxm_axpy(x, ltau, &cp->x1);
    cvxm_axpy(y, ltau, &cp->y1);
    cvxm_axpy(z, ltau, &cp->z1);

    cvxm_axpy(s, -1.0, z);
    *kappa = lkappa - ltau;

    //cvx_mat_printf(stdout, "%e", x, "f6_no_ir(x)");
    //cvx_mat_printf(stdout, "%e", z, "f6_no_ir(z)");
    //cvx_mat_printf(stdout, "%e", s, "f6_no_ir(s)");
    //printf("f6_no_ir:end tau=%e, kappa=%e\n", *tau, *kappa);
    return 0;
}


// f6(x, y, z, tau, s, kappa) solves the same system as f6_no_ir, 
// but applies iterative refinement.
static
int f6(cvx_conelp_problem_t *cp,
       cvx_matrix_t *x,
       cvx_matrix_t *y,
       cvx_matgrp_t *z_g,
       cvx_float_t *tau,
       cvx_matgrp_t *s_g,
       cvx_float_t *kappa,
       int refinement)
{
    cvx_matrix_t *s = s_g->mat, *z = z_g->mat;
    int /*lags = 0,*/ err = 0;
    double wtau, wkappa, wtau2, wkappa2;
    
    //printf("**f6:start tau=%e, kappa=%e, refinement=%d\n", *tau, *kappa, refinement);

    if (refinement > 0 /*|| (flags | CVX_DEBUG) != 0*/) {
        cvxm_copy(&cp->wx, x, 0);
        cvxm_copy(&cp->wy, y, 0);
        cvxm_copy(&cp->wz, z, 0);
        cvxm_copy(&cp->ws, s, 0);
        wtau   = *tau;
        wkappa = *kappa;
    }
    err = f6_no_ir(cp, x, y, z_g, tau, s_g, kappa);
    for (int i = 0; i < refinement; i++) {
        cvxm_copy(&cp->wx2, &cp->wx, 0);
        cvxm_copy(&cp->wy2, &cp->wy, 0);
        cvxm_copy(&cp->wz2, &cp->wz, 0);
        cvxm_copy(&cp->ws2, &cp->ws, 0);
        wtau2   = wtau;
        wkappa2 = wkappa;

        cvx_res(cp, x, y, z_g, *tau, s_g, *kappa, &cp->wx2, &cp->wy2, &cp->wz2_g, &wtau2,
                &cp->ws2_g, &wkappa2, &cp->W, cp->dg, &cp->lmbda_g);

        //cvx_mat_printf(stdout, "%e", &cp->wx2, "f6:wx2");
        //cvx_mat_printf(stdout, "%e", &cp->wy2, "f6:wy2");
        //cvx_mat_printf(stdout, "%e", &cp->wz2, "f6:wz2");
        //cvx_mat_printf(stdout, "%e", &cp->ws2, "f6:ws2");
        f6_no_ir(cp, &cp->wx2, &cp->wy2, &cp->wz2_g, &wtau2, &cp->ws2_g, &wkappa2);
        
        //cvx_mat_printf(stdout, "%e", &cp->wz2, "f6_no_ir:wz2");
        //cvx_mat_printf(stdout, "%e", &cp->ws2, "f6_no_ir:ws2");

        cvxm_axpy(x, 1.0, &cp->wx2);
        cvxm_axpy(y, 1.0, &cp->wy2);
        cvxm_axpy(s, 1.0, &cp->ws2);
        cvxm_axpy(z, 1.0, &cp->wz2);
        
        *tau   += wtau2;
        *kappa += wkappa2;
    }
    //printf("**f6:end tau=%e, kappa=%e, refinement=%d\n", *tau, *kappa, refinement);
    return err;
}


cvx_conelp_problem_t *cvx_conelp_setup(cvx_conelp_problem_t *prob,
                                       cvx_matrix_t *c,
                                       cvx_matrix_t *G,
                                       cvx_matrix_t *h,
                                       cvx_matrix_t *A,
                                       cvx_matrix_t *b,
                                       cvx_dimset_t *dims,
                                       cvx_kktsolver_t *kktsolver)
{

    cvx_size_t mc, nc, mG, nG, mA, nA, mh, nh, mb, nb;
    

    mc = nc = mG = nG = mA = nA = mh = nh = mb = nb = 0;
    if (! prob)
        return prob;

    if (!c) {
        prob->error = CVX_ERR_NULLCOST;
        return (cvx_conelp_problem_t *)0;
    }

    cvxm_size(&mc, &nc, c);
    if (G)
        cvxm_size(&mG, &nG, G);
    if (A)
        cvxm_size(&mA, &nA, A);
    if (h)
        cvxm_size(&mh, &nh, h);
    if (b)
        cvxm_size(&mb, &nb, b);
    
    if (nc > 1 || mc < 1) {
        prob->error = CVX_ERR_DIMC;
        return (cvx_conelp_problem_t *)0;
    }

    cvx_size_t cdim = cvx_dimset_sum(dims, CVXDIM_LINEAR) +
        cvx_dimset_sum(dims, CVXDIM_SOCP) +
        cvx_dimset_sum_squared(dims, CVXDIM_SDP);

    cvx_size_t cdim_pckd = cvx_dimset_sum(dims, CVXDIM_LINEAR) +
        cvx_dimset_sum(dims, CVXDIM_SOCP) +
        cvx_dimset_sum_packed(dims, CVXDIM_SDP);
    
    cvx_size_t cdim_diag = cvx_dimset_sum(dims, CVXDIM_LINEAR) +
        cvx_dimset_sum(dims, CVXDIM_SOCP) +
        cvx_dimset_sum(dims, CVXDIM_SDP);
    
    if (nh > 1 || mh != cdim) {
        prob->error = CVX_ERR_DIMH;
        return (cvx_conelp_problem_t *)0;
    }

    if (mG != cdim || nG != mc) {
        prob->error = CVX_ERR_DIMG;
        return (cvx_conelp_problem_t *)0;
    }
    if (nA != mc || mA != mb) {
        prob->error = CVX_ERR_DIMA;
        return (cvx_conelp_problem_t *)0;
    }
    if (nb != 1) {
        prob->error = CVX_ERR_DIMB;
        return (cvx_conelp_problem_t *)0;
    }
    if ( mb > mc || mb + cdim_pckd < mc) {
        prob->error = CVX_ERR_RANK;
        return (cvx_conelp_problem_t *)0;
    }

    // make space reservations; first init matrix to zero size; compute total space needed
    // and make one huge allocation; then divide space to matrices by mapping them over

    cvx_size_t sdim = cvx_dimset_sum(dims, CVXDIM_SDP);
    cvx_size_t total = 0;
    total += 7*mc*nc;  // for x, dx, rx, hrx, x1, wx, wx2
    total += 7*mb*nb;  // for y, dy, ry, hry, y1, wy, wy2
    total += 5*cdim;   // for s, ds, ws, ws2, ws3
    total += 8*cdim;   // for z, dz, rz, hrz, z1, wz, wz2, wz3
    total += 2*cdim_diag + 2; // for lmbda, lmbdasq
    total += 2*sdim;   // for sigs, sigz
    total += mh*nh;    // for th

    cvx_float_t *space = (cvx_float_t *)calloc(total, sizeof(cvx_float_t));
    if (!space) {
        prob->error = CVX_ERR_MEMORY;
        return (cvx_conelp_problem_t *)0;
    }
    prob->mlen = total;
    prob->memory = space;

    prob->c = c;
    prob->G = G;
    prob->h = h;
    prob->A = A;
    prob->b = b;
    prob->dims = dims;

    prob->primal_x = (cvx_matrix_t *)0;
    prob->primal_s = (cvx_matrix_t *)0;
    prob->dual_y = (cvx_matrix_t *)0;
    prob->dual_z = (cvx_matrix_t *)0;

    prob->cdim = cdim;
    prob->cdim_diag = cdim_diag;
    
    // init index sets for access to s, z and friends
    cvx_index_init(&prob->index_full, dims, 0);
    cvx_index_init(&prob->index_packed, dims, 1); // TODO: is this needed??
    cvx_index_init(&prob->index_diag, dims, 2);
    cvx_index_init(&prob->index_sig, dims, 3);

    // map matrices to data space;
    cvx_size_t offset = 0;

    // result matrices: x, y, s, z  (TODO: separate space for results?)
    cvxm_map_data(&prob->x, mc, nc, &space[offset]);
    offset += mc*nc;
    cvxm_map_data(&prob->y, mb, nb, &space[offset]);
    offset += mb*nb;
    cvxm_map_data(&prob->s, cdim, 1, &space[offset]);
    offset += cdim;
    cvxm_map_data(&prob->z, cdim, 1, &space[offset]);
    offset += cdim;

    // dx, dy, ds, dz
    cvxm_map_data(&prob->dx, mc, nc, &space[offset]);
    offset += mc*nc;
    cvxm_map_data(&prob->dy, mb, nb, &space[offset]);
    offset += mb*nb;
    cvxm_map_data(&prob->ds, cdim, 1, &space[offset]);
    offset += cdim;
    cvxm_map_data(&prob->dz, cdim, 1, &space[offset]);
    offset += cdim;

    // rx, ry, rz
    cvxm_map_data(&prob->rx, mc, nc, &space[offset]);
    offset += mc*nc;
    cvxm_map_data(&prob->ry, mb, nb, &space[offset]);
    offset += mb*nb;
    cvxm_map_data(&prob->rz, cdim, 1, &space[offset]);
    offset += cdim;
    
    // x1, y1, z1
    cvxm_map_data(&prob->x1, mc, nc, &space[offset]);
    offset += mc*nc;
    cvxm_map_data(&prob->y1, mb, nb, &space[offset]);
    offset += mb*nb;
    cvxm_map_data(&prob->z1, cdim, 1, &space[offset]);
    offset += cdim;

    // hrx, hry, hrz
    cvxm_map_data(&prob->hrx, mc, nc, &space[offset]);
    offset += mc*nc;
    cvxm_map_data(&prob->hry, mb, nb, &space[offset]);
    offset += mb*nb;
    cvxm_map_data(&prob->hrz, cdim, 1, &space[offset]);
    offset += cdim;

    // wx, wy, ws, wz
    cvxm_map_data(&prob->wx, mc, nc, &space[offset]);
    offset += mc*nc;
    cvxm_map_data(&prob->wy, mb, nb, &space[offset]);
    offset += mb*nb;
    cvxm_map_data(&prob->ws, cdim, 1, &space[offset]);
    offset += cdim;
    cvxm_map_data(&prob->wz, cdim, 1, &space[offset]);
    offset += cdim;

    // wx2, wy2, ws2, wz2
    cvxm_map_data(&prob->wx2, mc, nc, &space[offset]);
    offset += mc*nc;
    cvxm_map_data(&prob->wy2, mb, nb, &space[offset]);
    offset += mb*nb;
    cvxm_map_data(&prob->ws2, cdim, 1, &space[offset]);
    offset += cdim;
    cvxm_map_data(&prob->wz2, cdim, 1, &space[offset]);
    offset += cdim;
   
    // ws3, wz3
    cvxm_map_data(&prob->ws3, cdim, 1, &space[offset]);
    offset += cdim;
    cvxm_map_data(&prob->wz3, cdim, 1, &space[offset]);
    offset += cdim;

    // lmbda, lmbdasq (length == cdim_diag + 1)
    cvxm_map_data(&prob->lmbda, cdim_diag + 1, 1, &space[offset]);
    offset += cdim_diag + 1;
    cvxm_map_data(&prob->lmbdasq, cdim_diag + 1, 1, &space[offset]);
    offset += cdim_diag + 1;

    // sigs, sigz  (length == sum(dimset.S)) will hold eigenvalues of 'S' space
    cvxm_map_data(&prob->sigs, sdim, 1, &space[offset]);
    offset += sdim;
    cvxm_map_data(&prob->sigz, sdim, 1, &space[offset]);
    offset += sdim;

    // th
    cvxm_map_data(&prob->th, mh*nh, 1, &space[offset]);
    offset += mh*nh;

    // setup matrix group variables
    cvx_mgrp_init(&prob->h_g,   prob->h,    &prob->index_full);
    cvx_mgrp_init(&prob->s_g,   &prob->s,   &prob->index_full);
    cvx_mgrp_init(&prob->z_g,   &prob->z,   &prob->index_full);
    cvx_mgrp_init(&prob->z1_g,  &prob->z1,  &prob->index_full);
    cvx_mgrp_init(&prob->th_g,  &prob->th,  &prob->index_full);
    cvx_mgrp_init(&prob->ds_g,  &prob->ds,  &prob->index_full);
    cvx_mgrp_init(&prob->dz_g,  &prob->dz,  &prob->index_full);
    cvx_mgrp_init(&prob->rs_g,  &prob->rs,  &prob->index_full);
    cvx_mgrp_init(&prob->rz_g,  &prob->rz,  &prob->index_full);
    cvx_mgrp_init(&prob->hrs_g, &prob->hrs, &prob->index_full);
    cvx_mgrp_init(&prob->hrz_g, &prob->hrz, &prob->index_full);
    cvx_mgrp_init(&prob->ws_g,  &prob->ws,  &prob->index_full);
    cvx_mgrp_init(&prob->wz_g,  &prob->wz,  &prob->index_full);
    cvx_mgrp_init(&prob->ws2_g, &prob->ws2, &prob->index_full);
    cvx_mgrp_init(&prob->wz2_g, &prob->wz2, &prob->index_full);
    cvx_mgrp_init(&prob->ws3_g, &prob->ws3, &prob->index_full);
    cvx_mgrp_init(&prob->wz3_g, &prob->wz3, &prob->index_full);

    cvx_mgrp_init(&prob->sigs_g,   &prob->sigs,   &prob->index_sig);
    cvx_mgrp_init(&prob->sigz_g,   &prob->sigz,   &prob->index_sig);

    cvx_mgrp_init(&prob->lmbda_g,   &prob->lmbda,   &prob->index_diag);
    cvx_mgrp_init(&prob->lmbdasq_g, &prob->lmbdasq, &prob->index_diag);

    // allocate scaling
    cvx_scaling_init(&prob->W, dims);

    if (kktsolver) {
        prob->solver = kktsolver;
    } else {
        cvx_ldlsolver_init(&prob->__S, prob, dims, dims->mnl);
        prob->solver = &prob->__S;
    }
    
    // TODO: initialized workspace; needed if SDP constraints
    __mblk_empty(&prob->work);
    cvx_size_t maxsdp = cvx_dimset_max(dims, CVXDIM_SDP);
    if (maxsdp > 0)
        __mblk_init(&prob->work, 6*maxsdp*maxsdp);

    return prob;
}

void cvx_conelp_set_start(cvx_conelp_problem_t *prob,
                          cvx_matrix_t *primal_x,
                          cvx_matrix_t *primal_s,
                          cvx_matrix_t *dual_y,
                          cvx_matrix_t *dual_z)
{
    if (! prob)
        return;
    // TODO: check sizes; and must give primal_x, primal_s or neither (likewise for y,z)
    prob->primal_x = primal_x;
    prob->primal_s = primal_s;
    prob->dual_y = dual_y;
    prob->dual_z = dual_z;
}

void cvx_conelp_init_scaling(cvx_scaling_t *W)
{
    cvx_matrix_t x;
    // D = ones
    cvx_scaling_elem(&x, W, CVXWS_D, 0);
    cvxm_mkconst(&x, 1.0);
    // DI = ones
    cvx_scaling_elem(&x, W, CVXWS_DI, 0);
    cvxm_mkconst(&x, 1.0);
    if (W->vcount > 0) {
        // BETA = ones
        cvx_scaling_elem(&x, W, CVXWS_BETA, 0);
        cvxm_mkconst(&x, 1.0);
        // V is unit vector
        for (int k = 0; k < W->vcount; k++) {
            cvx_scaling_elem(&x, W, CVXWS_V, k);
            cvxm_mkident(&x);
        }
    }
    // R & RTI are identity
    for (int k = 0; k < W->rcount; k++) {
        cvx_scaling_elem(&x, W, CVXWS_R, k);
        cvxm_mkident(&x);
        cvx_scaling_elem(&x, W, CVXWS_RTI, k);
        cvxm_mkident(&x);
    }
}

int cvx_conelp_compute_start(cvx_conelp_problem_t *prob)
{
    cvx_stats_t *stats = &prob->stats;

    int primalstart = ! (prob->primal_x && prob->primal_s);
    int dualstart = ! (prob->dual_y && prob->dual_z);

    //printf("solver.mnl: %ld\n", prob->solver->u.ldl.mnl);
    if (primalstart || dualstart) {      
        cvx_conelp_init_scaling(&prob->W);       
        cvx_kktfactor(prob->solver, &prob->W, __cvxnil, __cvxnil);
    }
    if (primalstart) {
        // minimize    || G * x - h ||^2
        // subject to  A * x = b
        //
        // by solving
        //
        //     [ 0   A'  G' ]   [ x  ]   [ 0 ]
        //     [ A   0   0  ] * [ dy ] = [ b ].
        //     [ G   0  -I  ]   [ -s ]   [ h ]
        cvxm_scale(&prob->x, 0.0, 0);
        cvxm_copy(&prob->dy, prob->b, 0);
        cvxm_copy(&prob->s, prob->h, 0);
        //cvx_mgrp_printf(stdout, EFRM, &prob->s_g, "initial s");
        
        cvx_kktsolve(prob->solver, &prob->x, &prob->dy, &prob->s_g);
        cvxm_scale(&prob->s, -1.0, 0);
    } else {
        cvxm_copy(&prob->x, prob->primal_x, 0);
        cvxm_copy(&prob->s, prob->primal_s, 0);
    }
    //cvx_mat_printf(stdout, "%e", &prob->s, "primal s");
    //cvx_mat_printf(stdout, "%e", &prob->x, "primal x");

    prob->ts = cvx_max_step(&prob->s_g, __nilgrp, &prob->work);
    //fprintf(stdout, "initial ts=%f\n", prob->ts);

    if (prob->ts >= 0.0 && ! primalstart) {
        prob->error = CVX_ERR_NEG_INITIAL_S;
        return -1;
    }
    if (dualstart) {
        // minimize   || z ||^2
        // subject to G'*z + A'*y + c = 0
        //
        // by solving 
        //
        //     [ 0   A'  G' ] [ dx ]   [ -c ]
        //     [ A   0   0  ] [ y  ] = [  0 ].
        //     [ G   0  -I  ] [ z  ]   [  0 ]
        cvxm_scale(&prob->y, 0.0, 0);
        cvxm_copy(&prob->dx, prob->c, 0);
        cvxm_scale(&prob->dx, -1.0, 0);
        cvxm_scale(&prob->z, 0.0, 0);
        //cvx_mgrp_printf(stdout, EFRM, &prob->z_g, "initial z");
        cvx_kktsolve(prob->solver, &prob->dx, &prob->y, &prob->z_g);
    } else {
        cvxm_copy(&prob->y, prob->dual_y, 0);
        cvxm_copy(&prob->z, prob->dual_z, 0);
    }
    //cvx_mat_printf(stdout, "%e", &prob->z, "dual z");
    //cvx_mat_printf(stdout, "%e", &prob->y, "dual y");
        
    prob->tz = cvx_max_step(&prob->z_g,  __nilgrp, &prob->work);
    //fprintf(stdout, "initial tz=%f\n", prob->tz);
    if (prob->tz >= 0.0 && ! dualstart) {
        prob->error = CVX_ERR_NEG_INITIAL_Z;
        return -1;
    }

    stats->resx0 = __maxvec(2, (cvx_float_t []){1.0, cvxm_nrm2(prob->c)});
    stats->resy0 = __maxvec(2, (cvx_float_t []){1.0, cvxm_nrm2(prob->b)});
    stats->resz0 = __maxvec(2, (cvx_float_t []){1.0, cvx_snrm2(&prob->h_g)});

    //printf("initial: resz0=%.4f, resy0=%.4f, resz0=%.4f\n", stats->resx0, stats->resy0, stats->resz0);
    prob->nrms = cvx_snrm2(&prob->s_g);
    prob->nrmz = cvx_snrm2(&prob->z_g);

    //printf("initial nrms=%.5f, nrmz=%.5f\n", prob->nrms, prob->nrmz);
    
    stats->gap = stats->pcost = stats->dcost = stats->relgap = 0.0;
    
    return 0;
}



int cvx_conelp_ready(cvx_conelp_problem_t *prob, cvx_stats_t *stats, int iter, int stat)
{
    //cvx_stats_t *stats = &prob->stats;

    if (stat == CVX_STAT_OPTIMAL && iter == -1) {
        // constructed initial point is feasible and optimal
        cvx_mksymm(&prob->s_g);
        cvx_mksymm(&prob->z_g);

        // rx = A'*y + G'z + c
        cvxm_copy(&prob->rx, prob->c, CVX_ALL);
        cvxm_mvmult(1.0, &prob->rx, 1.0, prob->A, &prob->y, CVX_TRANS);
        cvxm_mvmult(1.0, &prob->rx, 1.0, prob->G, &prob->z, CVX_TRANS);
        stats->resx = cvxm_nrm2(&prob->rx);

        // ry = b - A*x  ; TODO - computes -b - A*x ;; check
        cvxm_copy(&prob->ry, prob->b, CVX_ALL);
        cvxm_mvmult(-1.0, &prob->ry, -1.0, prob->A, &prob->x, CVX_TRANS); 
        stats->resy = cvxm_nrm2(&prob->ry);
        
        // rz = s + G*x - h
        cvxm_copy(&prob->rz, &prob->s, CVX_ALL);
        cvxm_mvmult(1.0, &prob->rz, 1.0, prob->G, &prob->x, 0);
        cvxm_axpy(&prob->rz, 1.0, prob->h);
        stats->resz = cvx_snrm2(&prob->rz_g);

        stats->pres = __MAX2(stats->resy/stats->resy0, stats->resz/stats->resz0);
        stats->dres = stats->resx/stats->resx0;
        stats->cx = cvxm_dot(prob->c, &prob->x);
        stats->by = cvxm_dot(prob->b, &prob->y);
        stats->hz = cvx_sdot(&prob->h_g, &prob->z_g);
        
        prob->solution.x = &prob->x;
        prob->solution.s = &prob->s;
        prob->solution.y = &prob->y;
        prob->solution.z = &prob->z;

        prob->solution.status = stat;
        prob->solution.gap = stats->gap;
        prob->solution.relative_gap = stats->relgap;
        prob->solution.primal_objective = stats->cx;
        prob->solution.dual_objective = - (stats->by + stats->hz);
        prob->solution.primal_infeasibility = stats->pres;
        prob->solution.dual_infeasibility = stats->dres;
        prob->solution.primal_slack = - prob->ts;
        prob->solution.dual_slack = - prob->tz;
        prob->solution.primal_residual_cert = __NaN();
        prob->solution.dual_residual_cert = __NaN();
        prob->solution.iterations = 0;
    }
    else if (stat == CVX_STAT_UNKNOWN || stat == CVX_STAT_OPTIMAL) {
        cvxm_scale(&prob->x, 1.0/prob->tau, CVX_ALL);
        cvxm_scale(&prob->y, 1.0/prob->tau, CVX_ALL);
        cvxm_scale(&prob->s, 1.0/prob->tau, CVX_ALL);
        cvxm_scale(&prob->z, 1.0/prob->tau, CVX_ALL);
        cvx_mksymm(&prob->s_g);
        cvx_mksymm(&prob->z_g);
        prob->ts = cvx_max_step(&prob->s_g, __nilgrp, &prob->work);
        prob->tz = cvx_max_step(&prob->s_g, __nilgrp, &prob->work);
        
        prob->error = stat == CVX_STAT_UNKNOWN ? CVX_ERR_MAXITER : 0;
        prob->solution.x = &prob->x;
        prob->solution.s = &prob->s;
        prob->solution.y = &prob->y;
        prob->solution.z = &prob->z;

        prob->solution.status = stat;

        prob->solution.gap = stats->gap;
        prob->solution.relative_gap = stats->relgap;
        prob->solution.primal_objective = stats->cx;
        prob->solution.dual_objective = - (stats->by + stats->hz);
        prob->solution.primal_infeasibility = stats->pres;
        prob->solution.dual_infeasibility = stats->dres;
        prob->solution.primal_slack = - prob->ts;
        prob->solution.dual_slack = - prob->tz;
        if (stat == CVX_STAT_OPTIMAL) {
            prob->solution.primal_residual_cert = __NaN();
            prob->solution.dual_residual_cert = __NaN();
        } else {
            prob->solution.primal_residual_cert = stats->pinfres;
            prob->solution.dual_residual_cert = stats->dinfres;
        }
        prob->solution.iterations = iter;
    }
    else if (stat == CVX_STAT_PRIMAL_INFEASIBLE) {
        prob->solution.status = stat;
        prob->error = stat;
        cvxm_scale(&prob->y, 1.0/(-stats->hz - stats->by), CVX_ALL);
        cvxm_scale(&prob->z, 1.0/(-stats->hz - stats->by), CVX_ALL);
        cvx_mksymm(&prob->z_g);

        prob->solution.x = __cvxnil;
        prob->solution.s = __cvxnil;
        prob->solution.y = __cvxnil;
        prob->solution.z = __cvxnil;

        prob->solution.status = stat;

        prob->solution.gap = __NaN();
        prob->solution.relative_gap = __NaN();
        prob->solution.primal_objective = __NaN();
        prob->solution.dual_objective = 1.0;
        prob->solution.primal_infeasibility = __NaN();
        prob->solution.dual_infeasibility = __NaN();
        prob->solution.primal_slack = __NaN();
        prob->solution.dual_slack = - prob->tz;
        prob->solution.primal_residual_cert = stats->pinfres;
        prob->solution.dual_residual_cert = __NaN();
        prob->solution.iterations = iter;
    }
    else if (stat == CVX_STAT_DUAL_INFEASIBLE) {
        prob->error = stat;
        prob->solution.status = stat;
        cvxm_scale(&prob->x, 1.0/-stats->cx, CVX_ALL);
        cvxm_scale(&prob->s, 1.0/-stats->cx, CVX_ALL);
        cvx_mksymm(&prob->s_g);

        prob->solution.x = __cvxnil;
        prob->solution.s = __cvxnil;
        prob->solution.y = __cvxnil;
        prob->solution.z = __cvxnil;

        prob->solution.status = stat;

        prob->solution.gap = __NaN();
        prob->solution.relative_gap = __NaN();
        prob->solution.primal_objective = 1.0;
        prob->solution.dual_objective = __NaN();
        prob->solution.primal_infeasibility = __NaN();
        prob->solution.dual_infeasibility = __NaN();
        prob->solution.primal_slack = - prob->ts;
        prob->solution.dual_slack = __NaN();
        prob->solution.primal_residual_cert = __NaN();
        prob->solution.dual_residual_cert = stats->dinfres;
        prob->solution.iterations = iter;
    }
    return -stat;
}

int cvx_conelp_solve(cvx_conelp_problem_t *prob, cvx_solopts_t *opts)
{
    cvx_stats_t *stats = &prob->stats;
    
    cvx_float_t max_nrms, max_nrmz;
    cvx_matrix_t lk, sk, zk, ls, lz;
    
    cvx_float_t abstol = CVX_ABSTOL;
    cvx_float_t reltol = CVX_RELTOL;
    cvx_float_t feastol = CVX_FEASTOL;

    int maxiter = opts->max_iter > 0 ? opts->max_iter : CVX_MAXITER;
    int refinement = opts->refinement > 0 ? opts->refinement : 0;

    if (prob->dims->qlen > 0 || prob->dims->slen > 0)
        refinement = 1;   

    int primalstart = ! (prob->primal_x && prob->primal_s);
    int dualstart = ! (prob->dual_y && prob->dual_z);
    
    prob->error = 0;
    
    if (primalstart && dualstart) {
        // check for constructed initial point
        stats->gap = cvx_sdot(&prob->s_g, &prob->z_g);
        stats->pcost = cvxm_dot(prob->c, &prob->x);
        stats->dcost =
            - cvxm_dot(prob->b, &prob->y) - cvx_sdot(&prob->h_g, &prob->z_g);

        if (stats->pcost < 0.0) {
            stats->relgap = stats->gap / - stats->pcost;
        } else if (stats->dcost > 0.0) {
            stats->relgap = stats->gap / stats->dcost;
        } else {
            stats->relgap = __NaN();
        }

        //printf("gap: %e, pcost: %e, dcost: %e, relgap: %e\n", stats->gap, stats->pcost, stats->dcost, stats->relgap);

        if (prob->ts <= 0.0 && prob->tz < 0 &&
            (stats->gap <= abstol ||
             (!isnan(stats->relgap) && stats->relgap <= reltol))) {

            // initial point is feasible and optimal
            return cvx_conelp_ready(prob, stats, -1, CVX_STAT_OPTIMAL);
        }    

        max_nrms = __MAX2(1.0, prob->nrms);
        max_nrmz = __MAX2(1.0, prob->nrmz);
        if (prob->ts >= - 1e-8 * max_nrms) {
            cvx_mgrp_update_sz(&prob->s_g, 1.0+prob->ts, 0);
        }

        if (prob->tz >= -1e-8 * max_nrmz) {
            cvx_mgrp_update_sz(&prob->z_g, 1.0+prob->tz, 0);
        }
        //cvx_mat_printf(stdout, "%e", &prob->s, "s:");
        //cvx_mat_printf(stdout, "%e", &prob->z, "z:");
    } else if (primalstart && ! dualstart) {
        max_nrms = __MAX2(1.0, prob->nrms);
        if (prob->ts >= - 1e-8 * max_nrms) {
            cvx_mgrp_update_sz(&prob->s_g, 1.0+prob->ts, 0);
        }

    } else if (dualstart && ! primalstart) {
        max_nrmz = __MAX2(1.0, prob->nrmz);
        if (prob->tz >= -1e-8 * max_nrmz) {
            cvx_mgrp_update_sz(&prob->z_g, 1.0+prob->tz, 0);
        }
    }

    prob->tau = 1.0;
    prob->kappa = 1.0;
    prob->wkappa3 = 0.0;
    cvxm_copy(&prob->rx,  prob->c, 0);
    cvxm_copy(&prob->hrx, prob->c, 0);
    cvxm_copy(&prob->ry,  prob->b, 0);
    cvxm_copy(&prob->hry, prob->b, 0);

    stats->gap = cvx_sdot(&prob->s_g, &prob->z_g);
    
    //printf("preloop: resz0=%.4f, resy0=%.4f, resz0=%.4f\n", stats->resx0, stats->resy0, stats->resz0);
    // -----------------------------------------------------------------------------

    for (int iter = 0; iter <= maxiter; iter++) {
        
        // hrx = -A'*y - G'*z
        cvxm_scale(&prob->hrx, 0.0, CVX_ALL);
        cvxm_mvmult(0.0, &prob->hrx, -1.0, prob->A, &prob->y, CVX_TRANS);
        //cvx_mat_printf(stdout, "%e", &prob->z, "hresx z");
        cvx_sgemv(1.0, &prob->hrx, -1.0, prob->G, &prob->z_g, CVX_TRANS);
        //cvxm_mvmult(1.0, &prob->hrx, -1.0, prob->G, &prob->z, CVX_TRANS);
        //stats->hresx = cvxm_nrm2(&prob->hrx);
        stats->hresx = SQRT(cvxm_dot(&prob->hrx, &prob->hrx));
        //printf("hresx: %e\n", stats->hresx);
        //cvx_mat_printf(stdout, "%e", &prob->hrx, "hrx");
        
        // rx  = hrx - c*tau
        //     = -A'*y - G'*z - c*tau
        cvxm_copy(&prob->rx, &prob->hrx, 0);
        cvxm_axpy(&prob->rx, -prob->tau, prob->c);
        stats->resx = cvxm_nrm2(&prob->rx) / prob->tau;
        //printf("resx : %e\n", stats->resx);
        //cvx_mat_printf(stdout, "%e", &prob->rx, "rx");

        // hry = A*x
        cvxm_scale(&prob->hry, 0.0, CVX_ALL);
        cvxm_mvmult(0.0, &prob->hry, 1.0, prob->A, &prob->x, 0);
        stats->hresy = cvxm_nrm2(&prob->hry);
        //printf("hresy: %e\n", stats->hresy);

        // ry  = hry - b*tau
        //     = A*x - b*tau
        cvxm_copy(&prob->ry, &prob->hry, 0);
        cvxm_axpy(&prob->ry, -prob->tau, prob->b);
        stats->resy = cvxm_nrm2(&prob->ry) / prob->tau;
        //printf("resy : %e\n", stats->resy);

        // hrz = s + G*x
        cvxm_scale(&prob->hrz, 0.0, CVX_ALL);
        cvxm_mvmult(0.0, &prob->hrz, 1.0, prob->G, &prob->x, 0);
        cvxm_axpy(&prob->hrz, 1.0, &prob->s);
        stats->hresz = cvx_snrm2(&prob->hrz_g);
        //printf("hresz: %e\n", stats->hresz);

        //  rz = hrz = h*tau
        //     = s + G*x - h*tau
        cvxm_copy(&prob->rz, &prob->hrz, 0);
        cvxm_axpy(&prob->rz, -prob->tau, prob->h);

        if (iter > 0) {
            //cvx_mat_printf(stdout, "%e", &prob->hrz, "hrz (iter > 0)");
            //cvx_mat_printf(stdout, "%e", &prob->rz, "rz (iter > 0)");
        }

        //stats->resz = cvx_snrm2(&prob->rz_g) / prob->tau;
        cvx_float_t tmp = cvx_snrm2(&prob->rz_g);
        stats->resz = tmp / prob->tau;
        //printf("resz : %e [%e/%e]\n", stats->resz, tmp, prob->tau);

        // rt = kappa + c'*x + b'*y + h'*z '
        stats->cx = cvxm_dot(prob->c, &prob->x);
        stats->by = cvxm_dot(prob->b, &prob->y);
        stats->hz = cvx_sdot(&prob->h_g, &prob->z_g);
        stats->rt = prob->kappa + stats->cx + stats->by + stats->hz;
        //printf("rt   : %e\n", stats->rt);
        //printf("cx   : %e\n", stats->cx);
        //printf("by   : %e\n", stats->by);
        //printf("hz   : %e\n", stats->hz);

        // statistics for stopping
        stats->pcost = stats->cx / prob->tau;
        stats->dcost = -(stats->by + stats->hz) / prob->tau;
        if (stats->pcost < 0.0) {
            stats->relgap = stats->gap / -stats->pcost;
        } else if (stats->dcost > 0.0) {
            stats->relgap = stats->gap / stats->dcost;
        } else {
            stats->relgap = __NaN();
        }

        //printf("pres: max(%e/%e, %e/%e)\n", stats->resy, stats->resy0, stats->resz, stats->resz0);
        //printf("dres: %e/%e\n", stats->resx, stats->resx0);
        
        stats->pres = __MAX2((stats->resy/stats->resy0), (stats->resz/stats->resz0));
        stats->dres = stats->resx / stats->resx0;
        stats->pinfres = __NaN();
        //printf("pinfres.0 %e\n", stats->pinfres);
        if (stats->hz + stats->by < 0.0) {
            stats->pinfres = stats->hresx / stats->resx0 / (-stats->hz - stats->by);
        }
        stats->dinfres = __NaN();
        if (stats->cx < 0.0) {
            stats->dinfres = __MAX2((stats->hresy/stats->resy0),  (stats->hresz/stats->resz0)) / (-stats->cx);
        }
        //printf("pinfres: %e, dinfres: %e\n", stats->pinfres, stats->dinfres);
        if (opts->show_progress > 0) {
            if (iter == 0) {
                fprintf(stderr, "%10s %12s %10s %8s %7s %5s\n",
                    "pcost", "dcost", "gap", "pres", "dres", "k/t");
            }
            fprintf(stderr, "%2d: %11.4e %11.4e %4.0e %7.0e %7.0e %7.0e\n",
                    iter, stats->pcost, stats->dcost, stats->gap, stats->pres,
                    stats->dres, prob->kappa/prob->tau);
        }
        // ---------------------------------------------------------------------
        // test for stopping criteria

        if (iter == maxiter) {
            return cvx_conelp_ready(prob, stats, iter, CVX_STAT_UNKNOWN);

        }
        else  if (stats->pres <= feastol &&
                    stats->dres <= feastol &&
                    (stats->gap <= abstol ||
                     (!isnan(stats->relgap) && stats->relgap < reltol))) {

            return cvx_conelp_ready(prob, stats, iter, CVX_STAT_OPTIMAL);
        }
        else if (! isnan(stats->pinfres) && stats->pinfres <= feastol) {
            return cvx_conelp_ready(prob, stats, iter, CVX_STAT_PRIMAL_INFEASIBLE);
        }
        else if (! isnan(stats->dinfres) && stats->dinfres <= feastol) {
            return cvx_conelp_ready(prob, stats, iter, CVX_STAT_DUAL_INFEASIBLE);
        }

        // -----------------------------------------------------------------------
        // Compute initial scaling W:
        // 
        //     W * z = W^{-T} * s = lambda
        //     dg * tau = 1/dg * kappa = lambdag.
        if (iter == 0) {
            cvx_compute_scaling(&prob->W, &prob->s_g, &prob->z_g, &prob->lmbda_g, &prob->work);
            //     dg = sqrt( kappa / tau )
            //     dgi = sqrt( tau / kappa )
            //     lambda_g = sqrt( tau * kappa )  
            // 
            // lambda_g is stored in the last position of lmbda.
            prob->dg  = SQRT(prob->kappa / prob->tau);
            prob->dgi = SQRT(prob->tau / prob->kappa);
            cvxm_set(&prob->lmbda, prob->cdim_diag, 0, SQRT(prob->tau * prob->kappa));
            //cvx_mat_printf(stdout, "%e", &prob->lmbda, "lmbda");
            //printf("dg: %e, dgi: %e\n", prob->dg, prob->dgi);
        }
        
        // lmdasq = lmda o lmbda
        cvx_ssqr(&prob->lmbdasq_g, &prob->lmbda_g/*, 0*/);
        cvx_float_t lambda_g0 = cvxm_get(&prob->lmbda, prob->cdim_diag, 0);
        cvxm_set(&prob->lmbdasq, prob->cdim_diag, 0, lambda_g0 * lambda_g0);
        //cvx_mat_printf(stdout, "%e", &prob->lmbdasq, "lmbdasq");

        // Solve
        //
        //     [ 0   A'  G'    ] [ x1        ]          [ c ]
        //     [-A   0   0     ]*[ y1        ] = -dgi * [ b ].
        //     [-G   0   W'*W  ] [ W^{-1}*z1 ]          [ h ]

        cvx_kktfactor(prob->solver, &prob->W, __cvxnil, __cvxnil);

        cvxm_copy(&prob->x1, prob->c, CVX_ALL);
        cvxm_scale(&prob->x1, -1.0, CVX_ALL);
        cvxm_copy(&prob->y1, prob->b, CVX_ALL);
        cvxm_copy(&prob->z1, prob->h, CVX_ALL);
        
        int err = cvx_kktsolve(prob->solver, &prob->x1, &prob->y1, &prob->z1_g);
        //cvx_mat_printf(stdout, "%e", &prob->x1, "solved x");
        //cvx_mat_printf(stdout, "%e", &prob->z1, "solved z");

        cvxm_scale(&prob->x1, prob->dgi, CVX_ALL);
        cvxm_scale(&prob->y1, prob->dgi, CVX_ALL);
        cvxm_scale(&prob->z1, prob->dgi, CVX_ALL);

        if (err < 0) {
            if (iter == 0 && ! primalstart && ! dualstart) {
                prob->error = CVX_ERR_RANK;
                return -1;
            } else {
                prob->error = CVX_ERR_SINGULAR;
                return cvx_conelp_ready(prob, stats, iter, CVX_STAT_SINGULAR);
            }
        }

        cvxm_copy(&prob->th, prob->h, CVX_ALL);
        cvx_scale(&prob->th_g, &prob->W, CVX_TRANS|CVX_INV, &prob->work);
        //cvx_mat_printf(stdout, "%e", &prob->th, "scaled th");

        cvx_float_t nrm, mu, sigma, step, tt, tk;
        
        // TODO: check this!!!
        nrm = cvxm_nrm2(&prob->lmbda);
        mu  = (nrm * nrm) / (1.0 + (cvx_float_t)prob->cdim_diag);
        sigma = 0.0;
        
        for (int i = 0; i < 2; i++) {
            // Solve
            // 
            // [ 0         ]   [  0   A'  G'  c ] [ dx        ]
            // [ 0         ]   [ -A   0   0   b ] [ dy        ]
            // [ W'*ds     ] - [ -G   0   0   h ] [ W^{-1}*dz ]
            // [ dg*dkappa ]   [ -c' -b' -h'  0 ] [ dtau/dg   ]
            // 
            //               [ rx   ]
            //               [ ry   ]
            // = - (1-sigma) [ rz   ]
            //               [ rtau ]
            // 
            // lmbda o (dz + ds) = -lmbda o lmbda + sigma*mu*e
            // lmbdag * (dtau + dkappa) = - kappa * tau + sigma*mu
            // 
            // ds = -lmbdasq                          if i == 0
            //    = -lmbdasq - dsa o dza + sigma*mu*e if i == 1
            // dkappa = -lambdasq[-1]                            if i == 0 
            //        = -lambdasq[-1] - dkappaa*dtaua + sigma*mu if i == 1.
            
            cvx_mgrp_copy_lambda(&prob->ds_g, &prob->lmbdasq_g);
            prob->dkappa = cvxm_get(&prob->lmbdasq, prob->cdim_diag, 0);
            //cvx_mat_printf(stdout, "%e", &prob->ds, "ds");
            //printf("dkappa=%e, mu=%e, sigma=%e, refinement=%d, i=%d\n", prob->dkappa, mu, sigma, refinement, i);
            
            if (i == 1) {
                // ws3 holds ds o dz ; wkappa3 holds dtau*dkappa 
                cvxm_axpy(&prob->ds, 1.0, &prob->ws3);

                // update ds with -sigma*mu
                cvx_mgrp_update_sz(&prob->ds_g, -sigma*mu, 0);
                prob->dkappa = prob->dkappa + prob->wkappa3 - sigma*mu;
                //cvx_mat_printf(stdout, "%e", &prob->ds, "from saved ws3->ds (i=1)");
            }
            // (dx, dy, dz, dtau) = (1-sigma)*(rx, ry, rz, rt)
            cvxm_copy(&prob->dx, &prob->rx, CVX_ALL);
            cvxm_scale(&prob->dx, 1.0 - sigma, CVX_ALL);
            cvxm_copy(&prob->dy, &prob->ry, CVX_ALL);
            cvxm_scale(&prob->dy, 1.0 - sigma, CVX_ALL);
            cvxm_copy(&prob->dz, &prob->rz, CVX_ALL);
            cvxm_scale(&prob->dz, 1.0 - sigma, CVX_ALL);
            prob->dtau = (1.0 - sigma) * stats->rt;
            if (i == 1) {
                //cvx_mat_printf(stdout, "%e", &prob->rx, "rx (i=1)");
                //cvx_mat_printf(stdout, "%e", &prob->dx, "dx (i=1)");
                //cvx_mat_printf(stdout, "%e", &prob->rz, "rz (i=1)");
                //cvx_mat_printf(stdout, "%e", &prob->dz, "dz (i=1)");
                //cvx_mat_printf(stdout, "%e", &prob->ds, "ds (i=1)");
            }
            if (i == 0) {
                //cvx_mat_printf(stdout, "%e", &prob->dx, "dx (i=0)");
            }
            f6(prob, &prob->dx, &prob->dy,
               &prob->dz_g, &prob->dtau, &prob->ds_g, &prob->dkappa, refinement);
            //cvx_mat_printf(stdout, "%e", &prob->dx, "f6(dx)");
            //printf("f6: dtau=%e, dkappa=%e\n", prob->dtau, prob->dkappa);
            if (i == 0) {
                // Save ds o dz and dkappa * dtau for Mehrotra correction
                cvxm_copy(&prob->ws3, &prob->ds, CVX_ALL);
                cvx_sprod(&prob->ws3_g, &prob->dz_g, 0, &prob->work);
                prob->wkappa3 = prob->dtau * prob->dkappa;
                //cvx_mat_printf(stdout, "%e", &prob->ws3, "saved ds->ws3 (i=0)");
            }

            // Maximum step to boundary.
            // 
            // If i is 1, also compute eigenvalue decomposition of the 's' 
            // blocks in ds, dz.  The eigenvectors Qs, Qz are stored in 
            // ds_k, dz_k.  The eigenvalues are stored in sigs, sigz.

            cvx_scale2(&prob->ds_g, &prob->lmbda_g, 0, &prob->work);
            cvx_scale2(&prob->dz_g, &prob->lmbda_g, 0, &prob->work);
            //cvx_mat_printf(stdout, "%e", &prob->dz, "inloop scale2(dz)");
            //cvx_mat_printf(stdout, "%e", &prob->ds, "inloop scale2(ds)");
            //cvx_mat_printf(stdout, "%e", &prob->lmbda, "lmbda");

            if (i == 0) {
                prob->ts = cvx_max_step(&prob->ds_g, __nilgrp, &prob->work);
                prob->tz = cvx_max_step(&prob->dz_g, __nilgrp, &prob->work);
            } else {
                prob->ts = cvx_max_step(&prob->ds_g, &prob->sigs_g, &prob->work);
                prob->tz = cvx_max_step(&prob->dz_g, &prob->sigz_g, &prob->work);
            }
            
            tt = -prob->dtau / cvxm_get(&prob->lmbda, prob->cdim_diag, 0);
            tk = -prob->dkappa / cvxm_get(&prob->lmbda, prob->cdim_diag, 0);
            cvx_float_t t = __maxvec(5, (cvx_float_t[]){0.0, prob->ts, prob->tz, tt, tk});
            //printf("ts=%e, tz=%e, tt=%e, tk=%e, t=%e\n", prob->ts, prob->tz, tt, tk, t);
            //cvx_mat_printf(stdout, "%e", &prob->dz, "inloop maxstep dz");
            if (t == 0.0) {
                step = 1.0;
            } else {
                if (i == 0) {
                    step = __minvec(2, (cvx_float_t []){1.0, 1.0/t});
                } else {
                    step = __minvec(2, (cvx_float_t []){1.0, STEP/t});
                }
            }
            if (i == 0) {
                // sigma = (1 -step)^3
                sigma = (1.0 - step) * (1.0 - step) * (1.0 - step);
            }
            //printf("sigma=%e, step=%e\n", sigma, step);
            //cvx_mat_printf(stdout, "%e", &prob->dz, "loop end dz");
            //cvx_mat_printf(stdout, "%e", &prob->ds, "loop end ds");
        }
        //printf("** tau = %e, kappa = %e\n", prob->tau, prob->kappa);
        //printf("** step = %e, sigma = %e\n", step, sigma);
        //cvx_mat_printf(stdout, "%e", &prob->dz, "after dz");
        //cvx_mat_printf(stdout, "%e", &prob->ds, "after ds");

        // update x, y
        cvxm_axpy(&prob->x, step, &prob->dx);
        cvxm_axpy(&prob->y, step, &prob->dy);
        //cvx_mat_printf(stdout, "%e", &prob->x, "update x");

        // Replace 'l' and 'q' blocks of ds and dz with the updated 
        // variables in the current scaling.
        // Replace 's' blocks of ds and dz with the factors Ls, Lz in a 
        // factorization Ls*Ls', Lz*Lz' of the updated variables in the 
        // current scaling.
        //
        // ds := e + step*ds for 'l' and 'q' blocks.
        // dz := e + step*dz for 'l' and 'q' blocks.
        cvx_mgrp_scale_sz(&prob->ds_g, step, CVXDIM_LINEAR|CVXDIM_SOCP);
        cvx_mgrp_update_sz(&prob->ds_g, 1.0, CVXDIM_LINEAR|CVXDIM_SOCP);

        cvx_mgrp_scale_sz(&prob->dz_g, step, CVXDIM_LINEAR|CVXDIM_SOCP);
        cvx_mgrp_update_sz(&prob->dz_g, 1.0, CVXDIM_LINEAR|CVXDIM_SOCP);

        //cvx_mat_printf(stdout, "%e", &prob->dz, "update dz");
        //cvx_mat_printf(stdout, "%e", &prob->ds, "update ds");
        // ds := H(lambda)^{-1/2} * ds and dz := H(lambda)^{-1/2} * dz.
        // 
        // This replaces the 'l' and 'q' components of ds and dz with the
        // updated variables in the current scaling.  
        // The 's' components of ds and dz are replaced with 
        // 
        // diag(lmbda_k)^{1/2} * Qs * diag(lmbda_k)^{1/2} 
        // diag(lmbda_k)^{1/2} * Qz * diag(lmbda_k)^{1/2} 
        //cvx_mat_printf(stdout, "%e", &prob->dz, "pre-scale2 dz");

        cvx_scale2(&prob->ds_g, &prob->lmbda_g, CVX_INV, &prob->work);
        cvx_scale2(&prob->dz_g, &prob->lmbda_g, CVX_INV, &prob->work);

        //cvx_mat_printf(stdout, "%e", &prob->dz, "scale2(dz)");
        //cvx_mat_printf(stdout, "%e", &prob->ds, "scale2(ds)");

        // sigs := ( e + step*sigs ) ./ lambda for 's' blocks.
        // sigz := ( e + step*sigz ) ./ lambda for 's' blocks.
        cvxm_scale(&prob->sigs, step, CVX_ALL);
        cvxm_scale(&prob->sigz, step, CVX_ALL);
        cvxm_add(&prob->sigs, 1.0, CVX_ALL);
        cvxm_add(&prob->sigz, 1.0, CVX_ALL);
        //cvx_mat_printf(stdout, "%e", &prob->sigs, "scaled sigs");
        //cvx_mat_printf(stdout, "%e", &prob->sigz, "scaled sigz");

        for (int k = 0; k < cvx_mgrp_count(&prob->lmbda_g, CVXDIM_SDP); k++) {
            cvx_mgrp_elem(&lk, &prob->lmbda_g, CVXDIM_SDP, k);
            // sigs ./ lmbda
            cvx_mgrp_elem(&sk, &prob->sigs_g, CVXDIM_SDP, k);
            cvxm_solve_diag(&sk, 1.0, &lk, CVX_RIGHT);
            // sigz ./ lmbda
            cvx_mgrp_elem(&sk, &prob->sigz_g, CVXDIM_SDP, k);
            cvxm_solve_diag(&sk, 1.0, &lk, CVX_RIGHT);
        }
        //cvx_mat_printf(stdout, "%e", &prob->sigs, "scaled sigs");
        //cvx_mat_printf(stdout, "%e", &prob->sigz, "scaled sigz");
        
        // TODO: missing somethings..!!!!
        // divide ds, dz by lmbda S blocks;; 
        for (int k = 0; k < cvx_mgrp_count(&prob->ds_g, CVXDIM_SDP); k++) {
            cvx_float_t a;
            int m = cvx_mgrp_elem(&sk, &prob->ds_g, CVXDIM_SDP, k);
            cvx_mgrp_elem(&zk, &prob->dz_g, CVXDIM_SDP, k);
            cvx_mgrp_elem(&ls, &prob->sigs_g, CVXDIM_SDP, k);
            cvx_mgrp_elem(&lz, &prob->sigz_g, CVXDIM_SDP, k);
                      
            for (int i = 0; i < m; i++) {
                cvxm_view_map(&lk, &sk, 0, i, m, 1);
                a = SQRT(cvxm_get(&ls, i, 0));
                cvxm_scale(&lk, a, 0);

                cvxm_view_map(&lk, &zk, 0, i, m, 1);
                a = SQRT(cvxm_get(&lz, i, 0));
                cvxm_scale(&lk, a, 0);
            }
        }

        //cvx_mat_printf(stdout, "%e", &prob->lmbda, "lambda");
        //cvx_mat_printf(stdout, "%e", &prob->lmbda, "pre scaling lambda");
        //cvx_mat_printf(stdout, "%e", &prob->ds, "pre scaling ds");
        //cvx_mat_printf(stdout, "%e", &prob->dz, "pre scaling dz");
        cvx_update_scaling(&prob->W, &prob->lmbda_g, &prob->ds_g, &prob->dz_g, &prob->work);
        //cvx_scaling_printf(stdout, "%e", &prob->W, "W");
        //cvx_mat_printf(stdout, "%e", &prob->lmbda, "scaled lambda");
        
        // For kappa, tau block: 
        //
        //     dg := sqrt( (kappa + step*dkappa) / (tau + step*dtau) ) 
        //         = dg * sqrt( (1 - step*tk) / (1 - step*tt) )
        //
        //     lmbda[-1] := sqrt((tau + step*dtau) * (kappa + step*dkappa))
        //                = lmbda[-1] * sqrt(( 1 - step*tt) * (1 - step*tk))
        prob->dg *= SQRT(1.0 - step*tk)/SQRT(1.0-step*tt);
        prob->dgi = 1.0/prob->dg;
        cvx_float_t a = SQRT(1.0 - step*tk) * SQRT(1.0-step*tt);
        cvx_float_t g = cvxm_get(&prob->lmbda, prob->cdim_diag, 0);
        cvxm_set(&prob->lmbda, prob->cdim_diag, 0, a*g);
        //printf("dg=%e, dgi=%e\n", prob->dg, prob->dgi);
        //cvx_mat_printf(stdout, "%e", &prob->lmbda, "updated lambda");
        
        // Unscale s, z, tau, kappa (unscaled variables are used only to 
        // compute feasibility residuals).
        cvx_mgrp_copy_lambda(&prob->s_g, &prob->lmbda_g);
        cvx_scale(&prob->s_g, &prob->W, CVX_TRANS, &prob->work);
        //cvx_mat_printf(stdout, "%e", &prob->s, "unscaled s");

        cvx_mgrp_copy_lambda(&prob->z_g, &prob->lmbda_g);
        cvx_scale(&prob->z_g, &prob->W, CVX_INV, &prob->work);
        //cvx_mat_printf(stdout, "%e", &prob->z, "unscaled z");

        prob->kappa = cvxm_get(&prob->lmbda, prob->cdim_diag, 0) / prob->dgi;
        prob->tau   = cvxm_get(&prob->lmbda, prob->cdim_diag, 0) * prob->dgi;
        cvxm_map_data(&lk, prob->cdim_diag, 1, cvxm_data(&prob->lmbda, 0));
        g = cvxm_nrm2(&lk) / prob->tau;
        stats->gap = g*g;
        //printf(" ** kappa=%e, tau=%e, gap=%e, sigma=%e\n", prob->kappa, prob->tau, stats->gap, sigma);

    }

    // we never reach here; TODO: fail if we do
    return 0;
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
