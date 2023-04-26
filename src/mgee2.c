/* C source code for the corrected GEE2 for ordinal data with */
/* misclassification in both the response and the covariate. */
/* It treats the error parameters as known. */
/* The uncertainy in estimated error parameters are */
/* accounted for in the R function mgee2() */

/* --- date: 2010/03/06 */

#include "mgee2.h"



void Cgetmgee2v_i(
            double *S_DM_i,
            double *S_Ytilde_i,
            double *S_assocDM_i,
            double *S_Ztilde_i,
            double *S_Ytilde_i_Mat_ext,
            double *S_invPast_ij_all,
            double *S_P_ij_all,
            double *S_L_ij,
            double *S_wgt_i,
            double *S_X_i_ENUM,
            double *S_XDM_ENUM,
            double *S_Xtilde_i_Mat_ext,
            double *S_invGast_ij_all,
            double *S_G_ij_all,
            double *S_L_x_ij,
            double *S_delta_i,
            int *S_m,
            int *S_K,
            int *S_K_x,
            int *S_p_b,
            int *S_p_a,
            int *S_p_g,
            int *S_p_vp,
            double *S_beta,
            double *S_alpha,
            double *S_U_i,
            double *S_M_i,
            double *S_Lambda_g_i,
            double *S_Lambda_vp_i
            )
{
    /* ---------  MAIN DECLS ---------- */
    /* in maindecls.c */

    MATRIX *Ytilde_i, *Ytilde_i_Mat_ext, *DM_i, *Ztilde_i, *assocDM_i;
    MATRIX *invPast_ij_all, *P_ij_all, *L_ij;
    MATRIX *X_i_ENUM, *Xtilde_i_Mat_ext, *wgt_i, *XDM_ENUM;
    MATRIX *invGast_ij_all, *G_ij_all, *L_x_ij;
    MATRIX *delta_i;

    MATRIX *beta, *alpha, *U_i, *M_i, *Lambda_g_i, *Lambda_vp_i;
    MATRIX *dpXt_i_dvpT_l, *dYtil_i_dgT, *dZtil_i_dgT;


    int m, K, K_x, p_a, p_b, p_t, p_g, p_g_1, p_vp, p_vp_1;
    int nlevs_X, nlevs_Y;
    int l;
    int rr, cc;
    int one;
    int nENUM;
    int nrow_assocDM_i;
    int *one_ptr;

    MATRIX *lambda_i, *mu_i;
    MATRIX *zeta_i, *xi_i;
    MATRIX *dlambda_i_dbetaT, *dmu_i_dbetaT,
        *dzeta_i_dbetaT, *dzeta_i_dalphaT,
        *dxi_i_dalphaT, *dxi_i_dbetaT;

    MATRIX  *D_1i,*D_2i, *V_1i, *invV_1i,
        *invV_2i, *U_1i_l, *U_2i_l, *U_i_l,
        *U_1i,  *U_2i,
        *M_1i_l, *M_2i_l, *M_1i, *M_2i,
        *M_21i_l, *M_21i,
        *Lambda_g_i_l_1, *Lambda_g_i_l_2, *Lambda_g_i_l;

    /* MATRIX *printmat;
    int kk; */

    /* ---------  Initialize data -------- */
    m = *S_m;
    K = *S_K;
    K_x = *S_K_x;
    nlevs_X = K_x + 1;
    nlevs_Y = K + 1;
    p_a = *S_p_a;
    p_b = *S_p_b;
    p_t = p_a + p_b;
    p_g = *S_p_g;
    p_g_1 = (int) p_g/(K * (K + 1));
    p_vp = *S_p_vp;
    p_vp_1 = (int) p_vp/(K_x * (K_x + 1));
    one = 1;
    nENUM = pow((K_x + 1), m);
    nrow_assocDM_i = (int) (pow(K, 2) * (m * (m - 1))/2);

    one_ptr = &one;

    DM_i = form_matrix(S_DM_i, m * K, p_b, PERMANENT);
    Ytilde_i = form_matrix(S_Ytilde_i, m * K, 1, PERMANENT);
    Ztilde_i = form_matrix(S_Ztilde_i, nrow_assocDM_i, 1, PERMANENT);
    Ytilde_i_Mat_ext = form_matrix(S_Ytilde_i_Mat_ext, m, nlevs_Y, PERMANENT);
    assocDM_i = form_matrix(S_assocDM_i, nrow_assocDM_i, p_a, PERMANENT);
    invPast_ij_all = form_matrix(S_invPast_ij_all, K, K * m, PERMANENT);
    P_ij_all = form_matrix(S_P_ij_all, K + 1, (K + 1) * m, PERMANENT);
    L_ij = form_matrix(S_L_ij, p_g_1, 1, PERMANENT);
    wgt_i = form_matrix(S_wgt_i, 1, nENUM, PERMANENT);
    XDM_ENUM = form_matrix(S_XDM_ENUM, m * K, nENUM * K_x, PERMANENT);
    X_i_ENUM = form_matrix(S_X_i_ENUM, m, nENUM, PERMANENT);
    Xtilde_i_Mat_ext = form_matrix(S_Xtilde_i_Mat_ext, m, nlevs_X, PERMANENT);
    invGast_ij_all = form_matrix(S_invGast_ij_all, K_x, K_x * m, PERMANENT);
    G_ij_all = form_matrix(S_G_ij_all, K_x + 1, (K_x + 1) * m, PERMANENT);
    L_x_ij = form_matrix(S_L_x_ij, p_vp_1, 1, PERMANENT);
    delta_i = form_matrix(S_delta_i, m, 1, PERMANENT);
    beta = form_matrix(S_beta, p_b, 1, PERMANENT);
    alpha = form_matrix(S_alpha, p_a, 1, PERMANENT);
    U_i = form_matrix(S_U_i, p_t, 1, PERMANENT);
    M_i = form_matrix(S_M_i, p_t, p_t, PERMANENT);
    Lambda_g_i = form_matrix(S_Lambda_g_i, p_t, p_g, PERMANENT);
    Lambda_vp_i = form_matrix(S_Lambda_vp_i, p_t, p_vp, PERMANENT);

    lambda_i = VC_GEE_create_matrix(m * K, 1, PERMANENT);
    mu_i = VC_GEE_create_matrix(m * K, 1, PERMANENT);
    dlambda_i_dbetaT = VC_GEE_create_matrix(m * K,p_b,PERMANENT);
    dmu_i_dbetaT = VC_GEE_create_matrix(m * K,p_b,PERMANENT);
    zeta_i = VC_GEE_create_matrix(nrow_assocDM_i,1,PERMANENT);
    dzeta_i_dalphaT = VC_GEE_create_matrix(nrow_assocDM_i,p_a,PERMANENT);
    dzeta_i_dbetaT = VC_GEE_create_matrix(nrow_assocDM_i,p_b,PERMANENT);
    xi_i = VC_GEE_create_matrix(nrow_assocDM_i, 1, PERMANENT);
    dxi_i_dalphaT = VC_GEE_create_matrix(nrow_assocDM_i,p_a,PERMANENT);
    dxi_i_dbetaT = VC_GEE_create_matrix(nrow_assocDM_i,p_b,PERMANENT);

    D_1i = VC_GEE_create_matrix(p_b, m * K, PERMANENT);
    D_2i = VC_GEE_create_matrix(p_a, nrow_assocDM_i, PERMANENT);
    V_1i = VC_GEE_create_matrix(m * K, m * K, PERMANENT);
    invV_1i = VC_GEE_create_matrix(m * K, m * K, PERMANENT);
    invV_2i = VC_GEE_create_matrix(nrow_assocDM_i, nrow_assocDM_i, PERMANENT);
    U_1i_l = VC_GEE_create_matrix(p_b, 1, PERMANENT);
    U_2i_l = VC_GEE_create_matrix(p_a, 1, PERMANENT);
    U_1i = VC_GEE_create_matrix(p_b, 1, PERMANENT);
    U_2i = VC_GEE_create_matrix(p_a, 1, PERMANENT);
    U_i_l = VC_GEE_create_matrix(p_t, 1, PERMANENT);

    M_1i_l = VC_GEE_create_matrix(p_b, p_b, PERMANENT);
    M_2i_l = VC_GEE_create_matrix(p_a, p_a, PERMANENT);
    M_1i = VC_GEE_create_matrix(p_b, p_b, PERMANENT);
    M_2i = VC_GEE_create_matrix(p_a, p_a, PERMANENT);

    M_21i_l = VC_GEE_create_matrix(p_a, p_b, PERMANENT);
    M_21i = VC_GEE_create_matrix(p_a, p_b, PERMANENT);

    Lambda_g_i_l_1 = VC_GEE_create_matrix(p_b, p_g, PERMANENT);
    Lambda_g_i_l_2 = VC_GEE_create_matrix(p_a, p_g, PERMANENT);
    Lambda_g_i_l = VC_GEE_create_matrix(p_t, p_g, PERMANENT);
    dpXt_i_dvpT_l = VC_GEE_create_matrix(1, p_vp, PERMANENT);
    dYtil_i_dgT = get_dYtil_i_dgT(
        Ytilde_i_Mat_ext,
        invPast_ij_all,
        P_ij_all,
        L_ij,
        delta_i,
        m,
        K);
    dZtil_i_dgT = get_dZtil_i_dgT(
        dYtil_i_dgT,
        Ytilde_i,
        m,
        K);
    make_permanent(dYtil_i_dgT);
    make_permanent(dZtil_i_dgT);

    /*  ----------  core loop  ---------- */
        for (l = 0; l < nENUM; l++)
        {
            /* Rprintf("l = %d\n", l); */

                /* Plug each possible XDM_i into DM_i */
                cols_plug(XDM_ENUM,  l * K_x, (l + 1) * K_x - 1,
                    DM_i, K);

                /* First order quantities */
                /* ---------------------- */
                get_lambda_i(lambda_i, DM_i, beta);
                get_dlambda_i_dbetaT(dlambda_i_dbetaT, lambda_i, DM_i);
                get_mu_i(mu_i, lambda_i, K, m);
                get_dmu_i_dbetaT(dmu_i_dbetaT, dlambda_i_dbetaT, K, m);

                /* Quantities related to bivariate cumulatives */
                /*  zeta_i, dzeta_i_dalphaT, and dzeta_i_dbetaT */
                get_bivar_cumuls_i(zeta_i,
                    dzeta_i_dalphaT, dzeta_i_dbetaT,
                    lambda_i, DM_i, assocDM_i, alpha,
                    K, p_b, p_a, m);

                /* Get xi_i, dxi_i_dalphaT, dxi_i_dbetaT, V_1i, invV_2i   */
                /* ------------------------------------------- */

                get_bivar_marginals_i(xi_i, V_1i, invV_2i, dxi_i_dalphaT,
                    dxi_i_dbetaT, dzeta_i_dalphaT, dzeta_i_dbetaT,
                    lambda_i, mu_i, DM_i, zeta_i, assocDM_i, alpha,
                    K, p_b, p_a, m);

                get_mattransp(dmu_i_dbetaT, D_1i);

                cholinv(V_1i, invV_1i);

                get_estfun(D_1i, invV_1i, Ytilde_i, mu_i, U_1i_l);

                scalar_times_matrix(MEL(wgt_i, 0, l), U_1i_l);

                matrix_addto(U_1i_l, U_1i);

                get_mattransp(dxi_i_dalphaT, D_2i);

                get_estfun(D_2i, invV_2i, Ztilde_i, xi_i, U_2i_l);
                scalar_times_matrix(MEL(wgt_i, 0, l), U_2i_l);

                matrix_addto(U_2i_l, U_2i);

                get_dvd(D_1i, invV_1i, dmu_i_dbetaT, M_1i_l);
                scalar_times_matrix(MEL(wgt_i, 0, l), M_1i_l);
                matrix_addto(M_1i_l, M_1i);

                get_dvd(D_2i, invV_2i, dxi_i_dalphaT, M_2i_l);
                scalar_times_matrix(MEL(wgt_i, 0, l), M_2i_l);
                matrix_addto(M_2i_l, M_2i);

                get_dvd(D_2i, invV_2i, dxi_i_dbetaT, M_21i_l);
                scalar_times_matrix(MEL(wgt_i, 0, l), M_21i_l);
                matrix_addto(M_21i_l, M_21i);

                VC_GEE_plug(U_1i_l, U_i_l, 0, 0);
                VC_GEE_plug(U_2i_l, U_i_l, p_b, 0);

                matrix_addto(U_i_l, U_i, 0, 0);

                /* Account for estimated $\eta$ */
                if ( MEL(wgt_i, 0, l) != 0 ) {
                    /* Update Lambda_g_i */
                    get_dvd(D_1i, invV_1i, dYtil_i_dgT, Lambda_g_i_l_1);
                    get_dvd(D_2i, invV_2i, dZtil_i_dgT, Lambda_g_i_l_2);
                    VC_GEE_plug(Lambda_g_i_l_1, Lambda_g_i_l, 0, 0);
                    VC_GEE_plug(Lambda_g_i_l_2, Lambda_g_i_l, p_b, 0);
                    scalar_times_matrix(MEL(wgt_i, 0, l), Lambda_g_i_l);
                    matrix_addto(Lambda_g_i_l, Lambda_g_i);

                    /* Get dpXt_i_dvpT_l */
                    get_dpXt_i_dvpT_l(
                        dpXt_i_dvpT_l,
                        l,
                        X_i_ENUM,
                        wgt_i,
                        Xtilde_i_Mat_ext,
                        invGast_ij_all,
                        G_ij_all,
                        L_x_ij,
                        delta_i,
                        m,
                        K_x,
                        p_vp);

                    /* Update Lambda_vp_i */
                    for (rr = 0; rr < p_t; rr++)
                    {
                        for (cc = 0; cc < p_vp; cc++)
                        {
                            MEL(Lambda_vp_i, rr, cc) += (
                                MEL(U_i_l, rr, 0) *
                                MEL(dpXt_i_dvpT_l, 0, cc) );
                        }
                    }
                } /* end of if weight is not 0 */

        } /* end of for l */


        VC_GEE_plug(M_1i, M_i, 0, 0);
        VC_GEE_plug(M_21i, M_i, p_b, 0);
        VC_GEE_plug(M_2i, M_i, p_b, p_b);

    /* ---------  Return to R/S  --------- */
    to_S(U_i, S_U_i)
    to_S(M_i, S_M_i)
    to_S(Lambda_g_i, S_Lambda_g_i)
    to_S(Lambda_vp_i, S_Lambda_vp_i)

    /* free memory*/
    VC_GEE_destroy_matrix(DM_i);
    VC_GEE_destroy_matrix(Ytilde_i);
    VC_GEE_destroy_matrix(Ytilde_i_Mat_ext);
    VC_GEE_destroy_matrix(Ztilde_i);
    VC_GEE_destroy_matrix(assocDM_i);
    VC_GEE_destroy_matrix(invGast_ij_all);
    VC_GEE_destroy_matrix(G_ij_all);
    VC_GEE_destroy_matrix(L_x_ij);
    VC_GEE_destroy_matrix(X_i_ENUM);
    VC_GEE_destroy_matrix(Xtilde_i_Mat_ext);
    VC_GEE_destroy_matrix(wgt_i);
    VC_GEE_destroy_matrix(XDM_ENUM);
    VC_GEE_destroy_matrix(invGast_ij_all);
    VC_GEE_destroy_matrix(G_ij_all);
    VC_GEE_destroy_matrix(L_x_ij);
    VC_GEE_destroy_matrix(delta_i);
    VC_GEE_destroy_matrix(dpXt_i_dvpT_l);
    VC_GEE_destroy_matrix(dYtil_i_dgT);
    VC_GEE_destroy_matrix(dZtil_i_dgT);
    VC_GEE_destroy_matrix(Lambda_g_i_l_1) ;
    VC_GEE_destroy_matrix(Lambda_g_i_l_2);
    VC_GEE_destroy_matrix(Lambda_g_i_l);


    VC_GEE_destroy_matrix(lambda_i);
    VC_GEE_destroy_matrix(mu_i);
    VC_GEE_destroy_matrix(zeta_i);
    VC_GEE_destroy_matrix(xi_i);
    VC_GEE_destroy_matrix(dlambda_i_dbetaT);
    VC_GEE_destroy_matrix(dmu_i_dbetaT);
    VC_GEE_destroy_matrix(dzeta_i_dbetaT);
    VC_GEE_destroy_matrix(dzeta_i_dalphaT);
    VC_GEE_destroy_matrix(dxi_i_dalphaT);
    VC_GEE_destroy_matrix(dxi_i_dbetaT);

    VC_GEE_destroy_matrix(D_1i);
    VC_GEE_destroy_matrix(D_2i);
    VC_GEE_destroy_matrix(V_1i);
    VC_GEE_destroy_matrix(invV_1i);
    VC_GEE_destroy_matrix(invV_2i);
    VC_GEE_destroy_matrix(U_1i_l);
    VC_GEE_destroy_matrix(U_2i_l);
    VC_GEE_destroy_matrix(U_i_l);
    VC_GEE_destroy_matrix(U_1i);
    VC_GEE_destroy_matrix(U_2i);
    VC_GEE_destroy_matrix(M_1i_l);
    VC_GEE_destroy_matrix(M_2i_l);
    VC_GEE_destroy_matrix(M_1i);
    VC_GEE_destroy_matrix(M_2i);
    VC_GEE_destroy_matrix(M_21i_l);
    VC_GEE_destroy_matrix(M_21i);

    VC_GEE_destroy_matrix(beta);
    VC_GEE_destroy_matrix(alpha);
    VC_GEE_destroy_matrix(U_i);
    VC_GEE_destroy_matrix(M_i);
    VC_GEE_destroy_matrix(Lambda_g_i);
    VC_GEE_destroy_matrix(Lambda_vp_i);
}




void Cgetmgee2_i(double *S_DM_i,
            double *S_Ytilde_i,
            double *S_assocDM_i,
            double *S_Ztilde_i,
            double *S_wgt_i,
            double *S_XDM_ENUM,
            int *S_m,
            int *S_K,
            int *S_K_x,
            int *S_p_b,
            int *S_p_a,
            double *S_beta,
            double *S_alpha,
            double *S_U_i,
            double *S_M_i,
            double *S_Sigma_i
            )
{
    /* ---------  MAIN DECLS ---------- */
    /* in maindecls.c */

    MATRIX *Ytilde_i, *DM_i, *Ztilde_i, *assocDM_i, *wgt_i, *XDM_ENUM;

    MATRIX *beta, *alpha;

    int m, K, K_x, p_a, p_b, p_t;
    int l;
    int one;
    int nENUM,  nrow_DM_i, nXDM_ENUM;
    int nrow_assocDM_i;
    int *one_ptr;


    MATRIX *lambda_i, *mu_i;
    MATRIX *zeta_i, *xi_i;
    MATRIX *dlambda_i_dbetaT, *dmu_i_dbetaT,
        *dzeta_i_dbetaT, *dzeta_i_dalphaT,
        *dxi_i_dalphaT, *dxi_i_dbetaT;

    MATRIX  *D_1i,*D_2i, *V_1i, *invV_1i,
        *invV_2i, *U_1i_l, *U_2i_l,
        *U_1i,  *U_2i,
        *M_1i_l, *M_2i_l, *M_1i, *M_2i,
        *M_21i_l, *M_21i;

    MATRIX *U_i, *M_i, *Sigma_i;

    /* ---------  Initialize data -------- */
    m = *S_m;
    K = *S_K;
    K_x = *S_K_x;
    p_a = *S_p_a;
    p_b = *S_p_b;
    p_t = p_a + p_b;
    one = 1;
    nENUM = pow((K_x + 1), m);
    nXDM_ENUM = nENUM * K_x;
    nrow_DM_i = m * K;
    nrow_assocDM_i = (int) (pow(K, 2) * (m * (m - 1))/2.);

    one_ptr = &one;

    from_S(S_Ytilde_i, &nrow_DM_i, &one, Ytilde_i)
    from_S(S_DM_i, &nrow_DM_i, S_p_b, DM_i)
    from_S(S_Ztilde_i, &nrow_assocDM_i, &one, Ztilde_i)
    from_S(S_assocDM_i, &nrow_assocDM_i, S_p_a, assocDM_i)
    from_S(S_wgt_i, &one, &nENUM, wgt_i)
    from_S(S_XDM_ENUM, &nrow_DM_i, &nXDM_ENUM, XDM_ENUM)
    make_permanent(Ytilde_i);
    make_permanent(DM_i);
    make_permanent(Ztilde_i);
    make_permanent(assocDM_i);
    make_permanent(wgt_i);
    make_permanent(XDM_ENUM);

    from_S(S_beta, S_p_b, one_ptr, beta)
    from_S(S_alpha, S_p_a, one_ptr, alpha)
    from_S(S_U_i, &p_t, one_ptr, U_i)
    from_S(S_M_i, &p_t, &p_t, M_i)
    from_S(S_Sigma_i, &p_t, &p_t, Sigma_i)

    make_permanent(beta);
    make_permanent(alpha);
    make_permanent(U_i);
    make_permanent(M_i);
    make_permanent(Sigma_i);

    /*
    Ytilde_i->nrows = nrow_DM_i;
    Ytilde_i->ncols = 1;
    Ytilde_i->data = S_Ytilde_i;
    Ytilde_i->permanence = PERMANENT; */

    lambda_i = VC_GEE_create_matrix(
        nrow_DM_i, 1, PERMANENT);
    mu_i = VC_GEE_create_matrix(
        nrow_DM_i, 1, PERMANENT);
    dlambda_i_dbetaT = VC_GEE_create_matrix(
        nrow_DM_i,p_b,PERMANENT);
    dmu_i_dbetaT = VC_GEE_create_matrix(
        nrow_DM_i,p_b,PERMANENT);
    zeta_i = VC_GEE_create_matrix(
        nrow_assocDM_i,1,PERMANENT);
    dzeta_i_dalphaT = VC_GEE_create_matrix(
        nrow_assocDM_i,p_a,PERMANENT);
    dzeta_i_dbetaT = VC_GEE_create_matrix(
        nrow_assocDM_i,p_b,PERMANENT);
    xi_i = VC_GEE_create_matrix(
        nrow_assocDM_i, 1, PERMANENT);
    dxi_i_dalphaT = VC_GEE_create_matrix(
        nrow_assocDM_i,p_a,PERMANENT);
    dxi_i_dbetaT = VC_GEE_create_matrix(
        nrow_assocDM_i,p_b,PERMANENT);

    D_1i = VC_GEE_create_matrix(p_b, nrow_DM_i, PERMANENT);
    D_2i = VC_GEE_create_matrix(p_a, nrow_assocDM_i, PERMANENT);
    V_1i = VC_GEE_create_matrix(nrow_DM_i,
        nrow_DM_i, PERMANENT);
    invV_1i = VC_GEE_create_matrix(nrow_DM_i,
        nrow_DM_i, PERMANENT);
    invV_2i = VC_GEE_create_matrix(nrow_assocDM_i,
        nrow_assocDM_i, PERMANENT);
    U_1i_l = VC_GEE_create_matrix(p_b, 1, PERMANENT);
    U_2i_l = VC_GEE_create_matrix(p_a, 1, PERMANENT);
    U_1i = VC_GEE_create_matrix(p_b, 1, PERMANENT);
    U_2i = VC_GEE_create_matrix(p_a, 1, PERMANENT);

    M_1i_l = VC_GEE_create_matrix(p_b, p_b, PERMANENT);
    M_2i_l = VC_GEE_create_matrix(p_a, p_a, PERMANENT);
    M_1i = VC_GEE_create_matrix(p_b, p_b, PERMANENT);
    M_2i = VC_GEE_create_matrix(p_a, p_a, PERMANENT);

    M_21i_l = VC_GEE_create_matrix(p_a, p_b, PERMANENT);
    M_21i = VC_GEE_create_matrix(p_a, p_b, PERMANENT);
    M_i = VC_GEE_create_matrix(p_t, p_t, PERMANENT);

    /*  ----------  core loop  ---------- */
            for (l = 0; l < nENUM; l++)
            {
                /* Plug each possible XDM_i into DM_i */
                cols_plug(XDM_ENUM,  l * K_x, (l + 1) * K_x - 1,
                    DM_i, K);

                /* First order quantities */
                /* ---------------------- */
                get_lambda_i(lambda_i, DM_i, beta);
                get_dlambda_i_dbetaT(dlambda_i_dbetaT, lambda_i, DM_i);
                get_mu_i(mu_i, lambda_i, K, m);
                get_dmu_i_dbetaT(dmu_i_dbetaT, dlambda_i_dbetaT, K, m);

                /* Quantities related to bivariate cumulatives */
                /*  zeta_i, dzeta_i_dalphaT, and dzeta_i_dbetaT */
                get_bivar_cumuls_i(zeta_i,
                    dzeta_i_dalphaT, dzeta_i_dbetaT,
                    lambda_i, DM_i, assocDM_i, alpha,
                    K, p_b, p_a, m);

                /* Get xi_i, dxi_i_dalphaT, dxi_i_dbetaT, V_1i, invV_2i   */
                /* ------------------------------------------- */

                get_bivar_marginals_i(xi_i, V_1i, invV_2i, dxi_i_dalphaT,
                    dxi_i_dbetaT, dzeta_i_dalphaT, dzeta_i_dbetaT,
                    lambda_i, mu_i, DM_i, zeta_i, assocDM_i, alpha,
                    K, p_b, p_a, m);

                get_mattransp(dmu_i_dbetaT, D_1i);

                cholinv(V_1i, invV_1i);

                get_estfun(D_1i, invV_1i, Ytilde_i, mu_i, U_1i_l);

                scalar_times_matrix(MEL(wgt_i, 0, l), U_1i_l);

                matrix_addto(U_1i_l, U_1i);

                get_mattransp(dxi_i_dalphaT, D_2i);

                get_estfun(D_2i, invV_2i, Ztilde_i, xi_i, U_2i_l);
                scalar_times_matrix(MEL(wgt_i, 0, l), U_2i_l);

                matrix_addto(U_2i_l, U_2i);

                get_dvd(D_1i, invV_1i, dmu_i_dbetaT, M_1i_l);
                scalar_times_matrix(MEL(wgt_i, 0, l), M_1i_l);
                matrix_addto(M_1i_l, M_1i);

                get_dvd(D_2i, invV_2i, dxi_i_dalphaT, M_2i_l);
                scalar_times_matrix(MEL(wgt_i, 0, l), M_2i_l);
                matrix_addto(M_2i_l, M_2i);

                get_dvd(D_2i, invV_2i, dxi_i_dbetaT, M_21i_l);
                scalar_times_matrix(MEL(wgt_i, 0, l), M_21i_l);
                matrix_addto(M_21i_l, M_21i);

            } /* end of for l */

            /* Things to be returned to R */

            VC_GEE_plug(U_1i, U_i, 0, 0);
            VC_GEE_plug(U_2i, U_i, p_b, 0);
            VC_GEE_plug(M_1i, M_i, 0, 0);
            VC_GEE_plug(M_21i, M_i, p_b, 0);
            VC_GEE_plug(M_2i, M_i, p_b, p_b);
            Sigma_i = get_outer(U_i, U_i);


    /* free memory*/
    VC_GEE_destroy_matrix(DM_i);
    VC_GEE_destroy_matrix(Ytilde_i);
    VC_GEE_destroy_matrix(Ztilde_i);
    VC_GEE_destroy_matrix(assocDM_i);
    VC_GEE_destroy_matrix(wgt_i);
    VC_GEE_destroy_matrix(XDM_ENUM);


    VC_GEE_destroy_matrix(lambda_i);
    VC_GEE_destroy_matrix(mu_i);
    VC_GEE_destroy_matrix(zeta_i);
    VC_GEE_destroy_matrix(xi_i);
    VC_GEE_destroy_matrix(dlambda_i_dbetaT);
    VC_GEE_destroy_matrix(dmu_i_dbetaT);
    VC_GEE_destroy_matrix(dzeta_i_dbetaT);
    VC_GEE_destroy_matrix(dzeta_i_dalphaT);
    VC_GEE_destroy_matrix(dxi_i_dalphaT);
    VC_GEE_destroy_matrix(dxi_i_dbetaT);

    VC_GEE_destroy_matrix(D_1i);
    VC_GEE_destroy_matrix(D_2i);
    VC_GEE_destroy_matrix(V_1i);
    VC_GEE_destroy_matrix(invV_1i);
    VC_GEE_destroy_matrix(invV_2i);
    VC_GEE_destroy_matrix(U_1i_l);
    VC_GEE_destroy_matrix(U_2i_l);
    VC_GEE_destroy_matrix(U_1i);
    VC_GEE_destroy_matrix(U_2i);
    VC_GEE_destroy_matrix(M_1i_l);
    VC_GEE_destroy_matrix(M_2i_l);
    VC_GEE_destroy_matrix(M_1i);
    VC_GEE_destroy_matrix(M_2i);
    VC_GEE_destroy_matrix(M_21i_l);
    VC_GEE_destroy_matrix(M_21i);
    VC_GEE_destroy_matrix(beta);
    VC_GEE_destroy_matrix(alpha);

    /* ---------  Return to R/S  --------- */
    to_S(U_i, S_U_i)
    to_S(M_i, S_M_i)
    to_S(Sigma_i, S_Sigma_i)

    VC_GEE_destroy_matrix(U_i);
    VC_GEE_destroy_matrix(M_i);
    VC_GEE_destroy_matrix(Sigma_i);
}


void Cmgee2(double *S_DM,
            double *S_Ytilde,
            double *S_wgt,
            double *S_extID,
            double *S_XDM_ENUM,
            double *S_assocDM,
            double *S_Ztilde,
            double *S_assocID,
            int *S_n,
            int *S_m,
            int *S_K,
            int *S_K_x,
            int *S_p_b,
            int *S_p_a,
            double *S_beta,
            double *S_alpha,
            double *S_Gamma,
            double *S_Sigma,
            int *S_convergence,
            int *S_iter,
            double *tol
            )
{
    /* ---------  MAIN DECLS ---------- */
    /* in maindecls.c */

/* MATRIX *DMin, *Ytildein,
        *extIDin, *assocIDin, *assocDMin, *Ztildein,   */
    MATRIX *Ytilde_i, *DM_i, *Ztilde_i, *assocDM_i, *wgt, *XDM_ENUM;
    /* MATRIX **DM, **Ytilde, **Ztilde, **assocDM; */

    MATRIX *beta, *alpha, *theta, *theta_o;

    int n, m, K, K_x, p_a, p_b, p_t;
    int i, l, j, p;
    int one, maxit, iter, last_iter, convergence;
    int nENUM, n_extID, len_Z, nrow_DM_i, nXDM_ENUM;
    int nrow_assocDM_i;
    int *one_ptr;

    double dn, dif;

    MATRIX *lambda_i, *mu_i;
    MATRIX *zeta_i, *xi_i;
    MATRIX *dlambda_i_dbetaT, *dmu_i_dbetaT,
        *dzeta_i_dbetaT, *dzeta_i_dalphaT,
        *dxi_i_dalphaT, *dxi_i_dbetaT;

    MATRIX  *D_1i,*D_2i, *V_1i, *invV_1i,
        *invV_2i, *U_1i_l, *U_2i_l,
        *U_1i,  *U_2i, *U_beta, *U_alpha,
        *M_1i_l, *M_2i_l,
        *M_1, *invM_1, *M_2, *invM_2;

    MATRIX *U_i, *M_21i_l, *M_21,*Sigma, *Gamma;


    /* ---------  Initialize data -------- */
    n = *S_n;
    dn = (double)n;
    m = *S_m;
    K = *S_K;
    K_x = *S_K_x;
    p_a = *S_p_a;
    p_b = *S_p_b;
    p_t = p_a + p_b;
    iter = 0;
    last_iter = 0;
    maxit = *S_iter;
    convergence = *S_convergence;
    one = 1;
    dif = 1.;
    nENUM = pow((K_x + 1), m);
    nXDM_ENUM = nENUM * K_x;
    nrow_DM_i = m * K;
    nrow_assocDM_i = (int) (pow(K, 2) * (m * (m - 1))/2.);
    n_extID = (n * m * K);
    len_Z = (n * nrow_assocDM_i);

    one_ptr = &one;

    from_S(S_beta, S_p_b, one_ptr, beta)
    from_S(S_alpha, S_p_a, one_ptr, alpha)
    from_S(S_wgt, S_n, &nENUM, wgt)
    from_S(S_XDM_ENUM, &nrow_DM_i, &nXDM_ENUM, XDM_ENUM)
    make_permanent(beta);
    make_permanent(alpha);
    make_permanent(wgt);
    make_permanent(XDM_ENUM);

    Ytilde_i = VC_GEE_create_matrix(nrow_DM_i, 1, PERMANENT);
    DM_i = VC_GEE_create_matrix(nrow_DM_i, p_b, PERMANENT);
    Ztilde_i = VC_GEE_create_matrix(nrow_assocDM_i, 1, PERMANENT);
    assocDM_i = VC_GEE_create_matrix(nrow_assocDM_i, p_a, PERMANENT);

    /* Things in the loop */

    theta = get_rbind(beta, alpha);
    make_permanent(theta);
    theta_o = VC_GEE_create_matrix(p_t, 1, PERMANENT);

        U_beta = VC_GEE_create_matrix(p_b, 1, PERMANENT);
        U_alpha= VC_GEE_create_matrix(p_a, 1, PERMANENT);
        M_1 = VC_GEE_create_matrix(p_b,p_b,PERMANENT);
        M_2 = VC_GEE_create_matrix(p_a,p_a,PERMANENT);
            U_1i = VC_GEE_create_matrix(p_b, 1, PERMANENT);
            U_2i = VC_GEE_create_matrix(p_a, 1, PERMANENT);

            lambda_i = VC_GEE_create_matrix(
                    nrow_DM_i, 1, PERMANENT);
            mu_i = VC_GEE_create_matrix(
                    nrow_DM_i, 1, PERMANENT);
            dlambda_i_dbetaT = VC_GEE_create_matrix(
                    nrow_DM_i,p_b,PERMANENT);
            dmu_i_dbetaT = VC_GEE_create_matrix(
                    nrow_DM_i,p_b,PERMANENT);
               zeta_i = VC_GEE_create_matrix(
                    nrow_assocDM_i,1,PERMANENT);
                dzeta_i_dalphaT = VC_GEE_create_matrix(
                    nrow_assocDM_i,p_a,PERMANENT);
                dzeta_i_dbetaT = VC_GEE_create_matrix(
                    nrow_assocDM_i,p_b,PERMANENT);
                xi_i = VC_GEE_create_matrix(
                    nrow_assocDM_i, 1, PERMANENT);
                dxi_i_dalphaT = VC_GEE_create_matrix(
                    nrow_assocDM_i,p_a,PERMANENT);
                dxi_i_dbetaT = VC_GEE_create_matrix(
                    nrow_assocDM_i,p_b,PERMANENT);
                V_1i = VC_GEE_create_matrix(nrow_DM_i,
                       nrow_DM_i, PERMANENT);
                invV_1i = VC_GEE_create_matrix(nrow_DM_i,
                      nrow_DM_i, PERMANENT);
                invV_2i = VC_GEE_create_matrix(nrow_assocDM_i,
                       nrow_assocDM_i, PERMANENT);


    D_1i = VC_GEE_create_matrix(p_b, nrow_DM_i, PERMANENT);
    D_2i = VC_GEE_create_matrix(p_a, nrow_assocDM_i, PERMANENT);
    U_1i_l = VC_GEE_create_matrix(p_b, 1, PERMANENT);
    U_2i_l = VC_GEE_create_matrix(p_a, 1, PERMANENT);
    M_1i_l = VC_GEE_create_matrix(p_b, p_b, PERMANENT);
    M_2i_l = VC_GEE_create_matrix(p_a, p_a, PERMANENT);
    invM_1 = VC_GEE_create_matrix(p_b, p_b, PERMANENT);
    invM_2 = VC_GEE_create_matrix(p_a, p_a, PERMANENT);


    /*  ----------  core loop  ---------- */
    do {

        VC_GEE_plug(beta, theta_o, 0, 0);
        VC_GEE_plug(alpha, theta_o, p_b, 0);
        set_zero(U_beta);
        set_zero(U_alpha);
        set_zero(M_1);
        set_zero(M_2);

        for (i = 0; i < n; i++)
        {
            set_zero(U_1i);
            set_zero(U_2i);
            for (j = 0; j < Ytilde_i->nrows; j++)
            {
                MEL(Ytilde_i, j, 0) = *(S_Ytilde + i * nrow_DM_i + j);
                for (p = 0; p < p_b; p++)
                {
                    MEL(DM_i, j, p) = *(S_DM + (i * nrow_DM_i + j) * p_b + p);
                }
            }
            for (j = 0; j < Ztilde_i->nrows; j++)
            {
                MEL(Ztilde_i, j, 0) = *(S_Ztilde + i * nrow_assocDM_i + j);
                for (p = 0; p < p_b; p++)
                {
                    MEL(assocDM_i, j, p) = *(S_assocDM +
                        (i * nrow_assocDM_i + j) * p_a + p);
                }
            }

            for (l = 0; l < nENUM; l++)
            {
                /* Plug each possible XDM_i into DM_i */
                cols_plug(XDM_ENUM,  l * K_x, (l + 1) * K_x - 1,
                    DM_i, K);

                /* First order quantities */
                /* ---------------------- */
                get_lambda_i(lambda_i, DM_i, beta);
                get_dlambda_i_dbetaT(dlambda_i_dbetaT, lambda_i, DM_i);
                get_mu_i(mu_i, lambda_i, K, m);
                get_dmu_i_dbetaT(dmu_i_dbetaT, dlambda_i_dbetaT, K, m);

                /* Quantities related to bivariate cumulatives */
                /*  zeta_i, dzeta_i_dalphaT, and dzeta_i_dbetaT */
                get_bivar_cumuls_i(zeta_i,
                    dzeta_i_dalphaT, dzeta_i_dbetaT,
                    lambda_i, DM_i, assocDM_i, alpha,
                    K, p_b, p_a, m);

                /* Get xi_i, dxi_i_dalphaT, dxi_i_dbetaT, V_1i, invV_2i   */
                /* ------------------------------------------- */

                get_bivar_marginals_i(xi_i, V_1i, invV_2i, dxi_i_dalphaT,
                    dxi_i_dbetaT, dzeta_i_dalphaT, dzeta_i_dbetaT,
                    lambda_i, mu_i, DM_i, zeta_i, assocDM_i, alpha,
                    K, p_b, p_a, m);

                get_mattransp(dmu_i_dbetaT, D_1i);

                cholinv(V_1i, invV_1i);

                get_estfun(D_1i, invV_1i, Ytilde_i, mu_i, U_1i_l);

                scalar_times_matrix(MEL(wgt, i, l), U_1i_l);

                matrix_addto(U_1i_l, U_1i);

                get_mattransp(dxi_i_dalphaT, D_2i);

                get_estfun(D_2i, invV_2i, Ztilde_i, xi_i, U_2i_l);
                scalar_times_matrix(MEL(wgt, i, l), U_2i_l);

                matrix_addto(U_2i_l, U_2i);

                get_dvd(D_1i, invV_1i, dmu_i_dbetaT, M_1i_l);
                scalar_times_matrix(MEL(wgt, i, l), M_1i_l);
                matrix_addto(M_1i_l, M_1);

                get_dvd(D_2i, invV_2i, dxi_i_dalphaT, M_2i_l);
                scalar_times_matrix(MEL(wgt, i, l), M_2i_l);
                matrix_addto(M_2i_l, M_2);

            } /* end of for l */

            matrix_addto(U_1i, U_beta);
            matrix_addto(U_2i, U_alpha);

        } /* end of for i */

        scalar_times_matrix(1./dn, U_beta);
        scalar_times_matrix(1./dn, U_alpha);
        scalar_times_matrix(1./dn, M_1);
        scalar_times_matrix(1./dn, M_2);

        /* Fisher scoring algorithm */
        cholinv(M_1, invM_1);
        cholinv(M_2, invM_2);

        fisherscoring(0.8, invM_1, U_beta, beta);
        fisherscoring(0.8, invM_2, U_alpha, alpha);

        VC_GEE_plug(beta, theta, 0, 0);
        VC_GEE_plug(alpha, theta, p_b, 0);
        dif = get_max_reldif(theta, theta_o);
        iter++;

    } while ( (dif > *tol) && (iter < maxit) );


    /* Calculate the bread and meat of the sandwich */

    Gamma = VC_GEE_create_matrix(p_t,p_t,PERMANENT);
    M_21 = VC_GEE_create_matrix(p_a,p_b,PERMANENT);
    M_21i_l = VC_GEE_create_matrix(p_a,p_b,PERMANENT);
    U_i = VC_GEE_create_matrix(p_t,1,PERMANENT);
    Sigma = VC_GEE_create_matrix(p_t,p_t,PERMANENT);
        for (i = 0; i < n; i++)
        {
            set_zero(U_1i);
            set_zero(U_2i);
            set_zero(U_i);
            Ytilde_i->data = (S_Ytilde + i * nrow_DM_i);
            DM_i->data = (S_DM + i * nrow_DM_i * p_b);
            Ztilde_i->data = (S_Ztilde + i * nrow_assocDM_i);
            assocDM_i->data = (S_assocDM + i * nrow_assocDM_i * p_a);

            for (l = 0; l < nENUM; l++)
            {
                /* Plug each possible XDM_i into DM_i */
                cols_plug(XDM_ENUM,  l * K_x, (l + 1) * K_x - 1,
                    DM_i, K);

                /* First order quantities */
                /* ---------------------- */
                get_lambda_i(lambda_i, DM_i, beta);
                get_dlambda_i_dbetaT(dlambda_i_dbetaT, lambda_i, DM_i);
                get_mu_i(mu_i, lambda_i, K, m);
                get_dmu_i_dbetaT(dmu_i_dbetaT, dlambda_i_dbetaT, K, m);

                /* Quantities related to bivariate cumulatives */
                /*  zeta_i, dzeta_i_dalphaT, and dzeta_i_dbetaT */
                get_bivar_cumuls_i(zeta_i,
                    dzeta_i_dalphaT, dzeta_i_dbetaT,
                    lambda_i, DM_i, assocDM_i, alpha,
                    K, p_b, p_a, m);

                /* Get xi_i, dxi_i_dalphaT, dxi_i_dbetaT, V_1i, invV_2i   */
                /* ------------------------------------------- */

                get_bivar_marginals_i(xi_i, V_1i, invV_2i, dxi_i_dalphaT,
                    dxi_i_dbetaT, dzeta_i_dalphaT, dzeta_i_dbetaT,
                    lambda_i, mu_i, DM_i, zeta_i, assocDM_i, alpha,
                    K, p_b, p_a, m);

                get_mattransp(dmu_i_dbetaT, D_1i);
                cholinv(V_1i, invV_1i);

                get_estfun(D_1i, invV_1i, Ytilde_i, mu_i, U_1i_l);
                scalar_times_matrix(MEL(wgt, i, l), U_1i_l);
                matrix_addto(U_1i_l, U_1i);

                get_mattransp(dxi_i_dalphaT, D_2i);
                get_estfun(D_2i, invV_2i, Ztilde_i, xi_i, U_2i_l);
                scalar_times_matrix(MEL(wgt, i, l), U_2i_l);
                matrix_addto(U_2i_l, U_2i);

                get_dvd(D_1i, invV_1i, dmu_i_dbetaT, M_1i_l);
                scalar_times_matrix(MEL(wgt, i, l), M_1i_l);
                matrix_addto(M_1i_l, M_1);

                get_dvd(D_2i, invV_2i, dxi_i_dalphaT, M_2i_l);
                scalar_times_matrix(MEL(wgt, i, l), M_2i_l);
                matrix_addto(M_2i_l, M_2);

                get_dvd(D_2i, invV_2i, dxi_i_dbetaT, M_21i_l);
                scalar_times_matrix(MEL(wgt, i, l), M_21i_l);
                matrix_addto(M_21i_l, M_21);
            } /* end of for l */

            VC_GEE_plug(U_1i, U_i, 0, 0);
            VC_GEE_plug(U_2i, U_i, p_b, 0);
            add_outer_colvec_to(U_i, Sigma);

        } /* end of for i */

    VC_GEE_plug(M_1, Gamma, 0, 0);
    VC_GEE_plug(M_21, Gamma, p_b, 0);
    VC_GEE_plug(M_2, Gamma, p_b, p_b);
    scalar_times_matrix(1./dn, Gamma);
    scalar_times_matrix(1./dn, Sigma);

    /* ---------  Return to R/S  --------- */
    to_S(beta, S_beta)
    to_S(alpha, S_alpha)
    to_S(Gamma, S_Gamma)
    to_S(Sigma, S_Sigma)
    *S_convergence = (int) ((dif < *tol) && (iter < maxit));
    *S_iter = iter;

    /* free memory*/
    VC_GEE_destroy_matrix(DM_i);
    VC_GEE_destroy_matrix(Ytilde_i);
    VC_GEE_destroy_matrix(Ztilde_i);
    VC_GEE_destroy_matrix(assocDM_i);
    VC_GEE_destroy_matrix(wgt);
    VC_GEE_destroy_matrix(XDM_ENUM);


    VC_GEE_destroy_matrix(lambda_i);
    VC_GEE_destroy_matrix(mu_i);
    VC_GEE_destroy_matrix(zeta_i);
    VC_GEE_destroy_matrix(xi_i);
    VC_GEE_destroy_matrix(dlambda_i_dbetaT);
    VC_GEE_destroy_matrix(dmu_i_dbetaT);
    VC_GEE_destroy_matrix(dzeta_i_dbetaT);
    VC_GEE_destroy_matrix(dzeta_i_dalphaT);
    VC_GEE_destroy_matrix(dxi_i_dalphaT);
    VC_GEE_destroy_matrix(dxi_i_dbetaT);

    VC_GEE_destroy_matrix(D_1i);
    VC_GEE_destroy_matrix(D_2i);
    VC_GEE_destroy_matrix(V_1i);
    VC_GEE_destroy_matrix(invV_1i);
    VC_GEE_destroy_matrix(invV_2i);
    VC_GEE_destroy_matrix(U_1i_l);
    VC_GEE_destroy_matrix(U_2i_l);
    VC_GEE_destroy_matrix(U_1i);
    VC_GEE_destroy_matrix(U_2i);
    VC_GEE_destroy_matrix(U_beta);
    VC_GEE_destroy_matrix(U_alpha);
    VC_GEE_destroy_matrix(M_1i_l);
    VC_GEE_destroy_matrix(M_2i_l);
    VC_GEE_destroy_matrix(M_1);
    VC_GEE_destroy_matrix(M_2);
    VC_GEE_destroy_matrix(invM_1);
    VC_GEE_destroy_matrix(invM_2);

    VC_GEE_destroy_matrix(M_21i_l);
    VC_GEE_destroy_matrix(M_21);
    VC_GEE_destroy_matrix(U_i);

    VC_GEE_destroy_matrix(beta);
    VC_GEE_destroy_matrix(alpha);
    VC_GEE_destroy_matrix(theta);
    VC_GEE_destroy_matrix(theta_o);
    VC_GEE_destroy_matrix(Sigma);
    VC_GEE_destroy_matrix(Gamma);
} /* end of Cmgee2() */


/* printmat = VC_GEE_matcopy(theta);
Rprintf("nrows: %d\n", ((printmat->nrows)));
Rprintf("ncols: %d\n", ((printmat->ncols)));
for(kk = 0; kk < (printmat->nrows)*(printmat->ncols); kk++) {
Rprintf(" %f\n", (*(printmat->data + kk)));
} */


void Cgetordgee2_i(double *S_DM_i,
            double *S_Y_i,
            double *S_assocDM_i,
            double *S_Z_i,
            int *S_m,
            int *S_K,
            int *S_p_b,
            int *S_p_a,
            double *S_beta,
            double *S_alpha,
            double *S_U_i,
            double *S_M_i,
            double *S_Sigma_i
            )
{
    /* ---------  MAIN DECLS ---------- */
    /* in maindecls.c */

    MATRIX *Y_i, *DM_i, *Z_i, *assocDM_i;

    MATRIX *beta, *alpha;

    int m, K, p_a, p_b, p_t;
    int one;
    int nrow_DM_i;
    int nrow_assocDM_i;
    int *one_ptr;


    MATRIX *lambda_i, *mu_i;
    MATRIX *zeta_i, *xi_i;
    MATRIX *dlambda_i_dbetaT, *dmu_i_dbetaT,
        *dzeta_i_dbetaT, *dzeta_i_dalphaT,
        *dxi_i_dalphaT, *dxi_i_dbetaT;

    MATRIX  *D_1i,*D_2i, *V_1i, *invV_1i,
        *invV_2i,
        *U_1i,  *U_2i,
        *M_1i, *M_2i,
        *M_21i;

    MATRIX *U_i, *M_i, *Sigma_i;



    /* ---------  Initialize data -------- */
    m = *S_m;
    K = *S_K;
    p_a = *S_p_a;
    p_b = *S_p_b;
    p_t = p_a + p_b;
    one = 1;
    nrow_DM_i = m * K;
    //Yuliang: changed __int64 to long long
    nrow_assocDM_i = (int) (pow(K, 2) * ((long long)m * ((long long)m - 1))/2.);

    //Yuliang: check m, K, nrow_DM_i
    /*Rprintf("m: %d, K: %d, p_a: %d, p_b: %d, p_t: %d\n", m,K,p_a,p_b,p_t);
    Rprintf("nrow_DM_i: %d\n", nrow_DM_i);*/
    /*if ( m <0 )
        Rprintf("m<0, m=%d\n",m);
    else if ( K < 0 )
        Rprintf("K<0, K=%d\n",K);
    else if ( p_a < 0 )
        Rprintf("p_a<0, p_a=%d\n",p_a);
    else if (p_b<0)
        Rprintf("p_b<0, p_b=%d\n",p_b);
    else if (nrow_DM_i<0)
        Rprintf("nrow_DM_i<0, nrow_DM_i=%d\n",nrow_DM_i);
    else if (nrow_assocDM_i<0)
        Rprintf("nrow_assocDM_i<0, nrow_assocDM_i=%d\n",nrow_assocDM_i);
    else
        Rprintf("all checked positvie.\n");*/
    one_ptr = &one;

    from_S(S_Y_i, &nrow_DM_i, &one, Y_i)
    from_S(S_DM_i, &nrow_DM_i, S_p_b, DM_i)
    from_S(S_Z_i, &nrow_assocDM_i, &one, Z_i)
    from_S(S_assocDM_i, &nrow_assocDM_i, S_p_a, assocDM_i)
    make_permanent(Y_i);
    make_permanent(DM_i);
    make_permanent(Z_i);
    make_permanent(assocDM_i);

    from_S(S_beta, S_p_b, one_ptr, beta)
    from_S(S_alpha, S_p_a, one_ptr, alpha)
    from_S(S_U_i, &p_t, one_ptr, U_i)
    from_S(S_M_i, &p_t, &p_t, M_i)
    from_S(S_Sigma_i, &p_t, &p_t, Sigma_i)

    make_permanent(beta);
    make_permanent(alpha);
    make_permanent(U_i);
    make_permanent(M_i);
    make_permanent(Sigma_i);

    /*Yuliang : check the value of nrow_DM_i*/
    /*Rprintf("nrow_DM_i: %d\n", nrow_DM_i);
    Rprintf("nrow_assocDM_i: %d\n", nrow_assocDM_i);*/

    lambda_i = VC_GEE_create_matrix(
        nrow_DM_i, 1, PERMANENT);
    mu_i = VC_GEE_create_matrix(
        nrow_DM_i, 1, PERMANENT);
    dlambda_i_dbetaT = VC_GEE_create_matrix(
        nrow_DM_i,p_b,PERMANENT);
    dmu_i_dbetaT = VC_GEE_create_matrix(
        nrow_DM_i,p_b,PERMANENT);
    zeta_i = VC_GEE_create_matrix(
        nrow_assocDM_i,1,PERMANENT);
    dzeta_i_dalphaT = VC_GEE_create_matrix(
        nrow_assocDM_i,p_a,PERMANENT);
    dzeta_i_dbetaT = VC_GEE_create_matrix(
        nrow_assocDM_i,p_b,PERMANENT);
    xi_i = VC_GEE_create_matrix(
        nrow_assocDM_i, 1, PERMANENT);
    dxi_i_dalphaT = VC_GEE_create_matrix(
        nrow_assocDM_i,p_a,PERMANENT);
    dxi_i_dbetaT = VC_GEE_create_matrix(
        nrow_assocDM_i,p_b,PERMANENT);




    D_1i = VC_GEE_create_matrix(p_b, nrow_DM_i, PERMANENT);
    D_2i = VC_GEE_create_matrix(p_a, nrow_assocDM_i, PERMANENT);
    V_1i = VC_GEE_create_matrix(nrow_DM_i,
        nrow_DM_i, PERMANENT);
    invV_1i = VC_GEE_create_matrix(nrow_DM_i,
        nrow_DM_i, PERMANENT);
    invV_2i = VC_GEE_create_matrix(nrow_assocDM_i,
        nrow_assocDM_i, PERMANENT);
    U_1i = VC_GEE_create_matrix(p_b, 1, PERMANENT);
    U_2i = VC_GEE_create_matrix(p_a, 1, PERMANENT);



    M_1i = VC_GEE_create_matrix(p_b, p_b, PERMANENT);
    M_2i = VC_GEE_create_matrix(p_a, p_a, PERMANENT);



    M_21i = VC_GEE_create_matrix(p_a, p_b, PERMANENT);
    M_i = VC_GEE_create_matrix(p_t, p_t, PERMANENT);



    /*  ----------  core loop  ---------- */

                /* First order quantities */
                /* ---------------------- */
                get_lambda_i(lambda_i, DM_i, beta);
                get_dlambda_i_dbetaT(dlambda_i_dbetaT, lambda_i, DM_i);
                get_mu_i(mu_i, lambda_i, K, m);
                get_dmu_i_dbetaT(dmu_i_dbetaT, dlambda_i_dbetaT, K, m);

                /* Quantities related to bivariate cumulatives */
                /*  zeta_i, dzeta_i_dalphaT, and dzeta_i_dbetaT */
                get_bivar_cumuls_i(zeta_i,
                    dzeta_i_dalphaT, dzeta_i_dbetaT,
                    lambda_i, DM_i, assocDM_i, alpha,
                    K, p_b, p_a, m);

                /* Get xi_i, dxi_i_dalphaT, dxi_i_dbetaT, V_1i, invV_2i   */
                /* ------------------------------------------- */

                    //Yuliang: This function is in question in 64bit

                get_bivar_marginals_i(xi_i, V_1i, invV_2i, dxi_i_dalphaT,
                    dxi_i_dbetaT, dzeta_i_dalphaT, dzeta_i_dbetaT,
                    lambda_i, mu_i, DM_i, zeta_i, assocDM_i, alpha,
                    K, p_b, p_a, m);

                //Yuliang
                //Rprintf("get_bivar_marginals_i run successfully in Cgetordgee2_i breakpoint.\n");

                get_mattransp(dmu_i_dbetaT, D_1i);

                cholinv(V_1i, invV_1i);

                get_estfun(D_1i, invV_1i, Y_i, mu_i, U_1i);
                get_mattransp(dxi_i_dalphaT, D_2i);

                get_estfun(D_2i, invV_2i, Z_i, xi_i, U_2i);

                get_dvd(D_1i, invV_1i, dmu_i_dbetaT, M_1i);

                get_dvd(D_2i, invV_2i, dxi_i_dalphaT, M_2i);

                get_dvd(D_2i, invV_2i, dxi_i_dbetaT, M_21i);


            /* Things to be returned to R */

            VC_GEE_plug(U_1i, U_i, 0, 0);
            VC_GEE_plug(U_2i, U_i, p_b, 0);
            VC_GEE_plug(M_1i, M_i, 0, 0);
            VC_GEE_plug(M_21i, M_i, p_b, 0);
            VC_GEE_plug(M_2i, M_i, p_b, p_b);
            Sigma_i = get_outer(U_i, U_i);


    /* free memory*/
    VC_GEE_destroy_matrix(DM_i);
    VC_GEE_destroy_matrix(Y_i);
    VC_GEE_destroy_matrix(Z_i);
    VC_GEE_destroy_matrix(assocDM_i);

    VC_GEE_destroy_matrix(lambda_i);
    VC_GEE_destroy_matrix(mu_i);
    VC_GEE_destroy_matrix(zeta_i);
    VC_GEE_destroy_matrix(xi_i);
    VC_GEE_destroy_matrix(dlambda_i_dbetaT);
    VC_GEE_destroy_matrix(dmu_i_dbetaT);
    VC_GEE_destroy_matrix(dzeta_i_dbetaT);
    VC_GEE_destroy_matrix(dzeta_i_dalphaT);
    VC_GEE_destroy_matrix(dxi_i_dalphaT);
    VC_GEE_destroy_matrix(dxi_i_dbetaT);

    VC_GEE_destroy_matrix(D_1i);
    VC_GEE_destroy_matrix(D_2i);
    VC_GEE_destroy_matrix(V_1i);
    VC_GEE_destroy_matrix(invV_1i);
    VC_GEE_destroy_matrix(invV_2i);

    VC_GEE_destroy_matrix(U_1i);
    VC_GEE_destroy_matrix(U_2i);

    VC_GEE_destroy_matrix(M_1i);
    VC_GEE_destroy_matrix(M_2i);

    VC_GEE_destroy_matrix(M_21i);
    VC_GEE_destroy_matrix(beta);
    VC_GEE_destroy_matrix(alpha);

    /* ---------  Return to R/S  --------- */
    to_S(U_i, S_U_i)
    to_S(M_i, S_M_i)
    to_S(Sigma_i, S_Sigma_i)

    VC_GEE_destroy_matrix(U_i);
    VC_GEE_destroy_matrix(M_i);
    VC_GEE_destroy_matrix(Sigma_i);


}


/* -----------   Support functions  2010/03/08   ----------- */


static MATRIX *VC_GEE_create_matrix(nrows, ncols, permanence)
int nrows, ncols, permanence;
{
    MATRIX *tmp;
    double *head;
    int i;

    tmp = (MATRIX *) calloc (1, sizeof (struct matrix));

    if (tmp == NULL)
    {
	error("VC_GEE_create_matrix: malloc failed %d",
	      sizeof(struct matrix));
    }

    tmp->nrows = nrows;
    tmp->ncols = ncols;
    tmp->permanence = permanence;

    /*Yuliang: check if nrows and ncols are valid numbers*/
    //if (nrows<0 | ncols<0){
    //    error("invalid nrow=%d or ncol=%d",nrows,ncols);
    //}

    tmp->data = (double *) calloc (1,  nrows * ncols * sizeof (double)) ;
    /*tmp->data = Calloc (nrows * ncols, double) ;*/
    if (tmp->data == NULL)
    {
	error("VC_GEE_create_matrix: malloc failed, nrows=%d ncols=%d",
	      nrows, ncols);
    }

    head = tmp->data;
    for (i = 0 ; i < nrows*ncols ; i++)
    {
	*head = (double)0.;
	head++;
    }

    return tmp;
}

static void VC_GEE_destroy_matrix(mat)
MATRIX *mat;
{
    if (mat == (MATRIX *) NULL) return;

    mat->nrows = 0;
    mat->ncols = 0;
    /*if (mat->data != (double *)NULL) free((char *) mat->data);*/
    mat->data = (double *)NULL;
    /*if (mat != (MATRIX *)NULL) free((char *) mat);*/
    mat = (MATRIX *)NULL;
}



static MATRIX *VC_GEE_transp(mat)
MATRIX *mat;
{
    double *telem, *inelem, *tbase;
    int nelem;
    MATRIX *tmp;

    tmp = VC_GEE_create_matrix(mat->ncols, mat->nrows, EPHEMERAL);
    inelem = mat->data;
    tbase = tmp->data;
    telem = tbase;
    for (nelem = 0 ; nelem < (mat->ncols * mat->nrows) ; nelem++)
    {
	*telem = *(inelem++);
	telem += mat->nrows;
	if (nelem % mat->ncols == (mat->ncols)-1)
	    telem = ++tbase;
    }
    if (is_ephemeral(mat)) VC_GEE_destroy_matrix(mat);
    return tmp;

}


static MATRIX *VC_GEE_extract_rows(Source,VC_GEE_start,end)
MATRIX *Source;
int VC_GEE_start, end;
{
    MATRIX *temp;
    int rows_to_get, i, j;

    rows_to_get = end - VC_GEE_start + 1;

    //Yuliang
    //Rprintf("Inside VC_GEE_extract_rows: checkpoint 1. \n");
    //Rprintf("rows_to_get=%d, \n Source->ncols=%d\n",rows_to_get,Source->ncols);

    temp = VC_GEE_create_matrix(rows_to_get,Source->ncols,EPHEMERAL);

    for (i = 0 ; i < rows_to_get ; i++)
    {
	for (j = 0 ; j < Source->ncols ; j++)
	{
	    *(ELREF(temp,i,j)) = *(ELREF(Source,VC_GEE_start,j));
	}
	VC_GEE_start++;
    }
    return temp;
}

static MATRIX *VC_GEE_extract_cols(x, VC_GEE_start, end)
MATRIX *x;
int VC_GEE_start, end;
{
    MATRIX *tmp;
    tmp = VC_GEE_transp(x);
    tmp = VC_GEE_extract_rows(tmp, VC_GEE_start, end);
    tmp = VC_GEE_transp(tmp);
    free_if_ephemeral(x);
    return tmp;
}

static MATRIX *VC_GEE_matcopy(inmat)
MATRIX *inmat;
{
    int i, j;
    MATRIX *outmat;

    outmat = VC_GEE_create_matrix(inmat->nrows,inmat->ncols,EPHEMERAL);
    for (i = 0 ; i < inmat->nrows ; i++)
    {
	for (j = 0 ; j < inmat->ncols ; j++)
	{
	    *(ELREF(outmat,i,j)) = *(ELREF(inmat,i,j));
	}
    }
    return outmat;
}

static int VC_GEE_split(matptr, discptr, matarrptr)
MATRIX *matptr, *discptr, *matarrptr[];
{
    int i, iVC_GEE_start, k, VC_GEE_start, end;
    if (discptr->ncols != 1)
    {
	error("VC_GEE_split: discriminator must be column vec.\nVC_GEE_split: ncols = %d", discptr->ncols);
    }

    k = 0;

    iVC_GEE_start = (int)MEL(discptr, 0, 0);
    VC_GEE_start = 0;
    end = 0;
    for (i = 1 ; i <= discptr->nrows ; i++)
    {
	if (i == discptr->nrows || MEL(discptr, i, 0) != iVC_GEE_start)
	{
	    matarrptr[k] = VC_GEE_matcopy(VC_GEE_extract_rows(matptr, VC_GEE_start, end));
	    make_permanent(matarrptr[k]);
	    k++;
	    VC_GEE_start = end+1;
	    if (i < discptr->nrows)
		iVC_GEE_start = MEL(discptr, i, 0);
	}
	if (VC_GEE_start < discptr->nrows) end++ ;
    }
    return k;
}


static void VC_GEE_plug(VC_GEE_plugm, socket, row, col)
int row, col;
MATRIX *VC_GEE_plugm, *socket;
{
    int pcol, prow;
    double *sockload, *VC_GEE_plughead, *sockrow_VC_GEE_start;
    int i,j;

    pcol = VC_GEE_plugm->ncols;
    prow = VC_GEE_plugm->nrows;

    if (pcol+col > socket->ncols || prow+row > socket->nrows)
    {
	error("M+-: VC_GEE_plug: socket too small");
    }

    sockload = socket->data + col + row*(socket->ncols);
    VC_GEE_plughead = VC_GEE_plugm->data;
    sockrow_VC_GEE_start = sockload;

    for (i = 0 ; i < prow ; i++)
    {
	sockload = sockrow_VC_GEE_start;
	for (j = 0 ; j < pcol ; j++)
	{
	    *(sockload++) = *(VC_GEE_plughead++);
	}
	sockrow_VC_GEE_start += socket->ncols;
    }
    free_if_ephemeral(VC_GEE_plugm);
}

static MATRIX *VC_GEE_form_diag(vec)
MATRIX *vec;
{
    MATRIX *tmp;
    int i, ord;

    ord = vec->nrows;
    tmp = VC_GEE_create_matrix(ord,  ord, EPHEMERAL);
    for (i = 0 ; i < ord ; i++)
	*(ELREF(tmp,i,i)) = MEL(vec,i,0);
    free_if_ephemeral(vec);
    return tmp;
}


#define get_nelem(x) (((x)->nrows) * ((x)->ncols))

static double VC_GEE_elsum(x)
MATRIX *x;
{
    double t=0.;
    double *loc;
    int i, nelem;

    nelem = get_nelem(x);
    loc = x->data;
    for (i = 0 ; i < nelem ; i++)
	t += *(loc++);
    if (is_ephemeral(x)) VC_GEE_destroy_matrix(x);
    return t;
}

static MATRIX *VC_GEE_matabs(x)
MATRIX *x;
{
    double *load, *look;
    MATRIX *tmp;
    int nelem, i;

    nelem = get_nelem(x);
    tmp = VC_GEE_create_matrix(x->nrows, x->ncols, EPHEMERAL);
    load = tmp->data;
    look = x->data;
    for (i = 0 ; i < nelem ; i++)
	*(load++) = fabs(*look++);
    free_if_ephemeral(x);
    return tmp ;
}


static double VC_GEE_matmax(x)
MATRIX *x;
{
    double t;
    double *loc;
    int i, nelem;

    nelem = get_nelem(x);
    loc = x->data;
    t = MEL(x,0,0);
    for (i = 0 ; i < nelem ; i++)
    {
	if (*(loc) > t) t = *(loc);
	loc++;
    }
    free_if_ephemeral(x);
    return t;
}


static MATRIX *VC_GEE_matexp(x)
MATRIX *x;
{
    double *load, *look;
    MATRIX *tmp;
    int nelem, i;

    nelem = get_nelem(x);
    tmp = VC_GEE_create_matrix(x->nrows, x->ncols, EPHEMERAL);
    load = tmp->data;
    look = x->data;
    for (i = 0 ; i < nelem ; i++)
	*(load++) = exp(*look++);
    free_if_ephemeral(x);
    return tmp ;
}


static MATRIX *VC_GEE_matadd(mat1, mat2)
MATRIX *mat1, *mat2;
{
    MATRIX *result;
    double *mat1base, *mat2base, *resbase;
    int i, j;
    if ((mat1->ncols != mat2->ncols) || (mat1->nrows != mat2->nrows))
    {
	error("VC_GEE_matadd: args (%dx%d) + (%dx%d) don't conform.\nfatal error",
	      mat1->nrows, mat1->ncols, mat2->nrows, mat2->ncols);
    }
    result = VC_GEE_create_matrix(mat1->nrows, mat1->ncols, EPHEMERAL);
    resbase = result->data;
    mat1base = mat1->data;
    mat2base = mat2->data;
    for (j = 0 ; j < result->nrows ; j++)
    {
	for (i = 0 ; i < result->ncols ; i++)
	{
	    *resbase = *mat1base + *mat2base ;
	    resbase++ ; mat1base++ ; mat2base++ ;
	}
    }
    if (is_ephemeral(mat1)) VC_GEE_destroy_matrix(mat1);
    if (is_ephemeral(mat2)) VC_GEE_destroy_matrix(mat2);
    return result;
}


static MATRIX *VC_GEE_matsub(mat1, mat2)
MATRIX *mat1, *mat2;
{
    MATRIX *result;
    double *mat1base, *mat2base, *resbase;
    int i, j;
    if ((mat1->ncols != mat2->ncols) || (mat1->nrows != mat2->nrows))
    {
	error("VC_GEE_matsub: args (%dx%d) + (%dx%d) don't conform.\n",
	      mat1->nrows, mat1->ncols, mat2->nrows, mat2->ncols);
    }
    result = VC_GEE_create_matrix(mat1->nrows, mat1->ncols, EPHEMERAL);
    resbase = result->data;
    mat1base = mat1->data;
    mat2base = mat2->data;
    for (j = 0 ; j < result->nrows ; j++)
    {
	for (i = 0 ; i < result->ncols ; i++)
	{
	    *resbase = *mat1base - *mat2base ;
	    resbase++ ; mat1base++ ; mat2base++ ;
	}
    }
    if (is_ephemeral(mat1)) VC_GEE_destroy_matrix(mat1);
    if (is_ephemeral(mat2)) VC_GEE_destroy_matrix(mat2);
    return result;
}


static MATRIX *VC_GEE_matmult(mat1, mat2)
MATRIX *mat1, *mat2;
{
    double *mat1base, *mat1loc, *mat2base, *mat2loc, *resbase;
    MATRIX *result;
    int i, rows, j;

    if (mat1->ncols != mat2->nrows)
    {
	error("VC_GEE_matmult: args (%dx%d) * (%dx%d) don't conform.\n",
	      mat1->nrows, mat1->ncols, mat2->nrows, mat2->ncols);
    }

    result = VC_GEE_create_matrix(mat1->nrows, mat2->ncols, EPHEMERAL);

    resbase = result->data;
    mat1base = mat1->data;
    mat2base = mat2->data;

    for (j = 0 ; j < result->nrows ; j++)
    {
	for (i = 0 ; i < result->ncols ; i++)
	{
	    mat1loc = mat1base;
	    mat2loc = mat2base;
	    for (rows = 0 ; rows < mat2->nrows ; rows++)
	    {
		*resbase += *(mat1loc++) * *mat2loc;
		mat2loc += mat2->ncols;
	    }
	    ++resbase;
	    ++mat2base;
	}
	mat1base += mat1->ncols;
	mat2base = mat2->data;
    }
    if (is_ephemeral(mat1)) VC_GEE_destroy_matrix(mat1);
    if (is_ephemeral(mat2)) VC_GEE_destroy_matrix(mat2);
    return result;
}


static MATRIX *VC_GEE_px1_times_pxq(px1, pxq)
MATRIX *px1, *pxq;
{
    MATRIX *tmp;
    double *load, colel;
    int i, j;

    if (px1->ncols != 1)
    {
	error("M+-: VC_GEE_px1_times_pxq: arg1 not a col-vec");
    }
    if (px1->nrows != pxq->nrows)
    {
	error("M+-: VC_GEE_px1_times_pxq: args not conforming");
    }
    tmp = VC_GEE_matcopy(pxq);
    load = tmp->data;
    for (i = 0 ; i < tmp->nrows ; i++)
    {
	colel = MEL(px1, i, 0);
	for (j = 0 ; j < tmp->ncols ; j++)
	{
	    *load *= colel ;
	    load++ ;
	}
    }
    free_if_ephemeral(px1);
    free_if_ephemeral(pxq);
    return tmp;
}


static MATRIX *VC_GEE_pxq_divby_px1(pxq, px1)
MATRIX *px1, *pxq;
{
    MATRIX *tmp;
    double *load, colel;
    int i, j;
    if (px1->ncols != 1)
    {
	error("M+-: VC_GEE_pxq_divby_px1: arg2 not a col-vec");
    }
    if (px1->nrows != pxq->nrows)
    {
	error("M+-: VC_GEE_pxq_divby_px1: args not conforming");
    }

    tmp = VC_GEE_matcopy(pxq);
    load = tmp->data;
    for ( i = 0 ; i < tmp->nrows ; i++)
    {
	colel = MEL(px1, i, 0);
	for (j = 0 ; j < tmp->ncols ; j++)
	{
	    *load = (*load) / colel ;
	    load++ ;
	}
    }
    free_if_ephemeral(px1);
    free_if_ephemeral(pxq);
    return tmp;
}


static MATRIX *VC_GEE_scalar_times_matrix(a, X)
double a;
MATRIX *X;
{
    MATRIX *tmp;
    double *tbase;
    int i, nelem;
    tmp = VC_GEE_matcopy(X);
    nelem = get_nelem(tmp);
    tbase = tmp->data;
    for (i = 0 ; i < nelem ; i++) {
	*tbase *= a ;
	tbase++ ;
    }
    free_if_ephemeral(X);
    return tmp;
}


static MATRIX *VC_GEE_ident(ord)
int ord;
{
    MATRIX *I;
    int i;

    I = VC_GEE_create_matrix(ord, ord, EPHEMERAL);
    for (i = 0 ; i < ord ; i++)
	*(ELREF(I,i,i)) = (double)1.0;
    return I;
}


static MATRIX *VC_GEE_col_1s(k)
int k;
{
    MATRIX *tmp;
    int i;
    tmp = VC_GEE_create_matrix(k, 1, EPHEMERAL);
    for (i = 0 ; i < k ; i++)
    {
	MEL(tmp,i,0) = 1.;
    }
    return tmp;
}


static int VC_GEE_nchanges(X)
MATRIX *X;
{
    int tmp = 1, iVC_GEE_start, i;

    if (X->ncols != 1)
    {
	error("VC_GEE_nchanges:  must be column VC_GEE_vector; ncols = %d",
	      X->ncols);
    }

    iVC_GEE_start = MEL(X, 0, 0);

    for (i = 1 ; i < X->nrows ; i++)
    {
	if (MEL (X, i, 0) != iVC_GEE_start)
	{
	    tmp++;
	    iVC_GEE_start = MEL (X, i, 0);
	}
    }
    return tmp;
}


static MATRIX *VC_GEE_diag_as_vec(inmat)
MATRIX *inmat;
{
    int i;
    MATRIX *outmat;

    if(inmat->ncols!=inmat->nrows)
    {
	error("M+-: VC_GEE_diag_as_vec: arg is not a square matrix");
    }

    outmat= VC_GEE_create_matrix(inmat->nrows,1,EPHEMERAL);
    for(i= 0;i<inmat->nrows;i++)
    {
	*(ELREF(outmat,i,0))=  *(ELREF(inmat,i,i));
    }
    return outmat;
}


static MATRIX *VC_GEE_matsqrt(x)
MATRIX *x;
{
    int i,j;
    MATRIX *tmp;
    tmp= VC_GEE_matcopy(x);
    for(i= 0;i<x->ncols;i++)
    {
	for(j= 0;j<x->nrows;j++)
	{
	    MEL(tmp,i,j)= sqrt(MEL(x,i,j));
	}
    }
    if(is_ephemeral(x))VC_GEE_destroy_matrix(x);
    return tmp;
}


static MATRIX *VC_GEE_mat1over(x)
MATRIX *x;
{
    int i,j;
    MATRIX *tmp;
    tmp = VC_GEE_matcopy(x);
    for(i=0;i<x->ncols;i++)
    {
	for(j=0;j<x->nrows;j++)
	{
	    MEL(tmp,i,j)= 1./(MEL(x,i,j));
	}
    }
    if(is_ephemeral(x))VC_GEE_destroy_matrix(x);
    return tmp;
}

/* ---- Witten by Z. Chen ---- */

static MATRIX *get_seq1(el_start, end)
int el_start, end;
{
    int i;
    int nelem = end - el_start + 1;
    double *load;
    MATRIX *tmp;
    tmp = VC_GEE_create_matrix(nelem, 1, EPHEMERAL);
    load = tmp->data;
    for (i = 0; i < nelem; i++)
    {
        *load = (double) (el_start + i);
        load++;
    }
    return tmp;
}

static MATRIX *get_rep_scalar(a, nrep)
int a, nrep;
{
    int i;
    MATRIX *tmp;
    double *load;
    tmp = VC_GEE_create_matrix(nrep, 1, EPHEMERAL);
    load = tmp->data;
    for (i = 0; i < nrep; i++)
    {
        *load = a;
        load++;
    }
    return tmp;
}

static MATRIX *get_rep(x, nrep)
MATRIX *x;
int nrep;
{
    int i, j;
    double *load, *look;
    int nelem = get_nelem(x);
    MATRIX *tmp;
    tmp = VC_GEE_create_matrix((x->nrows)*nrep, x->ncols, EPHEMERAL);
    load = tmp->data;
    for (i = 0; i < nrep; i++)
    {
        look = x->data;
        for (j = 0; j < nelem; j++)
        {
            *load = *look;
            load++; look++;
        }
    }
    free_if_ephemeral(x);
    return tmp;
}

static MATRIX *get_kronecker(x, y)
MATRIX *x, *y;  /* x and y are colvecs */
{
    int i, j;
    double *load, *lookx, *looky;
    MATRIX *tmp;
    tmp = VC_GEE_create_matrix((x->nrows)*(y->nrows), 1, EPHEMERAL);
    load = tmp->data;
    lookx = x->data;
    for (i = 0; i < (x->nrows); i++)
    {
        looky = y->data;
        for (j = 0; j < (y->nrows); j++)
        {
            *load = (*lookx) * (*looky);
            load++; looky++;
        }
        lookx++;
    }
    free_if_ephemeral(x); free_if_ephemeral(y);
    return tmp;
}

/* Sum the elements in each row and get a column vector*/
static MATRIX *get_sum1row(inmat)
MATRIX *inmat;
{
    int i, j;
    double t;
    MATRIX *outmat;
    double *load, *look;

    outmat = VC_GEE_create_matrix((inmat->nrows), 1, EPHEMERAL);
    load = outmat->data;
    look = inmat->data;
    for (i = 0; i < (inmat->nrows); i++)
    {
        t = 0.;
        for (j = 0; j < (inmat->ncols); j++)
        {
            t += (*look);
            look++;
        }
        *load = t;
        load++;
    }
    free_if_ephemeral(inmat);
    return outmat;
}


/* Sum the elements in each column and get a row vector*/
static MATRIX *get_sum2col(inmat)
MATRIX *inmat;
{
    int i, j;
    double t;
    MATRIX *outmat;
    double *load, *look;

    outmat = VC_GEE_create_matrix(1, (inmat->ncols), EPHEMERAL);
    load = outmat->data;

    for (i = 0; i < (inmat->ncols); i++)
    {
        look = inmat->data + i;
        t = 0.;
        for (j = 0; j < (inmat->nrows); j++)
        {
            t += (*look);
            look += (inmat->ncols);
        }
        *load = t;
        load++;
    }
    free_if_ephemeral(inmat);
    return outmat;
}


static MATRIX *VC_GEE_matexpit(x)
MATRIX *x;
{
    double *load, *look;
    double exp_val;
    int nelem, i;
    MATRIX *tmp;

    nelem = get_nelem(x);
    tmp = VC_GEE_create_matrix(x->nrows, x->ncols, EPHEMERAL);
    load = tmp->data;
    look = x->data;
    for (i = 0 ; i < nelem ; i++)
    {
        exp_val = exp(*look);
        *(load) = exp_val/(1.0 + exp_val);
        load++; look++;
    }
    free_if_ephemeral(x);
    return tmp ;
}

static void get_dlambda_i_dbetaT(dlambda_i_dbetaT,
                                 lambda_i, DM_i)
MATRIX *dlambda_i_dbetaT, *lambda_i, *DM_i;
{
    double *load,  *look;
    double lam;
    int i, j;

    load = dlambda_i_dbetaT->data;
    look = DM_i->data;
    for (i = 0 ; i < (lambda_i->nrows); i++)
    {
        lam = MEL(lambda_i, i, 0);
        for (j = 0 ; j < (DM_i->ncols); j++)
        {
            *load = (lam * (1. - lam) * (*look));
            load++ ; look++;
        }
    }
}

static void get_dmu_i_dbetaT(dmu_i_dbetaT, dlambda_i_dbetaT, K, m)
MATRIX  *dmu_i_dbetaT, *dlambda_i_dbetaT;
int K;
int m;
{
    int i, j, k, p_b;
    double *load, *look, *look2;

    p_b = dlambda_i_dbetaT->ncols;

    load = dmu_i_dbetaT->data;
    look = dlambda_i_dbetaT->data;
    look2 = (look + p_b);

    for (i = 0; i < m; i++)
    {
        for (k = 0; k < (K - 1); k++)
        {

            for (j = 0; j < p_b; j++)
            {
                *load = ((*look) - (*look2));
                load++; look++; look2++;
            }
        }
        for (j = 0; j < p_b; j++) /* At level K */
        {
            *load = *look;
            load++; look++; look2++;
        }
    }

    free_if_ephemeral(dlambda_i_dbetaT);

}


static void get_lambda_i(lambda_i, DM_i, beta)
MATRIX *lambda_i, *DM_i, *beta;
{
    int i, j;
    double *load, *l_ptr, *r_ptr;
    double lin;

    load = lambda_i->data;
    l_ptr = DM_i->data;

    if ((DM_i->ncols) != (beta->nrows))
    {
        error("Matrix multiplications: dimensions do not match");
    }

    for (i = 0; i < (lambda_i->nrows); i++)
    {
        lin = 0.;
        r_ptr = beta->data;
        for (j = 0; j < (DM_i->ncols); j++)
        {
            lin = lin + ((*l_ptr) * (*r_ptr));
            l_ptr++; r_ptr++;
        }
        *load = (exp(lin)/(1. + exp(lin)));
        load++;
    }
}

static void get_mu_i(mu_i, lambda_i, K, m)
MATRIX *mu_i, *lambda_i;
int K, m;
{
    double *load, *look;
    int j, k;

    load = mu_i->data;
    look = lambda_i->data;
    for (j = 0; j < m; j++)
    {
        for (k = 0; k < (K - 1); k++)
        {
            *load = ((*look) - (*(look + 1)));
            load++; look++;
        }
        *load = *look; /* At level K, mu==lambda */
        load++; look++;
    }
}

static MATRIX *get_outer(vec1, vec2)
MATRIX *vec1, *vec2;
{
    MATRIX *tmp;
    double *load, *look1, *look2;
    int i, j;

    if ((vec1->ncols != 1) || (vec2->ncols != 1))
    {
        error("M+-: args not a col-vec");
    }

    tmp = VC_GEE_create_matrix(vec1->nrows, vec2->nrows, EPHEMERAL);
    load = tmp->data;
    look1 = vec1->data;
    for (i = 0; i < (vec1->nrows); i++)
    {
        look2 = vec2->data;
        for (j = 0; j < (vec2->nrows); j++)
        {
            *load = ((*look1) * (*look2));
            load++; look2++;
        }
        look1++;
    }

    free_if_ephemeral(vec1);  free_if_ephemeral(vec2);

    return tmp;
}

static MATRIX *get_rbind(mat1, mat2)
MATRIX *mat1, *mat2;
{
    MATRIX *outmat;
    int i;
    double *load, *look;

    if ((mat1->ncols) != (mat2->ncols))
    {
        error("Matrices row bind: columns do not match");
    }

    outmat = VC_GEE_create_matrix((mat1->nrows) + (mat2->nrows),
            mat1->ncols, EPHEMERAL);
    load = outmat->data;
    look = mat1->data;
    for (i = 0; i < ((mat1->nrows) * (mat1->ncols)); i++)
    {
        *load = *look;
        load++; look++;
    }
    look = mat2->data;
    for (i = 0; i < ((mat2->nrows) * (mat2->ncols)); i++)
    {
        *load = *look;
        load++; look++;
    }

    free_if_ephemeral(mat1); free_if_ephemeral(mat2);

    return outmat;
}

static MATRIX *get_cbind(mat1, mat2)
MATRIX *mat1, *mat2;
{
    MATRIX *outmat;
    int i, j;
    double *load, *look1, *look2;

    if ((mat1->nrows) != (mat2->nrows))
    {
        error("Matrices column bind: rows do not match");
    }

    outmat = VC_GEE_create_matrix(mat1->nrows,
            (mat1->ncols + mat2->ncols), EPHEMERAL);
    load = outmat->data;
    look1 = mat1->data;
    look2 = mat2->data;
    for (i = 0; i < mat1->nrows; i++)
    {
        for (j = 0; j < mat1->ncols; j++)
        {
            *load = *look1;
            load++; look1++;
        }
        for (j = 0; j < mat2->ncols; j++)
        {
            *load = *look2;
            load++; look2++;
        }
    }

    free_if_ephemeral(mat1); free_if_ephemeral(mat2);

    return outmat;
}


static int get_rowindex(m, K, j1, k1, j2, k2)
int m, K, j1, k1, j2, k2;
{
    int indx;
    indx = (int) (pow(K, 2) * (2*m - j1 - 1) * j1/2. +
                  pow(K, 2) * (j2 - j1 -1) + k1*K + k2);
    return indx;
}



/* Get zeta_i, dzeta_i_dalphaT, and dzeta_i_dbetaT */
static void get_bivar_cumuls_i(
     MATRIX *zeta_i,
     MATRIX *dzeta_i_dalphaT,
     MATRIX *dzeta_i_dbetaT,
     MATRIX *lambda_i,
     MATRIX *DM_i,
     MATRIX *assocDM_i,
     MATRIX *alpha,
     int K,
     int p_b,
     int p_a,
     int m
     )
{
    double *zeta_i_ptr;
    int j_a;
    double *dzda_ptr;
    int j_b;
    double *dzdb_ptr;

    double psi_i, lam1_i, lam2_i, a_i, d_i, f_i, g_i;
    double dada, ddda, dfda, dgda, dpda;
    double dadb, dddb, dfdb;

    int j1, j2, k1, k2, r;
    double *lam_i_ptr, *asDM_ptr;


    zeta_i_ptr = zeta_i->data;
    dzda_ptr = dzeta_i_dalphaT->data;
    dzdb_ptr = dzeta_i_dbetaT->data;
    asDM_ptr = assocDM_i->data;
    lam_i_ptr = lambda_i->data;

    r = 0;

    for (j1 = 0; j1 < (m - 1); j1++)
    {
        for (j2 = (j1 + 1); j2 < m; j2++)
        {
            for (k1 = 0; k1 < K; k1++)
            {
                for (k2 = 0; k2 < K; k2++)
                {
                    lam1_i = *(lam_i_ptr + j1 * K + k1);
                    lam2_i = *(lam_i_ptr + j2 * K + k2);
                    psi_i = exp(MEL(VC_GEE_matmult(
                        get_matrix_row(assocDM_i, r), alpha), 0, 0));
                    a_i = (1. - (1. - psi_i) *
                        (lam1_i + lam2_i));
                    d_i = (pow(a_i, 2) - 4. *
                        (psi_i - 1.) * psi_i *
                        lam1_i * lam2_i);
                    g_i = 2. * (psi_i - 1.);
                    f_i = (a_i- sqrt(d_i));
                    *zeta_i_ptr = f_i/g_i;

                    /*  derivatives */
                    for (j_a = 0; j_a < p_a; j_a++)
                    {
                        dpda = (psi_i * (*asDM_ptr));
                        dada = (dpda * (lam1_i + lam2_i));
                        ddda = (2. * a_i * dada -
                            4. * lam1_i * lam2_i *
                            (2. * psi_i - 1.) * dpda);
                        dfda = (dada - ddda/(2. * sqrt(d_i)));
                        dgda = (2. * dpda);
                        *dzda_ptr = (dfda/g_i -
                            f_i/pow(g_i, 2) * dgda);

                        dzda_ptr++; asDM_ptr++;
                    }

                    for (j_b = 0; j_b < p_b; j_b++)
                    {
                        dadb = ((psi_i - 1.) *
                            (lam2_i * (1. - lam2_i) *
                            MEL(DM_i, j2 * K + k2, j_b) +
                            lam1_i * (1. - lam1_i) *
                            MEL(DM_i, j1 * K + k1, j_b)));
                        dddb = (2. * a_i * dadb -
                            4. * (psi_i - 1.) * psi_i *
                            (lam1_i * (1. - lam2_i) *
                            lam2_i * MEL(DM_i, j2 * K + k2, j_b) +
                            lam2_i * (1. - lam1_i) * lam1_i *
                            MEL(DM_i, j1 * K + k1, j_b)));
                        dfdb = (dadb - dddb/(2. * sqrt(d_i)));
                        *dzdb_ptr = (dfdb/g_i);

                        dzdb_ptr++;
                    }
                    r++;
                    zeta_i_ptr++;
                } /* end of for k2 */
            } /* end of for k1 */
        } /* end of for j2 */
    } /* end of for j1 */
}


/* Get xi_i, dxi_i_dalphaT, dxi_i_dbetaT */
/*  V_1i, and diagonal blocks of invV_2i   */
/* -------------------------------------- */
static void get_bivar_marginals_i(
     MATRIX *xi_i,
     MATRIX *V_1i,
     MATRIX *invV_2i,
     MATRIX *dxi_i_dalphaT,
     MATRIX *dxi_i_dbetaT,

     MATRIX *dzeta_i_dalphaT,
     MATRIX *dzeta_i_dbetaT,
     MATRIX *lambda_i,
     MATRIX *mu_i,
     MATRIX *DM_i,
     MATRIX *zeta_i,
     MATRIX *assocDM_i,
     MATRIX *alpha,
     int K,
     int p_b,
     int p_a,
     int m
     )
{
    MATRIX *mu_i_j1, *tempmat1, *tempmat2;
    int j, j1, k1, j2, k2, k,  r, r1, r2, r3;

    r = 0; j = 0;
    for (j1 = 0; j1 < (m-1); j1++)
    {
        for (j2 = (j1 + 1); j2 < m; j2++)
        {
            for (k1 = 0; k1 < (K - 1); k1++)
            {
                /* Case 1): k1<(K-1), k2<(K-1) */
                /* --------------------------- */
                for (k2 = 0; k2 < (K - 1); k2++)
                {
                    r1 = get_rowindex(m, K, j1, k1 + 1, j2, k2);
                    r2 = get_rowindex(m, K, j1, k1, j2, k2 + 1);
                    r3 = get_rowindex(m, K, j1, k1 + 1, j2, k2 + 1);
                    MEL(xi_i, r, 0) = ((MEL(zeta_i, r, 0) -
                        MEL(zeta_i, r1, 0) -
                        MEL(zeta_i, r2, 0) + MEL(zeta_i, r3, 0)));
                    /* Off-diagonal blocks of V_1i */
                    MEL(V_1i, j1*K+k1, j2*K+k2) = MEL(xi_i, r, 0) -
                        MEL(mu_i, j1 * K + k1, 0) * MEL(mu_i, j2 * K + k2, 0);
                    MEL(V_1i, j2*K+k2, j1*K+k1) = MEL(V_1i, j1*K+k1, j2*K+k2);
                    /* dxi_i_dalphaT and dxi_i_dbetaT */
                    for (k = 0; k < p_a; k++)
                        MEL(dxi_i_dalphaT, r, k) = ((MEL(dzeta_i_dalphaT, r, k) -
                            MEL(dzeta_i_dalphaT, r1, k) -
                            MEL(dzeta_i_dalphaT, r2, k) +
                            MEL(dzeta_i_dalphaT, r3, k)));
                    for (k = 0; k < p_b; k++)
                        MEL(dxi_i_dbetaT, r, k) = ((MEL(dzeta_i_dbetaT, r, k) -
                            MEL(dzeta_i_dbetaT, r1, k) -
                            MEL(dzeta_i_dbetaT, r2, k) +
                            MEL(dzeta_i_dbetaT, r3, k)));
                    r++;
                }

                /* Case 2): k1<(K-1), k2==(K-1) */
                /* ---------------------------- */
                k2 = K - 1;
                r1 = get_rowindex(m, K, j1, k1 + 1, j2, k2);
                MEL(xi_i, r, 0) = (MEL(zeta_i, r, 0) - MEL(zeta_i, r1, 0));
                /* Off-diagonal blocks of V_1i */
                MEL(V_1i, j1*K+k1, j2*K+k2) = (MEL(xi_i, r, 0) -
                    MEL(mu_i, j1 * K + k1, 0) * MEL(mu_i, j2 * K + k2, 0));
                MEL(V_1i, j2*K+k2, j1*K+k1) = MEL(V_1i, j1*K+k1, j2*K+k2);
                /* dxi_i_dalphaT and dxi_i_dbetaT */
                for (k = 0; k < p_a; k++)
                    MEL(dxi_i_dalphaT, r, k) = (MEL(dzeta_i_dalphaT, r, k) -
                       MEL(dzeta_i_dalphaT, r1, k));
                for (k = 0; k < p_b; k++)
                    MEL(dxi_i_dbetaT, r, k) = (MEL(dzeta_i_dbetaT, r, k) -
                        MEL(dzeta_i_dbetaT, r1, k));
                r++;
            }
            /* Case 3): k1==(K-1), k2<(K-1) */
            /* ---------------------------- */
            k1 = K - 1;
            for (k2 = 0; k2 < (K - 1); k2++)
            {
               r2 = get_rowindex(m, K, j1, k1, j2, k2 + 1);
               /* xi_i */
               MEL(xi_i, r, 0) = (MEL(zeta_i, r, 0) -
                   MEL(zeta_i, r2, 0));
               /* Off-diagonal blocks of V_1i */
               MEL(V_1i, j1*K+k1, j2*K+k2) = (MEL(xi_i, r, 0) -
                   MEL(mu_i, j1 * K + k1, 0) * MEL(mu_i, j2 * K + k2, 0));
               MEL(V_1i, j2*K+k2, j1*K+k1) = MEL(V_1i, j1*K+k1, j2*K+k2);
               /* dzeta_i_dalphaT and dzeta_i_dbetaT */
               for (k = 0; k < p_a; k++)
                   MEL(dxi_i_dalphaT, r, k) = (MEL(dzeta_i_dalphaT, r, k) -
                       MEL(dzeta_i_dalphaT, r2, k));
               for (k = 0; k < p_b; k++)
                   MEL(dxi_i_dbetaT, r, k) = (MEL(dzeta_i_dbetaT, r, k) -
                       MEL(dzeta_i_dbetaT, r2, k));
               r++;
            }
            /* Case 4): k1==(K-1), k2==(K-1) */
            /* ----------------------------- */
            k2 = K - 1;
            MEL(xi_i, r, 0) = MEL(zeta_i, r, 0);
            /* Off-diagonal blocks of V_1i */
            MEL(V_1i, j1*K+k1, j2*K+k2) = (MEL(xi_i, r, 0) -
                MEL(mu_i, j1 * K + k1, 0) * MEL(mu_i, j2 * K + k2, 0));
            MEL(V_1i, j2*K+k2, j1*K+k1) = MEL(V_1i, j1*K+k1, j2*K+k2);

            //Yuliang
            /*Rprintf("Inside get_bivar_marginals_i: checkpoint 1. \n");
            Rprintf("j * pow(K, 2) = %d. \n (j + 1) * pow(K, 2) - 1 = %d\n",test_1, test_2);
            */

            //Yuliang: cannot do pow in 64bit. Use int temp_pow_K2
            int temp_pow_K2=pow(K, 2);


            /* Diagonal blocks in invV_2i */
            tempmat1 = VC_GEE_extract_rows(xi_i, j * temp_pow_K2,
                                           (j + 1) * temp_pow_K2 - 1);
            make_permanent(tempmat1);

            tempmat2 = VC_GEE_matsub(VC_GEE_form_diag(tempmat1),
                get_outer(tempmat1, tempmat1));


            VC_GEE_destroy_matrix(tempmat1);
            VC_GEE_plug(get_cholinv(tempmat2),
                invV_2i, j * temp_pow_K2, j * temp_pow_K2);// replaced by temp_pow_K2


            /* dxi_i_dalphaT and dxi_i_dbetaT */
            for (k = 0; k < p_a; k++)
                *(ELREF(dxi_i_dalphaT, r, k)) = MEL(dzeta_i_dalphaT, r, k);
            for (k = 0; k < p_b; k++)
                *(ELREF(dxi_i_dbetaT, r, k)) = MEL(dzeta_i_dbetaT, r, k);
            r++;
            j++;
        } /* end of for j2 */
    } /* end of for j1 */

    /* Diagonal blocks in V_1i  */
    /* ------------------------- */
    for (j1 = 0; j1 < m; j1++)
    {
        mu_i_j1 = VC_GEE_extract_rows(mu_i, j1*K, (j1+1)*K-1);
        make_permanent(mu_i_j1);
        tempmat1 = VC_GEE_matsub(VC_GEE_form_diag(mu_i_j1),
            get_outer(mu_i_j1, mu_i_j1));
        VC_GEE_plug(tempmat1, V_1i, j1*K, j1*K);
        VC_GEE_destroy_matrix(mu_i_j1);
    }
}


static MATRIX *get_cholinv(X)
MATRIX *X;
{
    MATRIX *Y, *tempmat;
    int nrows, ncols;
    double *y;
    int i, j;
    double *det = Calloc(2, double), *z = Calloc(X->nrows, double);
    /* double *det = (double *) calloc(2,  sizeof(double)),
        *z = (double *) calloc(X->nrows,  sizeof(double));  */
    double rcond;
    int job, info;

    nrows = X->nrows;
    ncols = X->ncols;

    tempmat = VC_GEE_matcopy(X);
    y = tempmat->data;

    F77_CALL(dpoco)(y, &nrows, &ncols, &rcond, z, &info);

    job = 01;

    /* Rprintf("pos-def : %d \n", 1 - info); */

    if (info == 0) {
        F77_CALL(dpodi)(y, &nrows, &ncols, det, &job);
    }
    from_S(y, &nrows, &ncols, Y)

    for (i = 1; i < nrows; i++)
    {
        for (j = 0; j < i; j++)
        {
            MEL(Y, i, j) = MEL(Y, j, i);
        }
    }

    Free(z); Free(det);
    /* free(z); free(det); */
    free_if_ephemeral(X);
    VC_GEE_destroy_matrix(tempmat);

    return Y;
}

static void matrix_copyto(inmat, outmat)
MATRIX *inmat, *outmat;
{
    int i, j;
    if (((inmat->nrows) != (outmat->nrows)) || ((inmat->ncols) != (outmat->ncols)))
    {
        error("Copy matrix to: dimensions of matrices do not match");
    }
    for (i = 0 ; i < (inmat->nrows); i++)
    {
        for (j = 0; j < (inmat->ncols); j++)
        {
            MEL(outmat, i, j) = MEL(inmat, i, j);
        }
    }
    free_if_ephemeral(inmat);
}

static void row_replace(inmat, i, outmat, j)
MATRIX *inmat, *outmat;
int i, j;
{
    double *load, *look;
    int col;

    if ((inmat->ncols) != (outmat->ncols))
    {
        error("Row replace: columns of matrices do not match");
    }

    look = (inmat->data + i * (inmat->ncols));
    load = (outmat->data + j * (outmat->ncols));

    for (col = 0; col < inmat->ncols; col++)
    {
        *load = *look;
        load++; look++;
    }
    free_if_ephemeral(inmat);
}

static void col_replace(inmat, i, outmat, j)
MATRIX *inmat, *outmat;
int i, j;
{
    double *load, *look;
    int row, in_ncols, out_ncols;

    if ((inmat->nrows) != (outmat->nrows))
    {
        error("Column replace: rows of matrices do not match");
    }

    look = (inmat->data + i);
    in_ncols = inmat->ncols;
    load = (outmat->data + j);
    out_ncols = outmat->ncols;

    for (row = 0; row < inmat->nrows; row++)
    {
        *load = *look;
        load += in_ncols;
        look += out_ncols;
    }
    free_if_ephemeral(inmat);
}

static void rows_plug(inmat, i1, i2, outmat, j)
MATRIX *inmat, *outmat;
int i1, i2, j;
{
    double *load, *look;
    int i, col;

    if ((inmat->ncols) != (outmat->ncols))
    {
        error("Row plug: columns do not match");
    }
    if ( (j + i2 - i1) > outmat->nrows )
    {
        error("Row plug: socket is too small");
    }

    look = (inmat->data + i1 * (inmat->ncols));
    load = (outmat->data + j * (outmat->ncols));

    for (i = 0; i < (i2 - i1 + 1); i++)
    {
        for (col = 0; col < inmat->ncols; col++)
        {
            *load = *look;
            load++; look++;
        }
    }

    free_if_ephemeral(inmat);
}

static void cols_plug(inmat, i1, i2, outmat, j)
MATRIX *inmat, *outmat;
int i1, i2, j;
{
    int i, k;
    if ((inmat->nrows) != (outmat->nrows))
    {
        error("Column plug: rows do not match");
    }
    if ( (j + i2 - i1) > outmat->ncols )
    {
        error("Column plug: socket is too small");
    }

    for (i = 0; i < (i2 - i1 + 1); i++)
    {
        for (k = 0; k < inmat->nrows; k++)
        {
            MEL(outmat, k, j + i) = MEL(inmat, k, i1 + i);
        }
    }
    free_if_ephemeral(inmat);
}


static void matrix_addto(inmat, outmat)
MATRIX *inmat, *outmat;
{
    double *load, *look;
    int i;

    if (((inmat->nrows)!=(outmat->nrows)) || ((inmat->ncols)!=(outmat->ncols)))
    {
        error("Matrix add to: dimensions do not match");
    }

    load = outmat->data;
    look = inmat->data;

    for (i = 0; i < ((inmat->nrows) * (inmat->ncols)); i++)
    {
        *load += *look;
        load++; look++;
    }
    free_if_ephemeral(inmat);
}


static MATRIX *matrix_subtract(mat1, mat2)
MATRIX *mat1, *mat2;
{
    MATRIX *outmat;
    double *load, *look1, *look2;
    int i;

    if (((mat1->nrows)!=(mat2->nrows)) || ((mat1->ncols)!=(mat2->ncols)))
    {
        error("Matrix subtract from: dimensions do not match");
    }

    outmat = VC_GEE_create_matrix(mat1->nrows, mat1->ncols, EPHEMERAL);
    load = outmat->data;
    look1 = mat1->data;
    look2 = mat2->data;
    for (i = 0 ; i < ((mat1->nrows) * (mat1->ncols)) ; i++)
    {
        *load = ((*look1) - (*look2));
        load++; look1++; look2++;
    }
    /* Should not free Ytilde[i] */
    /* free_if_ephemeral(mat1);  free_if_ephemeral(mat2); */
    return outmat;
}


static MATRIX *matrix_multiply(mat1, mat2)
MATRIX *mat1, *mat2;
{
    MATRIX *outmat;
    int i, j, k;
    double *load;

    if ((mat1->ncols) != (mat2->nrows))
    {
        error("Multiply two matrices: dimensions do not match");
    }

    outmat = VC_GEE_create_matrix(mat1->nrows, mat2->ncols, EPHEMERAL);

    load = outmat->data;

    for (i = 0; i < (mat1->nrows); i++)
    {
        for (j = 0; j < (mat2->ncols); j++)
        {
            for (k = 0; k < mat1->ncols; k++)
            {
                *load += (MEL(mat1, i, k) * MEL(mat2, k, j));
            }
            load++;
        }
    }
    /* Should not free DM[i] */
    /* free_if_ephemeral(mat1);  free_if_ephemeral(mat2); */
    return outmat;
}

static void matrix_elem_mult(inmat, outmat)
MATRIX *inmat, *outmat;
{
    double *load, *look;
    int i;

    if (((inmat->nrows)!=(outmat->nrows)) || ((inmat->ncols)!=(outmat->ncols)))
    {
        error("Element-wise matrices multiplication: dimensions do not match");
    }

    load = outmat->data;
    look = inmat->data;

    for (i = 0; i < ((inmat->nrows) * (inmat->ncols)); i++)
    {
        *load *= (*look);
        load++; look++;
    }
    free_if_ephemeral(inmat);
}


static void matrix_row_mult(rowvec, outmat) /* row vector times a matrix*/
MATRIX *rowvec, *outmat;
{
    double *load, *look;
    int i, j;

    if ((rowvec->ncols)!=(outmat->ncols))
    {
        error("Second dimension not the same");
    }

    load = outmat->data;

    for (i = 0; i < outmat->nrows; i++)
    {
        look = rowvec->data;
        for (j = 0; j < outmat->ncols; j++)
        {
            *load *= (*look);
            load++; look++;
        }
    }
    free_if_ephemeral(rowvec);
}


static void matrix_col_mult(colvec, outmat) /* col vector times a matrix*/
MATRIX *colvec, *outmat;
{
    double *load, *look;
    int i, j;

    if ((colvec->nrows)!=(outmat->nrows))
    {
        error("First dimension not the same");
    }
    load = outmat->data;
    look = colvec->data;
    for (i = 0; i < outmat->nrows; i++)
    {
        for (j = 0; j < outmat->ncols; j++)
        {
            *load *= (*look);
            load++;
        }
        look++;
    }
    free_if_ephemeral(colvec);
}


static void scalar_times_matrix(a, X)
double a;
MATRIX *X;
{
    double *load;
    int i, nelem;

    nelem = get_nelem(X);
    load = X->data;
    for (i = 0 ; i < nelem ; i++) {
        *load *= a ;
        load++ ;
    }
}


static MATRIX *get_matrix_row(X, i)
MATRIX *X;
int i;
{
    MATRIX *rowvec;
    int j;
    double *load;

    rowvec = VC_GEE_create_matrix(1, X->ncols, EPHEMERAL);
    load = rowvec->data;
    for (j = 0; j < rowvec->ncols; j++)
    {
        *load = MEL(X, i, j);
        load++;
    }
    return rowvec;
}

static void set_zero(X)
MATRIX *X;
{
    int i, nel ;
    double *load;

    nel = get_nelem(X);
    load = X->data;
    for (i = 0; i < nel; i++)
    {
        *load = 0.;
        load++;
    }
}



static void get_matsub(MATRIX *mat1,
              MATRIX *mat2, MATRIX *outmat)
{
    int i, j;

    if ((mat1->nrows != mat2->nrows) ||
        (mat1->nrows != outmat->nrows) ||
        (mat2->nrows != outmat->nrows) ||
        (mat1->ncols != mat2->ncols) ||
        (mat1->ncols != outmat->ncols) ||
        (mat2->ncols != outmat->ncols))
    {
        error("get_matsub(): Dimensions do not match");
    }
    for (i = 0; i < outmat->nrows; i++)
    {
        for (j = 0; j < outmat->ncols; j++)
        {
            MEL(outmat, i, j) = (MEL(mat1, i, j) -
                MEL(mat2, i, j));
        }
    }
}


static void get_matadd(MATRIX *mat1,
              MATRIX *mat2, MATRIX *outmat)
{
    int i, j;

    if ((mat1->nrows != mat2->nrows) ||
        (mat1->nrows != outmat->nrows) ||
        (mat2->nrows != outmat->nrows) ||
        (mat1->ncols != mat2->ncols) ||
        (mat1->ncols != outmat->ncols) ||
        (mat2->ncols != outmat->ncols))
    {
        error("get_matadd(): Dimensions do not match");
    }
    for (i = 0; i < outmat->nrows; i++)
    {
        for (j = 0; j < outmat->ncols; j++)
        {
            MEL(outmat, i, j) = (MEL(mat1, i, j) +
                MEL(mat2, i, j));
        }
    }
}

static void get_mattransp(MATRIX *X, MATRIX *Xtransp)
{
    int i, j;

    if ((X->nrows != Xtransp->ncols) ||
        (X->ncols != Xtransp->nrows))
    {
        error("get_mattransp(): Dimensions do not match");
    }
    for (i = 0; i < X->nrows; i++)
    {
        for (j = 0; j < X->ncols; j++)
        {
            MEL(Xtransp, j, i) = MEL(X, i, j);
        }
    }
}

static void get_matmult(MATRIX *mat1,
              MATRIX *mat2, MATRIX *outmat)
{
    int i, j, k;

    if ((mat1->nrows != outmat->nrows) ||
        (mat1->ncols != mat2->nrows) ||
        (mat2->ncols != outmat->ncols))
    {
        error("get_matmult(): Dimensions do not match");
    }
    for (i = 0; i < outmat->nrows; i++)
    {
        for (j = 0; j < outmat->ncols; j++)
        {
            for (k = 0; k < mat1->ncols; k++)
            {
                MEL(outmat, i, j) += (MEL(mat1, i, k) *
                MEL(mat2, k, j));
            }
        }
    }
}

static double get_max_reldif(MATRIX *theta,
                MATRIX *theta_o)
{
    int i;
    double maxreldif, reldif;

    maxreldif = fabs(MEL(theta, 0, 0)/MEL(theta_o, 0, 0) - 1.);
    for (i = 1; i < theta->nrows; i++)
    {
        reldif = fabs(MEL(theta, i, 0)/MEL(theta_o, i, 0) - 1.);
        if (maxreldif < reldif)
        {
             maxreldif = reldif;
        }
    }
    return maxreldif;
}

static void cholinv(MATRIX *X, MATRIX *Xinv)
{
    MATRIX *tempmat;
    int nrows, ncols;
    double *y;
    int i, j;
    double *det = Calloc(2, double), *z = Calloc(X->nrows, double);
    /* double *det = (double *) calloc(2,  sizeof(double)),
        *z = (double *) calloc(X->nrows,  sizeof(double));  */
    double rcond;
    int job, info;

    if ((X->nrows != Xinv->nrows) || (X->ncols != Xinv->ncols))
    {
        Free(z); Free(det);
        error("cholinv: Dimensions do not match");
    }

    nrows = X->nrows;
    ncols = X->ncols;

    matrix_copyto(X, Xinv);
    y = Xinv->data;

    F77_CALL(dpoco)(y, &nrows, &ncols, &rcond, z, &info);

    job = 01;

    /* Rprintf("pos-def : %d \n", 1 - info); */

    if (info == 0) {
        F77_CALL(dpodi)(y, &nrows, &ncols, det, &job);
    }

    from_S(y, &nrows, &ncols, tempmat)

    matrix_copyto(tempmat, Xinv);

    for (i = 1; i < nrows; i++)
    {
        for (j = 0; j < i; j++)
        {
            MEL(Xinv, i, j) = MEL(Xinv, j, i);
        }
    }

    Free(z); Free(det);
    /* free(z); free(det); */
    free_if_ephemeral(X);
    VC_GEE_destroy_matrix(tempmat);
}


static void outer_colvec_to(MATRIX *colvec,
                  MATRIX *outmat)
{
    int i, j;

    if((colvec->nrows != outmat->nrows) ||
        (colvec->nrows != outmat->ncols))
    {
        error("outer_colvec_to(): dimensions do not match");
    }
    for (i = 0; i < outmat->nrows; i++)
    {
        for (j = 0; j < outmat->ncols; j++)
        {
            MEL(outmat, i, j) = (
                MEL(colvec, i, 0) *
                MEL(colvec, j, 0) );
        }
    }
}


static void add_outer_colvec_to(MATRIX *colvec,
                  MATRIX *outmat)
{
    int i, j;

    if((colvec->nrows != outmat->nrows) ||
        (colvec->nrows != outmat->ncols))
    {
        error("add_outer_colvec_to(): dimensions do not match");
    }
    for (i = 0; i < outmat->nrows; i++)
    {
        for (j = 0; j < outmat->ncols; j++)
        {
            MEL(outmat, i, j) += (
                MEL(colvec, i, 0) *
                MEL(colvec, j, 0) );
        }
    }
}


static void get_estfun(MATRIX *D,
                MATRIX *invV,
                MATRIX *Y,
                MATRIX *mu,
                MATRIX *U)
{
    int i, j, k;
    double elem;

    if ( (D->ncols != invV->nrows) ||
        (invV->ncols != Y->nrows) ||
        (Y->nrows != mu->nrows) ||
        (U->nrows != D->nrows))
    {
        error("get_estfun(): Dimensions do not match");
    }
    set_zero(U);

    for (i = 0; i < D->nrows; i++)
    {
        for (j = 0; j < invV->ncols; j++)
        {
            elem = 0.;
            for (k = 0; k < D->ncols; k++)
            {
                elem += (MEL(D, i, k) * MEL(invV, k, j));
            }
            MEL(U, i, 0) += (elem * (MEL(Y, j, 0) - MEL(mu, j, 0)));
        }
    }
}

static void get_dvd(MATRIX *D,
                MATRIX *invV,
                MATRIX *dmu_dbetaT,
                MATRIX *M)
{
    int i, j, k, p;
    double elem;
    MATRIX *rowvec_DV;

    if ((D->ncols != invV->nrows) ||
        (invV->ncols != dmu_dbetaT->nrows) ||
        (M->nrows != D->nrows) ||
        (M->ncols != dmu_dbetaT->ncols))
    {
        error("get_dvd(): Dimensions do not match");
    }

    rowvec_DV = VC_GEE_create_matrix(1, invV->ncols, EPHEMERAL);

    set_zero(M);

    for (i = 0; i < D->nrows; i++)
    {
        for (p = 0; p < dmu_dbetaT->ncols; p++)
        {
            set_zero(rowvec_DV);
            elem = 0.;
            for (k = 0; k < invV->ncols; k++)
            {
                for (j = 0; j < D->ncols; j++)
                {
                    MEL(rowvec_DV, 0, k) += (MEL(D, i, j)
                        * MEL(invV, j, k));
                }
                MEL(M, i, p) += (MEL(rowvec_DV, 0, k) *
                    MEL(dmu_dbetaT, k, p));
            }
        }
    }
    VC_GEE_destroy_matrix(rowvec_DV);
}


static void fisherscoring(double scale,
                MATRIX *invM,
                MATRIX *U,
                MATRIX *par)
{
    int i, j;

    for (i = 0; i < par->nrows; i++)
    {
        for (j = 0; j < invM->ncols; j++)
        {
            MEL(par, i, 0) += ( scale *
                MEL(invM, i, j) * MEL(U, j, 0) );
        }
    }
}


static void get_dpXt_i_dvpT_l(
       MATRIX *dpXt_i_dvpT_l,
       int l,
       MATRIX *X_i_ENUM,
       MATRIX *wgt_i,
       MATRIX *Xtilde_i_Mat_ext,
       MATRIX *invGast_ij_all,
       MATRIX *G_ij_all,
       MATRIX *L_x_ij,
       MATRIX *delta_i,
       int m,
       int K_x)
{
    int j, q, k, r, p, col;

    set_zero(dpXt_i_dvpT_l);

    for (j = 0; j < m; j++)
    {
        if (MEL(delta_i, j, 0) != 1 )
        {
            q = (int) MEL(X_i_ENUM, j, l);
            /* q is the enumerated level for X_ij */
          if (q != 0) {
            col = 0;
            for (k = 0; k < (K_x + 1); k++)
            { /* indexing the true category in \pi_{ijkr} */
                for (r = 0; r < K_x; r++)
                { /* indexing the surrogate level in \pi_{ijkr} */
                    for (p = 0; p < L_x_ij->nrows; p++)
                    { /* indexing L_x_ij */

                        MEL(dpXt_i_dvpT_l, 0, col) += (
                            1./MEL(Xtilde_i_Mat_ext, j, q)
                            * (-1.) * MEL(invGast_ij_all, q - 1, j * K_x + r)
                            * MEL(G_ij_all, k, j * (K_x + 1) + r + 1)
                            * (1 - MEL(G_ij_all, k, j * (K_x + 1) + r + 1))
                            * MEL(Xtilde_i_Mat_ext, j, k)
                            * MEL(L_x_ij, p, 0) );
                        col++;
                    }
                }
            }
          } else {
            /* the case of q = 0 */
            col = 0;
            for (k = 0; k < (K_x + 1); k++)
            { /* k is indexing the true category in \pi_{ijkr} */
                for (r = 0; r < K_x; r++)
                { /* r is indexing the surrogate level in \pi_{ijkr} */
                    for (p = 0; p < L_x_ij->nrows; p++)
                    { /* p is indexing L_x_ij */
                        MEL(dpXt_i_dvpT_l, 0, col) += (
                            MEL(wgt_i, 0, l)/MEL(Xtilde_i_Mat_ext, j, q)
                            * get_1_colsum(invGast_ij_all, j * K_x + r)
                            * MEL(G_ij_all, k, j * (K_x + 1) + r + 1)
                            * (1 - MEL(G_ij_all, k, j * (K_x + 1) + r + 1))
                            * MEL(Xtilde_i_Mat_ext, j, k)
                            * MEL(L_x_ij, p, 0) );
                        col++;
                    }
                }
            }
          } /* end of else */
        } /* end of if delta_ij = 0*/
    }
    /* scalar_times_matrix(MEL(wgt_i, 0, l), dpXt_i_dvpT_l); */
}


static MATRIX *get_dYtil_i_dgT(
       MATRIX *Ytilde_i_Mat_ext,
       MATRIX *invPast_ij_all,
       MATRIX *P_ij_all,
       MATRIX *L_ij,
       MATRIX *delta_i,
       int m,
       int K)
{
    MATRIX *dYtil_i_dgT;
    int j, k, r, p, u, col;
    MATRIX *tempmat;

    dYtil_i_dgT = VC_GEE_create_matrix(
        m * K, (K + 1) * K * (L_ij->nrows), PERMANENT);
    tempmat = VC_GEE_create_matrix(
        K, (K + 1) * K * (L_ij->nrows), PERMANENT);

    for (j = 0; j < m; j++)
    {
        if (MEL(delta_i, j, 0) != 1.)
        {
            set_zero(tempmat);

            for (u = 0; u < K; u++)
            { /* u is indexing the level in dYtil_iju_dgT */
                col = 0;
                for (k = 0; k < (K + 1); k++)
                { /* k is indexing the true category in \tau_ijkr */
                    for (r = 0; r < K; r++)
                    { /* r is indexing the surrogate level in \tau_ijkr */
                        for (p = 0; p < L_ij->nrows; p++)
                        { /* p is indexing L_ij */
                            MEL(tempmat, u, col) += (
                                -1. * MEL(invPast_ij_all, u, j * K + r)
                                * MEL(P_ij_all, k, j * (K + 1) + r + 1)
                                * (1 - MEL(P_ij_all, k, j * (K + 1) + r + 1))
                                * MEL(Ytilde_i_Mat_ext, j, k)
                                * MEL(L_ij, p, 0) );
                            col++;
                        }
                    }
                }
            }
        }
        VC_GEE_plug(tempmat, dYtil_i_dgT, j * K, 0);
    }

    VC_GEE_destroy_matrix(tempmat);

    return dYtil_i_dgT;
}

static MATRIX *get_dZtil_i_dgT(
           MATRIX *dYtil_i_dgT,
           MATRIX *Ytilde_i,
           int m,
           int K)
{
    MATRIX *dZtil_i_dgT;
    int j1, j2, k1, k2, r, p;

    dZtil_i_dgT = VC_GEE_create_matrix(
        (int) (pow(K, 2) * ((m * (m - 1))/2)),
        dYtil_i_dgT->ncols,  PERMANENT);

    r = 0;

    for (j1 = 0; j1 < (m - 1); j1++)
    {
        for (j2 = (j1 + 1); j2 < m; j2++)
        {
            for (k1 = 0; k1 < K; k1++)
            {
                for (k2 = 0; k2 < K; k2++)
                {
                    for (p = 0; p < (dYtil_i_dgT->ncols); p++)
                    {
                        MEL(dZtil_i_dgT, r, p) =
                            ( MEL(dYtil_i_dgT, j1 * K + k1, p) *
                              MEL(Ytilde_i, j2 * K + k2, 0) +
                              MEL(dYtil_i_dgT, j2 * K + k2, p) *
                              MEL(Ytilde_i, j1 * K + k1, 0) );
                    }
                    r++;
                }
            }
        }
    }

    return dZtil_i_dgT;
}

static MATRIX *form_matrix(
                      double *dblptr ,
                      int nrows ,
                      int ncols,
                      int permanence)
{
    MATRIX *outmat;
    int i, j;
    double *load;

    outmat = VC_GEE_create_matrix(nrows, ncols, permanence);
    load = dblptr;
    for ( j = 0 ; j < ncols ; j++ )
    {
        for ( i = 0 ; i < nrows ; i++ )
        {
            MEL(outmat , i , j ) = (double) * (load ++);
        }
    }

    return outmat;
}


static double get_1_colsum(MATRIX *X, int i)
{
    double elsum;
    int j;

    elsum = 0.;
    for (j = 0; j < X->nrows; j++)
    {
        elsum += MEL(X, j, i);
    }

    return elsum;
}

static double get_1_rowsum(MATRIX *X, int i)
{
    double elsum;
    int j;

    elsum = 0.;
    for (j = 0; j < X->ncols; j++)
    {
        elsum += MEL(X, i, j);
    }

    return elsum;
}

/* End of file*/


