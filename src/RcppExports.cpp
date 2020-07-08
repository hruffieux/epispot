// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "epispot_types.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// coreDualLoop
void coreDualLoop(const MapMat X, const MapMat Y, MapArr2D gam_vb, const MapArr2D log_Phi_mu_theta_plus_zeta, const MapArr2D log_1_min_Phi_mu_theta_plus_zeta, const double log_sig2_inv_vb, const MapArr1D log_tau_vb, MapMat beta_vb, MapMat X_beta_vb, MapArr2D mu_beta_vb, const MapArr1D sig2_beta_vb, const MapArr1D tau_vb, const MapArr1D shuffled_ind, const double c);
RcppExport SEXP _epispot_coreDualLoop(SEXP XSEXP, SEXP YSEXP, SEXP gam_vbSEXP, SEXP log_Phi_mu_theta_plus_zetaSEXP, SEXP log_1_min_Phi_mu_theta_plus_zetaSEXP, SEXP log_sig2_inv_vbSEXP, SEXP log_tau_vbSEXP, SEXP beta_vbSEXP, SEXP X_beta_vbSEXP, SEXP mu_beta_vbSEXP, SEXP sig2_beta_vbSEXP, SEXP tau_vbSEXP, SEXP shuffled_indSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MapMat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const MapMat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< MapArr2D >::type gam_vb(gam_vbSEXP);
    Rcpp::traits::input_parameter< const MapArr2D >::type log_Phi_mu_theta_plus_zeta(log_Phi_mu_theta_plus_zetaSEXP);
    Rcpp::traits::input_parameter< const MapArr2D >::type log_1_min_Phi_mu_theta_plus_zeta(log_1_min_Phi_mu_theta_plus_zetaSEXP);
    Rcpp::traits::input_parameter< const double >::type log_sig2_inv_vb(log_sig2_inv_vbSEXP);
    Rcpp::traits::input_parameter< const MapArr1D >::type log_tau_vb(log_tau_vbSEXP);
    Rcpp::traits::input_parameter< MapMat >::type beta_vb(beta_vbSEXP);
    Rcpp::traits::input_parameter< MapMat >::type X_beta_vb(X_beta_vbSEXP);
    Rcpp::traits::input_parameter< MapArr2D >::type mu_beta_vb(mu_beta_vbSEXP);
    Rcpp::traits::input_parameter< const MapArr1D >::type sig2_beta_vb(sig2_beta_vbSEXP);
    Rcpp::traits::input_parameter< const MapArr1D >::type tau_vb(tau_vbSEXP);
    Rcpp::traits::input_parameter< const MapArr1D >::type shuffled_ind(shuffled_indSEXP);
    Rcpp::traits::input_parameter< const double >::type c(cSEXP);
    coreDualLoop(X, Y, gam_vb, log_Phi_mu_theta_plus_zeta, log_1_min_Phi_mu_theta_plus_zeta, log_sig2_inv_vb, log_tau_vb, beta_vb, X_beta_vb, mu_beta_vb, sig2_beta_vb, tau_vb, shuffled_ind, c);
    return R_NilValue;
END_RCPP
}
// coreDualInfoLoop
void coreDualInfoLoop(const MapMat V, const MapMat W, MapArr1D rho_vb, const MapArr1D log_om, const MapArr1D log_1_min_om, const double s2, MapVec xi_vb, MapMat mat_v_mu, MapArr1D mu_xi_vb, const double sig2_xi_vb, const MapArr1D shuffled_ind, const double c);
RcppExport SEXP _epispot_coreDualInfoLoop(SEXP VSEXP, SEXP WSEXP, SEXP rho_vbSEXP, SEXP log_omSEXP, SEXP log_1_min_omSEXP, SEXP s2SEXP, SEXP xi_vbSEXP, SEXP mat_v_muSEXP, SEXP mu_xi_vbSEXP, SEXP sig2_xi_vbSEXP, SEXP shuffled_indSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MapMat >::type V(VSEXP);
    Rcpp::traits::input_parameter< const MapMat >::type W(WSEXP);
    Rcpp::traits::input_parameter< MapArr1D >::type rho_vb(rho_vbSEXP);
    Rcpp::traits::input_parameter< const MapArr1D >::type log_om(log_omSEXP);
    Rcpp::traits::input_parameter< const MapArr1D >::type log_1_min_om(log_1_min_omSEXP);
    Rcpp::traits::input_parameter< const double >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< MapVec >::type xi_vb(xi_vbSEXP);
    Rcpp::traits::input_parameter< MapMat >::type mat_v_mu(mat_v_muSEXP);
    Rcpp::traits::input_parameter< MapArr1D >::type mu_xi_vb(mu_xi_vbSEXP);
    Rcpp::traits::input_parameter< const double >::type sig2_xi_vb(sig2_xi_vbSEXP);
    Rcpp::traits::input_parameter< const MapArr1D >::type shuffled_ind(shuffled_indSEXP);
    Rcpp::traits::input_parameter< const double >::type c(cSEXP);
    coreDualInfoLoop(V, W, rho_vb, log_om, log_1_min_om, s2, xi_vb, mat_v_mu, mu_xi_vb, sig2_xi_vb, shuffled_ind, c);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_epispot_coreDualLoop", (DL_FUNC) &_epispot_coreDualLoop, 14},
    {"_epispot_coreDualInfoLoop", (DL_FUNC) &_epispot_coreDualInfoLoop, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_epispot(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
