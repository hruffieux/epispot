/*
*
* This file is part of the `epispot` R package:
*     https://github.com/hruffieux/epispot
*
* Functions for computationally expensive multiconcave block coordinate ascent updates.
*
* These functions use Eigen::Map to pass large matrices by reference from R.
* Given dimensionalities involved in some applications, copying such matrices
* would imply a prohibitive RAM overconsumption.
*
*
*/

#include "utils.h"

// for epispot_core function
// [[Rcpp::export]]
void coreBatch(const MapMat X,
                   const MapMat Y,
                   MapArr2D gam_vb,
                   const MapArr1D log_om,
                   const MapArr1D log_1_min_om,
                   const double log_sig2_inv_vb,
                   const MapArr1D log_tau_vb,
                   MapMat m1_beta,
                   MapMat mat_x_m1,
                   MapArr2D mu_beta_vb,
                   const MapArr1D sig2_beta_vb,
                   const MapArr1D tau_vb) {

  mu_beta_vb = (X.transpose() * (Y - mat_x_m1) +
    (Y.rows() - 1) * m1_beta).array().rowwise() * (sig2_beta_vb * tau_vb).transpose();

  gam_vb = exp(-logOnePlusExpMat(((mu_beta_vb.square().array().rowwise() / (-2 * sig2_beta_vb.transpose())).rowwise() -
    (log_tau_vb.transpose() / 2 + log(sig2_beta_vb).transpose() / 2)).colwise()
                                   + (log_1_min_om - log_om) - log_sig2_inv_vb / 2));

  m1_beta = mu_beta_vb * gam_vb;
  mat_x_m1 = X * m1_beta;

}



// for epispot_z_core and epispot_mix_core function
// [[Rcpp::export]]
void coreZBatch(const MapMat X,
                    const MapMat Y,
                    MapArr2D gam_vb,
                    const MapArr1D log_om,
                    const MapArr1D log_1_min_om,
                    const double log_sig2_inv_vb,
                    const MapArr1D log_tau_vb,
                    MapMat m1_beta,
                    MapMat mat_x_m1,
                    MapMat mat_z_mu,
                    MapArr2D mu_beta_vb,
                    const MapArr1D sig2_beta_vb,
                    const MapArr1D tau_vb) {


  mu_beta_vb = (X.transpose() * (Y - mat_x_m1 - mat_z_mu) +
    (Y.rows() - 1) * m1_beta).array().rowwise() * (sig2_beta_vb * tau_vb).transpose();

  gam_vb = exp(-logOnePlusExpMat(((mu_beta_vb.square().array().rowwise() / (-2 * sig2_beta_vb.transpose())).rowwise() -
    (log_tau_vb.transpose() / 2 + log(sig2_beta_vb).transpose() / 2)).colwise()
                                   + (log_1_min_om - log_om) - log_sig2_inv_vb / 2));

  m1_beta = mu_beta_vb * gam_vb;
  mat_x_m1 = X * m1_beta;

}


// for epispot_probit_core function
// [[Rcpp::export]]
void coreProbitBatch(const MapMat X,
                     const MapMat W,
                     MapArr2D gam_vb,
                     const MapArr1D log_om,
                     const MapArr1D log_1_min_om,
                     const double log_sig2_inv_vb,
                     MapMat m1_beta,
                     MapMat mat_x_m1,
                     MapMat mat_z_mu,
                     MapArr2D mu_beta_vb,
                     const double sig2_beta_vb) {


  mu_beta_vb = sig2_beta_vb * (X.transpose() * (W - mat_x_m1 - mat_z_mu) +
    (W.rows() - 1) * m1_beta).array();

  gam_vb = exp(-logOnePlusExpMat((mu_beta_vb.square().array() / (-2 * sig2_beta_vb)).colwise() +
    (log_1_min_om - log_om) - log_sig2_inv_vb / 2 - log(sig2_beta_vb) / 2));

  m1_beta = mu_beta_vb * gam_vb;
  mat_x_m1 = X * m1_beta;

}

