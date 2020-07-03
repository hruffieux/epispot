# This file is part of the `epispot` R package:
#     https://github.com/hruffieux/epispot
#
# Internal functions gathering the ELBO terms common to core algorithms.
#

########################################################
## E log p(beta, gamma | rest) - E log q(beta, gamma) ##
########################################################

e_beta_gamma_dual_info_ <- function(V, gam_vb, log_sig2_inv_vb, log_tau_vb,
                                    log_1_min_Phi_mat_v_mu, log_Phi_mat_v_mu, 
                                    mu_c_vb, m2_beta, sig2_beta_vb, sig2_c_vb, 
                                    sig2_rho_vb, list_sig2_theta_vb, 
                                    sig2_inv_vb, tau_vb, zeta_vb, 
                                    resp_spec = FALSE, bool_blocks = FALSE,
                                    bool_modules = FALSE, vec_fac_bl_y = NULL,
                                    vec_fac_bl_theta = NULL) {
  
  stopifnot(!(bool_blocks & bool_modules))
  
  eps <- .Machine$double.eps^0.75 # to control the argument of the log when gamma is very small
  
  arg <- log_sig2_inv_vb * gam_vb / 2 +
    sweep(gam_vb, 2, log_tau_vb, `*`) / 2 -
    sweep(m2_beta, 2, tau_vb, `*`) * sig2_inv_vb / 2 +
    gam_vb * log_Phi_mat_v_mu +
    (1 - gam_vb) * log_1_min_Phi_mat_v_mu -
    sig2_rho_vb / 2 + 1 / 2 * sweep(gam_vb, 2, log(sig2_beta_vb) + 1, `*`) -
    gam_vb * log(gam_vb + eps) - (1 - gam_vb) * log(1 - gam_vb + eps)
  
  if (is.list(list_sig2_theta_vb)) {
    
    arg <- sweep(arg, 1, unlist(lapply(list_sig2_theta_vb, diag)) / 2, `-`)
    
  } else if (!is.null(vec_fac_bl_theta)) {
    
    stopifnot(bool_blocks | bool_modules)
    
    if (bool_blocks) {
      
      bl_ids <- as.numeric(levels(vec_fac_bl_theta))
      n_bl <- length(bl_ids)
      
      vec_p_bl <- table(vec_fac_bl_theta)
      
      arg <- sweep(arg, 1, unlist(sapply(1:n_bl, function(bl) list_sig2_theta_vb[vec_fac_bl_theta == bl_ids[bl]])), `-`)
      
    } else { # bool_modules
      
      stopifnot(!is.null(vec_fac_bl_y))
      
      bl_ids_x <- as.numeric(levels(vec_fac_bl_theta))
      n_bl_x <- length(bl_ids_x)
      
      vec_p_bl <- table(vec_fac_bl_theta)
      
      bl_ids_y <- as.numeric(levels(vec_fac_bl_y))
      n_bl_y <- length(bl_ids_y)
      
      for(bl_y in 1:n_bl_y){
        
          arg[, vec_fac_bl_y == bl_ids_y[bl_y]] <- sweep(arg[, vec_fac_bl_y == bl_ids_y[bl_y], drop = FALSE], 1,
                                                         as.vector(unlist(sapply(1:n_bl_x, function(bl_x) list_sig2_theta_vb[vec_fac_bl_theta == bl_ids_x[bl_x], bl_y]))), `-`)
        
      }
    }
    
  } else {
    arg <- arg - list_sig2_theta_vb / 2
  }
  
  if (bool_blocks) {
    
    n_bl <- length(V)
    
    vec_arg <- as.vector(unlist(sapply(1:n_bl, function(bl) {
      V[[bl]]^2 %*% (zeta_vb[[bl]] * (sig2_c_vb[bl] + mu_c_vb[[bl]]^2) - (mu_c_vb[[bl]] * zeta_vb[[bl]])^2) / 2})))
    arg <- sweep(arg, 1, vec_arg, `-`)
    
  } else if (bool_modules) {
    
    n_bl_x <- length(V)
    
    bl_ids_y <- as.numeric(levels(vec_fac_bl_y))
    n_bl_y <- length(bl_ids_y)
    
    for(bl_y in 1:n_bl_y){
      
      vec_arg <- as.vector(unlist(sapply(1:n_bl_x, function(bl_x) {
        
        V[[bl_x]]^2 %*% (zeta_vb[[bl_x]][, bl_y] * (sig2_c_vb[bl_x, bl_y] + mu_c_vb[[bl_x]][, bl_y]^2) - (mu_c_vb[[bl_x]][, bl_y] * zeta_vb[[bl_x]][, bl_y])^2) / 2
      })))
      
      arg[, vec_fac_bl_y == bl_ids_y[bl_y]] <- sweep(arg[, vec_fac_bl_y == bl_ids_y[bl_y], drop = FALSE], 1, vec_arg, `-`)
    }
    
    
  } else {
    
    m1_c <- mu_c_vb * zeta_vb
    m2_c <- zeta_vb * (sig2_c_vb + mu_c_vb^2)
    
    if (resp_spec) {
      arg <- arg - V^2 %*% (m2_c - m1_c^2) / 2
    } else {
      arg <- sweep(arg, 1, V^2 %*% (m2_c - m1_c^2) / 2, `-`)
    }
    
  }
  
  sum(arg)
  
}


######################################################
## E log p(c, zeta_vb | rest) - E log q(c, zeta_vb) ##
######################################################

e_c_zeta_ <- function(log_om, log_1_min_om, mu_c_vb, s2, sig2_c_vb,
                      zeta_vb, resp_spec = FALSE) {
  
  eps <- .Machine$double.eps^0.25 # to control the argument of the log when gamma is very small
  
  arg <- - log(s2) * zeta_vb / 2  - zeta_vb * (sig2_c_vb + mu_c_vb^2) /(2 * s2) +
    1 / 2 * zeta_vb * (log(sig2_c_vb) + 1) - zeta_vb * log(zeta_vb + eps) -
    (1 - zeta_vb) * log(1 - zeta_vb + eps)
  
  if (resp_spec) {
    arg <- arg + sweep(zeta_vb, 1, log_om, `*`) + sweep((1 - zeta_vb), 1, log_1_min_om, `*`)
  } else {
    arg <- arg + zeta_vb * log_om + (1 - zeta_vb) * log_1_min_om
  }
  
  sum(arg)
}

############################################
## E log p(rho | rest) - E log q(rho) ##
############################################

e_rho_ <- function(mu_rho_vb, n0, sig2_rho_vb, T0_inv, vec_sum_log_det_rho) {
  
  q <- length(mu_rho_vb)
  (vec_sum_log_det_rho - # vec_sum_log_det_rho = log(det(T0_inv)) + log(det(sig2_rho_vb))
      T0_inv * crossprod(mu_rho_vb - n0) -
      q * T0_inv * sig2_rho_vb + q) / 2 # trace of a product
  
}


##################################################
## E log p(sig2_inv | rest) - E log q(sig2_inv) ##
##################################################

e_sig2_inv_ <- function(lambda, lambda_vb, log_sig2_inv_vb, nu, nu_vb, sig2_inv_vb) {
  
  (lambda - lambda_vb) * log_sig2_inv_vb - (nu - nu_vb) * sig2_inv_vb +
    lambda * log(nu) - lambda_vb * log(nu_vb) - lgamma(lambda) + lgamma(lambda_vb)
  
}


########################################
## E log p(tau | rest) - E log q(tau) ##
########################################

e_tau_ <- function(eta, eta_vb, kappa, kappa_vb, log_tau_vb, tau_vb) {
  
  sum((eta - eta_vb) * log_tau_vb - (kappa - kappa_vb) * tau_vb +
        eta * log(kappa) - eta_vb * log(kappa_vb) - lgamma(eta) + lgamma(eta_vb))
  
}


############################################
## E log p(theta | rest) - E log q(theta) ##
############################################

# S0_inv is assumed to be block-diagonal
e_theta_ <- function(mu_theta_vb, list_S0_inv, list_sig2_theta_vb, 
                     vec_sum_log_det, vec_fac_bl = NULL, vec_fac_bl_y = NULL) {
  
    p <- length(mu_theta_vb)
    
    if (is.null(vec_fac_bl)) {

        arg <- vec_sum_log_det / 2 - # vec_sum_log_det[bl] = log(det(S0_inv_bl)) + log(det(sig2_theta_vb_bl))
          sum(list_S0_inv * mu_theta_vb^2 +
                list_S0_inv * list_sig2_theta_vb) / 2 + p / 2 # trace of a product
    
    } else if (!is.null(vec_fac_bl_y)){
      
      bl_ids_x <- as.numeric(levels(vec_fac_bl))
      n_bl_x <- length(bl_ids_x)
      
      vec_p_bl <- table(vec_fac_bl)
      
      bl_ids_y <- as.numeric(levels(vec_fac_bl_y))
      n_bl_y <- length(bl_ids_y)
      
        arg <-  sapply(1:n_bl_y, function(bl_y) {
          
          sapply(1:n_bl_x, function(bl_x) {
            
            vec_sum_log_det[bl_x, bl_y] / 2 - 
              sum(list_S0_inv[vec_fac_bl == bl_ids_x[bl_x], bl_y] * mu_theta_vb[vec_fac_bl == bl_ids_x[bl_x], bl_y]^2 +
                    list_S0_inv[vec_fac_bl == bl_ids_x[bl_x], bl_y] * list_sig2_theta_vb[vec_fac_bl == bl_ids_x[bl_x], bl_y] ) / 2 + vec_p_bl[bl_x] / 2 # trace of a product
            
          }) 
        })
        
    } else {
      
      bl_ids <- as.numeric(levels(vec_fac_bl))
      n_bl <- length(bl_ids)
      
      vec_p_bl <- table(vec_fac_bl)
      
      arg <- sapply(1:n_bl, function(bl) {
        vec_sum_log_det[bl] / 2 - # vec_sum_log_det[bl] = log(det(S0_inv_bl)) + log(det(sig2_theta_vb_bl))
          sum(list_S0_inv[vec_fac_bl == bl_ids[bl]] * mu_theta_vb[vec_fac_bl == bl_ids[bl]]^2 +
                list_S0_inv[vec_fac_bl == bl_ids[bl]] * list_sig2_theta_vb[vec_fac_bl == bl_ids[bl]]) / 2 + vec_p_bl[bl] / 2 # trace of a product
        
      })
    
    }
  
  sum(arg)
  
}


#######################
## E log p(y | rest) ##
#######################

e_y_ <- function(n, kappa, kappa_vb, log_tau_vb, m2_beta, sig2_inv_vb, tau_vb,
                 m2_alpha = NULL, zeta2_inv_vb = NULL) {
  
  arg <- -n / 2 * log(2 * pi) + n / 2 * log_tau_vb - tau_vb *
    (kappa_vb - colSums(m2_beta) * sig2_inv_vb / 2 - kappa)
  
  if (!is.null(m2_alpha))
    arg <- arg + tau_vb * crossprod(m2_alpha, zeta2_inv_vb) / 2
  
  sum(arg)
  
}
