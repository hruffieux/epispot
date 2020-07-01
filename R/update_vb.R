# This file is part of the `epispot` R package:
#     https://github.com/hruffieux/epispot
#
# Internal functions gathering the variational updates for the core algorithms.
# Besides improving code readability via modular programming, the main purpose
# is to avoid copy-and-paste programming, as most of these updates (or slightly
# modified versions) are used more than once in the different core algorithms.
# For this reason, we choose to create functions for most variational updates,
# even for those consisting in very basic operations.
# Note that we don't modularize the body of the core for loops for performance
# reasons.


#####################
## alpha's updates ##
#####################

update_m2_alpha_ <- function(mu_alpha_vb, sig2_alpha_vb, sweep = FALSE) {
  
  if(sweep) {
    
    sweep(mu_alpha_vb ^ 2, 1, sig2_alpha_vb, `+`)
    
  } else {
    
    sig2_alpha_vb + mu_alpha_vb ^ 2
  }
  
}


update_sig2_alpha_vb_ <- function(n, zeta2_inv_vb, tau_vb = NULL, intercept = FALSE, c = 1) {
  
  den <- n - 1 + zeta2_inv_vb
  
  if (intercept)
    den[1] <- den[1] + 1 # the first column of Z was not scaled, it is the intercept.
  
  if (is.null(tau_vb)) {
    
    1 / (c * den)
    
  } else {
    
    1 / (c * tcrossprod(den, as.matrix(tau_vb)))
    
  }
  
}


update_sig2_alpha_logit_vb_ <- function(Z, psi_vb, zeta2_inv_vb) {
  
  1 / sweep(2 * crossprod(Z ^ 2, psi_vb), 1, zeta2_inv_vb, `+`)
  
}


update_mat_z_mu_ <- function(Z, mu_alpha_vb) Z %*% mu_alpha_vb


####################
## beta's updates ##
####################

update_m1_beta_ <- function(gam_vb, mu_beta_vb) gam_vb * mu_beta_vb


update_g_m1_beta_ <- function(list_mu_beta_vb, gam_vb) {
  
  G <- length(list_mu_beta_vb)
  
  lapply(1:G, function(g) sweep(list_mu_beta_vb[[g]], 2, gam_vb[g, ], `*`))
  
}


update_m2_beta_ <- function(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = FALSE) {
  
  if(sweep) {
    
    sweep(mu_beta_vb ^ 2, 2, sig2_beta_vb, `+`) * gam_vb
    
  } else {
    
    (mu_beta_vb ^ 2 + sig2_beta_vb) * gam_vb
    
  }
  
}


update_sig2_beta_vb_ <- function(n, sig2_inv_vb, tau_vb = NULL, c = 1) {
  
  if(is.null(tau_vb)) {
    
    1 / (c * (n - 1 + sig2_inv_vb))
    
  } else {
    
    1 / (c * (n - 1 + sig2_inv_vb) * tau_vb)
    
  }
}


update_sig2_beta_logit_vb_ <- function(X, psi_vb, sig2_inv_vb) {
  
  1 / (2 * crossprod(X ^ 2, psi_vb) + sig2_inv_vb)
  
}


update_mat_x_m1_ <- function(X, m1_beta) X %*% m1_beta


update_g_mat_x_m1_ <- function(list_X, list_m1_beta) {
  
  G <- length(list_X)
  
  Reduce(`+`, lapply(1:G, function(g) list_X[[g]] %*% list_m1_beta[[g]]))
}


update_g_m1_btb_ <- function(gam_vb, list_mu_beta_vb, list_sig2_beta_star, tau_vb) { ## not list_sig2_beta_star_inv!
  
  d <- length(tau_vb)
  G <- length(list_mu_beta_vb)
  
  lapply(1:G, function(g) {
    gam_vb[g, ]^2 * colSums(list_mu_beta_vb[[g]]^2) + # colSums(A^2) = diag(crossprod(A))
      gam_vb[g, ] * (sum(diag(as.matrix(list_sig2_beta_star[[g]]))) / tau_vb +
                       sapply(1:d, function(k) (1-gam_vb[g, k]) * sum(list_mu_beta_vb[[g]][, k]^2))) # tr(mu_gt mu_gt^T) = sum(mu_gt^2)
  })
  
}


update_g_m1_btXtXb_ <- function(list_X, gam_vb, list_mu_beta_vb, list_sig2_beta_star, tau_vb) {
  
  d <- length(tau_vb)
  G <- length(list_mu_beta_vb)
  
  lapply(1:G, function(g) {
    gam_vb[g, ]^2 * colSums((list_X[[g]] %*% list_mu_beta_vb[[g]])^2) +
      gam_vb[g, ] * (sum(crossprod(list_X[[g]]) * list_sig2_beta_star[[g]]) / tau_vb +
                       sapply(1:d, function(k) (1-gam_vb[g, k]) * sum(crossprod(list_X[[g]]) * tcrossprod(list_mu_beta_vb[[g]][, k])))) # tr(AB^T) = sum_ij A_ij B_ij
  })
}



########################
## c0 and c's updates ##
########################

update_sig2_c0_vb_ <- function(d, s02, c = 1) 1 / (c * (d + (1/s02)))


update_sig2_c_vb_ <- function(p, s2, d = 1, c = 1) 1 / (c * (d * (p - 1) + (1/s2)))


update_mat_v_mu_ <- function(V, mu_0_s, mat_c, mu_0_t = NULL, resp_spec = FALSE) { # mu_0_s = mu_theta_vb, mu_0_t = mu_rho_vb, mat_c = zeta * mu_c_vb
  
  if (is.null(mu_0_t)) {
    sweep(V %*% mat_c, 1, mu_0_s, `+`)
  } else {
    
    if (resp_spec) {
      sweep(sweep(V %*% mat_c, 1, mu_0_s, `+`), 2, mu_0_t, `+`)
    } else {
      d <- length(mu_0_t)
      sweep(tcrossprod(mu_0_s + V %*% mat_c, rep(1, d)), 2, mu_0_t, `+`)
    }
    
  }
  
}


update_mat_v_mu_block_ <- function(list_V, mu_0_s, mu_0_t, list_mat_c, vec_fac_bl) {
  
  d <- length(mu_0_t)
  
  bl_ids <- as.numeric(levels(vec_fac_bl))
  n_bl <- length(bl_ids)
  
  list_bl <- lapply(1:n_bl, function(bl) {
    
    sweep(tcrossprod(as.vector(mu_0_s[vec_fac_bl == bl_ids[bl]] + list_V[[bl]] %*% list_mat_c[[bl]]), rep(1, d)), 2, mu_0_t, `+`)
    
  })
  
  as.matrix(plyr::rbind.fill.matrix(list_bl))
  
}


update_mat_v_mu_block_modules_ <- function(list_V, mu_0_s, mu_0_t, list_mat_c, vec_fac_bl_x, vec_fac_bl_y) {
  
  bl_ids_x <- as.numeric(levels(vec_fac_bl_x))
  n_bl_x <- length(bl_ids_x)
  
  bl_ids_y <- as.numeric(levels(vec_fac_bl_y))
  n_bl_y <- length(bl_ids_y)
  
  vec_d_bl <- table(vec_fac_bl_y)
  
  cbind_fill_matrix(lapply(1:n_bl_y, function(bl_y) {
    
    if (!is.null(dim(mu_0_s))) {
      mu_0_s <- mu_0_s[, bl_y]
    }
    
    plyr::rbind.fill.matrix(lapply(1:n_bl_x, function(bl_x) {
      
      sweep(tcrossprod(as.vector(list_V[[bl_x]] %*% list_mat_c[[bl_x]][, bl_y] + mu_0_s[vec_fac_bl_x == bl_ids_x[bl_x]]), rep(1, vec_d_bl[bl_y])),
            2, mu_0_t[vec_fac_bl_y == bl_ids_y[bl_y]], `+`)
    }))
    
  }))
  
  
}


###################
## chi's updates ##
###################

update_chi_vb_ <- function(X, Z, m1_beta, m2_beta, mat_x_m1, mat_z_mu, sig2_alpha_vb) {
  
  sqrt(X^2 %*% m2_beta + mat_x_m1^2 - X^2 %*% m1_beta^2 + Z^2 %*% sig2_alpha_vb +
         mat_z_mu^2 + 2 * mat_x_m1 * mat_z_mu)
}


update_psi_logit_vb_ <- function(chi_vb) {
  
  exp(log(exp(log_sigmoid_(chi_vb)) - 1 / 2) - log(2 * chi_vb))
  
}


#####################
## rho's updates ##
#####################


update_mu_rho_vb_ <- function(W, mat_add, n0, sig2_rho_vb, T0_inv, is_mat = FALSE, c = 1) {
  
  
  if (is_mat) {
    as.vector(c * sig2_rho_vb * (colSums(W) + T0_inv * n0 - colSums(mat_add))) # mat_add <- sweep(mat_v_mu, 1, mu_rho_vb, `-`)
  } else {
    # as.vector(sig2_rho_vb %*% (colSums(W) + T0_inv %*% n0 - sum(mu_theta_vb)))
    # sig2_rho_vb and T0_inv is stored as a scalar which represents the value on the diagonal of the corresponding diagonal matrix
    as.vector(c * sig2_rho_vb * (colSums(W) + T0_inv * n0 - sum(mat_add))) # mat_add = mu_theta_vb
  }
  
}


update_sig2_rho_vb_ <- function(p, T0_inv) {
  
  # sig2_rho_vb and T0_inv are stored as scalars which represent the value on the diagonal of the corresponding diagonal matrix
  1 / (T0_inv + p)
  # as.matrix(solve(T0_inv + diag(p, nrow(T0_inv))))
  
}


#####################
## sigma's updates ##
#####################

update_lambda_vb_ <- function(lambda, sum_gam, c = 1) c * (lambda + sum_gam / 2) - c + 1


update_g_lambda_vb_ <- function(lambda, g_sizes, rs_gam) lambda + sum(g_sizes * rs_gam) / 2


update_nu_vb_ <- function(nu, m2_beta, tau_vb, c = 1) c * as.numeric(nu + crossprod(tau_vb, colSums(m2_beta)) / 2)


update_g_nu_vb_ <- function(nu, list_m1_btb, tau_vb) nu + sum(tau_vb * Reduce(`+`, list_m1_btb))/2


update_nu_bin_vb_ <- function(nu, m2_beta) nu + sum(m2_beta) / 2


update_log_sig2_inv_vb_ <- function(lambda_vb, nu_vb) digamma(lambda_vb) - log(nu_vb)


###################
## tau's updates ##
###################

update_eta_vb_ <- function(n, eta, gam_vb, c = 1) c * (eta + n / 2 + colSums(gam_vb) / 2) - c + 1


update_g_eta_vb_ <- function(n, eta, g_sizes, gam_vb) eta + n / 2 + as.numeric(crossprod(gam_vb, g_sizes)) / 2


update_eta_z_vb_ <- function(n, q, eta, gam_vb, c = 1) c * (eta + n / 2 + colSums(gam_vb) / 2 + q / 2) - c + 1


update_kappa_vb_ <- function(Y, kappa, mat_x_m1, m1_beta, m2_beta, sig2_inv_vb, c = 1) {
  
  n <- nrow(Y)
  
  c * (kappa + (colSums(Y^2) - 2 * colSums(Y * mat_x_m1)  +
                  (n - 1 + sig2_inv_vb) * colSums(m2_beta) +
                  colSums(mat_x_m1^2) - (n - 1) * colSums(m1_beta^2))/ 2)
  
}


update_g_kappa_vb_ <- function(Y, list_X, kappa, list_m1_beta, list_m1_btb,
                               list_m1_btXtXb, mat_x_m1, sig2_inv_vb) {
  
  n <- nrow(Y)
  G <- length(list_m1_beta)
  
  # avoid using do.call() as can trigger node stack overflow
  kappa + (colSums(Y^2) - 2 * colSums(Y * mat_x_m1)  +
             Reduce(`+`, list_m1_btXtXb) +
             sig2_inv_vb * Reduce(`+`, list_m1_btb) +
             colSums(mat_x_m1^2) -
             Reduce(`+`, lapply(1:G, function(g) colSums((list_X[[g]] %*% list_m1_beta[[g]])^2) ))) / 2
}


update_kappa_z_vb_ <- function(Y, Z, kappa, mu_alpha_vb, m1_beta, m2_alpha,
                               m2_beta, mat_x_m1, mat_z_mu, sig2_inv_vb,
                               zeta2_inv_vb, intercept = FALSE, c = 1) {
  n <- nrow(Y)
  
  kappa_vb <- c * (kappa + (colSums(Y^2) - 2 * colSums(Y * (mat_x_m1 + mat_z_mu))  +
                              (n - 1 + sig2_inv_vb) * colSums(m2_beta) +
                              colSums(mat_x_m1^2) - (n - 1) * colSums(m1_beta^2) +
                              (n - 1) * colSums(m2_alpha) +
                              crossprod(m2_alpha, zeta2_inv_vb) +
                              colSums(mat_z_mu^2) - (n - 1) * colSums(mu_alpha_vb^2) +
                              2 * colSums(mat_x_m1 * mat_z_mu))/ 2)
  
  if (intercept)
    kappa_vb <- kappa_vb + c * (m2_alpha[1, ] - (mu_alpha_vb[1, ])^2) / 2
  
  kappa_vb
}


update_log_tau_vb_ <- function(eta_vb, kappa_vb) digamma(eta_vb) - log(kappa_vb)


#####################
## theta's updates ##
#####################

update_mu_theta_vb_ <- function(W, sig2_theta_vb,
                                mat_add = 0, is_mat = FALSE, c = 1, 
                                vec_fac_bl = NULL) {

    
    if (is.null(vec_fac_bl)) {
      
      if (is_mat) {
        
        mu_theta_vb <- c * sig2_theta_vb * (rowSums(W) - rowSums(mat_add)) # mat_add = sweep(mat_v_mu, 1, mu_theta_vb, `-`)
        
      } else {
        
        mu_theta_vb <- c * sig2_theta_vb * (rowSums(W) - sum(mat_add)) # mat_add = mu_rho_vb
        
      }
      
    } else { #if (is.null(vec_fac_bl_y)){
      
      bl_ids <- as.numeric(levels(vec_fac_bl))
      n_bl <- length(bl_ids)
      
      if (is_mat) {
        
        unlist(sapply(1:n_bl, function(bl) c * sig2_theta_vb[vec_fac_bl == bl_ids[bl]] * (rowSums(W[vec_fac_bl == bl_ids[bl], , drop = FALSE])
                                                                                          - rowSums(mat_add[vec_fac_bl == bl_ids[bl], , drop = FALSE])))) # mat_add = sweep(mat_v_mu, 1, mu_theta_vb, `-`)
        
      } else {
        
        unlist(sapply(1:n_bl, function(bl) c * sig2_theta_vb[vec_fac_bl == bl_ids[bl]] * (rowSums(W[vec_fac_bl == bl_ids[bl], , drop = FALSE]) 
                                                                                          - sum(mat_add)))) # mat_add = mu_rho_vb
        
      }
      
    }
  
}


update_sig2_theta_vb_ <- function(d, p, s02, c = 1) {
  
    S0_inv <- 1 / s02 # stands for a diagonal matrix of size p with this value on the (constant) diagonal
    sig2_theta_vb <- as.numeric(update_sig2_c0_vb_(d, s02, c = c)) # idem
    
    vec_sum_log_det_theta <- - sum(log(s02) + log(d + S0_inv))
  
  create_named_list_(S0_inv, sig2_theta_vb, vec_sum_log_det_theta)
}



#################
## W's updates ##
#################

update_W_info_ <- function(gam_vb, mat_v_mu, log_1_pnorm, log_pnorm, c = 1) {
  
  sqrt_c <- sqrt(c)
  
  if (!isTRUE(all.equal(c, 1))) {
    
    sqrt_c <- sqrt(c)
    
    log_pnorm <- pnorm(sqrt_c * mat_v_mu, log.p = TRUE)
    log_1_pnorm <- pnorm(sqrt_c * mat_v_mu, log.p = TRUE, lower.tail = FALSE)
    
  } else {
    
    sqrt_c <- 1
  }
  
  imr0 <- inv_mills_ratio_(0, sqrt_c * mat_v_mu, log_1_pnorm, log_pnorm)
  
  (gam_vb * (inv_mills_ratio_(1, sqrt_c * mat_v_mu, log_1_pnorm, log_pnorm) - imr0) + imr0) / sqrt_c + mat_v_mu
  
}

update_W_probit_ <- function(Y, mat_z_mu, mat_x_m1) {
  
  mat_z_mu + mat_x_m1 + inv_mills_ratio_matrix_(Y, mat_z_mu + mat_x_m1)
  
}


update_W_struct_ <- function(gam_vb, mu_theta_vb) {
  
  log_pnorm <- pnorm(mu_theta_vb, log.p = TRUE)
  log_1_pnorm <- pnorm(mu_theta_vb, log.p = TRUE, lower.tail = FALSE)
  
  imr0 <- inv_mills_ratio_(0, mu_theta_vb, log_1_pnorm, log_pnorm)
  
  sweep(sweep(gam_vb, 1, (inv_mills_ratio_(1, mu_theta_vb, log_1_pnorm, log_pnorm) - imr0), `*`),
        1,  mu_theta_vb + imr0, `+`)
  
}



####################
## zeta's updates ##
####################


update_phi_z_vb_ <- function(phi, d, c = 1) c * (phi + d / 2) - c + 1


update_xi_z_vb_ <- function(xi, tau_vb, m2_alpha, c = 1) c * (xi + m2_alpha %*% tau_vb / 2)


update_xi_bin_vb_ <- function(xi, m2_alpha) xi + rowSums(m2_alpha) / 2


update_log_zeta2_inv_vb_ <- function(phi_vb, xi_vb) digamma(phi_vb) - log(xi_vb)
