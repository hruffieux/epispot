# This file is part of the `epispot` R package:
#     https://github.com/hruffieux/epispot
#
# Internal core function to call the variational algorithm for dual propensity
# control with external information variables. Sparse regression with identity
# link, no fixed covariates. See help of `epispot` function for details.
#
epispot_core_ <- function(Y, X, V, list_hyper, gam_vb, mu_beta_vb,
                                  om, s02, s2, sig2_beta_vb, tau_vb, tol, 
                                  maxit, anneal_schedule, verbose, batch = "y",
                                  full_output = FALSE, debug = TRUE) {
  
  # Y centered, and X and V standardized.
  
  q <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  r <- ncol(V)
  
  
  with(list_hyper, { # list_init not used with the with() function to avoid
                     # copy-on-write for large objects
    
    # Preparing annealing if any
    #
    if (is.null(anneal_schedule)) {
      annealing <- FALSE
      c <- 1
    } else {
      annealing <- TRUE
      ladder <- get_annealing_ladder_(anneal_schedule, verbose)
      c <- ladder[1]
    }
    
    eps <- .Machine$double.eps^0.5
    
    # Parameter initialization here for the top level only
    #
    theta_vb <- rnorm(p, sd = 0.1) 
    zeta_vb <- rnorm(q, mean = n0, sd = sqrt(t02))
    mu_xi_vb <- rnorm(r, sd = 0.1) 
      
      rho_vb <- rbeta(r, shape1 = om + eps, shape2 = 1 - om + eps)
      
      log_om <- log(om + eps)
      log_1_min_om <- log(1 - om + eps)
    
    # Covariate-specific parameters: objects derived from s02
    #
    obj_theta_vb <- update_sig2_theta_vb_(q, p, s02, c = c)
    
    S0_inv <- obj_theta_vb$S0_inv
    sig2_theta_vb <- obj_theta_vb$sig2_theta_vb
    vec_sum_log_det_theta <- obj_theta_vb$vec_sum_log_det_theta
    
    # Response-specific parameters: objects derived from t02
    #
    T0_inv <- 1 / t02
    sig2_zeta_vb <- update_sig2_xi0_vb_(p, t02, c = c) # stands for a diagonal matrix of size q with this value on the (constant) diagonal
    vec_sum_log_det_zeta <- - q * (log(t02) + log(p + T0_inv))
    
    
    # External information effects
    #
    sig2_xi_vb <- update_sig2_xi_vb_(p, s2, q, c = c)
    
    # Stored/precomputed objects
    #
    beta_vb <- update_beta_vb_(gam_vb, mu_beta_vb)
    m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = TRUE)
    xi_vb <- update_beta_vb_(rho_vb, mu_xi_vb)
    
    X_beta_vb <- update_X_beta_vb_(X, beta_vb)
    mat_v_mu <- update_mat_v_mu_(V, theta_vb, xi_vb, zeta_vb)
    
    log_Phi_mat_v_mu <- pnorm(mat_v_mu, log.p = TRUE)
    log_1_min_Phi_mat_v_mu <- pnorm(mat_v_mu, lower.tail = FALSE, log.p = TRUE)
    
    converged <- FALSE
    lb_new <- -Inf
    it <- 0
    
    while ((!converged) & (it < maxit)) {
      
      lb_old <- lb_new
      it <- it + 1
      
      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste("Iteration ", format(it), "... \n", sep = ""))
      
      # % #
      lambda_vb <- update_lambda_vb_(lambda, sum(gam_vb), c = c)
      nu_vb <- update_nu_vb_(nu, m2_beta, tau_vb, c = c)
      
      sig2_inv_vb <- lambda_vb / nu_vb
      # % #
      
      # % #
      eta_vb <- update_eta_vb_(n, eta, gam_vb, c = c)
      kappa_vb <- update_kappa_vb_(Y, kappa, X_beta_vb, beta_vb, m2_beta, sig2_inv_vb, c = c)
      
      tau_vb <- eta_vb / kappa_vb
      # % #
      
      sig2_beta_vb <- update_sig2_beta_vb_(n, sig2_inv_vb, tau_vb, c = c)
      
      log_tau_vb <- update_log_tau_vb_(eta_vb, kappa_vb)
      log_sig2_inv_vb <- update_log_sig2_inv_vb_(lambda_vb, nu_vb)
      
      # different possible batch-coordinate ascent schemes:
      
      if (batch == "y") { # optimal scheme
        
        # log_Phi_mat_v_mu <- pnorm(mat_v_mu, log.p = TRUE)
        # log_1_min_Phi_mat_v_mu <- pnorm(mat_v_mu, lower.tail = FALSE, log.p = TRUE)
        
        # C++ Eigen call for expensive updates
        shuffled_ind <- as.numeric(sample(0:(p-1))) # Zero-based index in C++
        
        coreDualLoop(X, Y, gam_vb, log_Phi_mat_v_mu,
                     log_1_min_Phi_mat_v_mu, log_sig2_inv_vb,
                     log_tau_vb, beta_vb, X_beta_vb, mu_beta_vb,
                     sig2_beta_vb, tau_vb, shuffled_ind, c = c)
        
      } else if (batch == "0"){ # no batch, used only internally
        # schemes "x" of "x-y" are not batch concave
        # hence not implemented as they may diverge
        
        for (k in sample(1:q)) {
          
          for (j in sample(1:p)) {
            
            X_beta_vb[, k] <- X_beta_vb[, k] - X[, j] * beta_vb[j, k]
            
            mu_beta_vb[j, k] <- c * sig2_beta_vb[k] * tau_vb[k] * crossprod(Y[, k] - X_beta_vb[, k], X[, j])
            
            gam_vb[j, k] <- exp(-log_one_plus_exp_(c * (pnorm(mat_v_mu[j, k], lower.tail = FALSE, log.p = TRUE) -
                                                          pnorm(mat_v_mu[j, k], log.p = TRUE) -
                                                          log_tau_vb[k] / 2 - log_sig2_inv_vb / 2 -
                                                          mu_beta_vb[j, k] ^ 2 / (2 * sig2_beta_vb[k]) -
                                                          log(sig2_beta_vb[k]) / 2)))
            
            beta_vb[j, k] <- gam_vb[j, k] * mu_beta_vb[j, k]
            
            X_beta_vb[, k] <- X_beta_vb[, k] + X[, j] * beta_vb[j, k]
            
          }
        }
        
      } else {
        
        stop ("Batch scheme not defined. Exit.")
        
      }
      
      m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = TRUE)
      
      W <- update_W_info_(gam_vb, mat_v_mu, log_1_min_Phi_mat_v_mu, log_Phi_mat_v_mu, c = c) # we use info_ so that the second argument is a matrix
      
      mat_v_mu <- sweep(mat_v_mu, 1, theta_vb, `-`)
      theta_vb <- update_theta_vb_(W, sig2_theta_vb,
                                         mat_v_mu, is_mat = TRUE, 
                                         c = c)
      
      mat_v_mu <- sweep(sweep(mat_v_mu, 1, theta_vb, `+`), 2, zeta_vb, `-`)
      
      zeta_vb <- update_zeta_vb_(W, mat_v_mu, n0, sig2_zeta_vb, T0_inv, is_mat = TRUE, c = c)
      mat_v_mu <- sweep(mat_v_mu, 2, zeta_vb, `+`)
      
      
      if (batch == "y") { # optimal scheme
        
        # # C++ Eigen call for expensive updates
        shuffled_ind_info <- as.numeric(sample(0:(r-1))) # Zero-based index in C++
        
        coreDualInfoLoop(V, W, rho_vb, log_om, log_1_min_om, s2, xi_vb,
                         mat_v_mu, mu_xi_vb, sig2_xi_vb, shuffled_ind_info, c = c)
        
      } else {
        
        for (l in sample(1:r)) {
          
          mat_v_mu <- sweep(mat_v_mu, 1, V[, l] * xi_vb[l], `-`)
          
          mu_xi_vb[l] <- c * sig2_xi_vb * sum(crossprod(W - mat_v_mu, V[, l]))
          
          rho_vb[l] <- exp(-log_one_plus_exp_(c * (log_1_min_om[l] - log_om[l] +
                                                      log(s2) / 2 - log(sig2_xi_vb) / 2 -
                                                      mu_xi_vb[l] ^ 2 / (2 * sig2_xi_vb))))
          
          xi_vb[l] <- mu_xi_vb[l] * rho_vb[l]
          
          mat_v_mu <- sweep(mat_v_mu, 1, V[, l] * xi_vb[l], `+`)
          
        }
        
      }
      
      log_Phi_mat_v_mu <- pnorm(mat_v_mu, log.p = TRUE)
      log_1_min_Phi_mat_v_mu <- pnorm(mat_v_mu, lower.tail = FALSE, log.p = TRUE)
      
      if (annealing) {
        
        if (verbose & (it == 1 | it %% 5 == 0))
          cat(paste("Temperature = ", format(1 / c, digits = 4), "\n\n", sep = ""))
        
        sig2_theta_vb <- c * sig2_theta_vb
        sig2_zeta_vb <- c * sig2_zeta_vb
        sig2_xi_vb <- c * sig2_xi_vb
        
        c <- ifelse(it < length(ladder), ladder[it + 1], 1)
        
        sig2_theta_vb <- sig2_theta_vb / c
        sig2_zeta_vb <- sig2_zeta_vb / c
        sig2_xi_vb <- sig2_xi_vb / c
        
        if (isTRUE(all.equal(c, 1))) {
          
          annealing <- FALSE
          
          if (verbose)
            cat("** Exiting annealing mode. **\n\n")
        }
        
      } else {
        
        lb_new <- elbo_(Y, V, eta, eta_vb, gam_vb, 
                                  kappa, kappa_vb, lambda, lambda_vb, 
                                  log_1_min_Phi_mat_v_mu, log_Phi_mat_v_mu, n0, 
                                  mu_xi_vb, zeta_vb, theta_vb, nu, nu_vb, om,
                                  sig2_beta_vb, S0_inv, s2, sig2_xi_vb, sig2_theta_vb,
                                  sig2_inv_vb, sig2_zeta_vb, T0_inv, tau_vb, rho_vb, 
                                  beta_vb, m2_beta, X_beta_vb, 
                                  vec_sum_log_det_zeta, vec_sum_log_det_theta)
        
        if (verbose & (it == 1 | it %% 5 == 0))
          cat(paste("ELBO = ", format(lb_new), "\n\n", sep = ""))
        
        if (debug && lb_new + eps < lb_old)
          stop("ELBO not increasing monotonically. Exit. ")
        
        converged <- (abs(lb_new - lb_old) < tol)
        
      }
    }
    
    if (verbose) {
    
      if (converged) {
        cat(paste0("Convergence obtained after ", format(it), " iterations. \n",
                  "Optimal marginal log-likelihood variational lower bound ",
                  "(ELBO) = ", format(lb_new), ". \n\n"))
      } else {
        warning("Maximal number of iterations reached before convergence. Exit.")
      }
    }
    
    lb_opt <- lb_new
    
    if (full_output) { # for internal use only
      
      create_named_list_(eta, eta_vb, gam_vb, kappa, kappa_vb, 
                         lambda, lambda_vb, n0, mu_beta_vb, mu_xi_vb, 
                         zeta_vb, theta_vb, nu, nu_vb, om, sig2_beta_vb, 
                         S0_inv, s2, sig2_xi_vb, sig2_theta_vb, sig2_inv_vb, 
                         sig2_zeta_vb, T0_inv, tau_vb, rho_vb, beta_vb, m2_beta, 
                         X_beta_vb, mat_v_mu, vec_sum_log_det_zeta,
                         vec_sum_log_det_theta, lb_opt, it)
      
    } else {
      
      names_x <- colnames(X)
      names_y <- colnames(Y)
      names_v <- colnames(V)
      
      rownames(gam_vb) <- rownames(beta_vb) <- names_x
      colnames(gam_vb) <- colnames(beta_vb) <- names_y
      
      names(theta_vb) <- names_x
      names(zeta_vb) <- names_y
      
      names(rho_vb) <- names_v
      names(xi_vb) <- names_v
      
      rownames(mat_v_mu) <- names_x
      colnames(mat_v_mu) <- names_y
      
      diff_lb <- abs(lb_opt - lb_old)
      
      create_named_list_(beta_vb, gam_vb, xi_vb, rho_vb, theta_vb, zeta_vb,
                         converged, it, lb_opt, diff_lb)
      
    }
  })
  
}



# Internal function which implements the marginal log-likelihood variational
# lower bound (ELBO) corresponding to the `epispot_struct_core` algorithm.
#
elbo_ <- function(Y, V, eta, eta_vb, gam_vb, kappa, 
                            kappa_vb, lambda, lambda_vb, log_1_min_Phi_mat_v_mu, 
                            log_Phi_mat_v_mu, n0, mu_xi_vb, zeta_vb, 
                            theta_vb, nu, nu_vb, om, sig2_beta_vb, S0_inv, 
                            s2, sig2_xi_vb, sig2_theta_vb, sig2_inv_vb, 
                            sig2_zeta_vb, T0_inv, tau_vb, rho_vb, beta_vb,
                            m2_beta, X_beta_vb, vec_sum_log_det_zeta, vec_sum_log_det_theta) {
  
  n <- nrow(Y)
  
  # needed for monotonically increasing elbo.
  #
  eta_vb <- update_eta_vb_(n, eta, gam_vb)
  kappa_vb <- update_kappa_vb_(Y, kappa, X_beta_vb, beta_vb, m2_beta, sig2_inv_vb)
  
  lambda_vb <- update_lambda_vb_(lambda, sum(gam_vb))
  nu_vb <- update_nu_vb_(nu, m2_beta, tau_vb)
  
  log_tau_vb <- update_log_tau_vb_(eta_vb, kappa_vb)
  log_sig2_inv_vb <- update_log_sig2_inv_vb_(lambda_vb, nu_vb)
  
    eps <- .Machine$double.eps^0.5
    log_om <- log(om + eps)
    log_1_min_om <- log(1 - om + eps)
  
  elbo_A <- e_y_(n, kappa, kappa_vb, log_tau_vb, m2_beta, sig2_inv_vb, tau_vb)
  
  elbo_B <- e_beta_gamma_(V, gam_vb, log_sig2_inv_vb, log_tau_vb,
                                    log_1_min_Phi_mat_v_mu, log_Phi_mat_v_mu, 
                                    mu_xi_vb, m2_beta, sig2_beta_vb, sig2_xi_vb, 
                                    sig2_zeta_vb, sig2_theta_vb, sig2_inv_vb, 
                                    tau_vb, rho_vb)
  
  elbo_C <- e_theta_(theta_vb, S0_inv, sig2_theta_vb, vec_sum_log_det_theta)
  
  elbo_D <- e_zeta_(zeta_vb, n0, sig2_zeta_vb, T0_inv, vec_sum_log_det_zeta)
  
  elbo_E <- e_xi_rho_(log_om, log_1_min_om, mu_xi_vb, s2, sig2_xi_vb, rho_vb)
  
  elbo_F <- e_tau_(eta, eta_vb, kappa, kappa_vb, log_tau_vb, tau_vb)
  
  elbo_G <- e_sig2_inv_(lambda, lambda_vb, log_sig2_inv_vb, nu, nu_vb, sig2_inv_vb)
  
  as.numeric(elbo_A + elbo_B + elbo_C + elbo_D + elbo_E + elbo_F + elbo_G)
  
}

