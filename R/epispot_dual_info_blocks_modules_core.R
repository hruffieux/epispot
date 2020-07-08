# This file is part of the `epispot` R package:
#     https://github.com/hruffieux/epispot
#
# Internal core function to call the variational algorithm for dual propensity
# control with external information variables. Sparse regression with identity
# link, no fixed covariates. See help of `epispot` function for details.
#
epispot_dual_info_blocks_modules_core_ <- function(Y, X, list_V, vec_fac_bl_x,
                                                 vec_fac_bl_y, module_names, list_hyper, 
                                                 gam_vb, mu_beta_vb, om, s02, s2, 
                                                 sig2_beta_vb, 
                                                 tau_vb, tol, maxit, anneal_schedule, verbose, batch = "y", 
                                                 full_output = FALSE, debug = TRUE) {
  
  # EB s02 specific to each grid element (i.e., also modules), no choice because part of vbem procedure.
  
  # Y centered, and X and V standardized.
  
  q <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  
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
    
    
    bl_ids_x <- as.numeric(levels(vec_fac_bl_x))
    n_bl_x <- length(bl_ids_x)
    
    vec_p_bl <- table(vec_fac_bl_x)
    
    bl_ids_y <- as.numeric(levels(vec_fac_bl_y))
    n_bl_y <- length(bl_ids_y)
    
    vec_q_bl <- table(vec_fac_bl_y)
    
    vec_r_bl <- sapply(list_V, function(V) ncol(V))# size of the V_bl (after removal of cst and coll annotations in each)
    
    # sanity checks
    #
    stopifnot(length(om) == n_bl_x & nrow(s2) == n_bl_x)
    stopifnot(all(sapply(om, function(om) ncol(om) == n_bl_y)) & ncol(s2) == n_bl_y)
    stopifnot(all(sapply(1:n_bl_x, function(bl_x) nrow(om[[bl_x]]) == vec_r_bl[bl_x])))
    stopifnot(sum(sapply(list_V, function(V) nrow(V))) == p)
    
    # Parameter initialization here for the top level only
    #
    theta_vb <- matrix(rnorm(p * n_bl_y, sd = 0.1), nrow = p)
    zeta_vb <- rnorm(q, mean = n0, sd = sqrt(t02))
    mu_xi_vb <- lapply(vec_r_bl, function(r_bl) matrix(rnorm(r_bl * n_bl_y, sd = 0.1), nrow = r_bl)) 
    
    
    # om is a list of length n_bl_x
    rho_vb <- lapply(om, function(om_bl) matrix(rbeta(length(om_bl), shape1 = om_bl + eps, shape2 = 1 - om_bl + eps), ncol = n_bl_y))
    
    log_om <- lapply(om, function(om_bl) log(om_bl + eps))
    log_1_min_om <- lapply(om, function(om_bl) log(1 - om_bl + eps))
    
      # Covariate-specific parameters: objects derived from s02
      #
      obj_theta_vb <- lapply(1:(n_bl_x * n_bl_y), function(bl) {
        
        bl_x <- ceiling(bl / n_bl_y)
        bl_y <- bl %% n_bl_y
        if (bl_y == 0)
          bl_y <- n_bl_y
        
        update_sig2_theta_vb_(vec_q_bl[bl_y], vec_p_bl[bl_x], s02[vec_fac_bl_x == bl_ids_x[bl_x], bl_y], c = c)
        
      })

      S0_inv <- do.call(cbind, lapply(obj_theta_vb, `[[`, "S0_inv"))
      sig2_theta_vb <- do.call(cbind, lapply(obj_theta_vb, `[[`, "sig2_theta_vb"))
      

    vec_sum_log_det_theta <- matrix(unlist(lapply(obj_theta_vb, `[[`, "vec_sum_log_det_theta")), nrow = n_bl_x, byrow = TRUE)
    
    # Response-specific parameters: objects derived from t02
    #
    T0_inv <- 1 / t02
    sig2_zeta_vb <- update_sig2_xi0_vb_(p, t02, c = c) # stands for a diagonal matrix of size q with this value on the (constant) diagonal
    vec_sum_log_det_zeta <- - q * (log(t02) + log(p + T0_inv))
    
    
    # External information effects
    #
    sig2_xi_vb <- matrix(sapply(1:n_bl_y, function(bl_y) {
      sapply(1:n_bl_x, function(bl_x) update_sig2_xi_vb_(vec_p_bl[bl_x], s2[bl_x, bl_y], vec_q_bl[bl_y], c = c))
    }), nrow = n_bl_x)
    
  
    # Stored/precomputed objects
    #
    beta_vb <- update_beta_vb_(gam_vb, mu_beta_vb)
    m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = TRUE)
    xi_vb <- lapply(1:n_bl_x, function(bl) update_beta_vb_(rho_vb[[bl]], mu_xi_vb[[bl]])) # matrix of size r x n_bl_x
    
    X_beta_vb <- update_X_beta_vb_(X, beta_vb)
    mat_v_mu <- update_mat_v_mu_block_modules_(list_V, theta_vb, zeta_vb, xi_vb, vec_fac_bl_x, vec_fac_bl_y)
    
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
      
      for (bl_y in sample(1:n_bl_y)) {
        
        mat_v_mu[, vec_fac_bl_y == bl_ids_y[bl_y]] <- sweep(mat_v_mu[, vec_fac_bl_y == bl_ids_y[bl_y], drop = FALSE], 1, theta_vb[, bl_y], `-`)
        
        theta_vb[, bl_y] <- update_theta_vb_(W[, vec_fac_bl_y == bl_ids_y[bl_y], drop = FALSE], 
                                                 sig2_theta_vb[, bl_y],
                                                 mat_v_mu[, vec_fac_bl_y == bl_ids_y[bl_y], drop = FALSE], is_mat = TRUE, 
                                                 c = c, vec_fac_bl = vec_fac_bl_x)
      
      
        mat_v_mu[, vec_fac_bl_y == bl_ids_y[bl_y]] <- sweep(mat_v_mu[, vec_fac_bl_y == bl_ids_y[bl_y], drop = FALSE], 1, theta_vb[, bl_y], `+`)
      }
      
      mat_v_mu <- sweep(mat_v_mu, 2, zeta_vb, `-`)
      
      zeta_vb <- update_zeta_vb_(W, mat_v_mu, n0, sig2_zeta_vb, T0_inv, is_mat = TRUE, c = c)
      mat_v_mu <- sweep(mat_v_mu, 2, zeta_vb, `+`)
      
      
      if (batch == "y") {
        
        for (bl_x in sample(1:n_bl_x)) {
          
          for (bl_y in sample(1:n_bl_y)) {
              
            # # C++ Eigen call for expensive updates
            shuffled_ind_info_bl <- as.numeric(sample(0:(vec_r_bl[bl_x]-1))) # Zero-based index in C++
            
            mat_v_mu_bl <- mat_v_mu[vec_fac_bl_x == bl_ids_x[bl_x], vec_fac_bl_y == bl_ids_y[bl_y], drop = FALSE]
            rho_vb_bl <- rho_vb[[bl_x]][, bl_y]
            xi_vb_bl <- xi_vb[[bl_x]][, bl_y]
            mu_xi_vb_bl <- mu_xi_vb[[bl_x]][, bl_y]
            
            coreDualInfoLoop(list_V[[bl_x]],
                             W[vec_fac_bl_x == bl_ids_x[bl_x], vec_fac_bl_y == bl_ids_y[bl_y], drop = FALSE],
                             rho_vb_bl,
                             log_om[[bl_x]][, bl_y], log_1_min_om[[bl_x]][, bl_y], s2[bl_x, bl_y],
                             xi_vb_bl,
                             mat_v_mu_bl,
                             mu_xi_vb_bl, sig2_xi_vb[bl_x, bl_y],
                             shuffled_ind_info_bl, c = c)
            
            mat_v_mu[vec_fac_bl_x == bl_ids_x[bl_x], vec_fac_bl_y == bl_ids_y[bl_y]] <- mat_v_mu_bl
            rho_vb[[bl_x]][, bl_y] <- rho_vb_bl
            xi_vb[[bl_x]][, bl_y] <- xi_vb_bl
            mu_xi_vb[[bl_x]][, bl_y] <- mu_xi_vb_bl
            
          }
        }
        
      } else {
          
        for (bl_x in sample(1:n_bl_x)) {
          
          for (bl_y in sample(1:n_bl_y)) {
            
            for (l in sample(1:vec_r_bl[bl_x])) {
              
              mat_v_mu[vec_fac_bl_x == bl_ids_x[bl_x], vec_fac_bl_y == bl_ids_y[bl_y]] <- sweep(mat_v_mu[vec_fac_bl_x == bl_ids_x[bl_x], vec_fac_bl_y == bl_ids_y[bl_y], drop = FALSE], 1,
                                                                                                list_V[[bl_x]][, l] * xi_vb[[bl_x]][l, bl_y], `-`)
              
              mu_xi_vb[[bl_x]][l, bl_y] <- c * sig2_xi_vb[bl_x, bl_y] *
                sum(crossprod(W[vec_fac_bl_x == bl_ids_x[bl_x], vec_fac_bl_y == bl_ids_y[bl_y], drop = FALSE] - mat_v_mu[vec_fac_bl_x == bl_ids_x[bl_x], vec_fac_bl_y == bl_ids_y[bl_y], drop = FALSE], list_V[[bl_x]][, l]))
              
              rho_vb[[bl_x]][l, bl_y] <- exp(-log_one_plus_exp_(c * (log_1_min_om[[bl_x]][l, bl_y] - log_om[[bl_x]][l, bl_y] +
                                                                        log(s2[bl_x, bl_y]) / 2 - log(sig2_xi_vb[bl_x, bl_y]) / 2 -
                                                                        mu_xi_vb[[bl_x]][l, bl_y] ^ 2 / (2 * sig2_xi_vb[bl_x, bl_y]))))
              
              xi_vb[[bl_x]][l, bl_y] <- mu_xi_vb[[bl_x]][l, bl_y] * rho_vb[[bl_x]][l, bl_y]
              
              
              mat_v_mu[vec_fac_bl_x == bl_ids_x[bl_x], vec_fac_bl_y == bl_ids_y[bl_y]] <- sweep(mat_v_mu[vec_fac_bl_x == bl_ids_x[bl_x], vec_fac_bl_y == bl_ids_y[bl_y], drop = FALSE], 1,
                                                                                                list_V[[bl_x]][, l] * xi_vb[[bl_x]][l, bl_y], `+`)
              
            }
          }
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
        
        lb_new <- elbo_dual_info_blocks_modules_(Y, list_V, eta, eta_vb, gam_vb, kappa, 
                                                 kappa_vb, lambda, lambda_vb, 
                                                 log_1_min_om, log_om, 
                                                 log_1_min_Phi_mat_v_mu, log_Phi_mat_v_mu, 
                                                 n0, 
                                                 mu_xi_vb, zeta_vb, theta_vb, nu, 
                                                 nu_vb, sig2_beta_vb, S0_inv, s2, 
                                                 sig2_xi_vb, sig2_theta_vb, sig2_inv_vb, 
                                                 sig2_zeta_vb, T0_inv, tau_vb, rho_vb, 
                                                 beta_vb, m2_beta, X_beta_vb, 
                                                 vec_sum_log_det_zeta,
                                                 vec_sum_log_det_theta, 
                                                 vec_fac_bl_x, vec_fac_bl_y)
        
        if (verbose & (it == 1 | it %% 5 == 0))
          cat(paste("ELBO = ", format(lb_new), "\n\n", sep = ""))
        
        if (debug && lb_new + eps < lb_old)
          stop("ELBO not increasing monotonically. Exit. ")
        
        converged <- (abs(lb_new - lb_old) < tol)
        
      }
    }
    
    
    if (verbose) {
      if (converged) {
        cat(paste("Convergence obtained after ", format(it), " iterations. \n",
                  "Optimal marginal log-likelihood variational lower bound ",
                  "(ELBO) = ", format(lb_new), ". \n\n", sep = ""))
      } else {
        warning("Maximal number of iterations reached before convergence. Exit.")
      }
    }
    
    lb_opt <- lb_new
    
    if (full_output) { # for internal use only
      
      create_named_list_(Y, list_V, eta, eta_vb, gam_vb, kappa, kappa_vb, lambda,
                         lambda_vb, n0, mu_xi_vb, zeta_vb, theta_vb, nu, nu_vb, om,
                         sig2_beta_vb, S0_inv, s2, sig2_xi_vb, sig2_theta_vb,
                         sig2_inv_vb, sig2_zeta_vb, T0_inv, tau_vb, rho_vb, beta_vb,
                         m2_beta, X_beta_vb, mat_v_mu, vec_sum_log_det_zeta,
                         vec_sum_log_det_theta, vec_fac_bl_x, vec_fac_bl_y, lb_opt, it)
      
    } else {
      
      names_x <- colnames(X)
      names_y <- colnames(Y)
      names_v <- lapply(list_V, function(V) colnames(V))

      rownames(gam_vb) <- rownames(beta_vb) <- names_x
      colnames(gam_vb) <- colnames(beta_vb) <- names_y
      
      rownames(theta_vb) <- names_x
      colnames(theta_vb) <- colnames(s02)
      names(zeta_vb) <- names_y
      
      xi_vb <- lapply(1:n_bl_x, function(bl) {
        rownames(xi_vb[[bl]]) <- colnames(list_V[[bl]])
        colnames(xi_vb[[bl]]) <- module_names
        xi_vb[[bl]]})
      
      rho_vb <- lapply(1:n_bl_x, function(bl) {
        rownames(rho_vb[[bl]]) <- colnames(list_V[[bl]])
        colnames(rho_vb[[bl]]) <- module_names
        rho_vb[[bl]]})
      
      if (n_bl_x > 1) {
        names(rho_vb) <- names(xi_vb) <- paste0("bl_", 1:n_bl_x)
      } else {
        rho_vb <- rho_vb[[1]]
        xi_vb <- xi_vb[[1]]
      }

      diff_lb <- abs(lb_opt - lb_old)
      
      create_named_list_(beta_vb, gam_vb, xi_vb, rho_vb, theta_vb, zeta_vb, 
                         converged, it, lb_opt, diff_lb)
      
    }
  })
  
}



# Internal function which implements the marginal log-likelihood variational
# lower bound (ELBO) corresponding to the `epispot_struct_core` algorithm.
#
elbo_dual_info_blocks_modules_ <- function(Y, list_V, eta, eta_vb, gam_vb, kappa, 
                                           kappa_vb, lambda, lambda_vb, 
                                           log_1_min_om, log_om, 
                                           log_1_min_Phi_mat_v_mu, log_Phi_mat_v_mu, 
                                           n0, 
                                           mu_xi_vb, zeta_vb, theta_vb, nu, 
                                           nu_vb, sig2_beta_vb, S0_inv, s2, 
                                           sig2_xi_vb, sig2_theta_vb, sig2_inv_vb, 
                                           sig2_zeta_vb, T0_inv, tau_vb, rho_vb, 
                                           beta_vb, m2_beta, X_beta_vb, 
                                           vec_sum_log_det_zeta,
                                           vec_sum_log_det_theta, 
                                           vec_fac_bl_x, vec_fac_bl_y) {
  
  n <- nrow(Y)
  n_bl_x <- length(list_V)
  bl_ids_y <- as.numeric(levels(vec_fac_bl_y))
  n_bl_y <- length(bl_ids_y)
  
  q <- ncol(Y)
  
  # needed for monotonically increasing elbo.
  #
  eta_vb <- update_eta_vb_(n, eta, gam_vb)
  kappa_vb <- update_kappa_vb_(Y, kappa, X_beta_vb, beta_vb, m2_beta, sig2_inv_vb)
  
  lambda_vb <- update_lambda_vb_(lambda, sum(gam_vb))
  nu_vb <- update_nu_vb_(nu, m2_beta, tau_vb)
  
  log_tau_vb <- update_log_tau_vb_(eta_vb, kappa_vb)
  log_sig2_inv_vb <- update_log_sig2_inv_vb_(lambda_vb, nu_vb)
  
  elbo_A <- e_y_(n, kappa, kappa_vb, log_tau_vb, m2_beta, sig2_inv_vb, tau_vb)
  
  elbo_B <- e_beta_gamma_dual_info_(list_V, gam_vb, log_sig2_inv_vb, log_tau_vb,
                                    log_1_min_Phi_mat_v_mu, log_Phi_mat_v_mu, 
                                    mu_xi_vb, m2_beta, sig2_beta_vb,
                                    sig2_xi_vb, sig2_zeta_vb, sig2_theta_vb,
                                    sig2_inv_vb, tau_vb, rho_vb, 
                                    bool_modules = TRUE, 
                                    vec_fac_bl_y = vec_fac_bl_y,
                                    vec_fac_bl_theta = vec_fac_bl_x)
  
  
  elbo_C <- e_theta_(theta_vb, S0_inv, sig2_theta_vb,
                     vec_sum_log_det_theta, vec_fac_bl = vec_fac_bl_x, 
                     vec_fac_bl_y = vec_fac_bl_y)
  
  elbo_D <- e_zeta_(zeta_vb, n0, sig2_zeta_vb, T0_inv, vec_sum_log_det_zeta)
  
  elbo_E <- sum(sapply(1:n_bl_x, function(bl_x) {
    sapply(1:n_bl_y, function(bl_y) {
      e_xi_rho_(log_om[[bl_x]][, bl_y], log_1_min_om[[bl_x]][, bl_y], 
                mu_xi_vb[[bl_x]][, bl_y], s2[bl_x, bl_y], sig2_xi_vb[bl_x, bl_y], rho_vb[[bl_x]][, bl_y])
    })
  }))
  
  elbo_F <- e_tau_(eta, eta_vb, kappa, kappa_vb, log_tau_vb, tau_vb)
  
  elbo_G <- e_sig2_inv_(lambda, lambda_vb, log_sig2_inv_vb, nu, nu_vb, sig2_inv_vb)
  
  as.numeric(elbo_A + elbo_B + elbo_C + elbo_D + elbo_E + elbo_F + elbo_G) 
  
}

