# This file is part of the `epispot` R package:
#     https://github.com/hruffieux/epispot
#
# Internal core function to call the variational algorithm for dual propensity
# control with external information variables. Sparse regression with identity
# link, no fixed covariates. See help of `epispot` function for details.
#
epispot_dual_info_blocks_modules_core_ <- function(Y, X, list_V, vec_fac_bl_x,
                                                 vec_fac_bl_y, list_hyper, 
                                                 gam_vb, mu_beta_vb, sig2_beta_vb, 
                                                 tau_vb, tol, maxit, anneal, verbose, batch = "y", 
                                                 full_output = FALSE, debug = TRUE) {
  
  # EB s02 specific to each grid element (i.e., also modules), no choice because part of vbem procedure.
  
  # Y centered, and X and V standardized.
  
  d <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  
  with(list_hyper, { # list_init not used with the with() function to avoid
    # copy-on-write for large objects
    
    # Preparing annealing if any
    #
    if (is.null(anneal)) {
      annealing <- FALSE
      c <- 1
    } else {
      annealing <- TRUE
      ladder <- get_annealing_ladder_(anneal, verbose)
      c <- ladder[1]
    }
    
    eps <- .Machine$double.eps^0.5
    
    
    bl_ids_x <- as.numeric(levels(vec_fac_bl_x))
    n_bl_x <- length(bl_ids_x)
    
    vec_p_bl <- table(vec_fac_bl_x)
    
    bl_ids_y <- as.numeric(levels(vec_fac_bl_y))
    n_bl_y <- length(bl_ids_y)
    
    vec_d_bl <- table(vec_fac_bl_y)
    
    vec_r_bl <- sapply(list_V, function(V) ncol(V))# size of the V_bl (after removal of cst and coll annotations in each)
    
    # sanity checks
    #
    stopifnot(length(om_vb) == n_bl_x & nrow(s2) == n_bl_x)
    stopifnot(all(sapply(om_vb, function(om) ncol(om) == n_bl_y)) & ncol(s2) == n_bl_y)
    stopifnot(all(sapply(1:n_bl_x, function(bl_x) nrow(om_vb[[bl_x]]) == vec_r_bl[bl_x])))
    stopifnot(sum(sapply(list_V, function(V) nrow(V))) == p)
   
    # Choose m0 so that, `a priori' (i.e. before optimization), E_p_gam is as specified by the user. 
    # In fact, we assume that the variance of theta (s0^2 in the hyperparameter doc) 
    # is very small so that the shift is negligeable: we set m0 to 0.
    #
    m0 <- rep(0, p)
    
    # Parameter initialization here for the top level only
    #
    mu_theta_vb <- matrix(rnorm(p * n_bl_y, sd = 0.1), nrow = p)
    mu_rho_vb <- rnorm(d, mean = n0, sd = sqrt(t02))
    mu_c_vb <- lapply(vec_r_bl, function(r_bl) matrix(rnorm(r_bl * n_bl_y, sd = 0.1), nrow = r_bl)) 
    
    
    # om_vb is a list of length n_bl_x
    zeta_vb <- lapply(om_vb, function(om_bl) matrix(rbeta(length(om_bl), shape1 = om_bl + eps, shape2 = 1 - om_bl + eps), ncol = n_bl_y))
    
    log_om_vb <- lapply(om_vb, function(om_bl) log(om_bl + eps))
    log_1_min_om_vb <- lapply(om_vb, function(om_bl) log(1 - om_bl + eps))
    
      # Covariate-specific parameters: objects derived from s02
      #
      obj_theta_vb <- lapply(1:(n_bl_x * n_bl_y), function(bl) {
        
        bl_x <- ceiling(bl / n_bl_y)
        bl_y <- bl %% n_bl_y
        if (bl_y == 0)
          bl_y <- n_bl_y
        
        update_sig2_theta_vb_(vec_d_bl[bl_y], vec_p_bl[bl_x], s02[vec_fac_bl_x == bl_ids_x[bl_x], bl_y], c = c)
        
      })

      S0_inv <- do.call(cbind, lapply(obj_theta_vb, `[[`, "S0_inv"))
      sig2_theta_vb <- do.call(cbind, lapply(obj_theta_vb, `[[`, "sig2_theta_vb"))
      

    vec_sum_log_det_theta <- matrix(unlist(lapply(obj_theta_vb, `[[`, "vec_sum_log_det_theta")), nrow = n_bl_x, byrow = TRUE)
    
    # Response-specific parameters: objects derived from t02
    #
    T0_inv <- 1 / t02
    sig2_rho_vb <- update_sig2_c0_vb_(p, t02, c = c) # stands for a diagonal matrix of size d with this value on the (constant) diagonal
    vec_sum_log_det_rho <- - d * (log(t02) + log(p + T0_inv))
    
    
    # External information effects
    #
    sig2_c_vb <- matrix(sapply(1:n_bl_y, function(bl_y) {
      sapply(1:n_bl_x, function(bl_x) update_sig2_c_vb_(vec_p_bl[bl_x], s2[bl_x, bl_y], vec_d_bl[bl_y], c = c))
    }), nrow = n_bl_x)
    
  
    # Stored/precomputed objects
    #
    m1_beta <- update_m1_beta_(gam_vb, mu_beta_vb)
    m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = TRUE)
    m1_c <- lapply(1:n_bl_x, function(bl) update_m1_beta_(zeta_vb[[bl]], mu_c_vb[[bl]])) # matrix of size r x n_bl_x
    
    mat_x_m1 <- update_mat_x_m1_(X, m1_beta)
    mat_v_mu <- update_mat_v_mu_block_modules_(list_V, mu_theta_vb, mu_rho_vb, m1_c, vec_fac_bl_x, vec_fac_bl_y)
    
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
      kappa_vb <- update_kappa_vb_(Y, kappa, mat_x_m1, m1_beta, m2_beta, sig2_inv_vb, c = c)
      
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
                     log_tau_vb, m1_beta, mat_x_m1, mu_beta_vb,
                     sig2_beta_vb, tau_vb, shuffled_ind, c = c)
        
      } else if (batch == "0"){ # no batch, used only internally
        # schemes "x" of "x-y" are not batch concave
        # hence not implemented as they may diverge
        
        for (k in sample(1:d)) {
          
          for (j in sample(1:p)) {
            
            mat_x_m1[, k] <- mat_x_m1[, k] - X[, j] * m1_beta[j, k]
            
            mu_beta_vb[j, k] <- c * sig2_beta_vb[k] * tau_vb[k] * crossprod(Y[, k] - mat_x_m1[, k], X[, j])
            
            gam_vb[j, k] <- exp(-log_one_plus_exp_(c * (pnorm(mat_v_mu[j, k], lower.tail = FALSE, log.p = TRUE) -
                                                          pnorm(mat_v_mu[j, k], log.p = TRUE) -
                                                          log_tau_vb[k] / 2 - log_sig2_inv_vb / 2 -
                                                          mu_beta_vb[j, k] ^ 2 / (2 * sig2_beta_vb[k]) -
                                                          log(sig2_beta_vb[k]) / 2)))
            
            m1_beta[j, k] <- gam_vb[j, k] * mu_beta_vb[j, k]
            
            mat_x_m1[, k] <- mat_x_m1[, k] + X[, j] * m1_beta[j, k]
            
          }
        }
        
      } else {
        
        stop ("Batch scheme not defined. Exit.")
        
      }
      
      m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = TRUE)
      
      W <- update_W_info_(gam_vb, mat_v_mu, log_1_min_Phi_mat_v_mu, log_Phi_mat_v_mu, c = c) # we use info_ so that the second argument is a matrix
      
      for (bl_y in sample(1:n_bl_y)) {
        
        mat_v_mu[, vec_fac_bl_y == bl_ids_y[bl_y]] <- sweep(mat_v_mu[, vec_fac_bl_y == bl_ids_y[bl_y], drop = FALSE], 1, mu_theta_vb[, bl_y], `-`)
        
        mu_theta_vb[, bl_y] <- update_mu_theta_vb_(W[, vec_fac_bl_y == bl_ids_y[bl_y], drop = FALSE], 
                                                 m0, S0_inv[, bl_y], sig2_theta_vb[, bl_y],
                                                 mat_v_mu[, vec_fac_bl_y == bl_ids_y[bl_y], drop = FALSE], is_mat = TRUE, 
                                                 c = c, vec_fac_bl = vec_fac_bl_x)
      
      
        mat_v_mu[, vec_fac_bl_y == bl_ids_y[bl_y]] <- sweep(mat_v_mu[, vec_fac_bl_y == bl_ids_y[bl_y], drop = FALSE], 1, mu_theta_vb[, bl_y], `+`)
      }
      
      mat_v_mu <- sweep(mat_v_mu, 2, mu_rho_vb, `-`)
      
      mu_rho_vb <- update_mu_rho_vb_(W, mat_v_mu, n0, sig2_rho_vb, T0_inv, is_mat = TRUE, c = c)
      mat_v_mu <- sweep(mat_v_mu, 2, mu_rho_vb, `+`)
      
      
      if (batch == "y") {
        
        for (bl_x in sample(1:n_bl_x)) {
          
          for (bl_y in sample(1:n_bl_y)) {
              
            # # C++ Eigen call for expensive updates
            shuffled_ind_info_bl <- as.numeric(sample(0:(vec_r_bl[bl_x]-1))) # Zero-based index in C++
            
            mat_v_mu_bl <- mat_v_mu[vec_fac_bl_x == bl_ids_x[bl_x], vec_fac_bl_y == bl_ids_y[bl_y], drop = FALSE]
            zeta_vb_bl <- zeta_vb[[bl_x]][, bl_y]
            m1_c_bl <- m1_c[[bl_x]][, bl_y]
            mu_c_vb_bl <- mu_c_vb[[bl_x]][, bl_y]
            
            coreDualInfoLoop(list_V[[bl_x]],
                             W[vec_fac_bl_x == bl_ids_x[bl_x], vec_fac_bl_y == bl_ids_y[bl_y], drop = FALSE],
                             zeta_vb_bl,
                             log_om_vb[[bl_x]][, bl_y], log_1_min_om_vb[[bl_x]][, bl_y], s2[bl_x, bl_y],
                             m1_c_bl,
                             mat_v_mu_bl,
                             mu_c_vb_bl, sig2_c_vb[bl_x, bl_y],
                             shuffled_ind_info_bl, c = c)
            
            mat_v_mu[vec_fac_bl_x == bl_ids_x[bl_x], vec_fac_bl_y == bl_ids_y[bl_y]] <- mat_v_mu_bl
            zeta_vb[[bl_x]][, bl_y] <- zeta_vb_bl
            m1_c[[bl_x]][, bl_y] <- m1_c_bl
            mu_c_vb[[bl_x]][, bl_y] <- mu_c_vb_bl
            
          }
        }
        
      } else {
          
        for (bl_x in sample(1:n_bl_x)) {
          
          for (bl_y in sample(1:n_bl_y)) {
            
            for (l in sample(1:vec_r_bl[bl_x])) {
              
              mat_v_mu[vec_fac_bl_x == bl_ids_x[bl_x], vec_fac_bl_y == bl_ids_y[bl_y]] <- sweep(mat_v_mu[vec_fac_bl_x == bl_ids_x[bl_x], vec_fac_bl_y == bl_ids_y[bl_y], drop = FALSE], 1,
                                                                                                list_V[[bl_x]][, l] * m1_c[[bl_x]][l, bl_y], `-`)
              
              mu_c_vb[[bl_x]][l, bl_y] <- c * sig2_c_vb[bl_x, bl_y] *
                sum(crossprod(W[vec_fac_bl_x == bl_ids_x[bl_x], vec_fac_bl_y == bl_ids_y[bl_y], drop = FALSE] - mat_v_mu[vec_fac_bl_x == bl_ids_x[bl_x], vec_fac_bl_y == bl_ids_y[bl_y], drop = FALSE], list_V[[bl_x]][, l]))
              
              zeta_vb[[bl_x]][l, bl_y] <- exp(-log_one_plus_exp_(c * (log_1_min_om_vb[[bl_x]][l, bl_y] - log_om_vb[[bl_x]][l, bl_y] +
                                                                        log(s2[bl_x, bl_y]) / 2 - log(sig2_c_vb[bl_x, bl_y]) / 2 -
                                                                        mu_c_vb[[bl_x]][l, bl_y] ^ 2 / (2 * sig2_c_vb[bl_x, bl_y]))))
              
              m1_c[[bl_x]][l, bl_y] <- mu_c_vb[[bl_x]][l, bl_y] * zeta_vb[[bl_x]][l, bl_y]
              
              
              mat_v_mu[vec_fac_bl_x == bl_ids_x[bl_x], vec_fac_bl_y == bl_ids_y[bl_y]] <- sweep(mat_v_mu[vec_fac_bl_x == bl_ids_x[bl_x], vec_fac_bl_y == bl_ids_y[bl_y], drop = FALSE], 1,
                                                                                                list_V[[bl_x]][, l] * m1_c[[bl_x]][l, bl_y], `+`)
              
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
        sig2_rho_vb <- c * sig2_rho_vb
        sig2_c_vb <- c * sig2_c_vb
        
        c <- ifelse(it < length(ladder), ladder[it + 1], 1)
        
        sig2_theta_vb <- sig2_theta_vb / c
        sig2_rho_vb <- sig2_rho_vb / c
        sig2_c_vb <- sig2_c_vb / c
        
        if (isTRUE(all.equal(c, 1))) {
          
          annealing <- FALSE
          
          if (verbose)
            cat("** Exiting annealing mode. **\n\n")
        }
        
      } else {
        
        lb_new <- elbo_dual_info_blocks_modules_(Y, list_V, eta, eta_vb, gam_vb, kappa, 
                                                 kappa_vb, lambda, lambda_vb, 
                                                 log_1_min_om_vb, log_om_vb, 
                                                 log_1_min_Phi_mat_v_mu, log_Phi_mat_v_mu, 
                                                 m0, n0, 
                                                 mu_c_vb, mu_rho_vb, mu_theta_vb, nu, 
                                                 nu_vb, sig2_beta_vb, S0_inv, s2, 
                                                 sig2_c_vb, sig2_theta_vb, sig2_inv_vb, 
                                                 sig2_rho_vb, T0_inv, tau_vb, zeta_vb, 
                                                 m1_beta, m2_beta, mat_x_m1, 
                                                 vec_sum_log_det_rho,
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
                         lambda_vb, m0, n0, mu_c_vb, mu_rho_vb, mu_theta_vb, nu, nu_vb, om_vb,
                         sig2_beta_vb, S0_inv, s2, sig2_c_vb, sig2_theta_vb,
                         sig2_inv_vb, sig2_rho_vb, T0_inv, tau_vb, zeta_vb, m1_beta,
                         m2_beta, mat_x_m1, mat_v_mu, vec_sum_log_det_rho,
                         vec_sum_log_det_theta, vec_fac_bl_x, vec_fac_bl_y)
      
    } else {
      
      names_x <- colnames(X)
      names_y <- colnames(Y)
      names_v <- lapply(list_V, function(V) colnames(V))
      
      rownames(gam_vb) <- rownames(mu_beta_vb) <- names_x
      colnames(gam_vb) <- colnames(mu_beta_vb) <- names_y
      
      rownames(mu_theta_vb) <- names_x
      colnames(mu_theta_vb) <- paste0("module_", 1:n_bl_y)
      names(mu_rho_vb) <- names_y
      
      mu_c_vb <- lapply(1:n_bl_x, function(bl) {
        rownames(mu_c_vb[[bl]]) <- colnames(list_V[[bl]])
        colnames(mu_c_vb[[bl]]) <- paste0("module_", 1:n_bl_y)
        mu_c_vb[[bl]]})
      
      om_vb <- lapply(1:n_bl_x, function(bl) {
        rownames(om_vb[[bl]]) <- colnames(list_V[[bl]])
        colnames(om_vb[[bl]]) <- paste0("module_", 1:n_bl_y)
        om_vb[[bl]]})
      
      zeta_vb <- lapply(1:n_bl_x, function(bl) {
        rownames(zeta_vb[[bl]]) <- colnames(list_V[[bl]])
        colnames(zeta_vb[[bl]]) <- paste0("module_", 1:n_bl_y)
        zeta_vb[[bl]]})
      
      names(zeta_vb) <- names(mu_c_vb) <- names(om_vb) <- paste0("bl_", 1:n_bl_x)
      
      diff_lb <- abs(lb_opt - lb_old)
      
      create_named_list_(mu_beta_vb, mu_c_vb, om_vb, gam_vb, mu_theta_vb, mu_rho_vb, zeta_vb, 
                         converged, it, lb_opt, diff_lb)
      
    }
  })
  
}



# Internal function which implements the marginal log-likelihood variational
# lower bound (ELBO) corresponding to the `epispot_struct_core` algorithm.
#
elbo_dual_info_blocks_modules_ <- function(Y, list_V, eta, eta_vb, gam_vb, kappa, 
                                           kappa_vb, lambda, lambda_vb, 
                                           log_1_min_om_vb, log_om_vb, 
                                           log_1_min_Phi_mat_v_mu, log_Phi_mat_v_mu, 
                                           m0, n0, 
                                           mu_c_vb, mu_rho_vb, mu_theta_vb, nu, 
                                           nu_vb, sig2_beta_vb, S0_inv, s2, 
                                           sig2_c_vb, sig2_theta_vb, sig2_inv_vb, 
                                           sig2_rho_vb, T0_inv, tau_vb, zeta_vb, 
                                           m1_beta, m2_beta, mat_x_m1, 
                                           vec_sum_log_det_rho,
                                           vec_sum_log_det_theta, 
                                           vec_fac_bl_x, vec_fac_bl_y) {
  
  n <- nrow(Y)
  n_bl_x <- length(list_V)
  bl_ids_y <- as.numeric(levels(vec_fac_bl_y))
  n_bl_y <- length(bl_ids_y)
  
  d <- ncol(Y)
  
  # needed for monotonically increasing elbo.
  #
  eta_vb <- update_eta_vb_(n, eta, gam_vb)
  kappa_vb <- update_kappa_vb_(Y, kappa, mat_x_m1, m1_beta, m2_beta, sig2_inv_vb)
  
  lambda_vb <- update_lambda_vb_(lambda, sum(gam_vb))
  nu_vb <- update_nu_vb_(nu, m2_beta, tau_vb)
  
  log_tau_vb <- update_log_tau_vb_(eta_vb, kappa_vb)
  log_sig2_inv_vb <- update_log_sig2_inv_vb_(lambda_vb, nu_vb)
  
  elbo_A <- e_y_(n, kappa, kappa_vb, log_tau_vb, m2_beta, sig2_inv_vb, tau_vb)
  
  elbo_B <- e_beta_gamma_dual_info_(list_V, gam_vb, log_sig2_inv_vb, log_tau_vb,
                                    log_1_min_Phi_mat_v_mu, log_Phi_mat_v_mu, 
                                    mu_c_vb, m2_beta, sig2_beta_vb,
                                    sig2_c_vb, sig2_rho_vb, sig2_theta_vb,
                                    sig2_inv_vb, tau_vb, zeta_vb, 
                                    bool_modules = TRUE, 
                                    vec_fac_bl_y = vec_fac_bl_y,
                                    vec_fac_bl_theta = vec_fac_bl_x)
  
  
  elbo_C <- e_theta_(m0, mu_theta_vb, S0_inv, sig2_theta_vb,
                     vec_sum_log_det_theta, vec_fac_bl = vec_fac_bl_x, 
                     vec_fac_bl_y = vec_fac_bl_y)
  
  elbo_D <- e_rho_(mu_rho_vb, n0, sig2_rho_vb, T0_inv, vec_sum_log_det_rho)
  
  elbo_E <- sum(sapply(1:n_bl_x, function(bl_x) {
    sapply(1:n_bl_y, function(bl_y) {
      e_c_zeta_(log_om_vb[[bl_x]][, bl_y], log_1_min_om_vb[[bl_x]][, bl_y], 
                mu_c_vb[[bl_x]][, bl_y], s2[bl_x, bl_y], sig2_c_vb[bl_x, bl_y], zeta_vb[[bl_x]][, bl_y])
    })
  }))
  
  elbo_F <- e_tau_(eta, eta_vb, kappa, kappa_vb, log_tau_vb, tau_vb)
  
  elbo_G <- e_sig2_inv_(lambda, lambda_vb, log_sig2_inv_vb, nu, nu_vb, sig2_inv_vb)
  
  as.numeric(elbo_A + elbo_B + elbo_C + elbo_D + elbo_E + elbo_F + elbo_G) 
  
}

