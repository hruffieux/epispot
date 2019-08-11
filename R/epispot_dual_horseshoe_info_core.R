# This file is part of the `epispot` R package:
#     https://github.com/hruffieux/epispot
#
# Internal core function to call the variational algorithm for dual propensity
# control. Sparse regression with identity link, no fixed covariates.
# See help of `epispot` function for details.
#
epispot_dual_horseshoe_info_core_ <- function(Y, X, V, list_hyper, gam_vb, mu_beta_vb, 
                                            sig2_beta_vb, tau_vb, df, list_struct, eb,
                                            tol, maxit, anneal, verbose, batch = "y", 
                                            full_output = FALSE, debug = TRUE,
                                            trace_path = NULL) {
 
  # Y must have been centered, and X standardized.
  
  d <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  r <- ncol(V)

  # Preparing trace saving if any
  #
  trace_ind_max <- trace_var_max <- NULL
  
  with(list_hyper, { # list_init not used with the with() function to avoid
                     # copy-on-write for large objects

    shr_fac_inv <- d # = 1 / shrinkage_factor for global variance
    
    # Preparing annealing if any
    #
    if (df == 1) {
      anneal_scale <- TRUE # if TRUE, scale parameters s02 and bhs_vb also annealed.
    } else {
      anneal_scale <- FALSE # annealed bhs_vb updates not sable for df > 1
    }
    
    if (df > 1 && !is.null(anneal) && anneal_scale) {
      stop("Annealing for scale parameters may not be numerically stable with df > 1.")
    }
    
    if (is.null(anneal)) {
      annealing <- FALSE
      c <- c_s <- 1 # c_s for scale parameters
    } else {
      annealing <- TRUE
      ladder <- get_annealing_ladder_(anneal, verbose)
      c <- ladder[1]
      c_s <- ifelse(anneal_scale, c, 1)
    }
    
    eps <- .Machine$double.eps^0.5
    
    # Covariate-specific parameters: objects derived from s02, list_struct (possible block-wise in parallel)
    #
    if (is.null(list_struct)) { # /! here list_struct is used to define block-specific s02, not to inject structure!
      n_bl <- 1
      bl_lgths <- p
    } else {
      vec_fac_bl <- list_struct$vec_fac_st
      bl_ids <- list_struct$bl_ids <- unique(vec_fac_bl)
      n_bl <- list_struct$n_bl <- length(bl_ids) 
      bl_lgths <- list_struct$bl_lgths <- table(vec_fac_bl)
    }
    
    # Variance initialization
    #
    S0_inv_vb <- rgamma(n_bl, shape = sapply(bl_lgths, function(lgth) max(lgth, d)), rate = 1) 
    
    # Some hyperparameters
    #
    A2_inv <- 1 #  hyperparameter # TODO: see how to fix, sensitivity analysis
    
    # Choose m0 so that, `a priori' (i.e. before optimization), E_p_gam is as specified by the user. 
    # In fact, we assume that the variance of theta (s0^2 in the hyperparameter doc) 
    # is very small so that the shift is negligeable: we set m0 to 0.
    #
    m0 <- rep(0, p)
    
    # Parameter initialization here for the top level 
    #
    mu_theta_vb <- rnorm(p, sd = 1 / sqrt(S0_inv_vb[1] * shr_fac_inv)) 
    sig2_theta_vb <- 1 / (d + rgamma(p, shape = S0_inv_vb[1] * shr_fac_inv, rate = 1)) # initial guess assuming b_vb = 1
    
    mu_rho_vb <- rnorm(d, mean = n0, sd = sqrt(t02))
    mu_c_vb <- rnorm(r, sd = 0.1) 
    
    
    if (eb) {
      a <- b <- a_vb <- b_vb <- NULL
      
      zeta_vb <- rbeta(r, shape1 = om_vb + eps, shape2 = 1 - om_vb + eps)
      
      log_om_vb <- log(om_vb + eps)
      log_1_min_om_vb <- log(1 - om_vb + eps)
      
    } else {
      zeta_vb <- rbeta(r, shape1 = a, shape2 = b)
    }
    
    # Response-specific parameters: objects derived from t02
    #
    T0_inv <- 1 / t02
    sig2_rho_vb <- update_sig2_c0_vb_(p, t02, c = c) # stands for a diagonal matrix of size d with this value on the (constant) diagonal
    
    vec_sum_log_det_rho <- - d * (log(t02) + log(p + T0_inv))
    
    
    # External information effects
    #
    sig2_c_vb <- update_sig2_c_vb_(p, s2, d, c = c)
    
    
    # Stored/precomputed objects
    #
    m1_beta <- update_m1_beta_(gam_vb, mu_beta_vb)
    m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = TRUE)
    m1_c <- update_m1_beta_(zeta_vb, mu_c_vb)
    
    mat_x_m1 <- update_mat_x_m1_(X, m1_beta)
    mat_v_mu <- update_mat_v_mu_(V, mu_theta_vb, m1_c, mu_rho_vb)
    
    # Fixed VB parameter
    #
    lambda_a_inv_vb <- 1 # no change with annealing 
    
    converged <- FALSE
    lb_new <- -Inf
    it <- 0
    
    
    while ((!converged) & (it < maxit)) {
      
      lb_old <- lb_new
      it <- it + 1
      
      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste("Iteration ", format(it), "... \n", sep = ""))
      
      if (!eb) digam_sum <- digamma(c * (a + b + 1) - 2 * c + 2)
      
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
      #
      if (batch == "y") { # optimal scheme
        
        log_Phi_mat_v_mu <- pnorm(mat_v_mu, log.p = TRUE)
        
        log_1_min_Phi_mat_v_mu <- pnorm(mat_v_mu, lower.tail = FALSE, log.p = TRUE)
        
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
      
      W <- update_W_info_(gam_vb, mat_v_mu, c = c) # we use info_ so that the second argument is a matrix
      
      # keep this order!
      #
      
      if (is.null(list_struct)) {
        
        G_vb <- c_s * S0_inv_vb * shr_fac_inv * (mu_theta_vb^2 + sig2_theta_vb - 2 * mu_theta_vb * m0 + m0^2) / 2 / df 
        nu_a_inv_vb <- c_s * (A2_inv + S0_inv_vb) 
        
      } else {
        
        G_vb <- unlist(lapply(1:n_bl, function(bl) { c_s * S0_inv_vb[bl] * shr_fac_inv * 
            (mu_theta_vb[vec_fac_bl == bl_ids[bl]]^2 + sig2_theta_vb[vec_fac_bl == bl_ids[bl]] - 
               2 * mu_theta_vb[vec_fac_bl == bl_ids[bl]] * m0[vec_fac_bl == bl_ids[bl]] + 
               m0[vec_fac_bl == bl_ids[bl]]^2) / 2 })) / df
        
        nu_a_inv_vb <- sapply(1:n_bl, function(bl) c_s * (A2_inv + S0_inv_vb[bl]))
        
      } 
      
      
      if (annealing & anneal_scale) {
        
        bhs_vb <- update_annealed_b_vb_(G_vb, c_s, df)
        
      } else {
        
        if (df == 1) {
          
          Q_app <- sapply(G_vb, function(G_vb_s) Q_approx(G_vb_s))  # TODO implement a Q_approx for vectors
          
          bhs_vb <- 1 / (Q_app * G_vb) - 1
          
        } else if (df == 3) {
          
          Q_app <- sapply(G_vb, function(G_vb_s) Q_approx(G_vb_s))
          
          bhs_vb <- exp(-log(3) - log(G_vb) + log(1 - G_vb * Q_app) - log(Q_app * (1 + G_vb) - 1)) - 1 / 3
          
        } else {
          # also works for df = 3 but might be slightly less efficient than the above
          
          Q_app <- sapply(G_vb, function(G_vb_s) Q_approx(G_vb_s))
          
          exponent <- (df + 1) / 2
          
          bhs_vb <- sapply(1:p, function(j) {
            
            exp(log(compute_integral_hs_(df, G_vb[j] * df, m = exponent, n = exponent, Q_ab = Q_app[j])) -
                  log(compute_integral_hs_(df, G_vb[j] * df, m = exponent, n = exponent - 1, Q_ab = Q_app[j])))
            
          })
          
        }
        
      }
      
      a_inv_vb <- lambda_a_inv_vb / nu_a_inv_vb
      
      
      if (is.null(list_struct)) {
        
        sig2_theta_vb <- update_sig2_c0_vb_(d, 1 / (S0_inv_vb * bhs_vb * shr_fac_inv), c = c)
        
        
        
        mat_v_mu <- sweep(mat_v_mu, 1, mu_theta_vb, `-`)
        mu_theta_vb <- update_mu_theta_vb_(W, m0, S0_inv_vb * bhs_vb * shr_fac_inv, sig2_theta_vb,
                                           vec_fac_st = NULL, mat_v_mu, is_mat = TRUE, c = c)
        
        lambda_s0_vb <- update_lambda_vb_(1 / 2, p, c = c_s)
        
        nu_s0_vb <- c_s * (a_inv_vb + 
                             sum(bhs_vb * shr_fac_inv * (mu_theta_vb^2 + sig2_theta_vb - 2 * mu_theta_vb * m0 + m0^2)) / 2) 
        
      } else {
        
        sig2_theta_vb <- unlist(lapply(1:n_bl, function(bl) { 
          update_sig2_c0_vb_(d, 1 / (S0_inv_vb[bl] * bhs_vb[vec_fac_bl == bl_ids[bl]] * shr_fac_inv), c = c) }))
        
        
        mat_v_mu <- sweep(mat_v_mu, 1, mu_theta_vb, `-`)
        
        mu_theta_vb <- unlist(lapply(1:n_bl, function(bl) {
          update_mu_theta_vb_(W[vec_fac_bl == bl_ids[bl], , drop = FALSE], m0[vec_fac_bl == bl_ids[bl]], 
                              S0_inv_vb[bl] * bhs_vb[vec_fac_bl == bl_ids[bl]] * shr_fac_inv, sig2_theta_vb[vec_fac_bl == bl_ids[bl]],
                                           vec_fac_st = NULL, mat_v_mu[vec_fac_bl == bl_ids[bl], , drop = FALSE], is_mat = TRUE, c = c)}))
        
        lambda_s0_vb <- sapply(1:n_bl, function(bl) update_lambda_vb_(1 / 2, bl_lgths[bl], c = c_s))
        
        nu_s0_vb <- sapply(1:n_bl, function(bl) { c_s * (a_inv_vb[bl] + 
                                                           sum(bhs_vb[vec_fac_bl == bl_ids[bl]] * shr_fac_inv * 
                                                                 (mu_theta_vb[vec_fac_bl == bl_ids[bl]]^2 + 
                                                                    sig2_theta_vb[vec_fac_bl == bl_ids[bl]] - 
                                                                    2 * mu_theta_vb[vec_fac_bl == bl_ids[bl]] * m0[vec_fac_bl == bl_ids[bl]] + 
                                                                    m0[vec_fac_bl == bl_ids[bl]]^2)) / 2)}) 
        
      }
      
      S0_inv_vb <- as.numeric(lambda_s0_vb / nu_s0_vb)
      
      
      mat_v_mu <- sweep(sweep(mat_v_mu, 1, mu_theta_vb, `+`), 2, mu_rho_vb, `-`)
      
      mu_rho_vb <- update_mu_rho_vb_(W, mat_v_mu, n0, sig2_rho_vb, T0_inv, 
                                     is_mat = TRUE, c = c)
      mat_v_mu <- sweep(mat_v_mu, 2, mu_rho_vb, `+`)
      
      
      if (batch == "y") { # optimal scheme
        
        if (!eb) {
          log_om_vb <- update_log_om_vb(a, digam_sum, zeta_vb, c = c)
          log_1_min_om_vb <- update_log_1_min_om_vb(b, 1, digam_sum, zeta_vb, c = c)
        }
        
        # # C++ Eigen call for expensive updates
        shuffled_ind_info <- as.numeric(sample(0:(r-1))) # Zero-based index in C++
        
        coreDualInfoLoop(V, W, zeta_vb, log_om_vb, log_1_min_om_vb, s2, m1_c,
                         mat_v_mu, mu_c_vb, sig2_c_vb, shuffled_ind_info, c = c)
        
      } else {
        
        for (l in sample(1:r)) {
          
          if (!eb) {
            log_om_vb <- update_log_om_vb(a, digam_sum, zeta_vb, c = c)
            log_1_min_om_vb <- update_log_1_min_om_vb(b, 1, digam_sum, zeta_vb, c = c)
          }
          
          mat_v_mu <- sweep(mat_v_mu, 1, V[, l] * m1_c[l], `-`)
          
          mu_c_vb[l] <- c * sig2_c_vb * sum(crossprod(W - mat_v_mu, V[, l]))
          
          zeta_vb[l] <- exp(-log_one_plus_exp_(c * (log_1_min_om_vb[l] - log_om_vb[l] +
                                                      log(s2) / 2 - log(sig2_c_vb) / 2 -
                                                      mu_c_vb[l] ^ 2 / (2 * sig2_c_vb))))
          
          m1_c[l] <- mu_c_vb[l] * zeta_vb[l]
          
          mat_v_mu <- sweep(mat_v_mu, 1, V[, l] * m1_c[l], `+`)
          
        }
        
      }
      
      if (!eb) {
        a_vb <- update_a_vb(a, zeta_vb, c = c)
        b_vb <- update_b_vb(b, 1, zeta_vb, c = c)
        om_vb <- a_vb / (a_vb + b_vb)
      } 
      
      
      if (verbose & (it == 1 | it %% 5 == 0)) {
        
        if (is.null(list_struct)) {
          cat(paste0("Updated global variance: ", format(nu_s0_vb / (lambda_s0_vb - 1) / shr_fac_inv, digits = 4), ".\n"))
          cat("Updated local variational parameter 1 / mu_bhs_vb for local variances: \n")
          print(summary(1 / bhs_vb))
          cat("\n")
        } else {
          cat("Updated block-specific global variances: \n")
          print(summary(nu_s0_vb / (lambda_s0_vb - 1) / shr_fac_inv))
          cat("\n")
          cat("Updated local variational parameter 1 / mu_bhs_vb for local variances: \n")
          print(summary(1 / bhs_vb))
          cat("\n")
        }
        
      }
      
      if (!is.null(trace_path) && (it == 1 | it %% 25 == 0)) {
        
        list_traces <- plot_trace_var_hs_(b_vb, S0_inv_vb, d, it, trace_ind_max, trace_var_max, trace_path)
        trace_ind_max <- list_traces$trace_ind_max
        trace_var_max <- list_traces$trace_var_max
        
      }
      
      if (annealing) {
        
        if (verbose & (it == 1 | it %% 5 == 0))
          cat(paste("Temperature = ", format(1 / c, digits = 4), "\n\n", sep = ""))
        
        sig2_rho_vb <- c * sig2_rho_vb
        sig2_c_vb <- c * sig2_c_vb
        
        c <- ifelse(it < length(ladder), ladder[it + 1], 1)
        c_s <- ifelse(anneal_scale, c, 1)
        
        sig2_rho_vb <- sig2_rho_vb / c
        sig2_c_vb <- sig2_c_vb / c
        
        if (isTRUE(all.equal(c, 1))) {
          
          annealing <- FALSE
          
          if (verbose)
            cat("** Exiting annealing mode. **\n\n")
        }
        
        
      } else {
        
        lb_new <- elbo_dual_horseshoe_info_(Y, V, a, a_vb, a_inv_vb, A2_inv, b, 
                                            bhs_vb, b_vb, eta, eta_vb, G_vb, 
                                            gam_vb, kappa, kappa_vb, lambda,
                                            lambda_vb, lambda_a_inv_vb, 
                                            lambda_s0_vb, m0, n0, mu_c_vb, 
                                            mu_rho_vb, mu_theta_vb, nu, nu_vb, 
                                            nu_a_inv_vb, nu_s0_vb, om_vb, Q_app, 
                                            sig2_beta_vb, S0_inv_vb, s2, sig2_c_vb, 
                                            sig2_theta_vb, sig2_inv_vb, sig2_rho_vb,
                                            T0_inv, tau_vb, zeta_vb, m1_beta, 
                                            m2_beta, mat_x_m1, mat_v_mu,
                                            vec_sum_log_det_rho, list_struct, df, 
                                            eb, shr_fac_inv)
        
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
    
    s02_vb <- nu_s0_vb / (lambda_s0_vb - 1) / shr_fac_inv
    
    
    if (full_output) { # for internal use only
      
      create_named_list_(a, a_vb, a_inv_vb, A2_inv, b, bhs_vb, b_vb, eta, eta_vb, 
                         G_vb, gam_vb, kappa, kappa_vb, lambda, lambda_vb, 
                         lambda_a_inv_vb, lambda_s0_vb, m0, n0, mu_beta_vb, mu_c_vb, 
                         mu_rho_vb, mu_theta_vb, nu, nu_vb, nu_a_inv_vb, nu_s0_vb, 
                         om_vb, Q_app, sig2_beta_vb, S0_inv_vb, s02_vb, s2, sig2_c_vb, 
                         sig2_theta_vb, sig2_inv_vb, sig2_rho_vb, T0_inv, tau_vb, 
                         zeta_vb, m1_beta, m2_beta, mat_x_m1, mat_v_mu, 
                         vec_sum_log_det_rho, list_struct, df, eb, lb_opt, shr_fac_inv)
      
    } else {
      
      names_x <- colnames(X)
      names_y <- colnames(Y)
      names_v <- colnames(V)
      
      rownames(gam_vb) <- names_x
      colnames(gam_vb) <- names_y
      
      names(mu_theta_vb) <- names_x
      names(mu_rho_vb) <- names_y
      names(bhs_vb) <- names_x
      
      names(zeta_vb) <- names_v
      names(mu_c_vb) <- names_v
      names(om_vb) <- names_v
      
      rownames(mat_v_mu) <- names_x
      colnames(mat_v_mu) <- names_y
      
      diff_lb <- abs(lb_opt - lb_old)
      
      create_named_list_(bhs_vb, gam_vb, mat_v_mu, mu_c_vb, mu_theta_vb, mu_rho_vb, 
                         om_vb, S0_inv_vb, s02_vb, zeta_vb, converged, it, lb_opt, 
                         diff_lb, df,  trace_ind_max, trace_var_max, shr_fac_inv)
      
    }
  })
  
}



# Internal function which implements the marginal log-likelihood variational
# lower bound (ELBO) corresponding to the `epispot_struct_core` algorithm.
#
elbo_dual_horseshoe_info_ <- function(Y, V, a, a_vb, a_inv_vb, A2_inv, b, 
                                      bhs_vb, b_vb, eta, eta_vb, G_vb, 
                                      gam_vb, kappa, kappa_vb, lambda,
                                      lambda_vb, lambda_a_inv_vb, 
                                      lambda_s0_vb, m0, n0, mu_c_vb, 
                                      mu_rho_vb, mu_theta_vb, nu, nu_vb, 
                                      nu_a_inv_vb, nu_s0_vb, om_vb, Q_app, 
                                      sig2_beta_vb, S0_inv_vb, s2, sig2_c_vb, 
                                      sig2_theta_vb, sig2_inv_vb, sig2_rho_vb,
                                      T0_inv, tau_vb, zeta_vb, m1_beta, 
                                      m2_beta, mat_x_m1, mat_v_mu, 
                                      vec_sum_log_det_rho, 
                                      list_struct, df, eb, shr_fac_inv) {

  
  n <- nrow(Y)
  d <- ncol(Y)

  # needed for monotonically increasing elbo.
  #
  eta_vb <- update_eta_vb_(n, eta, gam_vb)
  kappa_vb <- update_kappa_vb_(Y, kappa, mat_x_m1, m1_beta, m2_beta, sig2_inv_vb)
  
  lambda_vb <- update_lambda_vb_(lambda, sum(gam_vb))
  nu_vb <- update_nu_vb_(nu, m2_beta, tau_vb)
  
  log_tau_vb <- update_log_tau_vb_(eta_vb, kappa_vb)
  log_sig2_inv_vb <- update_log_sig2_inv_vb_(lambda_vb, nu_vb)
  
  log_S0_inv_vb <- update_log_sig2_inv_vb_(lambda_s0_vb, nu_s0_vb)
  log_a_inv_vb <- update_log_sig2_inv_vb_(lambda_a_inv_vb, nu_a_inv_vb)
  
  
  if (eb) {
    eps <- .Machine$double.eps^0.5
    log_om_vb <- log(om_vb + eps)
    log_1_min_om_vb <- log(1 - om_vb + eps)
  } else {
    log_om_vb <- digamma(a_vb) - digamma(a_vb + b_vb)
    log_1_min_om_vb <- digamma(b_vb) - digamma(a_vb + b_vb)
  }
  
  
  if (!is.null(list_struct)) {
    n_bl <- list_struct$n_bl
    bl_ids <- list_struct$bl_ids
    bl_lgths <- list_struct$bl_lgths
    vec_fac_bl <- list_struct$vec_fac_st
  }
  
  
  elbo_A <- e_y_(n, kappa, kappa_vb, log_tau_vb, m2_beta, sig2_inv_vb, tau_vb)
  
  
  if (is.null(list_struct)) {
    
    elbo_B <- e_beta_gamma_dual_info_(V, gam_vb, log_sig2_inv_vb, log_tau_vb,
                                      mat_v_mu, mu_c_vb, m2_beta,
                                      sig2_beta_vb, sig2_c_vb, sig2_rho_vb,
                                      sig2_theta_vb, sig2_inv_vb, tau_vb, zeta_vb)
    
    
    elbo_C <- e_theta_hs_(bhs_vb, G_vb, log_S0_inv_vb + log(shr_fac_inv), m0, mu_theta_vb, 
                          Q_app, S0_inv_vb * shr_fac_inv, sig2_theta_vb, df)
    
  } else {
    elbo_B <- sum(sapply(1:n_bl, function(bl) {
      e_beta_gamma_dual_info_(V[vec_fac_bl == bl_ids[bl], , drop = FALSE], 
                              gam_vb[vec_fac_bl == bl_ids[bl], , drop = FALSE], 
                              log_sig2_inv_vb, log_tau_vb,
                              mat_v_mu[vec_fac_bl == bl_ids[bl], , drop = FALSE], 
                              mu_c_vb,
                              m2_beta[vec_fac_bl == bl_ids[bl], , drop = FALSE],
                              sig2_beta_vb, sig2_rho_vb,
                              sig2_theta_vb[vec_fac_bl == bl_ids[bl]], 
                              sig2_inv_vb, tau_vb, zeta_vb)}))
    
    elbo_C <- sum(sapply(1:n_bl, function(bl) { 
      e_theta_hs_(bhs_vb[vec_fac_bl == bl_ids[bl]], G_vb[vec_fac_bl == bl_ids[bl]], 
                  log_S0_inv_vb[bl] + log(shr_fac_inv), m0[vec_fac_bl == bl_ids[bl]],
                  mu_theta_vb[vec_fac_bl == bl_ids[bl]], 
                  Q_app[vec_fac_bl == bl_ids[bl]], S0_inv_vb[bl] * shr_fac_inv,
                  sig2_theta_vb[vec_fac_bl == bl_ids[bl]], df)}))
    
  }
  
  elbo_D <- e_rho_(mu_rho_vb, n0, sig2_rho_vb, T0_inv, vec_sum_log_det_rho)
  
  elbo_E <- e_c_zeta_(log_om_vb, log_1_min_om_vb, mu_c_vb, s2, sig2_c_vb, zeta_vb)
  
  elbo_F <- e_tau_(eta, eta_vb, kappa, kappa_vb, log_tau_vb, tau_vb)
  
  elbo_G <- sum(e_sig2_inv_hs_(a_inv_vb, lambda_s0_vb, log_a_inv_vb, log_S0_inv_vb, nu_s0_vb, S0_inv_vb)) # S0_inv_vb
  
  elbo_H <- sum(e_sig2_inv_(1 / 2, lambda_a_inv_vb, log_a_inv_vb, A2_inv, nu_a_inv_vb, a_inv_vb)) # a_inv_vb
  
  elbo_I <- e_sig2_inv_(lambda, lambda_vb, log_sig2_inv_vb, nu, nu_vb, sig2_inv_vb)
  
  if (eb) {
    elbo_J <- 0
  } else {
    elbo_J <- e_omega_(a, a_vb, b, b_vb, log_om_vb, log_1_min_om_vb)
  }
  
  elbo_A + elbo_B + elbo_C + elbo_D + elbo_E + elbo_F + elbo_G + elbo_H + elbo_I + elbo_J
  
}

