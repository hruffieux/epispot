epispot_dual_info_vbem_core_ <- function(Y, X, V, list_hyper, gam_vb, mu_beta_vb,
                                         sig2_beta_vb, tau_vb, om, list_struct, 
                                         bool_blocks, hs, df, eb_local_scale,
                                         tol, maxit, anneal, anneal_vb_em, verbose,
                                         adaptive_tol_em = FALSE) {

  if (eb_local_scale & !is.null(list_struct))
    stop("EB local scale not implemented with structured sparsity")
  
  r <- ncol(V)
  
  maxit_em <- ceiling(maxit / 5)
  converged_em <- FALSE
  it_em <- 0
  tol_em <- min(tol * 10, 1) # maximum within-EM VB tol = 1
  s2_min <- 1e-12
  lb_old <- -Inf
  
  if (adaptive_tol_em) {
    # maximum within-EM VB tol = 10
    vec_tol_vb_within_em <- unique(sapply(tol_em * c(1:5, 10, 20, 50, 100, 500, 1000), 
                                          function(tol_adapt) min(tol_adapt, 10))) 
    ind_status_conv <- 1
  } 
  tol_vb_within_em <- tol_em
  
  list_hyper$om_vb <- om #rep(1 / 2, r) # prior proportion of active annotations 
  
  vb <- create_named_list_(gam_vb, mu_beta_vb, sig2_beta_vb, tau_vb)
  
  
  while ((!converged_em) & list_hyper$s2 > s2_min & (it_em < maxit_em)) {
    
    it_em <- it_em + 1
    
    if (verbose) {
      cat("---------- VB updates ----------\n")
      cat(paste0(ifelse(adaptive_tol_em, "Adaptive within", "Within"), 
                 "-EM VB tolerance: ", format(tol_vb_within_em, digits = 3), "\n"))
    }

    if (hs) {
      vb <- epispot_dual_horseshoe_info_core_(Y, X, V, list_hyper, vb$gam_vb, 
                                            vb$mu_beta_vb, vb$sig2_beta_vb, 
                                            vb$tau_vb, df, list_struct, eb = TRUE, 
                                            tol_vb_within_em, maxit, anneal_vb_em, verbose = FALSE, 
                                            full_output = TRUE)
      list_hyper$s02 <- NULL
      
    } else {
      vb <- epispot_dual_info_core_(Y, X, V, list_hyper, vb$gam_vb, vb$mu_beta_vb,
                                  vb$sig2_beta_vb, vb$tau_vb, list_struct, 
                                  eb = TRUE, eb_local_scale, tol_vb_within_em, maxit, anneal_vb_em,
                                  verbose = FALSE, full_output = TRUE)
      
      if (eb_local_scale) {
        list_hyper$s02 <- vb$sig2_theta_vb + vb$mu_theta_vb^2 # vector of size lenght(vb$mu_theta_vb)
      } else {
        list_hyper$s02 <- sum(vb$sig2_theta_vb + vb$mu_theta_vb^2) / ncol(X) # scalaer
      }

    }
   
    
    if (verbose)
      cat("--- EM hyperparameter updates ---\n")
    
    list_hyper$om_vb <- vb$zeta_vb
    list_hyper$s2 <- sum(vb$zeta_vb * (vb$sig2_c_vb + vb$mu_c_vb^2)) / sum(vb$zeta_vb)
  
    if (verbose) {
      
      cat(paste0("EM iteration ", it_em, ". \n"))
      
      if (!is.null(list_hyper$s02)) {
        
        if (eb_local_scale) {
          cat("New value for hyperparameter s02 : ")
          print(summary(list_hyper$s02))
        } else {
          cat(paste0("New value for hyperparameter s02 : ", format(list_hyper$s02, digits = 4),". \n"))
        }
      }
      cat(paste0("New value for hyperparameter s2 : ", format(list_hyper$s2, digits = 4), ". \n"))
      cat("New values for hyperparameter omega : \n")
      print(summary(list_hyper$om_vb))
      
      cat("\n\n")
    }
    
    diff_lb <- abs(vb$lb_opt - lb_old)
    
    if (adaptive_tol_em) {
    
      sum_exceed <- sum(diff_lb > vec_tol_vb_within_em)
      
      if (sum_exceed == 0) {
        
        converged_em <- TRUE
        
      } else if (it_em == 1 | ind_status_conv > sum_exceed) {
        
        ind_status_conv <- sum_exceed
        tol_vb_within_em <- vec_tol_vb_within_em[ind_status_conv]
        
      }
    
    } else {
      
      converged_em <- (diff_lb < tol_em)
      
    }
    
    lb_old <- vb$lb_opt
  
  }
  
  if (bool_blocks) {
    
    out <- list("s2" = list_hyper$s2, "om" = list_hyper$om_vb, "s02" = list_hyper$s02)
    
  } else {
    
    if (verbose) {
      if (converged_em) {
        cat(paste0("Convergence of the EM hyperparameter optimization run obtained after ", format(it_em), " EM iterations. \n\n"))
      } else if (list_hyper$s2 <= s2_min) {
        cat(paste0("EM hyperparameter optimization run stopped after s2 getting below ", s2_min, ". \n\n"))
      } else if (!is.null(list_hyper$s02) && list_hyper$s02 <= s2_min) {
        cat(paste0("EM hyperparameter optimization run stopped after s02 getting below ", s2_min, ". \n\n"))
      } else {
        warning("Maximal number of EM iterations reached before convergence. Exit EM run. \n\n")
      }
      
      cat("======= Final VB run =======\n") 
      cat(paste0("Empirical-Bayes hyperparameters, s2 : ", format(list_hyper$s2, digits = 4), ", omega :\n"))
      print(summary(list_hyper$om_vb))
      
      if (!is.null(list_hyper$s02)) {
        
        if (eb_local_scale) {
          cat("s02 : ")
          print(summary(list_hyper$s02))
        } else {
          cat(paste0("s02 : ", format(list_hyper$s02, digits = 4),". \n"))
        }
        cat("\n\n")
      }

    }
  
    if (hs) {
      out <- epispot_dual_horseshoe_info_core_(Y, X, V, list_hyper, vb$gam_vb, 
                                             vb$mu_beta_vb, vb$sig2_beta_vb, 
                                             vb$tau_vb, df, list_struct, 
                                             eb = TRUE, tol, maxit, anneal, verbose)
    } else {
      out <- epispot_dual_info_core_(Y, X, V, list_hyper, vb$gam_vb, vb$mu_beta_vb,
                                   vb$sig2_beta_vb, vb$tau_vb, list_struct, 
                                   eb = TRUE, eb_local_scale, tol, maxit, anneal, verbose)
    }
   
  
    if (!is.null(list_hyper$s02)) {
      out$s02 <- list_hyper$s02
    }
    
    out$s2 <- list_hyper$s2
    
  }
  
  out
  
}

