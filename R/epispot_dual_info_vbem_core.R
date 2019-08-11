epispot_dual_info_vbem_core_ <- function(Y, X, V, list_hyper, gam_vb, mu_beta_vb,
                                       sig2_beta_vb, tau_vb, list_struct, 
                                       bool_blocks, hs, df, tol, maxit, anneal, verbose) {
  
  r <- ncol(V)
  
  maxit_em <- ceiling(maxit / 5)
  converged_em <- FALSE
  it_em <- 0
  tol_em <- min(tol * 10, 1) # 1e-3
  s2_min <- 1e-6
  lb_old <- -Inf
  list_hyper$om_vb <- rep(1 / 2, r) # prior proportion of active annotations 
  
  vb <- create_named_list_(gam_vb, mu_beta_vb, sig2_beta_vb, tau_vb)
  
  
  while ((!converged_em) & list_hyper$s2 > s2_min & (it_em < maxit_em)) {
    
    it_em <- it_em + 1
    
    if (verbose)
      cat("---------- VB updates ----------\n")
    
    if (hs) {
      vb <- epispot_dual_horseshoe_info_core_(Y, X, V, list_hyper, vb$gam_vb, 
                                            vb$mu_beta_vb, vb$sig2_beta_vb, 
                                            vb$tau_vb, df, list_struct, eb = TRUE, 
                                            tol_em, maxit, anneal, verbose = FALSE, 
                                            full_output = TRUE)
      list_hyper$s02 <- NULL
      
    } else {
      vb <- epispot_dual_info_core_(Y, X, V, list_hyper, vb$gam_vb, vb$mu_beta_vb,
                                  vb$sig2_beta_vb, vb$tau_vb, list_struct, 
                                  eb = TRUE, tol_em, maxit, anneal,
                                  verbose = FALSE, full_output = TRUE)
      
      list_hyper$s02 <- sum(vb$sig2_theta_vb + vb$mu_theta_vb^2) / ncol(X)
      
    }
   
    
    if (verbose)
      cat("--- EM hyperparameter updates ---\n")
    
    list_hyper$om_vb <- vb$zeta_vb
    list_hyper$s2 <- sum(vb$zeta_vb * (vb$sig2_c_vb + vb$mu_c_vb^2)) / sum(vb$zeta_vb)
  
    if (verbose) {
      
      cat(paste0("EM iteration ", it_em, ". \n"))
      
      if (!is.null(list_hyper$s02))
        cat(paste0("New value for hyperparameter s02 : ", format(list_hyper$s02, digits = 4),". \n"))
      
      cat(paste0("New value for hyperparameter s2 : ", format(list_hyper$s2, digits = 4), ". \n"))
      cat("New values for hyperparameter omega : \n")
      print(summary(list_hyper$om_vb))
      
      cat("\n\n")
    }
     
    converged_em <- (abs(vb$lb_opt - lb_old) < tol_em)
    
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
      
      if (!is.null(list_hyper$s02))
        cat(paste0("s02 : ", format(list_hyper$s02, digits = 4)," \n"))
      cat("\n\n")
    }
  
    if (hs) {
      out <- epispot_dual_horseshoe_info_core_(Y, X, V, list_hyper, vb$gam_vb, 
                                             vb$mu_beta_vb, vb$sig2_beta_vb, 
                                             vb$tau_vb, df, list_struct, 
                                             eb = TRUE, tol, maxit, anneal, verbose)
    } else {
      out <- epispot_dual_info_core_(Y, X, V, list_hyper, vb$gam_vb, vb$mu_beta_vb,
                                   vb$sig2_beta_vb, vb$tau_vb, list_struct, 
                                   eb = TRUE, tol, maxit, anneal, verbose)
    }
   
  
    if (!is.null(list_hyper$s02)) {
      out$s02 <- list_hyper$s02
    }
    
    out$s2 <- list_hyper$s2
    
  }
  
  out
  
}

