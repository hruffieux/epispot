# This file is part of the `epispot` R package:
#     https://github.com/hruffieux/epispot
#

#' Fit EPISPOT: annotation-driven approach for large-scale joint regression with
#' multiple responses.
#'
#' @param Y Response data matrix of dimension n x q, where n is the number of
#'   samples and q is the number of response variables; Y is centred prior to 
#'   the run.
#' @param X Input matrix of dimension n x p, where p is the number of candidate
#'   predictors. No intercept must be supplied; X is scaled prior to the run 
#'   - beware the interpretation of the regression estimates.
#' @param V Annotation matrix of dimension p x r, where r is the number of
#'   variables representing external information on the candidate predictors
#'   which may make their selection more or less likely. V is standardised prior
#'   to the run - beware the interpretation of the annotation effect estimates.
#'   Annotations with either concern or do not concern given predictors should 
#'   be coded as binary.
#' @param p0 Vector of size 2 whose entries are the prior expectation and 
#'   variance of the number of predictors associated with each response.
#'   Must be \code{NULL} if \code{list_init} and \code{list_hyper} are both 
#'   non-\code{NULL}.
#' @param anneal_schedule Parameters for annealing scheme. Must be a vector 
#'   whose first entry is the type of schedule: 1 = geometric spacing (default), 
#'   2 = harmonic spacing or 3 = linear spacing, the second entry is the initial 
#'   temperature (default is 2), and the third entry is the temperature grid 
#'   size (default is 10). If \code{NULL}, no annealing is performed.
#' @param list_modules An object of class "\code{modules}" containing settings 
#'   for parallel inference using modules of responses. Must be filled using
#'   the \code{\link{set_modules}} function or must be \code{NULL} for no
#'   partitioning into modules.
#' @param bin_annot_freq Minimal frequency for binary annotations (if any), 
#'                       i.e., annotation variables which concern less than 
#'                       \code{bin_annot_freq} x 100 % or more than 
#'                       (1 - \code{bin_annot_freq}) x 100 candidate predictors 
#'                       are removed prior to the analysis, as they may not be 
#'                       sufficiently informative. Default is 0.05. If set to 
#'                       NULL, no filter is applied.
#' @param list_hyper An object of class "\code{hyper}" containing the model
#'   hyperparameters. Must be filled using the \code{\link{set_hyper}}
#'   function or must be \code{NULL} for default hyperparameters.
#' @param list_init An object of class "\code{init}" containing the initial
#'   variational and EM parameters. Must be filled using the 
#'   \code{\link{set_init}} function or be \code{NULL} for a default 
#'   initialisation.
#' @param user_seed Seed set for reproducible default choices of hyperparameters
#'   (if \code{list_hyper} is \code{NULL}) and initial variational parameters
#'   (if \code{list_init} is \code{NULL}). Default is \code{NULL}, no
#'   seed set.
#' @param tol Tolerance for the stopping criterion.
#' @param adaptive_tol_em Boolean indicating whether the tolerance for the 
#'   within-EM variational runs should be controlled adaptively depending on the
#'   EM-convergence status in order to save computational time. Default is 
#'   \code{TRUE}.
#' @param maxit Maximum number of iterations allowed.
#' @param anneal_vbem Parameters for annealing scheme for the internal runs of 
#'   the variational EM algorithm. Default is geometric spacing, initial 
#'   temperature is 2 and grid size is 10. See \code{anneal_schedule}.
#' @param save_hyper If \code{TRUE}, the hyperparameters used for the model are
#'   saved as output.
#' @param save_init If \code{TRUE}, the initial variational parameters used for
#'   the inference are saved as output.
#' @param verbose If \code{TRUE}, messages are displayed during execution.
#'
#' @return An object of class "\code{epispot}" containing the following output:
#'  \item{beta_vb}{Estimated effect size matrix of dimension p x q. Entry (s, t) 
#'                 corresponds to the variational posterior mean 
#'                 (mu_beta_vb_st x gam_vb_st) of the regression effect between 
#'                 candidate predictor s and response t.}
#'  \item{gam_vb}{Posterior inclusion probability matrix of dimension p x q.
#'                Entry (s, t) corresponds to the posterior probability of
#'                association between candidate predictor s and response t.}
#'  \item{xi_vb}{Vector of size r, where entry l contains the estimated marginal
#'               effect of annotation l on the probabilities of association.}
#'  \item{rho_vb}{Posterior inclusion probability vector of size r for the
#'                annotation variables.}
#'  \item{theta_vb}{Vector of length p containing the posterior mean of theta. 
#'                  Entry s corresponds to the propensity of candidate predictor 
#'                  s to be included in the model.}
#'  \item{zeta_vb}{Vector of length q containing the posterior mean of zeta.
#'                 Entry t controls the proportion of predictors associated
#'                 with response t.}
#'  \item{converged}{A boolean indicating whether the algorithm has converged
#'                   before reaching \code{maxit} iterations.}
#'  \item{it}{Final number of iterations.}
#'  \item{lb_opt}{Optimized variational lower bound for the marginal
#'                log-likelihood (ELBO).}
#'  \item{diff_lb}{Difference in ELBO between the last and penultimate
#'                 iterations. This may be a useful diagnostic information when
#'                 convergence has not been reached before \code{maxit}.}
#'  \item{p0}{Vector of length two defining the applied sparsity control.}
#'  \item{rmvd_cst_x, rmvd_cst_v}{Vectors containing the indices of constant
#'                                variables in \code{X}, resp. \code{V}, removed
#'                                prior to the analysis. \code{NULL} if no 
#'                                constant variable removed.}
#'  \item{rmvd_coll_x, rmvd_cst_v}{Vectors containing the indices of variables
#'                                 in \code{X}, resp. \code{V}, removed prior
#'                                 to the analysis because collinear to other
#'                                 variables. The entry name indicates the
#'                                 corresponding variable kept in the analysis
#'                                 (i.e., that causing the collinearity for the
#'                                 entry in question). \code{NULL} if no 
#'                                 collinear variable removed.}
#'  \item{rmvd_bin_annot_freq_v}{If \code{bin_annot_freq} is non-NULL, vector
#'                               containing the indices of variables in \code{V}
#'                               removed prior to the analysis because they 
#'                               their frequency is lower than 
#'                               \code{bin_annot_freq} or higher than 
#'                               1-\code{bin_annot_freq}.}
#'  \item{list_hyper, list_init}{If \code{save_hyper}, resp. \code{save_init},
#'                               \code{TRUE}, hyperparameters, resp. initial
#'                               variational parameters, used for inference are
#'                               saved as output.}
#'  
#' @examples
#' seed <- 123; set.seed(seed)
#'
#' ###################
#' ## Simulate data ##
#' ###################
#'
#' ## Example using small problem sizes:
#' ##
#' n <- 50; p <- 60; p_act <- 10; q <- 25; q_act <- 15; r <- 10
#'
#' ## Candidate predictors (subject to selection)
#' ##
#' # Here example with common genetic variants under Hardy-Weinberg equilibrium
#' #
#' X_act <- matrix(rbinom(n * p_act, size = 2, p = 0.25), nrow = n)
#' X_inact <- matrix(rbinom(n * (p - p_act), size = 2, p = 0.25), nrow = n)
#'
#' # shuffle indices 
#' shuff_x_ind <- sample(p)
#' shuff_y_ind <- sample(q)
#' 
#' X <- cbind(X_act, X_inact)[, shuff_x_ind]
#'
#' # Association pattern and effect sizes
#' #
#' pat <- matrix(FALSE, ncol = q, nrow = p)
#' bool_x <- shuff_x_ind <= p_act
#' bool_y <- shuff_y_ind <= q_act
#' 
#' pat_act <- beta_act <- matrix(0, nrow = p_act, ncol = q_act)
#' pat_act[sample(p_act * q_act, floor(p_act * q_act / 5))] <- 1
#' beta_act[as.logical(pat_act)] <-  rnorm(sum(pat_act))
#' 
#' pat[bool_x, bool_y] <- pat_act
#' 
#' # Gaussian responses
#' #
#' Y_act <- matrix(rnorm(n * q_act, mean = X_act %*% beta_act), nrow = n)
#' Y_inact <- matrix(rnorm(n * (q - q_act)), nrow = n)
#'
#' Y <- cbind(Y_act, Y_inact)[, shuff_y_ind]
#'
#' # Annotation variables
#' #
#' V <- matrix(rnorm(p * r), nrow = p)
#' V[bool_x, ] <- rnorm(p_act * r, mean = 2)
#'
#'
#' ########################
#' ## Infer associations ##
#' ########################
#'
#' # Expectation and variance for the prior number of predictors associated with
#' # each response
#' #
#' p0 <- c(mean(colSums(pat)), 10)
#' 
#' # Inference without module
#' #
#' res_epispot <- epispot(Y = Y, X = X, V = V, p0 = p0, user_seed = seed)
#'
#' # Inference with modules
#' #
#' module_ids <- sample(c(rep(1, floor(q/2)), rep(2, q - floor(q/2)))) # 2 modules
#' list_modules <- set_modules(module_ids, n_cpus = 1)
#' 
#' res_epispot_modules <- epispot(Y = Y, X = X, V = V, p0 = p0, 
#'                                list_modules = list_modules, user_seed = seed)
#'
#' @references
#' H. Ruffieux, B. P. Fairfax, I. Nassiri, E. Vigorito, C. Walllace, S. 
#' Richardson, L. Bottolo. EPISPOT: an epigenome-driven approach for detecting 
#' and interpreting hotspots in molecular QTL studies, arXiv, 2020.
#'
#' @seealso \code{\link{set_hyper}}, \code{\link{set_init}}, 
#' \code{\link{set_modules}}.
#'
#' @export
#'
epispot <- function(Y, X, V, p0, anneal_schedule = c(1, 2, 10), 
                    list_modules = NULL, bin_annot_freq = 0.05, 
                    list_hyper = NULL, list_init = NULL, user_seed = NULL, 
                    tol = 0.1, adaptive_tol_em = TRUE, maxit = 1000, 
                    anneal_vbem = TRUE, save_hyper = FALSE, 
                    save_init = FALSE, verbose = TRUE) {
  
  if (verbose) cat("== Preparing the data ... \n")
  
  dat <- prepare_data_(Y, X, V, bin_annot_freq, user_seed, tol, maxit, verbose)
  
  X <- dat$X
  Y <- dat$Y
  V <- dat$V
  
  bool_rmvd_x <- dat$bool_rmvd_x
  bool_rmvd_v <- dat$bool_rmvd_v
  
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  r <- ncol(V)
  
  names_x <- colnames(X)
  names_y <- colnames(Y)
  names_v <- colnames(V)
  
  check_annealing_(anneal_schedule, maxit)
  
  if (anneal_vbem) {
    
    if (!is.null(anneal_schedule)) {
      anneal_schedule_vbem <- anneal_schedule
    } else {
      stop(paste0("anneal_vbem must be set to FALSE if no annealing is used ", 
                  "for the VB algorithm (i.e., anneal_schedule is NULL)"))
    }

  } else {
    
    anneal_schedule_vbem <- NULL
    
  }

  check_annealing_(anneal_schedule_vbem, maxit)
  
  if (verbose) cat("... done. == \n\n")
  
  if (!is.null(list_modules)) {
    
    list_modules <- prepare_blocks_(list_modules, q, bool_rmvd_x)
    
    if (!is.null(list_modules$order_y_ids)) { # for the case where modules are provided and responses are not grouped per module
      # we group them prior to the analysis and ungroup them after the run in epispot.R
      Y <- Y[, list_modules$order_y_ids, drop = FALSE]
      
    }
    
    n_bl_x <- list_modules$n_bl_x
    n_bl_y <- list_modules$n_bl_y
    
    n_cpus <- list_modules$n_cpus
    
    vec_fac_bl_x <- list_modules$vec_fac_bl_x
    vec_fac_bl_y <- list_modules$vec_fac_bl_y
    
  }
  
  if (is.null(list_hyper) | is.null(list_init)) {
    
    check_structure_(p0, "vector", "numeric", 2)
    check_positive_(p0) 
    
    if (p0[1] > p) {
      
      stop(paste0("The expected number of predictors associated with each ", 
                  "response must be lower than the total number of predictors. ",
                  "Please change p0[1]. Exit."))
      
    } else if (p0[1] > floor(p/2)) {
      
      warning(paste0("The expected number of predictors associated with each ", 
                     "response is large, which corresponds to a non-sparse ", 
                     "setting. Please consider lowering p0[1]."))
      
    }
    
    if (p0[2] > floor(p/2)) 
      warning(paste0("The prior variance for the number of predictors ", 
                     "associated with each response is large. Please consider ", 
                     "lowering p0[2]."))
    
    
  } else {
    
    if (!is.null(p0))
      warning(paste0("Provided argument p0 not used, as both list_hyper ",
                     "and list_init were provided."))
    
  }
  
  
  if (verbose) cat("== Preparing the hyperparameters ... \n\n")
  
  list_hyper <- prepare_list_hyper_(list_hyper, Y, p, p0, list_modules,
                                    bool_rmvd_x, names_x, names_y, verbose)
  
  if (verbose) cat("... done. == \n\n")
  
  if (verbose) cat("== Preparing the parameter initialization ... \n\n")
  
  list_init <- prepare_list_init_(list_init, Y, p, p0, r, list_modules,
                                  bool_rmvd_x, bool_rmvd_v, 
                                  user_seed, verbose)
  
  if (verbose) cat("... done. == \n\n")
  
  if (verbose){
    cat(paste0("===================================================== \n",
               "== EPISPOT: annotation-driven association mapping  == \n",
               "===================================================== \n\n"))
  }
  
  
  if (is.null(list_modules)) {
    
    vb <- epispot_vbem_core_(Y, X, V, list_hyper, list_init$gam_vb,
                                       list_init$mu_beta_vb, list_init$om,
                                       list_init$s02, list_init$s2,
                                       list_init$sig2_beta_vb, list_init$tau_vb,
                                       bool_blocks = FALSE, tol, maxit, anneal_schedule, 
                                       anneal_schedule_vbem, verbose,
                                       adaptive_tol_em = adaptive_tol_em)
    
  } else {
    
    list_pos_bl_x <- split(1:p, vec_fac_bl_x)
    
    split_bl_mat_x_ <- function(bl_x) {
      
      X_bl <- X[, list_pos_bl_x[[bl_x]], drop = FALSE]
      
      V_bl <- V[list_pos_bl_x[[bl_x]],, drop = FALSE]
      
      list_V_bl_bin_annot_freq <- rm_bin_annot_freq_(V_bl, bin_annot_freq, verbose = FALSE)
      V_bl <- list_V_bl_bin_annot_freq$mat
      bool_bin_annot_freq_v_bl <- list_V_bl_bin_annot_freq$bool_bin_annot_freq
      rmvd_bin_annot_freq_v_bl <- list_V_bl_bin_annot_freq$rmvd_bin_annot_freq
      
      V_bl <- scale(V_bl) # the VB algorithm assumes a scaled V
      
      list_V_bl_cst <- rm_constant_(V_bl, verbose = FALSE)
      V_bl <- list_V_bl_cst$mat
      bool_cst_v_bl <- list_V_bl_cst$bool_cst
      rmvd_cst_v_bl <- list_V_bl_cst$rmvd_cst
      
      list_V_bl_coll <- rm_collinear_(V_bl, verbose = FALSE)
      V_bl <- list_V_bl_coll$mat
      r <- ncol(V_bl)
      bool_coll_v_bl <- list_V_bl_coll$bool_coll
      rmvd_coll_v_bl <- list_V_bl_coll$rmvd_coll
      
      bool_rmvd_v_bl <- bool_bin_annot_freq_v_bl
      bool_rmvd_v_bl[!bool_bin_annot_freq_v_bl] <- bool_cst_v_bl
      bool_rmvd_v_bl[!bool_bin_annot_freq_v_bl][!bool_cst_v_bl] <- bool_coll_v_bl
      

      if (sum(!bool_rmvd_v_bl) == 0)
        stop(paste0("There exist one or more blocks for which no non-constant ",
                    "annotation variables remain. Try to use less blocks."))
      
      
      create_named_list_(X_bl, V_bl, bool_rmvd_v_bl, rmvd_cst_v_bl, rmvd_coll_v_bl, rmvd_bin_annot_freq_v_bl)
    }
    
    list_bl_mat_x <- parallel::mclapply(1:n_bl_x, function(bl_x) split_bl_mat_x_(bl_x), mc.cores = n_cpus)
    
    
    if (n_bl_y > 1) {
      
      list_pos_bl_y <- split(1:q, vec_fac_bl_y)
      
      split_bl_mat_y_ <- function(bl_y) {
        
        Y_bl <- Y[, list_pos_bl_y[[bl_y]], drop = FALSE]
        Y_bl
        
      }
      
      list_bl_mat_y <- parallel::mclapply(1:n_bl_y, 
                                          function(bl_y) split_bl_mat_y_(bl_y), 
                                          mc.cores = n_cpus)
    }
    
    epispot_bl_ <- function(bl) {
      
      if (n_bl_y > 1) {
        
        bl_x <- ceiling(bl / n_bl_y)
        bl_y <- bl %% n_bl_y
        if (bl_y == 0)
          bl_y <- n_bl_y
        
        # recover split Y matrix
        #
        Y_bl <- list_bl_mat_y[[bl_y]]
        
      } else {
        
        bl_x <- bl
        bl_y <- NULL
        
        Y_bl <- Y
        
      }
      
      # recover split X and V matrices 
      #
      list_bl_x <- list_bl_mat_x[[bl_x]]
      X_bl <- list_bl_x$X_bl
      V_bl <- list_bl_x$V_bl 
      bool_rmvd_v_bl <- list_bl_x$bool_rmvd_v_bl
      rmvd_cst_v_bl <- list_bl_x$rmvd_cst_v 
      rmvd_coll_v_bl <- list_bl_x$rmvd_coll_v
      rmvd_bin_annot_freq_v_bl <- list_bl_x$rmvd_bin_annot_freq_v
      
      
      # split hyperparameters and initial parameters
      #
      # mat_x:
      #
      list_hyper_bl <- list_hyper
      
      pos_x <- list_pos_bl_x[[bl_x]]
      list_hyper_bl$p_hyper <- length(pos_x)
      
      list_init_bl <- list_init
      
      list_init_bl$p_init <- length(pos_x)
      list_init_bl$gam_vb <- list_init_bl$gam_vb[pos_x,, drop = FALSE]
      list_init_bl$mu_beta_vb <- list_init_bl$mu_beta_vb[pos_x,, drop = FALSE]
      
      if (n_bl_y > 1) {
        
        pos_y <- list_pos_bl_y[[bl_y]]
        
        list_hyper_bl$q_hyper <- length(pos_y)
        
        list_hyper_bl$eta <- list_hyper_bl$eta[pos_y]
        list_hyper_bl$kappa <- list_hyper_bl$kappa[pos_y]
        list_hyper_bl$n0 <- list_hyper_bl$n0[pos_y]
        
        list_init_bl$q_init <- length(pos_y)
        list_init_bl$gam_vb <- list_init_bl$gam_vb[, pos_y, drop = FALSE]
        list_init_bl$mu_beta_vb <- list_init_bl$mu_beta_vb[, pos_y, drop = FALSE]
        list_init_bl$sig2_beta_vb <- list_init_bl$sig2_beta_vb[pos_y]
        
        list_init_bl$tau_vb <- list_init_bl$tau_vb[pos_y]
        
      }
      
      list_init_bl$mu_xi_vb <- list_init_bl$mu_xi_vb[!bool_rmvd_v_bl,, drop = FALSE]
      
      # adjust the sparsity level w.r.t. the blocks size
      
      p_bl <- ncol(X_bl) # block size
      q_bl <- ncol(Y_bl)
      
      p0_bl <- p0
      p0_bl[1] <- p0[1] * p_bl / p
      adj_hyper <- get_n0_t02(q_bl, p_bl, p0_bl)
      
      list_hyper_bl$n0 <- adj_hyper$n0
      list_hyper_bl$t02 <- adj_hyper$t02
      
      vb_bl <- epispot_vbem_core_(Y_bl, X_bl, V_bl, list_hyper_bl, 
                                            list_init_bl$gam_vb,
                                            list_init_bl$mu_beta_vb, 
                                            list_init$om, list_init$s02, 
                                            list_init$s2, 
                                            list_init_bl$sig2_beta_vb,
                                            list_init_bl$tau_vb, 
                                            bool_blocks = TRUE, tol, maxit, 
                                            anneal_schedule, anneal_schedule_vbem, verbose = TRUE,
                                            adaptive_tol_em = adaptive_tol_em)
      
      
      vb_bl$rmvd_cst_v <- rmvd_cst_v_bl
      vb_bl$rmvd_coll_v <- rmvd_coll_v_bl
      vb_bl$rmvd_bin_annot_freq_v <- rmvd_bin_annot_freq_v_bl
      vb_bl$V_bl <- V_bl
      
      
      vb_bl
    }
    
    list_vb <- parallel::mclapply(1:(n_bl_x * n_bl_y), 
                                  function(bl) epispot_bl_(bl), 
                                  mc.cores = n_cpus)
    
    if (n_bl_y > 1) {
      
      if(any(sapply(list_vb, function(vb) class(vb) == "try-error"))) {
        stop(paste0("For at least one of the block, no hyperparameter values ", 
                    "matching the expectation and variance of the number of ", 
                    "active predictors per responses supplied in p0. Please ", 
                    "change p0."))
      }
      
      if (n_bl_x > 1)
        stop("EB local scales not implemented for n_bl_x > 1. Exit")
      
      list_init$s02 <- sapply(list_vb, `[[`, "s02")
    
      list_init$s2 <- matrix(unlist(lapply(list_vb, `[[`, "s2")), 
                             nrow = n_bl_x, byrow = TRUE)  
      
      tmp_om <- lapply(list_vb, `[[`, "om")
      
      list_init$om <- parallel::mclapply(1:n_bl_x, function(bl_x) {
        cbind_fill_matrix(tmp_om[((bl_x - 1) * n_bl_y + 1) : (bl_x * n_bl_y)])
      }, mc.cores = n_cpus)
      
      # om list of length n_bl_x 
      # containing matrices of size r' x n_bl_y (sizes of om r' can be different due to cst or coll columns in V_bl removed)
      
    } else {
      
      if(any(sapply(list_vb, function(vb) class(vb) == "try-error"))) {
        
        stop(paste0("For at least one of the block, no hyperparameter values ", 
                    "matching the expectation and variance of the number of ", 
                    "active predictors per responses supplied in p0. Please ",
                    "change p0."))
        
      }
      
      list_init$s02 <- unlist(lapply(list_vb, `[[`, "s02"))
      list_init$s2 <- unlist(lapply(list_vb, `[[`, "s2")) # now it is a vector with the s2 corresponding to each predictor
      list_init$om <- lapply(list_vb, `[[`, "om") # om list of length n_bl_x (sizes of om can be different due to cst or coll columns in V_bl removed)
      
    }
    list_V <- lapply(list_bl_mat_x, `[[`, "V_bl") # V_bl without cst and coll and standardized in each block
    
    
    if (n_bl_y > 1) {
      
      vb <- epispot_blocks_modules_core_(Y, X, list_V, vec_fac_bl_x,
                                                   vec_fac_bl_y, 
                                                   list_modules$module_names,
                                                   list_hyper, 
                                                   list_init$gam_vb, 
                                                   list_init$mu_beta_vb, 
                                                   list_init$om,
                                                   list_init$s02, list_init$s2,
                                                   list_init$sig2_beta_vb, 
                                                   list_init$tau_vb,
                                                   tol, maxit, anneal_schedule, verbose)
      
      
    } else {
      
      vb <- epispot_blocks_core_(Y, X, list_V, vec_fac_bl_x, 
                                           list_hyper, list_init$gam_vb, 
                                           list_init$mu_beta_vb, 
                                           list_init$om, list_init$s02, 
                                           list_init$s2, list_init$sig2_beta_vb, 
                                           list_init$tau_vb, tol, maxit, anneal_schedule, 
                                           verbose)
      
    }
    
  }
  
  vb$p0 <- p0
  
  vb$rmvd_cst_x <- dat$rmvd_cst_x
  vb$rmvd_coll_x <- dat$rmvd_coll_x
  
  vb$rmvd_cst_v <- dat$rmvd_cst_v
  vb$rmvd_coll_v <- dat$rmvd_coll_v
  vb$rmvd_bin_annot_freq_v <- dat$rmvd_bin_annot_freq_v
  
  if (!is.null(list_modules$undo_order_y_ids)) {# for the case where modules are provided and responses are not grouped per module
    # we have grouped them prior to the analysis and ungroup them here after the run in epispot.R
    
    Y <- Y[, list_modules$undo_order_y_ids, drop = FALSE]
    
    vb$beta_vb <- vb$beta_vb[, list_modules$undo_order_y_ids, drop = FALSE] 
    vb$gam_vb <- vb$gam_vb[, list_modules$undo_order_y_ids, drop = FALSE] 
    vb$zeta_vb <- vb$zeta_vb[list_modules$undo_order_y_ids]
    
  }
  
  if (save_hyper) {
    
    if (!is.null(list_modules$undo_order_y_ids)) { 
      list_hyper$eta <- list_hyper$eta[list_modules$undo_order_y_ids]
      list_hyper$kappa <- list_hyper$kappa[list_modules$undo_order_y_ids]
      list_hyper$n0 <- list_hyper$n0[list_modules$undo_order_y_ids]
    }
    
    vb$list_hyper <- list_hyper
  }
  if (save_init) {
    
    if (!is.null(list_modules$undo_order_y_ids)) { 
      
      list_init$gam_vb <- list_init$gam_vb[, list_modules$undo_order_y_ids, drop = FALSE]
      list_init$mu_beta_vb <- list_init$mu_beta_vb[, list_modules$undo_order_y_ids, drop = FALSE]
      list_init$sig2_beta_vb <- list_init$sig2_beta_vb[list_modules$undo_order_y_ids]
      
    }
    
    vb$list_init <- list_init
    
  }
  
  class(vb) <- "epispot"
  
  
  if (verbose) cat("... done. == \n\n")
  
  vb
  
}
