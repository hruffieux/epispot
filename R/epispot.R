# This file is part of the `epispot` R package:
#     https://github.com/hruffieux/epispot
#

#' Fit sparse multivariate regression models using variational inference.
#'
#' Variational approximation procedure fitting sparse multivariate regression
#' models for combined selection of predictors and associated responses in
#' high-dimensional set-ups. Dependence across responses linked to the same
#' predictors is captured through the model hierarchical structure.
#'
#'
#' The optimization is made using efficient block coordinate ascent schemes, for
#' which convergence is ensured as the objective (elbo) is multiconcave
#' for the selected blocks, i.e., it is concave in each block of parameters
#' whose updates are made simultaneously, see Wu et al. (reference Section
#' below).
#'
#' The continuous response variables in \code{Y} (if any) will be centered
#' before application of the variational algorithm, and the candidate predictors
#' in \code{X} will be standardized.
#'
#'
#' @param Y Response data matrix of dimension n x d, where n is the number of
#'   samples and d is the number of response variables.
#' @param X Input matrix of dimension n x p, where p is the number of candidate
#'   predictors. \code{X} cannot contain NAs. No intercept must be supplied.
#' @param p0 Vector of size 2 whose entries are the prior expectation and 
#'   variance of the number of predictors associated with each response.
#'   Must be \code{NULL} if \code{list_init} and \code{list_hyper} are both 
#'   non-\code{NULL}.
#' @param V Annotation matrix of dimension p x r, where r is the number of
#'   variables representing external information on the candidate predictors
#'   which may make their selection more or less likely. \code{NULL} if no such
#'   information.
#' @param list_blocks An object of class "\code{blocks}" containing settings for
#'   parallel inference on a partitioned predictor space. Must be filled using
#'   the \code{\link{set_blocks}} function or must be \code{NULL} for no
#'   partitioning.
#' @param list_hyper An object of class "\code{hyper}" containing the model
#'   hyperparameters. Must be filled using the \code{\link{set_hyper}}
#'   function or must be \code{NULL} for default hyperparameters.
#' @param list_init An object of class "\code{init}" containing the initial
#'   variational parameters. Must be filled using the \code{\link{set_init}}
#'   function or be \code{NULL} for a default initialization.
#' @param user_seed Seed set for reproducible default choices of hyperparameters
#'   (if \code{list_hyper} is \code{NULL}) and initial variational parameters
#'   (if \code{list_init} is \code{NULL}). Default is \code{NULL}, no
#'   seed set.
#' @param tol Tolerance for the stopping criterion.
#' @param adaptive_tol_em Boolean indicating whether the tolerance for the 
#'   within-EM variational runs should be controlled adaptively depending on the
#'   EM-convergence status.
#' @param maxit Maximum number of iterations allowed.
#' @param anneal Parameters for annealing scheme. Must be a vector whose first
#'   entry is sets the type of ladder: 1 = geometric spacing, 2 = harmonic
#'   spacing or 3 = linear spacing, the second entry is the initial temperature,
#'   and the third entry is the ladder size. If \code{NULL} (default), no
#'   annealing is performed.
#' @param anneal_vb_em Parameters for annealing scheme for the internal runs of the 
#'   variational EM algorithm. 
#' @param save_hyper If \code{TRUE}, the hyperparameters used for the model are
#'   saved as output.
#' @param save_init If \code{TRUE}, the initial variational parameters used for
#'   the inference are saved as output.
#' @param verbose If \code{TRUE}, messages are displayed during execution.
#'
#' @return An object of class "\code{epispot}" containing the following output:
#'  \item{m1_beta}{Estimated effect size matrix of dimension p x q. Entry (s, t) 
#'                 corresponds to the variational posterior mean 
#'                 (mu_beta_vb_st x gam_vb_st) of the regression effect between 
#'                 candidate predictor s and response t.}
#'  \item{gam_vb}{Posterior inclusion probability matrix of dimension p x d.
#'                Entry (s, t) corresponds to the posterior probability of
#'                association between candidate predictor s and response t.}
#'  \item{mu_c_vb}{Vector of size r, where entry l contains the overall effect 
#'                 of annotation l on the
#'                 probabilities of associations.\code{NULL} if \code{V} is
#'                 \code{NULL}.}
#'  \item{mu_rho_vb}{Vector of length d containing the posterior mean of rho.
#'                   Entry t controls the proportion of predictors associated
#'                   with response t.}
#'  \item{mu_theta_vb}{Vector of length p containing the posterior mean of
#'                     theta. Entry s corresponds to the propensity of candidate
#'                     predictor s to be included in the model.}
#'  \item{om}{Vector of length p containing the posterior mean of omega.
#'               Entry s controls the proportion of responses associated with
#'               candidate predictor s. NULL if \code{V} is non-\code{NULL}.}
#'  \item{zeta_vb}{Posterior inclusion probability vector of size r for the
#'                 annotation variables. \code{NULL} if \code{V} is \code{NULL}.}
#'  \item{converged}{A boolean indicating whether the algorithm has converged
#'                   before reaching \code{maxit} iterations.}
#'  \item{it}{Final number of iterations.}
#'  \item{lb_opt}{Optimized variational lower bound for the marginal
#'                log-likelihood (ELBO).}
#'  \item{diff_lb}{Difference in ELBO between the last and penultimate
#'                 iterations. This may be a useful diagnostic information when
#'                 convergence has not been reached before \code{maxit}.}
#'  \item{p0}{Vector of length 2 defining the applied sparsity control.}
#'  \item{rmvd_cst_x}{Vectors containing the indices of constant
#'                                variables in \code{X} removed
#'                                prior to the analysis.}
#'  \item{rmvd_coll_x}{Vectors containing the indices of variables
#'                                  in \code{X} removed prior
#'                                  to the analysis because collinear to other
#'                                  variables. The entry name indicates the
#'                                  corresponding variable kept in the analysis
#'                                  (i.e., that causing the collinearity for the
#'                                  entry in question).}
#'  \item{list_hyper, list_init}{If \code{save_hyper}, resp. \code{save_init},
#'                               \code{TRUE}, hyperparameters, resp. initial
#'                               variational parameters, used for inference are
#'                               saved as output.}
#'  \item{...}{Other specific outputs are possible depending on the model used.}
#'  
#' @examples
#' seed <- 123; set.seed(seed)
#'
#' ###################
#' ## Simulate data ##
#' ###################
#'
#' ## Examples using small problem sizes:
#' ##
#' n <- 50; p <- 60; p_act <- 10; d <- 25; d_act <- 15; r <- 10
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
#' shuff_y_ind <- sample(d)
#' 
#' X <- cbind(X_act, X_inact)[, shuff_x_ind]
#'
#' # Association pattern and effect sizes
#' #
#' pat <- matrix(FALSE, ncol = d, nrow = p)
#' bool_x <- shuff_x_ind <= p_act
#' bool_y <- shuff_y_ind <= d_act
#' 
#' pat_act <- beta_act <- matrix(0, nrow = p_act, ncol = d_act)
#' pat_act[sample(p_act * d_act, floor(p_act * d_act / 5))] <- 1
#' beta_act[as.logical(pat_act)] <-  rnorm(sum(pat_act))
#' 
#' pat[bool_x, bool_y] <- pat_act
#' 
#' # Gaussian responses
#' #
#' Y_act <- matrix(rnorm(n * d_act, mean = X_act %*% beta_act), nrow = n)
#' Y_inact <- matrix(rnorm(n * (q - d_act)), nrow = n)
#'
#' Y <- cbind(Y_act, Y_inact)[, shuff_y_ind]
#'
#' # Annotation variables
#' #
#' V <- matrix(rnorm(p * r), nrow = p)
#' V[bool_x_act, ] <- rnorm(p_act * r, mean = 2)
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
#' # Inference without modules
#' #
#' res_epispot <- epispot(Y = Y, X = X, p0 = p0, V = V, user_seed = seed)
#'
#' # Inference with modules
#' #
#' module_ids <- sample(c(rep(1, floor(d/2)), rep(2, d - floor(d/2)))) # 2 modules
#' list_modules <- set_modules(module_ids, n_cpus = 1)
#' 
#' vb_m <- epispot(Y = Y, X = X, p0 = p0, V = V, list_blocks = list_modules,
#'                 user_seed = seed)
#'
#' @references
#' H. Ruffieux, B. P. Fairfax, E. Vigorito, C. Walllace, S. Richardson, 
#' L. Bottolo. EPISPOT: an epigenome-driven approach for detecting and 
#' interpreting hotspots in molecular QTL studies, arXiv, 2020.
#'
#' @seealso \code{\link{set_hyper}}, \code{\link{set_init}}, 
#' \code{\link{set_modules}}.
#'
#' @export
#'
epispot <- function(Y, X, p0, V, list_blocks = NULL, list_hyper = NULL, 
                    list_init = NULL, user_seed = NULL, tol = 1e-3, 
                    adaptive_tol_em = FALSE, maxit = 1000, anneal = NULL, 
                    anneal_vb_em = NULL, save_hyper = FALSE, save_init = FALSE, 
                    verbose = TRUE) {
  
  
  if (verbose) cat("== Preparing the data ... \n")
  
  dat <- prepare_data_(Y, X, V, user_seed, tol, maxit, verbose)
  
  X <- dat$X
  Y <- dat$Y
  V <- dat$V
  
  bool_rmvd_x <- dat$bool_rmvd_x
  bool_rmvd_v <- dat$bool_rmvd_v
  
  n <- nrow(X)
  p <- ncol(X)
  d <- ncol(Y)
  r <- ncol(V)
  
  names_x <- colnames(X)
  names_y <- colnames(Y)
  names_v <- colnames(V)
  
  check_annealing_(anneal, maxit)
  check_annealing_(anneal_vb_em, maxit)
  
  if (verbose) cat("... done. == \n\n")
  
  if (!is.null(list_blocks)) {
    
    list_blocks <- prepare_blocks_(list_blocks, d, bool_rmvd_x)
    
    if (!is.null(list_blocks$order_y_ids)) { # for the case where modules are provided and responses are not grouped per module
      # we group them prior to the analysis and ungroup them after the run in epispot.R
      Y <- Y[, list_blocks$order_y_ids, drop = FALSE]
      
    }
    
    n_bl_x <- list_blocks$n_bl_x
    n_bl_y <- list_blocks$n_bl_y
    
    n_cpus <- list_blocks$n_cpus
    
    vec_fac_bl_x <- list_blocks$vec_fac_bl_x
    vec_fac_bl_y <- list_blocks$vec_fac_bl_y
    
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
  
  list_hyper <- prepare_list_hyper_(list_hyper, Y, p, p0, list_blocks,
                                    bool_rmvd_x, names_x, names_y, verbose)
  
  if (verbose) cat("... done. == \n\n")
  
  if (verbose) cat("== Preparing the parameter initialization ... \n\n")
  
  list_init <- prepare_list_init_(list_init, Y, p, p0, r, list_blocks,
                                  bool_rmvd_x, bool_rmvd_v, 
                                  user_seed, verbose)
  
  if (verbose) cat("... done. == \n\n")
  
  if (verbose){
    cat(paste0("===================================================== \n",
               "== EPISPOT: annotation-driven association mapping  == \n",
               "===================================================== \n\n"))
  }
  
  
  if (is.null(list_blocks)) {
    
    vb <- epispot_dual_info_vbem_core_(Y, X, V, list_hyper, list_init$gam_vb,
                                       list_init$mu_beta_vb, list_init$om,
                                       list_init$s02, list_init$s2,
                                       list_init$sig2_beta_vb, list_init$tau_vb,
                                       bool_blocks = FALSE, tol, maxit, anneal, 
                                       anneal_vb_em, verbose,
                                       adaptive_tol_em = adaptive_tol_em)
    
  } else {
    
    list_pos_bl_x <- split(1:p, vec_fac_bl_x)
    
    split_bl_mat_x_ <- function(bl_x) {
      
      X_bl <- X[, list_pos_bl_x[[bl_x]], drop = FALSE]
      
      V_bl <- scale(V[list_pos_bl_x[[bl_x]],, drop = FALSE]) # the VB algorithm assumes a scaled V
      
      list_V_bl_cst <- rm_constant_(V_bl, verbose = FALSE)
      V_bl <- list_V_bl_cst$mat
      bool_cst_v_bl <- list_V_bl_cst$bool_cst
      rmvd_cst_v_bl <- list_V_bl_cst$rmvd_cst
      
      list_V_bl_coll <- rm_collinear_(V_bl, verbose = FALSE)
      V_bl <- list_V_bl_coll$mat
      r <- ncol(V_bl)
      bool_coll_v_bl <- list_V_bl_coll$bool_coll
      rmvd_coll_v_bl <- list_V_bl_coll$rmvd_coll
      
      bool_rmvd_v_bl <- bool_cst_v_bl
      bool_rmvd_v_bl[!bool_cst_v_bl] <- bool_coll_v_bl
      
      if (sum(!bool_rmvd_v_bl) == 0)
        stop(paste0("There exist one or more blocks for which no non-constant ",
                    "annotation variables remain. Try to use less blocks."))
      
      
      create_named_list_(X_bl, V_bl, bool_rmvd_v_bl, rmvd_cst_v_bl, rmvd_coll_v_bl)
    }
    
    list_bl_mat_x <- parallel::mclapply(1:n_bl_x, function(bl_x) split_bl_mat_x_(bl_x), mc.cores = n_cpus)
    
    
    if (n_bl_y > 1) {
      
      list_pos_bl_y <- split(1:d, vec_fac_bl_y)
      
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
        
        list_hyper_bl$d_hyper <- length(pos_y)
        
        list_hyper_bl$eta <- list_hyper_bl$eta[pos_y]
        list_hyper_bl$kappa <- list_hyper_bl$kappa[pos_y]
        list_hyper_bl$n0 <- list_hyper_bl$n0[pos_y]
        
        list_init_bl$d_init <- length(pos_y)
        list_init_bl$gam_vb <- list_init_bl$gam_vb[, pos_y, drop = FALSE]
        list_init_bl$mu_beta_vb <- list_init_bl$mu_beta_vb[, pos_y, drop = FALSE]
        list_init_bl$sig2_beta_vb <- list_init_bl$sig2_beta_vb[pos_y]
        
        list_init_bl$tau_vb <- list_init_bl$tau_vb[pos_y]
        
      }
      
      list_init_bl$mu_c_vb <- list_init_bl$mu_c_vb[!bool_rmvd_v_bl,, drop = FALSE]
      
      # adjust the sparsity level w.r.t. the blocks size
      
      p_bl <- ncol(X_bl) # block size
      d_bl <- ncol(Y_bl)
      
      p0_bl <- p0
      p0_bl[1] <- p0[1] * p_bl / p
      adj_hyper <- get_n0_t02(d_bl, p_bl, p0_bl)
      
      list_hyper_bl$n0 <- adj_hyper$n0
      list_hyper_bl$t02 <- adj_hyper$t02
      
      vb_bl <- epispot_dual_info_vbem_core_(Y_bl, X_bl, V_bl, list_hyper_bl, 
                                            list_init_bl$gam_vb,
                                            list_init_bl$mu_beta_vb, 
                                            list_init$om, list_init$s02, 
                                            list_init$s2, 
                                            list_init_bl$sig2_beta_vb,
                                            list_init_bl$tau_vb, 
                                            bool_blocks = TRUE, tol, maxit, 
                                            anneal, anneal_vb_em, verbose = TRUE,
                                            adaptive_tol_em = adaptive_tol_em)
      
      
      vb_bl$rmvd_cst_v <- rmvd_cst_v_bl
      vb_bl$rmvd_coll_v <- rmvd_coll_v_bl
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
      
      rownames(list_init$s02) <- colnames(X)
      colnames(list_init$s02) <- list_blocks$module_names
      
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
      
      vb <- epispot_dual_info_blocks_modules_core_(Y, X, list_V, vec_fac_bl_x,
                                                   vec_fac_bl_y, list_hyper, 
                                                   list_init$gam_vb, 
                                                   list_init$mu_beta_vb, 
                                                   list_init$om,
                                                   list_init$s02, list_init$s2,
                                                   list_init$sig2_beta_vb, 
                                                   list_init$tau_vb,
                                                   tol, maxit, anneal, verbose)
      
      
    } else {
      
      vb <- epispot_dual_info_blocks_core_(Y, X, list_V, vec_fac_bl_x, 
                                           list_hyper, list_init$gam_vb, 
                                           list_init$mu_beta_vb, 
                                           list_init$om, list_init$s02, 
                                           list_init$s2, list_init$sig2_beta_vb, 
                                           list_init$tau_vb, tol, maxit, anneal, 
                                           verbose)
      
    }
    
    vb$s02 <- list_init$s02
    
  }
  
  vb$p0 <- p0
  
  vb$rmvd_cst_x <- dat$rmvd_cst_x
  vb$rmvd_coll_x <- dat$rmvd_coll_x
  
  vb$rmvd_cst_v <- dat$rmvd_cst_v
  vb$rmvd_coll_v <- dat$rmvd_coll_v
  
  if (!is.null(list_blocks$undo_order_y_ids)) {# for the case where modules are provided and responses are not grouped per module
    # we have grouped them prior to the analysis and ungroup them here after the run in epispot.R
    
    Y <- Y[, list_blocks$undo_order_y_ids, drop = FALSE]
    
    vb$m1_beta <- vb$m1_beta[, list_blocks$undo_order_y_ids, drop = FALSE] 
    vb$gam_vb <- vb$gam_vb[, list_blocks$undo_order_y_ids, drop = FALSE] 
    vb$mu_rho_vb <- vb$mu_rho_vb[list_blocks$undo_order_y_ids]
    
  }
  
  if (save_hyper) {
    
    if (!is.null(list_blocks$undo_order_y_ids)) { 
      list_hyper$eta <- list_hyper$eta[list_blocks$undo_order_y_ids]
      list_hyper$kappa <- list_hyper$kappa[list_blocks$undo_order_y_ids]
      list_hyper$n0 <- list_hyper$n0[list_blocks$undo_order_y_ids]
    }
    
    vb$list_hyper <- list_hyper
  }
  if (save_init) {
    
    if (!is.null(list_blocks$undo_order_y_ids)) { 
      
      list_init$gam_vb <- list_init$gam_vb[, list_blocks$undo_order_y_ids, drop = FALSE]
      list_init$mu_beta_vb <- list_init$mu_beta_vb[, list_blocks$undo_order_y_ids, drop = FALSE]
      list_init$sig2_beta_vb <- list_init$sig2_beta_vb[list_blocks$undo_order_y_ids]
      
    }
    
    vb$list_init <- list_init
    
  }
  
  class(vb) <- "epispot"
  
  
  if (verbose) cat("... done. == \n\n")
  
  vb
  
}
