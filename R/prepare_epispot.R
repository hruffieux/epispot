# This file is part of the `epispot` R package:
#     https://github.com/hruffieux/epispot
#

# Internal function implementing sanity checks and needed preprocessing before
# the application of the different `epispot_*_core` algorithms.
#
prepare_data_ <- function(Y, X, V, bin_annot_freq, user_seed, tol, maxit, verbose) {
  
  check_structure_(bin_annot_freq, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(bin_annot_freq)) {
    check_zero_one_(bin_annot_freq) 
  }
  
  check_structure_(user_seed, "vector", "numeric", 1, null_ok = TRUE)

  if (!is.null(user_seed)) {
    
    if (verbose) cat(paste0("Seed set to user_seed ", user_seed, ". \n"))
    
    set.seed(user_seed)
    
  }
  
  check_structure_(tol, "vector", "numeric", 1)
  check_positive_(tol, eps=.Machine$double.eps)
  
  check_structure_(maxit, "vector", "numeric", 1)
  check_natural_(maxit)
  
  check_structure_(verbose, "vector", "logical", 1)
  
  check_structure_(X, "matrix", "numeric")
  
  n <- nrow(X)
  p <- ncol(X)
  
  check_structure_(Y, "matrix", "numeric")
  q <- ncol(Y)
  
  if (nrow(Y) != n) stop("X and Y must have the same number of samples.")
  
  if (is.null(rownames(X)) & is.null(rownames(Y)))
    rownames(X) <- rownames(Y) <- paste("Ind_", 1:n, sep="")
  else if (is.null(rownames(X))) rownames(X) <- rownames(Y)
  else if (is.null(rownames(Y))) rownames(Y) <- rownames(X)
  else if (any(rownames(X) != rownames(Y)))
    stop("The provided rownames of X and Y must be the same.")
  
  if (is.null(colnames(X))) colnames(X) <- paste("Cov_x_", 1:p, sep="")
  if (is.null(colnames(Y))) colnames(Y) <- paste("Resp_", 1:q, sep="")
  
  X <- scale(X)
  
  list_X_cst <- rm_constant_(X, verbose)
  X <- list_X_cst$mat
  bool_cst_x <- list_X_cst$bool_cst
  rmvd_cst_x <- list_X_cst$rmvd_cst
  
  list_X_coll <- rm_collinear_(X, verbose)
  X <- list_X_coll$mat
  
  bool_coll_x <- list_X_coll$bool_coll
  rmvd_coll_x <- list_X_coll$rmvd_coll
  
  bool_rmvd_x <- bool_cst_x
  bool_rmvd_x[!bool_cst_x] <- bool_coll_x
  
  check_structure_(V, "matrix", "numeric")
  
  r <- ncol(V)
  if (nrow(V) != p) stop("The number of rows of V must match the number of candidate predictors in X.")
  
  V <- V[!bool_rmvd_x, , drop = FALSE] # remove the rows corresponding to the removed candidate predictors
  
  if (is.null(rownames(V))) rownames(V) <- colnames(X)
  else if(any(rownames(V) != colnames(X)))
    stop("The provided rownames of V must be the same than those of X and Y or NULL.")
  
  if (is.null(colnames(V))) colnames(V) <- paste("Annot_", 1:r, sep="")
  

  list_V_bin_annot_freq <- rm_bin_annot_freq_(V, bin_annot_freq, verbose)
  V <- list_V_bin_annot_freq$mat
  bool_bin_annot_freq_v <- list_V_bin_annot_freq$bool_bin_annot_freq
  rmvd_bin_annot_freq_v <- list_V_bin_annot_freq$rmvd_bin_annot_freq
    
  V <- scale(V)
  
  list_V_cst <- rm_constant_(V, verbose)
  V <- list_V_cst$mat
  bool_cst_v <- list_V_cst$bool_cst
  rmvd_cst_v <- list_V_cst$rmvd_cst
  
  list_V_coll <- rm_collinear_(V, verbose)
  V <- list_V_coll$mat
  r <- ncol(V)
  bool_coll_v <- list_V_coll$bool_coll
  rmvd_coll_v <- list_V_coll$rmvd_coll
  
  bool_rmvd_v <- bool_bin_annot_freq_v
  bool_rmvd_v[!bool_bin_annot_freq_v] <- bool_cst_v
  bool_rmvd_v[!bool_bin_annot_freq_v][!bool_cst_v] <- bool_coll_v
  
  if (sum(!bool_rmvd_v) == 0)
    stop("All variables provided in V are removed by the QC filters. Exit")
  
  p <- ncol(X)
  if (p < 1) stop(paste("There must be at least 1 non-constant candidate ", 
                        "predictor stored in X."))
  
  Y <- scale(Y, center = TRUE, scale = FALSE)
  
  if (is.null(r) || r < 1) V <- NULL # in principle useless given the above assert.
  
  create_named_list_(Y, X, V,
                     bool_rmvd_x, bool_rmvd_v,
                     rmvd_cst_x, rmvd_cst_v,
                     rmvd_coll_x,  rmvd_coll_v,
                     rmvd_bin_annot_freq_v)
  
}


check_annealing_ <- function(anneal_schedule, maxit) {
  
  check_structure_(anneal_schedule, "vector", "numeric", 3, null_ok = TRUE)
  
  if (!is.null(anneal_schedule)) {
    
    check_natural_(anneal_schedule[c(1, 3)])
    check_positive_(anneal_schedule[2])
    
    stopifnot(anneal_schedule[1] %in% 1:3)
    
    if (anneal_schedule[2] < 1.5)
      stop(paste0("Initial temperature very small. May not be large enough ",
                  "for a successful exploration. Please increase it or select no annealing."))
    
    if (anneal_schedule[3] > 1000)
      stop(paste0("Temperature ladder size very large. This may be unnecessarily ",
                  "computationally demanding. Please decrease it."))
    
    if (maxit <= anneal_schedule[3])
      stop("The maximum number of iterations (maxit) must be strictly larger than the ladder size.")
    
  }
  
}


# Internal function implementing sanity checks and needed preprocessing for the
# model hyperparameters before the application of the different `epispot_*_core`
# algorithms.
#
prepare_list_hyper_ <- function(list_hyper, Y, p, p0, list_blocks, bool_rmvd_x, 
                                names_x, names_y, verbose) {
  
  q <- ncol(Y)
  
  if (is.null(list_hyper)) {
    
    if (verbose) cat("list_hyper set automatically. \n")
    
    list_hyper <- auto_set_hyper_(Y, p, p0)
    
  } else {
    
    if (!inherits(list_hyper, c("hyper", "out_hyper")))
      stop(paste("The provided list_hyper must be an object of class ``hyper'' ",
                 "or ``out_hyper''. \n",
                 "*** you must either use the function set_hyper to ",
                 "set your own hyperparameters or use list_hyper from a ``vb'' ",
                 "object or set the argument list_hyper to NULL for automatic choice. ***",
                 sep=""))
    
    if (inherits(list_hyper, "hyper")) {
      p_hyper_match <- length(bool_rmvd_x)
    } else {
      p_hyper_match <- p
    }
    
    
    if (list_hyper$q_hyper != q)
      stop(paste("The dimensions (q) of the provided hyperparameters ",
                 "(list_hyper) are not consistent with that of Y.\n", sep=""))
    
    if (list_hyper$p_hyper != p_hyper_match)
      stop(paste("The dimensions (p) of the provided hyperparameters ",
                 "(list_hyper) are not consistent with that of X.\n", sep=""))
    
    if (!is.null(list_blocks$order_y_ids)) { # for the case where modules are provided and responses are not grouped per module
      # we group them prior to the analysis and ungroup them after the run in epispot.R
      
      list_hyper$eta <- list_hyper$eta[list_blocks$order_y_ids]
      list_hyper$kappa <- list_hyper$kappa[list_blocks$order_y_ids]
      list_hyper$n0 <- list_hyper$n0[list_blocks$order_y_ids]
      
    }
    
    if (!is.null(names(list_hyper$eta)) && names(list_hyper$eta) != names_y)
      stop("Provided names for the entries of eta do not match the colnames of the continuous variables in Y")
    
    if (!is.null(names(list_hyper$kappa)) && names(list_hyper$kappa) != names_y)
      stop("Provided names for the entries of kappa do not match the colnames of the continuous variables in Y")
    
    if (!is.null(names(list_hyper$n0)) && names(list_hyper$n0) != names_y)
      stop("Provided names for the entries of n0 do not match the colnames of the continuous variables in Y")
    
  }
  
  class(list_hyper) <- "out_hyper"
  
  list_hyper
}


# Internal function implementing sanity checks and needed preprocessing for the
# starting values before the application of the different `epispot_*_core`
# algorithms.
#
prepare_list_init_ <- function(list_init, Y, p, p0, r, list_blocks,
                               bool_rmvd_x, bool_rmvd_v, user_seed, verbose) {
  
  q <- ncol(Y)
  n <- nrow(Y)
  
  if (is.null(list_init)) {
    
    if (verbose) cat(paste("list_init set automatically. \n", sep=""))
    
    list_init <- auto_set_init_(Y, p, p0, r, user_seed)
    
  } else {
    
    if (!inherits(list_init, c("init", "out_init")))
      stop(paste0("The provided list_init must be an object of class ``init'' ",
                  "or `` out_init''. \n *** you must either use the function ", 
                  "set_init to set your own initialization or use list_init ", 
                  "from a ``vb'' object or  set the argument list_init to NULL ", 
                  "for automatic initialization. ***"))
    
    if (inherits(list_init, "init")) {
      p_init_match <- length(bool_rmvd_x)
      r_init_match <- length(bool_rmvd_v)
    } else {
      p_init_match <- p
      r_init_match <- r
    }
    
    if (list_init$q_init != q)
      stop(paste("The dimensions (q) of the provided initial parameters ",
                 "(list_init) are not consistent with that of Y.\n", sep=""))
    
    if (list_init$p_init != p_init_match)
      stop(paste("The dimensions (p) of the provided initial parameters ",
                 "(list_init) are not consistent with that of X.\n", sep=""))
    
    if (list_init$r_init != r_init_match)
      stop(paste("The dimensions (r) of the provided initial parameters ",
                 "(list_init) are not consistent with that of V.\n", sep=""))
    
    
    if (!is.null(list_blocks$order_y_ids)) { # for the case where modules are provided and responses are not grouped per module
      # we group them prior to the analysis and ungroup them after the run in epispot.R
      
      list_init$gam_vb <- list_init$gam_vb[, list_blocks$order_y_ids, drop = FALSE]
      list_init$mu_beta_vb <- list_init$mu_beta_vb[, list_blocks$order_y_ids, drop = FALSE]
      list_init$sig2_beta_vb <- list_init$sig2_beta_vb[list_blocks$order_y_ids]
      
    }
    
    if (inherits(list_init, "init")) {
      
      list_init$gam_vb <- list_init$gam_vb[!bool_rmvd_x,, drop = FALSE]
      
      list_init$mu_beta_vb <- list_init$mu_beta_vb[!bool_rmvd_x,, drop = FALSE]
      
      list_init$om <- list_init$om[!bool_rmvd_v]
      
    }
    
  }
  
  class(list_init) <- "out_init"
  
  list_init
}


# Internal function implementing sanity checks and needed preprocessing to the
# settings provided by the user for block-wise parallel inference.
#
prepare_blocks_ <- function(list_blocks, q, bool_rmvd_x) {
  
  p <- length(bool_rmvd_x)
  
  if (inherits(list_blocks, "modules")) {
    
    if (list_blocks$q_modules != q) {
      stop("The number of responses provided to set_modules does not match that in Y. Exit.")
    }
    
    order_y_ids <- list_blocks$order_y_ids
    undo_order_y_ids <- list_blocks$undo_order_y_ids
    module_names <- list_blocks$module_names
    
    list_blocks <- set_blocks(c(length(bool_rmvd_x), q), 
                              list(1, list_blocks$pos_modules), 
                              n_cpus = list_blocks$n_cpus) 
    
  } else if (inherits(list_blocks, "blocks")) {
    
    order_y_ids <- undo_order_y_ids <- NULL
    
    if (list_blocks$bl_y$n_bl > 1) {
      module_names <- paste0("Module_", 1:list_blocks$bl_y$n_bl)
    } else {
      module_names <- NULL
    }
    
    
  } else {
    
    stop(paste0("The provided list_blocks must be an object of class ``blocks'' or ``modules''. \n",
                "*** you must either use the function set_blocks to give the settings ",
                "for parallels applications of epispot on blocks of candidate ",
                "predictors or set list_blocks to NULL to apply epispot jointly on ",
                "all the candidate predictors (sufficient RAM required). ***"))
    
  }
  
  if (!is.null(list_blocks$bl_y)) {
    
    if (list_blocks$bl_x$n_var_blocks != length(bool_rmvd_x))
      stop(paste0("The number of candidate predictors p provided to the function ", 
                  "set_modules is not consistent with X.\n"))
    
    vec_fac_bl_x <- list_blocks$bl_x$vec_fac_bl[!bool_rmvd_x]
    
    if (list_blocks$bl_y$n_var_blocks != q)
      stop(paste("The number of responses q provided to the function set_blocks ",
                 "is not consistent with Y.\n", sep=""))
    
    vec_fac_bl_y <- list_blocks$bl_y$vec_fac_bl
    
    tab_bl_y <- table(vec_fac_bl_y)
    pres_bl_y <- tab_bl_y > 0
    
    # in case a block was removed due to the above because of bool_rmvd_y
    n_bl_y  <- sum(pres_bl_y)
    
  } else {
    
    if (list_blocks$n_var_blocks != length(bool_rmvd_x))
      stop(paste("The number of candidate predictors p provided to the function set_blocks ",
                 "is not consistent with X.\n", sep="")) 
    
    vec_fac_bl_x <- list_blocks$vec_fac_bl[!bool_rmvd_x]
    vec_fac_bl_y <- NULL
    n_bl_y <- 1
    
  }
  
  tab_bl_x <- table(vec_fac_bl_x)
  pres_bl_x <- tab_bl_x > 0
  
  # in case a block was removed due to the above because of bool_rmvd_x
  n_bl_x  <- sum(pres_bl_x)
  if (is.null(list_blocks$bl_y) && list_blocks$n_cpus > n_bl_x) {
    n_cpus <- n_bl_x
  } else if (!is.null(list_blocks$bl_y) && list_blocks$n_cpus > (n_bl_x * n_bl_y)){
    n_cpus <- n_bl_x * n_bl_y
  } else {
    n_cpus <- list_blocks$n_cpus
  }
  
  create_named_list_(n_bl_x, n_bl_y, n_cpus, vec_fac_bl_x, vec_fac_bl_y,
                     order_y_ids, undo_order_y_ids, module_names)
  
}

# Internal function gathering settings for blockwise inference. The function 
# \code{set_modules} is a wrapper of this function.  \code{set_blocks} also 
# allows inference on partitioned predictor spaces, although this functionality
# is not accessible to the package user.
#
set_blocks <- function(tot, pos_bl, n_cpus, verbose = TRUE) {
  
  check_structure_(n_cpus, "vector", "numeric", 1)
  check_natural_(n_cpus)
  
  check_structure_(verbose, "vector", "logical", 1)
  
  list_blocks <- lapply(seq_along(tot), function(ii) {
    
    tt <- tot[ii]
    
    if (is.list(pos_bl)) {
      pb <- pos_bl[[ii]]  
    } else {
      pb <- pos_bl
    }
    
    check_structure_(tt, "vector", "numeric", 1)
    check_natural_(tt)
    
    check_structure_(pb, "vector", "numeric")
    check_natural_(pb)
    
    if (length(pb) > 25)
      warning(paste("The provided number of blocks may be too large for accurate ",
                    "inference. If possible, use less blocks.", sep = ""))
    
    if (any(pb < 1) | any(pb > tt))
      stop("The positions provided in pos_bl must range between 1 and total number of variables given in tot.")
    
    if (any(duplicated(pb)))
      stop("The positions provided in pos_bl must be unique.")
    
    if (any(pb != cummax(pb)))
      stop("The positions provided in pos_bl must be monotonically increasing.")
    
    vec_fac_bl <- as.factor(cumsum(seq_along(1:tt) %in% pb))
    
    n_bl <- length(unique(vec_fac_bl))
    
    n_var_blocks <- tt
    
    create_named_list_(n_var_blocks, n_bl, vec_fac_bl)
    
  })
  
  if (length(list_blocks) > 1 && list_blocks[[2]]$n_bl != 1) {
    names(list_blocks) <- c("bl_x", "bl_y")
    tot_n_bl <- list_blocks$bl_x$n_bl * list_blocks$bl_y$n_bl
  } else {
    list_blocks <- list_blocks[[1]]
    tot_n_bl <- list_blocks$n_bl
  }
  
  if (n_cpus > 1) {
    
    n_cpus_avail <- parallel::detectCores()
    if (n_cpus > n_cpus_avail) {
      n_cpus <- n_cpus_avail
      warning(paste("The number of CPUs specified exceeds the number of CPUs ",
                    "available on the machine. The latter has been used instead.",
                    sep=""))
    }
    
    if (n_cpus > tot_n_bl){
      message <- paste("The number of cpus in use is at most equal to the number of blocks.",
                       "n_cpus is therefore set to ", tot_n_bl, ". \n", sep ="")
      if(verbose) cat(message)
      else warning(message)
      n_cpus <- tot_n_bl
    }
    
    if (verbose) cat(paste("epispot applied in parallel on ", tot_n_bl,
                           " blocks or modules, using ", n_cpus, " CPUs.\n",
                           "Please make sure that enough RAM is available. \n", sep=""))
  }
  
  list_blocks$n_cpus <- n_cpus
  
  class(list_blocks) <- "blocks"
  
  list_blocks
}

#' Wrapper of set_blocks to specifically define module partitions.
#'
#' Parallel applications of the method on blocks of candidate predictors for
#' large datasets allows faster and less RAM-greedy executions.
#'
#' @param module_ids Numeric vector of size the number of response variables 
#'   which assigns each response to its corresponding module. See example below.
#'   Each module must have size greater than 10 to avoid unstable inference. 
#' @param module_map Vector of size the total number of modules where each value 
#'   the provides the ids of the module as in module_ids and the names of the 
#'   vector correspond to the desired module names. See example below. Default 
#'   is \code{NULL} for default module names.
#' @param n_cpus Number of CPUs to be used. If large, one should ensure that
#'   enough RAM will be available for parallel execution. Set to 1 for serial
#'   execution.
#'
#' @return An object of class "\code{modules}" preparing the settings for a
#'   parallel module-based EPISPOT run a form that can be passed to the 
#'   \code{\link{epispot}} function.
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
#' ########################
#' ## Infer associations ##
#' ########################
#'
#' module_ids <- sample(c(rep(1, floor(q/2)), rep(2, q - floor(q/2)))) # 2 modules
#' 
#' module_names <- sort(unique(module_ids))
#' names(module_names) <- paste0("m_", 1:length(unique(module_ids)))
#' 
#' list_modules <- set_modules(module_ids, module_names, n_cpus = 1)
#' 
#' # Expectation and variance for the prior number of predictors associated with
#' # each response
#' #
#' p0 <- c(mean(colSums(pat)), 10)
#' 
#' vb_m <- epispot(Y = Y, X = X, V = V, p0 = p0, list_modules = list_modules,
#'                 user_seed = seed)
#'
#'@seealso \code{\link{epispot}}.
#'
#' @export
set_modules <- function(module_ids, module_map = NULL, n_cpus = 1) {
  
  check_structure_(module_ids, "vector", "numeric")
  
  tb_module_ids <- table(module_ids)
  if(any(tb_module_ids < 10)) {
    stop("Modules with size < 10 are not accepted as may induce unstable inference. Exit.")
  }
  
  n_modules <- length(tb_module_ids)
  q_modules <- length(module_ids)
  
  if (is.null(module_map)) {
    
    module_names <- paste0("Module_", 1:n_modules)
    
  } else {
    
    if (length(module_map) != n_modules) {
      stop("The number of module names provided does not match the number of modules. Exit.")
    }
    
    sorted_module_map <- sort(module_map)
    
    if (!all(sorted_module_map == sort(unique(module_ids)))) {
      stop("The module map is inconsistent with the modules provided. Exit.")
    }
    
    module_names <- names(sorted_module_map)
    
  }
  
  order_y_ids <- order(module_ids)
  if (all(order_y_ids == 1:q_modules)) { # already in the correct order.
    order_y_ids <- undo_order_y_ids <- NULL
  } else {
    undo_order_y_ids <- order(order_y_ids) 
  }
  
  pos_modules <- c(1, cumsum(tb_module_ids[-n_modules]) + 1)
  names(pos_modules) <- module_names
  
  list_modules <- create_named_list_(order_y_ids, undo_order_y_ids, q_modules, pos_modules, 
                                     module_names, n_modules, n_cpus)
  class(list_modules) <- "modules"
  
  list_modules
  
}


