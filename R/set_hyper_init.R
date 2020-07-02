# This file is part of the `epispot` R package:
#     https://github.com/hruffieux/epispot
#

#' Gather model hyperparameters provided by the user.
#'
#' This function must be used to provide hyperparameter values for the model
#' used in \code{\link{epispot}}.
#'
#' The \code{\link{epispot}} function can also be used with default 
#' hyperparameter choices (without using \code{\link{set_hyper}}) by setting the 
#' argument \code{list_hyper} to \code{NULL}.
#' 
#' @param d Number of responses.
#' @param p Number of candidate predictors.
#' @param lambda Vector of length 1 providing the values of hyperparameter
#'   \eqn{\lambda} for the prior distribution of \eqn{\sigma^{-2}}. \eqn{\sigma}
#'   represents the typical size of nonzero effects.
#' @param nu Vector of length 1 providing the values of hyperparameter \eqn{\nu}
#'   for the prior distribution of \eqn{\sigma^{-2}}. \eqn{\sigma} represents
#'   the typical size of nonzero effects.
#' @param eta Vector of length 1 or d. Provides the values of
#'   hyperparameter \eqn{\eta} for the prior distributions of the continuous
#'   response residual precisions, \eqn{\tau}. If of length 1, the provided
#'   value is repeated d times.
#' @param kappa Vector of length 1 or d. Provides the values of hyperparameter
#'   \eqn{\kappa} for the prior distributions of the response residual
#'   precisions, \eqn{\tau}. If of length 1, the provided value is repeated d 
#'   times.
#' @param n0 Vector of length 1 or d providing the prior mean for the 
#'   response-specific effect \eqn{\rho}.
#' @param t02 Vector of length 1 providing the prior variance for the 
#'   response-specific effect \eqn{\rho}.
#'
#' @return An object of class "\code{hyper}" preparing user hyperparameter in a
#'   form that can be passed to the \code{\link{epispot}} function.
#'
#' @examples
#' 
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
#' Y_inact <- matrix(rnorm(n * (d - d_act)), nrow = n)
#'
#' Y <- cbind(Y_act, Y_inact)[, shuff_y_ind]
#'
#' # Annotation variables
#' #
#' V <- matrix(rnorm(p * r), nrow = p)
#' V[bool_x, ] <- rnorm(p_act * r, mean = 2)
#'
#' #############################
#' ## Specify hyperparameters ##
#' #############################
#' 
#' list_hyper <- set_hyper(d, p, lambda = 1, nu = 1, eta = 1, kappa = 1, 
#'                         n0 = -2, t02 = 1)
#'                         
#' ########################
#' ## Infer associations ##
#' ########################
#'
#' res_epispot <- epispot(Y = Y, X = X, V = V, p0 = p0, list_hyper = list_hyper, 
#'                        user_seed = seed)
#'
#' @seealso  \code{\link{set_init}}, \code{\link{epispot}}
#'
#' @export
#'
set_hyper <- function(d, p, lambda, nu, eta, kappa, n0, t02) {
  
  check_structure_(d, "vector", "numeric", 1)
  check_natural_(d)
  
  check_structure_(p, "vector", "numeric", 1)
  check_natural_(p)
  
  check_structure_(n0, "vector", "double", c(1, d))
  if (length(n0) == 1) n0 <- rep(n0, d)
  
  check_structure_(t02, "vector", "double", 1)
  check_positive_(t02)
  
  check_structure_(lambda, "vector", "double", 1)
  check_positive_(lambda)
  
  check_structure_(nu, "vector", "double", 1)
  check_positive_(nu)
  
  check_structure_(eta, "vector", "double", c(1, d))
  check_positive_(eta)
  if (length(eta) == 1) eta <- rep(eta, d)
  
  check_structure_(kappa, "vector", "double", c(1, d))
  check_positive_(kappa)
  if (length(kappa) == 1) kappa <- rep(kappa, d)
  
  d_hyper <- d
  p_hyper <- p
  
  list_hyper <- create_named_list_(d_hyper, p_hyper,
                                   eta, kappa, lambda, nu, n0, t02)
  
  class(list_hyper) <- "hyper"
  
  list_hyper
  
}


# Internal function setting default model hyperparameters when not provided by
# the user.
#
auto_set_hyper_ <- function(Y, p, p0) {
  
  d <- ncol(Y)
  
  lambda <- 1e-2
  nu <- 1
  
  check_positive_(lambda)
  check_positive_(nu)
  
  # hyperparameter set using the data Y
  eta <- 1 / median(apply(Y, 2, var)) #median to be consistent when doing permutations
  if (!is.finite(eta)) eta <- 1e3
  check_positive_(eta)
  
  eta <- rep(eta, d)
  kappa <- rep(1, d)
  check_positive_(kappa)
  
  E_p_t <- p0[1]
  V_p_t <- p0[2]
  
  dn <- 1e-6
  up <- 1e5
  
  # Get n0 and t02 similarly as for a_omega_t and b_omega_t in HESS
  # (specify expectation and variance of number of active predictors per response)
  #
  # Look at : gam_st | theta_s = 0
  #
  tryCatch(t02 <- uniroot(function(x)
    get_V_p_t(get_mu(E_p_t, x, p), x, p) - V_p_t,
    interval = c(dn, up))$root,
    error = function(e) {
      stop(paste0("No hyperparameter values matching the expectation and variance ",
                  "of the number of active predictors per responses supplied in p0.",
                  "Please change p0."))
    })
  check_positive_(t02)
  
  # n0 sets the level of sparsity.
  # n0 <- - get_mu(E_p_t, t02, p)
  n0 <- get_mu(E_p_t, t02, p)
  n0 <- rep(n0, d)
  
  d_hyper <- d
  p_hyper <- p
  
  list_hyper <- create_named_list_(d_hyper, p_hyper, 
                                   eta, kappa, lambda, nu, n0, t02)
  
  class(list_hyper) <- "out_hyper"
  
  list_hyper
  
}

#' Gather initial variational parameters provided by the user.
#'
#' This function must be used to provide initial values for the variational
#' parameters used in \code{\link{epispot}}.
#'
#' The \code{\link{epispot}} function can also be used with default initial
#' parameter choices (without using \code{\link{set_init}}) by setting
#' its argument \code{list_init} to \code{NULL}.
#' 
#' @param d Number of responses.
#' @param p Number of candidate predictors.
#' @param r Number of candidate annotations.
#' @param gam_vb Matrix of size p x d with initial values for the variational
#'   parameter yielding posterior probabilities of inclusion.
#' @param mu_beta_vb Matrix of size p x d with initial values for the
#'   variational parameter yielding regression coefficient estimates for
#'   predictor-response pairs included in the model.
#' @param om Vector of length r with initial values for the EM hyperparameter 
#'   yielding the probability parameters for the annotation variables.
#' @param s02 Vector of length 1 with initial value for the EM hyperparameter 
#'   yielding the variance of the predictor-specific modulations of the 
#'   probability parameter for the predictor-response associations (these inital 
#'   values are the same for all predictors, for simplicity).
#' @param s2 Vector of length 1 with initial value for the EM hyperparameter 
#'   yielding the variance of the annotation effects on the probability 
#'   parameter for the predictor-response associations.
#' @param sig2_beta_vb Vector of length d with initial values for
#'   the variational parameter yielding estimates of effect variances for
#'   predictor-response pairs included in the model. These values are the same
#'   for all the predictors (as a result of the predictor variables being
#'   standardized before the variational algorithm).
#' @param tau_vb Vector of length d with initial values for the variational parameter
#'   yielding estimates for the continuous response residual precisions. 
#'
#' @return An object of class "\code{init}" preparing user initial values for
#'   the variational parameters in a form that can be passed to the
#'   \code{\link{epispot}} function.
#'
#' @examples
#' 
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
#' Y_inact <- matrix(rnorm(n * (d - d_act)), nrow = n)
#'
#' Y <- cbind(Y_act, Y_inact)[, shuff_y_ind]
#'
#' # Annotation variables
#' #
#' V <- matrix(rnorm(p * r), nrow = p)
#' V[bool_x, ] <- rnorm(p_act * r, mean = 2)
#' 
#' ################################
#' ## Specify initial parameters ##
#' ################################
#' 
#' tau_vb <- rep(1, d)
#' 
#' gam_vb <- matrix(rbeta(p * d, shape1 = 1, shape2 = 4 * d - 1), nrow = p)
#' 
#' mu_beta_vb <- matrix(rnorm(p * d), nrow = p)
#' sig2_beta_vb <- 1 / rgamma(d, shape = 2, rate = 1)
#' 
#' om <- rep(0.5, r)
#' s02 <- 1 / d
#' s2 <- 0.01
#' 
#' list_init <- set_init(d, p, r, gam_vb, mu_beta_vb, om, s02, s2, sig2_beta_vb, 
#'                       tau_vb)
#' 
#' ########################
#' ## Infer associations ##
#' ########################
#'
#' res_epispot <- epispot(Y = Y, X = X, V = V, p0 = p0, list_init = list_init, 
#'                        user_seed = seed)
#'                        
#' @seealso  \code{\link{set_hyper}}, \code{\link{epispot}}
#'
#' @export
#'
set_init <- function(d, p, r, gam_vb, mu_beta_vb, om, s02, s2, sig2_beta_vb, 
                     tau_vb) {
  
  check_structure_(d, "vector", "numeric", 1)
  check_natural_(d)
  
  check_structure_(p, "vector", "numeric", 1)
  check_natural_(p)
  
  check_structure_(r, "vector", "numeric", 1)
  check_natural_(r)
  
  check_structure_(gam_vb, "matrix", "double", c(p, d))
  check_zero_one_(gam_vb)
  
  check_structure_(mu_beta_vb, "matrix", "double", c(p, d))
  
  check_structure_(om, "vector", "double", r)
  check_zero_one_(om)
  
  check_structure_(s02, "vector", "numeric", 1)
  check_positive_(s02)
  
  check_structure_(s2, "vector", "numeric", 1)
  check_positive_(s2)
  
  check_structure_(sig2_beta_vb, "vector", "double", d)
  check_positive_(sig2_beta_vb)
  
  check_structure_(tau_vb, "vector", "double", d)
  check_positive_(tau_vb)
  
  d_init <- d
  p_init <- p
  r_init <- r
  
  list_init <- create_named_list_(d_init, p_init, r_init, gam_vb, mu_beta_vb, om,
                                  s02, s2, sig2_beta_vb, tau_vb)
  
  class(list_init) <- "init"
  
  list_init
}


# Internal function setting default starting values when not provided by the user.
#
auto_set_init_ <- function(Y, p, p0, r, user_seed) {
  
  if (!is.null(user_seed)) set.seed(user_seed)
  
  # Initialisation not modified for dual = TRUE (should not matter, but maybe change this) ### TODO
  
  d <- ncol(Y)
  
  E_p_t <- p0[1]
  V_p_t <- p0[2]
  
  dn <- 1e-6
  up <- 1e5
  
  # Get n0 and t02 similarly as for a_omega_t and b_omega_t in HESS
  # (specify expectation and variance of number of active predictors per response)
  #
  # Look at : gam_st | theta_s = 0
  #
  tryCatch(t02 <- uniroot(function(x)
    get_V_p_t(get_mu(E_p_t, x, p), x, p) - V_p_t,
    interval = c(dn, up))$root,
    error = function(e) {
      stop(paste0("No hyperparameter values matching the expectation and variance ",
                  "of the number of active predictors per responses supplied in p0.",
                  "Please change p0."))
    })
  
  check_positive_(t02)
  
  # n0 sets the level of sparsity.
  # n0 <- - get_mu(E_p_t, t02, p)
  n0 <- get_mu(E_p_t, t02, p)
  
  s02 <- 1 / d
  check_positive_(s02)
  
  s2 <- 0.1
  check_positive_(s2)
  
  om <- rep(1 / r, r) 
  check_zero_one_(om)
  
  gam_vb <- matrix(pnorm(rnorm(p * d, mean = n0, sd = s02 + t02)), 
                   nrow = p)                                       
  check_zero_one_(gam_vb)
  
  mu_beta_vb <- matrix(rnorm(p * d), nrow = p)
  
  tau_vb <- 1 / median(apply(Y, 2, var))
  if (!is.finite(tau_vb)) tau_vb <- 1e3
  tau_vb <- rep(tau_vb, d)
  check_positive_(tau_vb)
  
  sig2_inv_vb <- 1e-2
  sig2_beta_vb <- 1 / rgamma(d, shape = 2, rate = 1 / (sig2_inv_vb * tau_vb))
  check_positive_(sig2_beta_vb)
  
  d_init <- d
  p_init <- p
  r_init <- r
  
  list_init <- create_named_list_(d_init, p_init, r_init, gam_vb, mu_beta_vb, om, 
                                  s02, s2, sig2_beta_vb, tau_vb)
  
  class(list_init) <- "out_init"
  
  list_init
}
