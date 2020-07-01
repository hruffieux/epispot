# This file is part of the `epispot` R package:
#     https://github.com/hruffieux/epispot
#

#' Gather model hyperparameters provided by the user.
#'
#' This function must be used to provide hyperparameter values for the model
#' used in \code{\link{epispot}}.
#'
#' The \code{\link{epispot}} function can also be used with default hyperparameter
#' choices (without using \code{\link{set_hyper}}) by setting its argument
#' \code{list_hyper} to \code{NULL}.
#'
#' @param d Number of responses.
#' @param p Number of candidate predictors.
#' @param lambda Vector of length 1 providing the values of hyperparameter
#'   \eqn{\lambda} for the prior distribution of \eqn{\sigma^{-2}}. \eqn{\sigma}
#'   represents the typical size of nonzero effects.
#' @param nu Vector of length 1 providing the values of hyperparameter \eqn{\nu}
#'   for the prior distribution of \eqn{\sigma^{-2}}. \eqn{\sigma} represents
#'   the typical size of nonzero effects.
#' @param a Vector of length 1 or p providing the values of hyperparameter
#'   \eqn{a} for the prior distributions for the proportion of responses
#'   associated with each candidate predictor, \eqn{\omega} (vector of length p).
#'   If of length 1, the provided value is repeated p times.
#' @param b Vector of length 1 or p providing the values of hyperparameter
#'   \eqn{b} for the prior distributions for the proportion of responses
#'   associated with each candidate predictor, \eqn{\omega} (vector of length p).
#'   If of length 1, the provided value is repeated p times.
#' @param eta Vector of length 1 or d. Provides the values of
#'   hyperparameter \eqn{\eta} for the prior distributions of the continuous
#'   response residual precisions, \eqn{\tau}. If of length 1, the provided
#'   value is repeated d times.
#' @param kappa Vector of length 1 or d. Provides the values of hyperparameter
#'   \eqn{\kappa} for the prior distributions of the response residual
#'   precisions, \eqn{\tau}. If of length 1, the provided value is repeated d 
#'   times.
#' @param r Number of variables representing external information on the
#'   candidate predictors. Default is \code{NULL}, for \code{V} \code{NULL}.
#' @param n0 Vector of length 1 or d. 
#' @param s02 Variance hyperparameter when \code{V} is 
#' non-\code{NULL} non-\code{NULL}. Default is \code{NULL}.
#' @param s2 Variance hyperparameter when \code{V} is non-\code{NULL}
#'   non-\code{NULL}. Default is \code{NULL}.
#' @param t02 Variance hyperparameter. 
#'
#' @return An object of class "\code{hyper}" preparing user hyperparameter in a
#'   form that can be passed to the \code{\link{epispot}} function.
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
#' n <- 200; p <- 200; p0 <- 20; d <- 20; d0 <- 15; q <- 2; r <- 3
#'
#' ## Candidate predictors (subject to selection)
#' ##
#' # Here we simulate common genetic variants (but any type of candidate
#' # predictors can be supplied).
#' # 0 = homozygous, major allele, 1 = heterozygous, 2 = homozygous, minor allele
#'
#' X_act <- matrix(rbinom(n * p0, size = 2, p = 0.25), nrow = n)
#' X_inact <- matrix(rbinom(n * (p - p0), size = 2, p = 0.25), nrow = n)
#'
#' shuff_x_ind <- sample(p)
#' X <- cbind(X_act, X_inact)[, shuff_x_ind]
#'
#' bool_x_act <- shuff_x_ind <= p0
#'
#' pat_act <- beta <- matrix(0, nrow = p0, ncol = d0)
#' pat_act[sample(p0*d0, floor(p0*d0/5))] <- 1
#' beta[as.logical(pat_act)] <-  rnorm(sum(pat_act))
#'
#' ## Covariates (not subject to selection)
#' ##
#' Z <- matrix(rnorm(n * q), nrow = n)
#'
#' alpha <-  matrix(rnorm(q * d), nrow = q)
#'
#' ## Gaussian responses
#' ##
#' Y_act <- matrix(rnorm(n * d0, mean = X_act %*% beta, sd = 0.5), nrow = n)
#' Y_inact <- matrix(rnorm(n * (d - d0), sd = 0.5), nrow = n)
#' shuff_y_ind <- sample(d)
#' Y <- cbind(Y_act, Y_inact)[, shuff_y_ind] + Z %*% alpha
#'
#' ## Binary responses
#' ##
#' Y_bin <- ifelse(Y > 0, 1, 0)
#' ## Informative annotation variables
#' ##
#' V <- matrix(rnorm(p * r), nrow = p)
#' V[bool_x_act, ] <- rnorm(p0 * r, mean = 2)
#'
#' ########################
#' ## Infer associations ##
#' ########################
#'
#' ## Continuous responses
#' ##
#'
#' # No covariate
#' #
#' # a and b chosen so that the prior mean number of responses associated with
#' # each candidate predictor is 1/4.
#' list_hyper_g <- set_hyper(d, p, lambda = 1, nu = 1, a = 1, b = 4*d-1,
#'                           eta = 1, kappa = apply(Y, 2, var),
#'                           link = "identity")
#'
#' # We take p0 = p0 (known here); this choice may result in variable
#' # selections that are (too) conservative in some cases. In practice, it is
#' # advised to set p0 as a slightly overestimated guess of p0, or perform
#' # cross-validation using function `set_cv'.
#'
#' vb_g <- epispot(Y = Y, X = X, p0 = p0, link = "identity",
#'               list_hyper = list_hyper_g, user_seed = seed)
#'
#' # With covariates
#' #
#' list_hyper_g_z <- set_hyper(d, p, lambda = 1, nu = 1, a = 1, b = 4*d-1,
#'                             eta = 1, kappa = apply(Y, 2, var),
#'                             link = "identity", q = q, phi = 1, xi = 1)
#'
#' vb_g_z <- epispot(Y = Y, X = X, p0 = p0, Z = Z, link = "identity",
#'                 list_hyper = list_hyper_g_z, user_seed = seed)
#'
#'
#' # With external annotation variables
#' #
#' list_hyper_g_v <- set_hyper(d, p, lambda = 1, nu = 1, a = NULL, b = NULL,
#'                             eta = 1, kappa = apply(Y, 2, var),
#'                             link = "identity", r = r, m0 = 0, s02 = 0.1,
#'                             s2 = 0.001)
#'
#' vb_g_v <- epispot(Y = Y, X = X, p0 = p0,  V = V, link = "identity",
#'                 list_hyper = list_hyper_g_v, user_seed = seed)
#'
#' ## Binary responses
#' ##
#' list_hyper_logit <- set_hyper(d, p, lambda = 1, nu = 1, a = 1, b = 4*d-1,
#'                               eta = NULL, kappa = NULL, link = "logit",
#'                               q = q, phi = 1, xi = 1)
#'
#' vb_logit <- epispot(Y = Y_bin, X = X, p0 = p0, Z = Z, link = "logit",
#'                   list_hyper = list_hyper_logit, user_seed = seed)
#'
#' list_hyper_probit <- set_hyper(d, p, lambda = 1, nu = 1, a = 1, b = 4*d-1,
#'                                eta = NULL, kappa = NULL, link = "probit",
#'                                q = q, phi = 1, xi = 1)
#'
#' vb_probit <- epispot(Y = Y_bin, X = X, p0 = p0, Z = Z, link = "probit",
#'                    list_hyper = list_hyper_probit, user_seed = seed)
#'
#'
#' ## Mix of continuous and binary responses
#' ##
#' Y_mix <- cbind(Y, Y_bin)
#' ind_bin <- (d+1):(2*d)
#'
#' list_hyper_mix <- set_hyper(2*d, p, lambda = 1, nu = 1, a = 1, b = 8*d-1,
#'                             eta = 1, kappa = apply(Y, 2, var), link = "mix",
#'                             ind_bin = ind_bin, q = q, phi = 1, xi = 1)
#'
#' vb_mix <- epispot(Y = Y_mix, X = X, p0 = p0, Z = Z, link = "mix",
#'                 ind_bin = ind_bin, list_hyper = list_hyper_mix,
#'                 user_seed = seed)
#'
#' @seealso  \code{\link{set_init}}, \code{\link{epispot}}
#'
#' @export
#'
set_hyper <- function(d, p, r, lambda, nu, eta, kappa, n0, t02) {
  
  check_structure_(d, "vector", "numeric", 1)
  check_natural_(d)
  
  check_structure_(p, "vector", "numeric", 1)
  check_natural_(p)
  
  check_structure_(r, "vector", "numeric", 1)
  check_natural_(r)
  
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
  r_hyper <- r
  
  list_hyper <- create_named_list_(d_hyper, p_hyper, r_hyper,
                                   eta, kappa, lambda, nu, n0, t02)
  
  class(list_hyper) <- "hyper"
  
  list_hyper
  
}


# Internal function setting default model hyperparameters when not provided by
# the user.
#
auto_set_hyper_ <- function(Y, p, p0, r) {
  
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
  r_hyper <- r
  
  list_hyper <- create_named_list_(d_hyper, p_hyper, r_hyper,
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
#' @param gam_vb Matrix of size p x d with initial values for the variational
#'   parameter yielding posterior probabilities of inclusion.
#' @param mu_beta_vb Matrix of size p x d with initial values for the
#'   variational parameter yielding regression coefficient estimates for
#'   predictor-response pairs included in the model.
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
#' seed <- 123; set.seed(seed)
#'
#' ###################
#' ## Simulate data ##
#' ###################
#'
#' ## Examples using small problem sizes:
#' ##
#' n <- 200; p <- 200; p0 <- 20; d <- 20; d0 <- 15; q <- 2
#'
#' ## Candidate predictors (subject to selection)
#' ##
#' # Here we simulate common genetic variants (but any type of candidate
#' # predictors can be supplied).
#' # 0 = homozygous, major allele, 1 = heterozygous, 2 = homozygous, minor allele
#'
#' X_act <- matrix(rbinom(n * p0, size = 2, p = 0.25), nrow = n)
#' X_inact <- matrix(rbinom(n * (p - p0), size = 2, p = 0.25), nrow = n)
#'
#' shuff_x_ind <- sample(p)
#' X <- cbind(X_act, X_inact)[, shuff_x_ind]
#'
#' bool_x_act <- shuff_x_ind <= p0
#'
#' pat_act <- beta <- matrix(0, nrow = p0, ncol = d0)
#' pat_act[sample(p0*d0, floor(p0*d0/5))] <- 1
#' beta[as.logical(pat_act)] <-  rnorm(sum(pat_act))
#'
#' ## Covariates (not subject to selection)
#' ##
#' Z <- matrix(rnorm(n * q), nrow = n)
#'
#' alpha <-  matrix(rnorm(q * d), nrow = q)
#'
#' ## Gaussian responses
#' ##
#' Y_act <- matrix(rnorm(n * d0, mean = X_act %*% beta, sd = 0.5), nrow = n)
#' Y_inact <- matrix(rnorm(n * (d - d0), sd = 0.5), nrow = n)
#' shuff_y_ind <- sample(d)
#' Y <- cbind(Y_act, Y_inact)[, shuff_y_ind] + Z %*% alpha
#'
#' ## Binary responses
#' ##
#' Y_bin <- ifelse(Y > 0, 1, 0)
#'
#' ########################
#' ## Infer associations ##
#' ########################
#'
#' ## Continuous responses
#' ##
#'
#' # No covariate
#' #
#' # gam_vb chosen so that the prior mean number of responses associated with
#' # each candidate predictor is 1/4.
#' gam_vb <- matrix(rbeta(p * d, shape1 = 1, shape2 = 4*d-1), nrow = p)
#' mu_beta_vb <- matrix(rnorm(p * d), nrow = p)
#' tau_vb <- 1 / apply(Y, 2, var)
#' sig2_beta_vb <- 1 / rgamma(d, shape = 2, rate = 1 / tau_vb)
#'
#' list_init_g <- set_init(d, p, gam_vb, mu_beta_vb, sig2_beta_vb, tau_vb,
#'                         link = "identity")
#'
#' # We take p0 = p0 (known here); this choice may result in variable
#' # selections that are (too) conservative in some cases. In practice, it is
#' # advised to set p0 as a slightly overestimated guess of p0, or perform
#' # cross-validation using function `set_cv'.
#'
#' vb_g <- epispot(Y = Y, X = X, p0 = p0, link = "identity",
#'               list_init = list_init_g)
#'
#' # With covariates
#' #
#' mu_alpha_vb <- matrix(rnorm(q * d), nrow = q)
#' sig2_alpha_vb <- 1 / matrix(rgamma(q * d, shape = 2, rate = 1), nrow = q)
#'
#' list_init_g_z <- set_init(d, p, gam_vb, mu_beta_vb, sig2_beta_vb, tau_vb,
#'                           link = "identity", q = q,
#'                           mu_alpha_vb = mu_alpha_vb,
#'                           sig2_alpha_vb = sig2_alpha_vb)
#'
#' vb_g_z <- epispot(Y = Y, X = X, p0 = p0, Z = Z, link = "identity",
#'                 list_init = list_init_g_z)
#'
#' ## Binary responses
#' ##
#' # gam_vb chosen so that the prior mean number of responses associated with
#' # each candidate predictor is 1/4.
#' sig2_beta_vb_logit <- 1 / t(replicate(p, rgamma(d, shape = 2, rate = 1)))
#'
#' list_init_logit <- set_init(d, p, gam_vb, mu_beta_vb, sig2_beta_vb_logit,
#'                             tau_vb = NULL, link = "logit", q = q,
#'                             mu_alpha_vb = mu_alpha_vb,
#'                             sig2_alpha_vb = sig2_alpha_vb)
#'
#' vb_logit <- epispot(Y = Y_bin, X = X, p0 = p0, Z = Z, link = "logit",
#'                   list_init = list_init_logit)
#'
#' sig2_alpha_vb_probit <- sig2_alpha_vb[, 1]
#' sig2_beta_vb_probit <- sig2_beta_vb[1]
#' list_init_probit <- set_init(d, p, gam_vb, mu_beta_vb, sig2_beta_vb_probit,
#'                              tau_vb = NULL, link = "probit", q = q,
#'                              mu_alpha_vb = mu_alpha_vb,
#'                              sig2_alpha_vb = sig2_alpha_vb_probit)
#'
#' vb_probit <- epispot(Y = Y_bin, X = X, p0 = p0, Z = Z, link = "probit",
#'                    list_init = list_init_probit)
#'
#' ## Mix of continuous and binary responses
#' ##
#' Y_mix <- cbind(Y, Y_bin)
#' ind_bin <- (d+1):(2*d)
#'
#' # gam_vb chosen so that the prior mean number of responses associated with
#' # each candidate predictor is 1/4.
#' gam_vb_mix <- matrix(rbeta(p * 2*d, shape1 = 1, shape2 = 8*d-1), nrow = p)
#' mu_beta_vb_mix <- matrix(rnorm(p * 2*d), nrow = p)
#' sig2_beta_vb_mix <- 1 / c(rgamma(d, shape = 2, rate = 1 / tau_vb),
#'                           rgamma(d, shape = 2, rate = 1))
#' mu_alpha_vb_mix <- matrix(rnorm(q * 2*d), nrow = q)
#' sig2_alpha_vb_mix <- 1 / matrix(rgamma(q * 2*d, shape = 2, rate = 1), nrow = q)
#'
#' list_init_mix <- set_init(2*d, p, gam_vb_mix, mu_beta_vb_mix,
#'                           sig2_beta_vb_mix, tau_vb, link = "mix",
#'                           ind_bin = ind_bin, q = q,
#'                           mu_alpha_vb = mu_alpha_vb_mix,
#'                           sig2_alpha_vb = sig2_alpha_vb_mix)
#'
#' vb_mix <- epispot(Y = Y_mix, X = X, p0 = p0, Z = Z, link = "mix",
#'                 ind_bin = ind_bin, list_init = list_init_mix)
#'
#' @seealso  \code{\link{set_hyper}}, \code{\link{epispot}}
#'
#' @export
#'
set_init <- function(d, p, r, gam_vb, mu_beta_vb, om, s02, s2, sig2_beta_vb, tau_vb) {
  
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
  
  # Initialisation not modified for dual = TRUE (should not matter, but maybe change this) ### TODO
  
  d <- ncol(Y)
  
  if (!is.null(user_seed)) set.seed(user_seed)
  
  
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

  gam_vb <- matrix(pnorm(rnorm(p * d, mean = n0, sd = s02 + t02)), # Phi(theta + chi), and not 1 - Phi(theta + chi)
                   nrow = p)                                            # as reparametrisation theta* = - theta, chi* = - chi
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
