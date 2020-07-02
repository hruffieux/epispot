rm(list = ls())

seed <- 123
set.seed(seed)

############################
## simulate basic dataset ##
############################

n <- 100; p <- 75; d <- 20; p_act <- 10; r <- 10

# candidate predictors (subject to selection)
X_act <- matrix(rbinom(n * p_act, size = 2, p = 0.2), nrow = n)
X_inact <- matrix(rbinom(n * (p - p_act), size = 2, p = 0.2), nrow = n)
X <- cbind(X_act, X_inact)[, sample(p)]

beta <-  matrix(rnorm(p_act * d), nrow = p_act)

# Gaussian outcomes
Y <- matrix(rnorm(n * d, mean = X_act %*% beta, sd = 1), nrow = n)

# remove constant variables (needed for checking dimension consistency)
X <- scale(X)
rm_cst <- function(mat_sc) mat_sc[, !is.nan(colSums(mat_sc))]
rm_coll <- function(mat_sc) mat_sc[, !duplicated(mat_sc, MARGIN = 2)]

X <- rm_cst(X)
X <- rm_coll(X)

p <- ncol(X)

V <- matrix(rnorm(p * r), nrow = p)


#############
## EPISPOT ##
#############

p0 <- c(5, 25)

# Continuous outcomes, no covariates
#
res_epispot <- epispot(Y = Y, X = X, V = V, p0 = p0, user_seed = seed)