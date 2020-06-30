#' epispot: a package for combined predictor and outcome selection in
#' high-dimensional set-ups using variational inference
#'
#' The epispot package provides a parallel variational expectation-maximisation 
#' algorithm for sparse regression with hierarchically-related responses and 
#' predictor-level information. This software has been used in the context of 
#' molecular quantitative trait locus mapping, with several thousand molecular 
#' levels (responses), genetic markers (candidate predictors) and epigenomic 
#' annotation marks (predictor-level covariates) and individuals (samples).
#'
#' @section epispot functions: set_hyper, set_init, epispot, set_blocks.
#'
#' @docType package
#' @name epispot-package
#' @useDynLib epispot, .registration = TRUE
#' @import RcppEigen
#' @importFrom Rcpp evalCpp
#' @importFrom stats cor dnorm median pnorm qnorm rbeta rbinom rgamma rnorm setNames uniroot var
#' @importFrom grDevices dev.off png
#' @importFrom graphics abline legend matplot points
NULL
