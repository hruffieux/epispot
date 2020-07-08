# This file is part of the `epispot` R package:
#     https://github.com/hruffieux/epispot
#

# Diverse utility functions implementing sanity checks, basic preprocessing,
# and ticks to prevent overflow/underflow.
#

check_natural_ <- function(x, eps = .Machine$double.eps^0.75){
  if (any(x < eps | abs(x - round(x)) > eps)) {
    stop(paste(deparse(substitute(x)),
               " must be natural.", sep=""))
  }
}

check_positive_ <- function(x, eps = .Machine$double.eps^0.75){
  if (any(x < eps)) {
    err_mess <- paste(deparse(substitute(x)), " must be positive, greater than ",
                      format(eps, digits = 3), ".", sep="")
    if (length(x) > 1) err_mess <- paste("All entries of ", err_mess, sep="")
    stop(err_mess)
  }
}

check_zero_one_ <- function(x){
  if (any(x < 0) | any(x > 1)) {
    err_mess <- paste(deparse(substitute(x)), " must lie between 0 and 1.", sep="")
    if (length(x) > 1) err_mess <- paste("All entries of ", err_mess, sep="")
    stop(err_mess)
  }
}

check_binary_ <-function(x) {
  identical(as.vector(x), as.numeric(as.logical(x)))
}


check_structure_ <- function(x, struct, type, size = NULL,
                             null_ok = FALSE,  inf_ok = FALSE, na_ok = FALSE) {
  if (type == "double") {
    bool_type <-  is.double(x)
    type_mess <- "a double-precision "
  } else if (type == "integer") {
    bool_type <- is.integer(x)
    type_mess <- "an integer "
  } else if (type == "numeric") {
    bool_type <- is.numeric(x)
    type_mess <- "a numeric "
  } else if (type == "logical") {
    bool_type <- is.logical(x)
    type_mess <- "a boolean "
  } else if (type == "string") {
    bool_type <- is.character(x)
    type_mess <- "string "
  }

  bool_size <- TRUE # for case size = NULL (no assertion on the size/dimension)
  size_mess <- ""
  if (struct == "vector") {
    bool_struct <- is.vector(x) & (length(x) > 0) # not an empty vector
    if (!is.null(size)) {
      bool_size <- length(x) %in% size
      size_mess <- paste(" of length ", paste(size, collapse=" or "), sep = "")
    }
  } else if (struct == "matrix") {
    bool_struct <- is.matrix(x) & (length(x) > 0) # not an empty matrix
    if (!is.null(size)) {
      bool_size <- all(dim(x) == size)
      size_mess <- paste(" of dimension ", size[1], " x ", size[2], sep = "")
    }
  }

  correct_obj <- bool_struct & bool_type & bool_size

  bool_null <- is.null(x)

  if (!is.list(x) & type != "string") {
    na_mess <- ""
    if (!na_ok) {
      if (!bool_null) correct_obj <- correct_obj & !any(is.na(x))
      na_mess <- " without missing value"
    }

    inf_mess <- ""
    if (!inf_ok) {
      if (!bool_null) correct_obj <- correct_obj & all(is.finite(x[!is.na(x)]))
      inf_mess <- ", finite"
    }
  } else {
    na_mess <- ""
    inf_mess <- ""
  }

  null_mess <- ""
  if (null_ok) {
    correct_obj <- correct_obj | bool_null
    null_mess <- " or must be NULL"
  }

  if(!(correct_obj)) {
    stop(paste(deparse(substitute(x)), " must be a non-empty ", type_mess, struct,
               size_mess, inf_mess, na_mess, null_mess, ".", sep = ""))
  }
}


create_named_list_ <- function(...) {
  setNames(list(...), as.character(match.call()[-1]))
}



get_annealing_ladder_ <- function(anneal_schedule, verbose) {

  # ladder set following:
  # Importance Tempering, Robert B. Gramacy & Richard J. Samworth, pp.9-10, arxiv v4

  k_m <- 1 / anneal_schedule[2]
  m <- anneal_schedule[3]

  if(anneal_schedule[1] == 1) {

    type <- "geometric"

    delta_k <- k_m^(1 / (1 - m)) - 1

    ladder <- (1 + delta_k)^(1 - m:1)

  } else if (anneal_schedule[1] == 2) { # harmonic spacing

    type <- "harmonic"

    delta_k <- ( 1 / k_m - 1) / (m - 1)

    ladder <- 1 / (1 + delta_k * (m:1 - 1))

  } else { # linear spacing

    type <- "linear"

    delta_k <- (1 - k_m) / (m - 1)

    ladder <- k_m + delta_k * (1:m - 1)
  }

  if (verbose)
    cat(paste0("** Annealing with ", type," spacing ** \n\n"))

  ladder

}


log_one_plus_exp_ <- function(x) { # computes log(1 + exp(x)) avoiding
                                   # numerical overflow
  m <- x
  m[x < 0] <- 0

  log(exp(x - m) + exp(- m)) + m
}


log_sigmoid_ <- function(chi) {

  - log(1 + exp(- chi)) # chi is always positive so no overflow possible (underflow neither, thanks to the "+1")

}

log_det <- function(list_mat) {

  if (is.list(list_mat)) {
    sapply(list_mat, function(mat) {
      log_det <- determinant(mat, logarithm = TRUE)
      log_det$modulus * log_det$sign
    })
  } else {
    log_det <- determinant(list_mat, logarithm = TRUE)
    log_det$modulus * log_det$sign
  }

}


inv_mills_ratio_ <- function(y, U, log_1_pnorm_U, log_pnorm_U) {

  stopifnot(y %in% c(0, 1))

  # writing explicitely the formula for pnorm(, log = TRUE) is faster...
  if (y == 1) {

    m <- exp(-U^2/2 - log(sqrt(2*pi)) - log_pnorm_U)
    m[m < -U] <- -U

  } else {

    m <- - exp(-U^2/2 - log(sqrt(2*pi)) - log_1_pnorm_U)
    m[m > -U] <- -U

  }

  m

}


inv_mills_ratio_matrix_ <- function(Y, U) {

  if (is.matrix(U)) m <- matrix(NA, nrow = nrow(U), ncol = ncol(U))
  else m <- rep(NA, length(U))

  U_1 <- U[Y==1]
  m_1 <- exp(dnorm(U_1, log = TRUE) - pnorm(U_1, log.p = TRUE))
  m_1[m_1 < -U_1] <- -U_1

  m[Y==1] <- m_1


  U_0 <- U[Y==0]
  m_0 <- - exp(dnorm(U_0, log = TRUE) - pnorm(U_0, lower.tail = FALSE, log.p = TRUE))
  m_0[m_0 > -U_0] <- -U_0

  m[Y==0] <- m_0

  m

}


log_sum_exp_ <- function(x) {
  # Computes log(sum(exp(x))

  if ( max(abs(x)) > max(x) )
    offset <- min(x)
  else
    offset <- max(x)
  log(sum(exp(x - offset))) + offset

}


# Functions for hyperparameter settings
#
E_Phi_X <- function(mu, s2, lower_tail = TRUE) {

  pnorm(mu / sqrt(1 + s2), lower.tail = lower_tail)

}

E_Phi_X_2 <- function(mu, s2) {

  pnorm(mu / sqrt(1 + s2)) -
    2 * PowerTOST::OwensT(mu / sqrt(1 + s2), 1 / sqrt(1 + 2 * s2))

}

get_V_p_t <- function(mu, s2, p) {
  p * (p - 1) * E_Phi_X_2(mu, s2) -
    p^2 * E_Phi_X(mu, s2)^2 +
    p * E_Phi_X(mu, s2)
}


get_mu <- function(E_p_t, s2, p) {

  sqrt(1 + s2) * qnorm(E_p_t / p)
  
}


get_n0_t02 <- function(q, p, p0) {
  
  E_p_t <- p0[1]
  V_p_t <- min(p0[2], floor(2 * p / 3))

  dn <- 1e-6
  up <- 1e5
  
  # Get n0 and t02
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
  
  # n0 sets the level of sparsity.
  n0 <- get_mu(E_p_t, t02, p)
  n0 <- rep(n0, q)
  
  create_named_list_(n0, t02)
}

rm_bin_annot_freq_ <- function(mat, bin_annot_freq, verbose) {
  
  if (!is.null(bin_annot_freq)) {

    bool_bin_annot_freq <- apply(mat, 2, function(x) {
      if (check_binary_(x) & 
          (sum(x) < bin_annot_freq * length(x) | sum(1-x) < bin_annot_freq * length(x))) {
        TRUE
      } else {
        FALSE
      }
    })
    
    if (any(bool_bin_annot_freq)) {
      
      rmvd_bin_annot_freq <- colnames(mat)[bool_bin_annot_freq]
      
      if (verbose) {
        if (sum(bool_bin_annot_freq) < 50) {
          cat(paste0("Annotation variable(s) ", 
                     paste(rmvd_bin_annot_freq, collapse=", "),
                     " concern(s) less than ", bin_annot_freq * 100, " or more than ",
                     (1-bin_annot_freq) * 100," candidate predictors. \n",
                     "Removing corresponding column(s) and saving its/their id(s) ",
                     "in the function output ... \n\n"))
        } else {
          cat(paste0(sum(bool_bin_annot_freq), " annotation variable(s) ", 
                     " concern less than ", bin_annot_freq * 100, " or more than ",
                     (1-bin_annot_freq) * 100," candidate predictors. \n",
                     "Removing corresponding column(s) and saving its/their id(s) ",
                     "in the function output ... \n\n"))
        }
      } 
      
      mat <- mat[, !bool_bin_annot_freq, drop = FALSE]
      
    } else {
      
      rmvd_bin_annot_freq  <- NULL
      
    }
    
  } else {
    
    bool_bin_annot_freq <- rep(FALSE, ncol(mat))
    rmvd_bin_annot_freq <- NULL
    
  }
  

  create_named_list_(mat, bool_bin_annot_freq, rmvd_bin_annot_freq)
  
}



rm_constant_ <- function(mat, verbose) {

  bool_cst <- is.nan(colSums(mat))

  if (any(bool_cst)) {

    rmvd_cst <- colnames(mat)[bool_cst]

    if (verbose) {
      if (sum(bool_cst) < 50) {
        cat(paste0("Variable(s) ", paste(rmvd_cst, collapse=", "), " constant. \n",
                   "Removing corresponding column(s) and saving its/their id(s) ",
                   "in the function output ... \n\n"))
      } else {
        cat(paste0(sum(bool_cst), " variables constant. \n Removing corresponding ", 
                   "column(s) and saving their ids in the function output ... \n\n"))
      }
    }

    mat <- mat[, !bool_cst, drop = FALSE]
  } else {
    rmvd_cst <- NULL
  }

  create_named_list_(mat, bool_cst, rmvd_cst)
}

rm_collinear_ <- function(mat, verbose) {

  bool_coll <- duplicated(mat, MARGIN = 2)

  if (any(bool_coll)) {

    mat_coll <- mat[, bool_coll, drop = FALSE]
    rmvd_coll <- colnames(mat_coll)

    if (verbose) {
      if (length(rmvd_coll) < 50) {
        cat(paste0("Presence of collinear variable(s). ",
                   paste(rmvd_coll, collapse=", "), " redundant. \n",
                   "Removing corresponding column(s) and saving its/their id(s) ",
                   "in the function output ... \n"))
      } else {
        cat(paste0("Presence of collinear variables. ", length(rmvd_coll),
                   " redundant.\n", "Removing corresponding columns and saving ",
                   "their ids in the function output ... \n"))
      }
    }

    # associate to each removed replicate the name of the covariate with which
    # it is duplicated and that is kept in the dataset
    bool_with_coll <- duplicated(mat, MARGIN = 2, fromLast = TRUE) & !bool_coll
    mat_with_coll <- mat[, bool_with_coll, drop = FALSE]

    assoc_coll <- colnames(mat_with_coll)[match(data.frame(mat_coll),
                                                data.frame(mat_with_coll))]
    names(rmvd_coll) <- assoc_coll

    mat <- mat[, !bool_coll, drop = FALSE]

  } else {
    rmvd_coll <- NULL
  }

  create_named_list_(mat, bool_coll, rmvd_coll)
}


make_xihunks_ <- function(x, n_g) split(x, factor(sort(rank(x) %% n_g)))


cbind_fill_matrix <- function(...) { # more efficient than do.call cbind
  tr <- lapply(..., as.matrix)
  tr <- lapply(..., t)
  t(as.matrix(plyr::rbind.fill.matrix(tr)))
}
