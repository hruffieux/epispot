#ifndef EPISPOT_UTILS_H_
#define EPISPOT_UTILS_H_

#include <RcppEigen.h>
#include "epispot_types.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

Arr1D logOnePlusExp(const  Arr1D& x);

Arr2D logOnePlusExpMat(const  Arr2D& x);

#endif // EPISPOT_UTILS_H_
