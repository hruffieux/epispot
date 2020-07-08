## EPISPOT – annotation-driven approach for large-scale joint regression with multiple responses

## Overview

**epispot** is an R package for parallel variational expectation-maximisation in 
sparse regression with hierarchically-related responses and predictor-level 
information. A notable application is in the context of molecular quantitative 
trait locus mapping, with several thousand molecular levels (responses), genetic 
markers (candidate predictors) and epigenetic annotation marks (predictor-level 
covariates) and individuals (samples).

## Warning

**This is a development branch**, it is not guaranteed to be stable at any given 
time and features are subject to change. Please use the 
[stable version](https://github.com/hruffieux/epispot), unless you want to test 
and report issues.

## Installation

To install, run the following commands in R:

``` r
# after having installed devtools (install.packages("devtools"))
devtools::install_github("hruffieux/epispot", ref = "devel")
```

## License and authors

This software uses the GPL v3 license, see [LICENSE](LICENSE).
Authors and copyright are provided in [DESCRIPTION](DESCRIPTION). 

## Issues

To report an issue, please use the 
[epispot issue tracker](https://github.com/hruffieux/epispot/issues) at 
github.com.
