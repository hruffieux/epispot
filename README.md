<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- First time: run usethis::use_readme_rmd() to create a pre-commit hook that 
prevents from committing if the README.Rmd has changed, but has not been 
re-knitted to generate an updated README.md -->

## EPISPOT – epigenome-informed approach for detecting and interpreting QTL effects <img src="man/figures/epispot_logo.png" align="right" height="150"/>

<!-- Run for the R CMD checks, run usethis::use_github_actions() to set up the pipeline, possibly modify the .yaml file and then: -->

<!-- [![R build status](https://github.com/hruffieux/epispot/workflows/R-CMD-check/badge.svg)](https://github.com/hruffieux/epispot/actions)  -->

[![License: GPL
v3](https://img.shields.io/badge/license-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![](https://img.shields.io/badge/devel%20version-0.1.3-blue.svg)](https://github.com/hruffieux/epispot)
[![](https://img.shields.io/github/languages/code-size/hruffieux/epispot.svg)](https://github.com/hruffieux/epispot)
[![](https://img.shields.io/badge/doi-10.1101/2020.09.21.305789-yellow.svg)](https://doi.org/10.1101/2020.09.21.305789)

## Overview

**epispot** is an R package for molecular QTL mapping that integrates
large-scale epigenetic annotations to enhance the detection and
interpretation of *cis-* and *trans-* acting QTLs, including hotspots.
It employs a Bayesian variational expectation-maximisation approach to
jointly analyse thousands of molecular phenotypes (e.g., gene, protein
or metabolite levels) and genetic variants while incorporating
epigenetic marks as candidate auxiliary data on the potential of genetic
variants to control traits.

H. Ruffieux, B. Fairfax, I. Nassiri, E. Vigorito, C. Wallace, S.
Richardson, L. Bottolo. EPISPOT : An epigenome-driven approach for
detecting and interpreting hotspots in molecular QTL studies, *The
American Journal of Human Genetics*, 108:1-18,
10.1016/j.ajhg.2021.04.010, 2021.

## Warning

**This is a development branch**, it is not guaranteed to be stable at
any given time and features are subject to change. Please use the
[stable version](https://github.com/hruffieux/epispot), unless you want
to test and report issues.

## Installation

To install, run the following commands in R:

``` r
if(!require(remotes)) install.packages("remotes")
remotes::install_github("hruffieux/epispot")
```

## License and authors

This software uses the GPL v3 license, see [LICENSE](LICENSE). Authors
and copyright are provided in [DESCRIPTION](DESCRIPTION).

## Issues

To report an issue, please use the [epispot issue
tracker](https://github.com/hruffieux/epispot/issues) at github.com.
