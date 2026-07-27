# libcmaesr

[![R-CMD-check](https://github.com/mlr-org/libcmaesr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mlr-org/libcmaesr/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/libcmaesr)](https://cran.r-project.org/package=libcmaesr)

A lightweight R interface to the [libcmaes](https://github.com/CMA-ES/libcmaes) C++ library for Covariance Matrix Adaptation Evolution Strategy (CMA-ES). It allows for the optimization of black-box functions using the CMA-ES algorithm and its variants.

A patched copy of libcmaes is bundled with the package, so no system dependencies are needed beyond a standard C++ toolchain (Eigen headers are provided by the RcppEigen package at build time).

## Installation

Install the released version from CRAN with:

```r
install.packages("libcmaesr")
```

Or install the development version from GitHub with:

```r
# install.packages("pak")
pak::pak("mlr-org/libcmaesr")
```

## Example

This is a basic example which shows you how to solve a common test problem, the sphere function:

```r
library(libcmaesr)

# define objective function
dim = 3
lambda = 5

fn = function(x) {
  sum(x^2)
}

x0 = rep(0.5, dim)
lower = rep(-1, dim)
upper = rep(1, dim)
fevals = 5000 * dim
algo = "bipop"
ctrl = cmaes_control(algo = algo, max_fevals = fevals, seed = 123, lambda = lambda)
res = cmaes(fn, x0, lower, upper, ctrl)
print(res)
```
