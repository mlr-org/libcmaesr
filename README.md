
<!-- README.md is generated from README.Rmd. Please edit that file -->

# libcmaesr

<!-- badges: start -->

[![R-CMD-check](https://github.com/mlr-org/libcmaesr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mlr-org/libcmaesr/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/libcmaesr)](https://cran.r-project.org/package=libcmaesr)
<!-- badges: end -->

A lightweight R interface to the
[libcmaes](https://github.com/CMA-ES/libcmaes) C++ library for
covariance matrix adaptation evolution strategy (CMA-ES). It allows for
the optimization of black-box functions using the CMA-ES algorithm and
its variants.

A patched copy of libcmaes is bundled with the package, so no system
dependencies are needed beyond a standard C++ toolchain (Eigen headers
are provided by the RcppEigen package at build time).

## Installation

Install the released version from CRAN with:

``` r
install.packages("libcmaesr")
```

Or install the development version from GitHub with:

``` r
# install.packages("pak")
pak::pak("mlr-org/libcmaesr")
```

## Example

This is a basic example which shows you how to solve a common test
problem, the sphere function. The objective function receives a numeric
vector and returns a single number. `cmaes_control()` collects the
settings of the algorithm, and `cmaes()` runs the optimization.

``` r
library(libcmaesr)

fn = function(x) sum(x^2)

dim = 3
x0 = rep(0.5, dim)
lower = rep(-1, dim)
upper = rep(1, dim)

control = cmaes_control(algo = "bipop", max_fevals = 5000 * dim, seed = 123, lambda = 5)
res = cmaes(fn, x0, lower, upper, control)

res$x
#> [1] -1.102e-09 -4.948e-10  1.933e-09
res$y
#> [1] 5.195e-18
res$status_msg
#> [1] "[Success] The optimization has converged"
```

## Batch evaluation

Set `batch = TRUE` to evaluate a whole generation at once. The objective
function then receives a matrix with `lambda` rows and one column per
dimension, and returns one objective value per row. This is useful when
the evaluation can be vectorized or parallelized.

``` r
fn_batch = function(x) rowSums(x^2)

control = cmaes_control(algo = "acmaes", max_fevals = 1000, seed = 123, lambda = 10)
res = cmaes(fn_batch, x0, lower, upper, control, batch = TRUE)

res$y
#> [1] 2.185e-16
```

## Algorithms

The variant is selected with the `algo` argument of `cmaes_control()`.
The default is `"acmaes"`, as recommended by the [libcmaes practical
hints](https://github.com/CMA-ES/libcmaes/wiki/Practical-hints). For
multimodal problems, `"ipop"` or `"bipop"` are usually the better
choice.

``` r
cmaes_algos
#>  [1] "cmaes"      "ipop"       "bipop"      "acmaes"     "aipop"      "abipop"     "sepcmaes"  
#>  [8] "sepipop"    "sepbipop"   "sepacmaes"  "sepaipop"   "sepabipop"  "vdcma"      "vdipopcma" 
#> [15] "vdbipopcma"
```
