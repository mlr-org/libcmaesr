# libcmaesr

A lightweight R interface to the [libcmaes](https://github.com/CMA-ES/libcmaes) C++ library for Covariance Matrix Adaptation Evolution Strategy (CMA-ES). It allows for the optimization of black-box functions using the CMA-ES algorithm and its variants.

## Installation

You can install the development version from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("mlr-org/libcmaesr")
```

## Example

This is a basic example which shows you how to solve a common test problem, the sphere function:

```r
library(libcmaesr)

# define objective function
dim = 3
fn = function(x) {
  apply(x, 1, function(row) sum(row^2))
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
