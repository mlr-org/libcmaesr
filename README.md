# libcmaesr

A lightweight R interface to the [libcmaes](https://github.com/CMA-ES/libcmaes) C++ library for Covariance Matrix Adaptation Evolution Strategy (CMA-ES). It allows for the optimization of black-box functions using the CMA-ES algorithm and its variants.

## Installation

You can install the development version from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("berndbischl/libcmaesr")
```

## Example

This is a basic example which shows you how to solve a common test problem, the sphere function:

```r
library(libcmaesr)

# define objective function
sphere = function(x) {
  sum(x^2)
}

# define search space
lower = rep(-5, 2)
upper = rep(5, 2)
x0 = rep(1, 2)

# run optimization
res = cmaes(sphere, x0, lower, upper, cmaes_control())
print(res)
```
