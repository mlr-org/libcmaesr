# libcmaesr

A lightweight R interface to the [libcmaes](https://github.com/CMA-ES/libcmaes) C++ library for Covariance Matrix Adaptation Evolution Strategy (CMA-ES). It allows for the optimization of black-box functions using the CMA-ES algorithm and its variants.

## Installation

You can install the development version from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("mlr-org/libcmaesr")
```

## System requirements

This package compiles native code and requires:

- C++11 compiler toolchain
- CMake (>= 3.15)
- Eigen3 (header-only)

Install hints:

- Ubuntu:
```bash
sudo apt-get update && sudo apt-get install -y cmake libeigen3-dev
```

- macOS (Homebrew):
```bash
brew install cmake eigen
```

- Windows:
  - Install Rtools (UCRT) matching your R version.
  - Optionally install Eigen and CMake in MSYS2 for local builds:
```bash
pacman -Sy --noconfirm
pacman -S --noconfirm mingw-w64-ucrt-x86_64-eigen3 mingw-w64-ucrt-x86_64-cmake
```

For the authoritative and up-to-date list of CI system dependencies, see the workflow: [`.github/workflows/R-CMD-check.yaml`](.github/workflows/R-CMD-check.yaml).

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
