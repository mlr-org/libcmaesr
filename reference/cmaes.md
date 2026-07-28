# Covariance Matrix Adaptation Evolution Strategy

Implements the CMA-ES variants provided by libcmaes, see here:
<https://github.com/CMA-ES/libcmaes/> via a very light-weight C wrapper.

The control structure allows access to most control params of the ES,
but CMAES is supposed to handle most of them internally. Quoting Niko
Hansen from here: <https://cma-es.github.io/>:

“The CMA-ES does not require a tedious parameter tuning for its
application. In fact, the choice of strategy internal parameters is not
left to the user (arguably with the exception of population size
\\\lambda\\). Finding good (default) strategy parameters is considered
as part of the algorithm design, and not part of its application — the
aim is to have a well-performing algorithm as is. The default population
size \\\lambda\\ is comparatively small to allow for fast convergence.
Restarts with increasing population size (Auger & Hansen 2005) improve
the global search performance. For the application of the CMA-ES, an
initial solution, an initial standard deviation (step-size, variables
should be defined such that the same standard deviations can be
reasonably applied to all variables, see also here) and, possibly, the
termination criteria (e.g. a function tolerance) need to be set by the
user. The most common applications are model calibration (e.g. curve
fitting) and shape optimisation.”

Whether you believe in this completely for any problem is up to you, but
the general idea is to run it in its defaults, and only change them if
you know what you are doing.

1.  libcmaes could handle unbounded search spaces, but this is currently
    not supported, you need to set lower and upper bounds.

2.  Noisy functions / noisy handling CMAES is currently not supported.

3.  Surrogate variants are currently not supported.

4.  Geno-Pheno transformation is automatically applied, in the sense
    that we use the linear scaling to handle the bounds.

5.  Setting gradients is currently not supported.

6.  OpenMP is currently not supported. libcmaes uses OpenMP for the
    population evaluation mainly, but also for some Eigen stuff. Calling
    into the R Api via threading is not allowed, which would happen in
    the former.

In general, more details can be found here:
<https://github.com/CMA-ES/libcmaes/wiki/>.

## Usage

``` r
cmaes(objective, x0, lower, upper, control = cmaes_control(), batch = FALSE)
```

## Arguments

- objective:

  (`function(x)`)  
  Objective function, to minimize. If `batch` is `FALSE`, `x` is a
  numeric vector of length `n` and the function must return a scalar
  `numeric`. If `batch` is `TRUE`, `x` is a numeric matrix with one row
  per candidate and `n` columns. The number of rows can vary between
  iterations, e.g., restart strategies like IPOP and BIPOP increase the
  population size over time, so do not assume it is always `lambda`. The
  function must return a numeric vector with one element per row of `x`.
  The latter usually reduces overhead and allows you to orchestrate
  parallelization yourself if you need it because the objective function
  is more expensive.

- x0:

  (`numeric(n)`)  
  Initial point. NB: This point is IGNORED if you also set `x0_lower`
  and `x0_upper` in the `control` object, but it must still be a vector
  of length `n`!

- lower:

  (`numeric(n)`)  
  Lower bounds of search space.

- upper:

  (`numeric(n)`)  
  Upper bounds of search space.

- control:

  (`cmaes_control`)  
  A control object created by
  [`cmaes_control()`](https://libcmaesr.mlr-org.com/reference/cmaes_control.md).
  Default is a control object with all parameters set to their default
  values.

- batch:

  (`logical(1)`)  
  Whether the objective function evaluates a batch of points at once.
  Default is `FALSE`.

## Value

(named `list`). List with elements:

- 'x': (`numeric(n)`)  
  The best point found, length corresponds to x0, lower and upper.

- 'y': (`numeric(1)`)  
  The objective value of the best point.

- 'edm': (`numeric(1)`)  
  Expected distance to the minimum.

- 'time': (`numeric(1)`)  
  The time taken to find the solution in seconds.

- 'status_code': (`integer(1)`)  
  The status code, indicating success, failure, or the reason for
  stopping. See here:
  <https://github.com/CMA-ES/libcmaes/wiki/Optimizing-a-function>

- 'status_msg': (`character(1)`)  
  A human-readable status message from libcmaes.

## Examples

``` r
# minimize a simple quadratic function
objective = function(x) sum(x^2)
control = cmaes_control(seed = 1, max_fevals = 500)
res = cmaes(objective, x0 = c(0.5, 0.5), lower = c(-5, -5), upper = c(5, 5), control = control)
res$x
#> [1] -5.644196e-11 -7.141949e-10
res$y
#> [1] 5.132601e-19
```
