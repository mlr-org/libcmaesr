# CMA-ES Control Object

Create a control object for the CMA-ES algorithm. For more information
on the parameters, see here:
<https://cma-es.github.io/libcmaes/doc/html/classlibcmaes_1_1CMAParameters.html>.

## Usage

``` r
cmaes_control(
  maximize = FALSE,
  algo = "acmaes",
  max_fevals = 100,
  max_iter = NA_integer_,
  ftarget = NA_real_,
  f_tolerance = NA_real_,
  x_tolerance = NA_real_,
  lambda = NA_integer_,
  sigma = NA_real_,
  max_restarts = NA_integer_,
  elitism = NA_integer_,
  tpa = NA_integer_,
  tpa_dsigma = NA_real_,
  seed = NA_integer_,
  quiet = TRUE,
  x0_lower = NULL,
  x0_upper = NULL
)
```

## Arguments

- maximize:

  (`logical(1)`)  
  Whether to maximize the objective function. Default is `FALSE`.

- algo:

  (`character(1)`)  
  The CMAES variant to use. Possible values are:
  [cmaes_algos](https://libcmaesr.mlr-org.com/reference/cmaes_algos.md).
  Default is "acmaes", as recommended by
  <https://github.com/CMA-ES/libcmaes/wiki/Practical-hints>. For
  multimodal problems, you likely want to use "ipop" or "bipop".

- max_fevals:

  (`integer(1)`)  
  The maximum number of function evaluations. NA to disable. Default is
  100.

- max_iter:

  (`integer(1)`)  
  The maximum number of ES iterations. NA to disable (default).

- ftarget:

  (`numeric(1)`)  
  Stop when this target function value is reached. NA to disable
  (default).

- f_tolerance:

  (`numeric(1)`)  
  Sets function tolerance as stopping criterion; monitors the (absolute)
  difference in function value over iterations and stops optimization
  when below tolerance. NA to disable (default).

- x_tolerance:

  (`numeric(1)`)  
  Sets parameter (absolute) tolerance as stopping criterion. This checks
  entries of the covariance matrix, only touch when you know what you
  are doing. NA to disable (default).

- lambda:

  (`integer(1)`)  
  Number of generated descendants per iteration. Must be at least 2; NA
  for default handling by libcmaes.

- sigma:

  (`numeric(1)`)  
  Initial sigma for covariance. NA for default handling by libcmaes.

- max_restarts:

  (`integer(1)`)  
  The maximum number of restarts, for IPOP and BIPOP. NA for default
  handling by libcmaes.

- elitism:

  (`integer(1)`)  
  Sets elitism:  

  0

  :   no elitism

  1

  :   elitism: reinjects the best-ever seen solution

  2

  :   initial elitism: reinject x0 as long as it is not improved upon

  3

  :   initial elitism on restart: restart if best encountered solution
      is not the the final solution and reinjects the best solution
      until the population has better fitness, in its majority

  NA for default handling by libcmaes.

- tpa:

  (`integer(1)`)  
  Activates / deactivates two-point adaptation step-size mechanism. 0:
  no, 1: auto, 2: yes. NA for default handling by libcmaes.

- tpa_dsigma:

  (`numeric(1)`)  
  Sets tpa dsigma value, use with care. NA for default handling by
  libcmaes.

- seed:

  (`integer(1)`)  
  The seed for the random number generator. If `NA` (default), the seed
  is generated randomly by R and thereby coupled to the RNG-state of R.
  Otherwise, the RNG of the libcmaes is different to the one in R and is
  hence not subject to R's seeding. Special value `0` is used for
  handling by libcmaes, where system time is used in libcmaes to seed.

- quiet:

  (`logical(1)`)  
  Whether to suppress libcmaes output. Internal logging of libcmaes is
  rerouted to Rprintf, so things like capture.output() will work. Useful
  for debugging. Default is `TRUE`.

- x0_lower:

  (`numeric`)  
  Optional lower bounds for randomizing the initial mean `x0`, also
  after restarts. Use `NULL` to disable. If this is non-`NULL`,
  `x0_upper` must also be set and have the same length as `x0_lower`.

- x0_upper:

  (`numeric`)  
  Optional upper bounds for randomizing the initial mean `x0`, also
  after restarts. Use `NULL` to disable.

## Value

A cmaes_control S3 object, which is a list with the passed arguments.

## Examples

``` r
control = cmaes_control(algo = "bipop", max_fevals = 1000, seed = 42)
print(control)
#> CMA-ES control object:
#> cmaes_control(algo = "bipop", max_fevals = 1000L, seed = 42L)
```
