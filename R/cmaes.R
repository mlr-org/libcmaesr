#' @title CMA-ES Control Object
#'
#' @description Create a control object for the CMA-ES algorithm.
#' For more information on the parameters, see here:
#' \url{https://cma-es.github.io/libcmaes/doc/html/classlibcmaes_1_1CMAParameters.html}.
#'
#' @param maximize (`logical(1)`)\cr
#'   Whether to maximize the objective function.
#'   Default is `FALSE`.
#' @param algo (`character(1)`)\cr
#'   The CMAES variant to use.
#'   Possible values are: [cmaes_algos].
#'   Default is "acmaes", as recommended by
#'   \url{https://github.com/CMA-ES/libcmaes/wiki/Practical-hints}.
#'   For multimodal problems, you likely want to use "ipop" or "bipop".
#' @param max_fevals (`integer(1)`)\cr
#'   The maximum number of function evaluations.
#'   NA to disable.
#'   Default is 100.
#' @param max_iter (`integer(1)`)\cr
#'   The maximum number of ES iterations.
#'   NA to disable (default).
#' @param ftarget (`numeric(1)`)\cr
#'   Stop when this target function value is reached.
#'   NA to disable (default).
#' @param f_tolerance (`numeric(1)`)\cr
#'   Sets function tolerance as stopping criterion;
#'   monitors the (absolute) difference in function value over iterations and
#'   stops optimization when below tolerance.
#'   NA to disable (default).
#' @param x_tolerance (`numeric(1)`)\cr
#'   Sets parameter (absolute) tolerance as stopping criterion.
#'   This checks entries of the covariance matrix, only touch when you know what you are doing.
#'   NA to disable (default).
#' @param lambda (`integer(1)`)\cr
#'   Number of generated descendants per iteration.
#'   Must be at least 2; NA for default handling by libcmaes.
#' @param sigma (`numeric(1)`)\cr
#'   Initial sigma for covariance.
#'   NA for default handling by libcmaes.
#' @param max_restarts (`integer(1)`)\cr
#'   The maximum number of restarts, for IPOP and BIPOP.
#'   NA for default handling by libcmaes.
#' @param elitism (`integer(1)`)\cr
#'   Sets elitism:\cr
#'   \describe{
#'   \item{0}{no elitism}
#'   \item{1}{elitism: reinjects the best-ever seen solution}
#'   \item{2}{initial elitism: reinject x0 as long as it is not improved upon}
#'   \item{3}{initial elitism on restart: restart if best encountered solution is not
#'     the the final solution and reinjects the best solution until
#'     the population has better fitness, in its majority}
#'   }
#'   NA for default handling by libcmaes.
#' @param tpa (`integer(1)`)\cr
#'   Activates / deactivates two-point adaptation step-size mechanism.
#'   0: no, 1: auto, 2: yes.
#'   NA for default handling by libcmaes.
#' @param tpa_dsigma (`numeric(1)`)\cr
#'   Sets tpa dsigma value, use with care.
#'   NA for default handling by libcmaes.
#' @param seed (`integer(1)`)\cr
#'   The seed for the random number generator.
#'   If `NA` (default), the seed is generated randomly by R and thereby coupled to the RNG-state of R.
#'   Otherwise, the RNG of the libcmaes is different to the one in R and is hence not subject to R's seeding.
#'   Special value `0` is used for handling by libcmaes, where system time is used in libcmaes to seed.
#' @param quiet (`logical(1)`)\cr
#'   Whether to suppress libcmaes output.
#'   Internal logging of libcmaes is rerouted to Rprintf, so things like
#'   capture.output() will work. Useful for debugging.
#'   Default is `TRUE`.
#' @param x0_lower (`numeric`)\cr
#'   Optional lower bounds for randomizing the initial mean `x0`,
#'   also after restarts.
#'   Use `NULL` to disable.
#'   If this is non-`NULL`, `x0_upper` must also be set and have the same length as `x0_lower`.
#' @param x0_upper (`numeric`)\cr
#'   Optional upper bounds for randomizing the initial mean `x0`,
#'   also after restarts.
#'   Use `NULL` to disable.
#' @return A cmaes_control S3 object, which is a list with the passed arguments.
#'
#' @export
cmaes_control = function(
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
) {
  assert_flag(maximize)
  assert_choice(algo, cmaes_algos)
  max_fevals = asInt(max_fevals, lower = 1, na.ok = TRUE)
  max_iter = asInt(max_iter, lower = 1, na.ok = TRUE)
  assert_number(ftarget, na.ok = TRUE)
  assert_number(f_tolerance, lower = 0, na.ok = TRUE)
  assert_number(x_tolerance, lower = 0, na.ok = TRUE)
  lambda = asInt(lambda, lower = 2, na.ok = TRUE)
  assert_number(sigma, lower = 0, na.ok = TRUE)
  max_restarts = asInt(max_restarts, lower = 0, na.ok = TRUE)
  elitism = asInt(elitism, lower = 0, upper = 3, na.ok = TRUE)
  tpa = asInt(tpa, lower = 0, upper = 2, na.ok = TRUE)
  assert_number(tpa_dsigma, lower = 0, na.ok = TRUE)
  seed = asInt(seed, na.ok = TRUE)
  assert_flag(quiet)
  if (is.null(x0_lower) != is.null(x0_upper)) {
    stop("'x0_lower' and 'x0_upper' must both be NULL or not NULL!")
  }
  assert_numeric(x0_lower, min.len = 1, any.missing = FALSE, finite = TRUE, null.ok = TRUE)
  assert_numeric(x0_upper, len = length(x0_lower), any.missing = FALSE, finite = TRUE, null.ok = TRUE)
  if (!is.null(x0_lower) && !is.null(x0_upper) && !all(x0_lower < x0_upper)) {
    stop("'x0_lower' must be strictly smaller than 'x0_upper'!")
  }

  res = list(
    maximize = maximize,
    algo = algo,
    max_fevals = max_fevals,
    max_iter = max_iter,
    ftarget = ftarget,
    f_tolerance = f_tolerance,
    x_tolerance = x_tolerance,
    lambda = lambda,
    sigma = sigma,
    max_restarts = max_restarts,
    elitism = elitism,
    tpa = tpa,
    tpa_dsigma = tpa_dsigma,
    seed = seed,
    quiet = quiet,
    x0_lower = x0_lower,
    x0_upper = x0_upper
  )
  set_class(res, "cmaes_control")
}

#' @export
print.cmaes_control = function(x, ...) {
  cc_not_default = Filter(function(x) !is.null(x) && (length(x) > 1 || !is.na(x)), x)
  cat("CMA-ES control object:\n")
  print(as.call(c(alist(cmaes_control), cc_not_default)), ...)
  invisible(x)
}

#' @title Covariance Matrix Adaptation Evolution Strategy
#' @description
#' Implements the CMA-ES variants provided by libcmaes, see here: \url{https://github.com/CMA-ES/libcmaes/} via
#' a very light-weight C wrapper.
#'
#' 2.The control structure allows access to most control params of the ES, but CMAES is supposed to handle most of them internally.
#' Quoting Niko Hansen from here: \url{https://cma-es.github.io/}:
#'
#' \dQuote{The CMA-ES does not require a tedious parameter tuning for its application. In fact, the choice of strategy internal parameters
#' is not left to the user (arguably with the exception of population size λ). Finding good (default) strategy parameters is considered
#' as part of the algorithm design, and not part of its application—the aim is to have a well-performing algorithm as is.
#' The default population size λ is comparatively small to allow for fast convergence. Restarts with increasing population size
#' (Auger & Hansen 2005) improve the global search performance. For the application of the CMA-ES, an initial solution,
#' an initial standard deviation (step-size, variables should be defined such that the same standard deviations can be
#' reasonably applied to all variables, see also here) and, possibly, the termination criteria (e.g. a function tolerance)
#' need to be set by the user. The most common applications are model calibration (e.g. curve fitting) and shape optimisation.}
#'
#' Whether you believe in this completely for any problem is up to you, but the general idea is to run it in its defaults,
#' and only change them if you know what you are doing.
#'
#' 1. libcmaes could handle unbounded search spaces, but this is currently not supported, you need to set
#' lower and upper bounds.
#'
#' 2. Noisy functions / noisy handling CMAES is currently not supported.
#'
#' 3. Surrogate variants are currently not supported.
#'
#' 4. Geno-Pheno transformation is automatically applied, in the sense that we use the linear scaling
#' to handle the bounds.
#'
#' 5. Setting gradients is currently not supported.
#'
#' 6. OpenMP is currently not supported. libcmaes uses OpenMP for the population evaluation mainly, but also for some Eigen stuff.
#' Calling into the R Api via threading is not allowed, which would happen in the former.
#'
#' In general, more details can be found here: \url{https://github.com/CMA-ES/libcmaes/wiki/}.
#'
#' @param objective (`function(x)`)\cr
#'   Objective function, to minimize.
#'   If `batch` is `FALSE`, `x` is a numeric vector of length `n`
#'   and the function must return a scalar `numeric`.
#'   If `batch` is `TRUE`, `x` is a numeric matrix with `lambda` rows and `n` columns.
#'   The function must return a numeric vector of length `lambda`.
#'   The latter usually reduces overhead and allows you to orchestrate parallelization
#'   yourself if you need it because the objective function is more expensive.
#' @param x0 (`numeric(n)`)\cr
#'   Initial point.
#'   NB: This point is IGNORED if you also set `x0_lower` and `x0_upper` in the `control` object,
#'   but it must still be a vector of length `n`!
#' @param lower (`numeric(n)`)\cr
#'   Lower bounds of search space.
#' @param upper (`numeric(n)`)\cr
#'   Upper bounds of search space.
#' @param control (`cmaes_control`)\cr
#'   A control object created by [cmaes_control()].
#'   Default is a control object with all parameters set to their default values.
#' @param batch (`logical(1)`)\cr
#'   Whether the objective function evaluates a batch of points at once.
#'   Default is `FALSE`.
#' @return (name `list`). List with elements:
#'   - 'x': (`numeric(n)`)\cr
#'     The best point found, length corresponds to x0, lower and upper.
#'   - 'y': (`numeric(1)`)\cr
#'     The objective value of the best point.
#'   - 'edm': (`numeric(1)`)\cr
#'     Expected distance to the minimum.
#'   - 'time': (`numeric(1)`)\cr
#'     The time taken to find the solution in seconds.
#'   - 'status_code': (`integer(1)`)\cr
#'     The status code, indicating success, failure, or the reason for stopping.
#'   - 'status_msg': (`character(1)`)\cr
#'     A human-readable status message from libcmaes.
#'     See here: \url{https://github.com/CMA-ES/libcmaes/wiki/Optimizing-a-function}
#' @useDynLib libcmaesr, .registration = TRUE
#' @export
cmaes = function(objective, x0, lower, upper, control = cmaes_control(), batch = FALSE) {
  assert_function(objective)
  assert_numeric(x0, min.len = 1, any.missing = FALSE, finite = TRUE)
  assert_numeric(lower, min.len = 1, any.missing = FALSE, finite = TRUE)
  assert_numeric(upper, min.len = 1, any.missing = FALSE, finite = TRUE)
  assert_class(control, "cmaes_control")
  assert_flag(batch)

  if (length(x0) != length(lower) || length(x0) != length(upper)) {
    stop("'x0', 'lower' and 'upper' must have the same length!")
  }
  if (!all(lower < upper)) {
    stop("'lower' must be strictly smaller than 'upper'!")
  }
  if (!all(x0 > lower) || !all(x0 < upper)) {
    stop("'x0' must be strictly between 'lower' and 'upper'!")
  }

  if (!is.null(control$x0_lower) && length(control$x0_lower) != length(lower)) {
    stop("'x0_lower' must have the same length as 'lower'!")
  }
  if (!is.null(control$x0_upper) && length(control$x0_upper) != length(upper)) {
    stop("'x0_upper' must have the same length as 'upper'!")
  }
  if (!all(control$x0_upper < upper)) {
    stop("'x0_upper' must be strictly smaller than 'upper'!")
  }
  if (!all(control$x0_lower > lower)) {
    stop("'x0_lower' must be strictly larger than 'lower'!")
  }

  if (is.na(control$seed)) {
    control$seed = as.integer(runif(1, 1, .Machine$integer.max))
  }

  if (batch) {
    .Call("c_cmaes_wrap_batch", objective, x0, lower, upper, control, PACKAGE = "libcmaesr")
  } else {
    .Call("c_cmaes_wrap_single", objective, x0, lower, upper, control, PACKAGE = "libcmaesr")
  }
}
