#' @title CMA-ES Control Object
#' @description Create a control object for the CMA-ES algorithm.
#' @param algo The algorithm to use.
#' @param max_fevals (`integer(1)`)\cr 
#'   The maximum number of function evaluations.
#'   NA to disable.
#' @param max_iter (`integer(1)`)\cr 
#'   The maximum number of ES iterations.
#'   NA to disable.
#' @param ftarget (`numeric(1)`)\cr 
#'   Stop when this target function value is reached.
#'   NA to disable.
#' @param lambda (`integer(1)`)\cr 
#'   Number of generated descendents per iteration.
#'   NA to for default handling by libcmaes.
#' @param sigma (`numeric(1)`)\cr 
#'   Initial sigma for covariance.
#'   NA to for default handling by libcmaes.
#' @param seed (`integer(1)`)\cr 
#'   The seed for the random number generator. If `NA`, the seed is set to 0.
#'   NB: The RNG of the libcmaes if different to the one in R and is hence not subject to R's seeding.
#'   NA to disable, time is used in libcmaes to seed.
#' @return A cmaes_control object.
#' @export
cmaes_control = function(algo = "aCMAES", max_fevals = 100, max_iter = NA_integer_, ftarget = NA_real_,
  lambda = NA_integer_, sigma = NA_real_,
  seed = NA_integer_) 
{
  assert_choice(algo, c("aCMAES", "BIPOP_CMAES"))
  max_fevals = asInt(max_fevals, lower = 1, na.ok = TRUE)
  max_iter = asInt(max_iter, lower = 1, na.ok = TRUE)
  assert_number(ftarget, lower = 0, na.ok = TRUE)
  lambda = asInt(lambda, lower = 1, na.ok = TRUE)
  assert_number(sigma, lower = 0, na.ok = TRUE)
  seed = asInt(seed, na.ok = TRUE)
  res = list(
    algo = algo,
    max_fevals = max_fevals,
    max_iter = max_iter,
    ftarget = ftarget,
    lambda = lambda,
    sigma = sigma,
    seed = seed
  )
  set_class(res, "cmaes_control")
}

#' @title Covariance Matrix Adaptation Evolution Strategy
#' @description The main function for the CMA-ES algorithm.
#'
#' @param objective The objective function to minimize.
#' @param lower The lower bounds of the search space.
#' @param upper The upper bounds of the search space.
#' @param control A control object created by `cmaes_control`.
#' @return `NULL` for now.
#' @useDynLib libcmaesr, .registration = TRUE
#' @export
cmaes = function(objective, x0, lower, upper, control) {
  assert_function(objective)
  assert_numeric(x0, min.len = 1, any.missing = FALSE, finite = TRUE)
  assert_numeric(lower, min.len = 1, any.missing = FALSE, finite = TRUE)
  assert_numeric(upper, min.len = 1, any.missing = FALSE, finite = TRUE)
  assert_numeric(x0, min.len = 1)
  if (length(lower) != length(x0) || length(upper) != length(x0)) {
    stop("'x0', 'lower', 'upper' must all have the same length!")
  }
  if (!all(lower < upper)) {
    stop("'lower' must be strictly smaller than 'upper'!")
  }
  if (!all(x0 > lower) || !all(x0 < upper)) {
    stop("'x0' must be strictly between 'lower' and 'upper'!")
  }
  assert_class(control, "cmaes_control")
  .Call("c_cmaes_wrap", objective, x0, lower, upper, control, PACKAGE = "libcmaesr")
}