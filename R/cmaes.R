#' @title CMA-ES Control Object
#' @description Create a control object for the CMA-ES algorithm.
#' @return A cmaes_control object.
#' @export
cmaes_control = function(algo = "aCMAES", max_fevals = 10000, max_iter = 100000, ftarget = 1e-8) {
  assert_choice(algo, c("aCMAES", "BIPOP_CMAES"))
  assert_int(max_fevals, lower = 1)
  assert_int(max_iter, lower = 1)
  assert_number(ftarget, lower = 0)
  res = list(
    algo = algo,
    max_fevals = max_fevals,
    max_iter = max_iter,
    ftarget = ftarget
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