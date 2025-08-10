f_sphere = function(x) {
  apply(x, 1, function(row) sum(row^2))
}


# helper to create a logged sphere objective with lambda-aware argument checks
make_logged_sphere_objective = function(dim, lambda,
  lower = rep(-1, dim), upper = rep(1, dim), x0 = rep(0.5, dim)) {

  fn = function(x) {
    if (!is.na(lambda)) {
      assert_matrix(x, nrows = lambda, ncols = dim)
    } else {
      assert_matrix(x, ncols = dim)
    }
    ys = apply(x, 1, function(row) sum(row^2))
    df = as.data.frame(x)
    names(df) = sprintf("x%d", seq_len(dim))
    df$y = ys
    eval_log <<- rbind(eval_log, df) # nolint
    ys
  }
  list(fn = fn, x0 = x0, lower = lower, upper = upper)
}


test_that("cmaes finds minimum of sphere function", {
  seed = 123
  # make sure to test the special case dim=1
  for (dim in 1:3) {
    fevals = 200 * dim
    for (lambda in c(3, 5, NA)) {
      for (algo in cmaes_algos) {
        # global log for all objective evaluations in this iteration
        eval_log <<- data.frame()

        obj = make_logged_sphere_objective(dim, lambda)
        ctrl = cmaes_control(algo = algo, max_fevals = fevals, lambda = lambda, seed = seed)
        res = cmaes(obj$fn, obj$x0, obj$lower, obj$upper, ctrl)

        # basic structure checks
        expect_list(res)
        expect_named(res, c("x", "y", "edm", "time", "status"), ignore.order = TRUE)
        expect_int(res$status, lower = 0)
        expect_true(all(abs(res$x) < 1e-3))
        # solution quality (should get very close to the optimum)
        expect_lt(res$y, 1e-6)

        # the log must have at least one evaluation and the right shape
        expect_gt(nrow(eval_log), 0)
        expect_equal(ncol(eval_log), dim + 1)  # dim cols for x + 1 for y
        if (is.na(lambda)) {
          expect_lte(nrow(eval_log), fevals + 10)
        } else {
          expect_lte(nrow(eval_log), fevals + lambda) # last pop could be over fevals
        }

        eval_x = eval_log
        eval_x$y = NULL

        # all x-evaluations must be within [lower, upper]
        for (k in seq_len(dim)) {
          expect_true(all(eval_x[, k] >= obj$lower[k]))
          expect_true(all(eval_x[, k] <= obj$upper[k]))
        }
      }
    }
  }
})


test_that("cmaes works with exception in objective", {
  dim = 2
  fn = function(x) stop("foo")

  x0 = rep(0.5, dim)
  lower = rep(-1, dim)
  upper = rep(1, dim)
  ctrl = cmaes_control(max_fevals = 10)
  msg = capture.output(type = "message", {
    x = capture.output(type = "output", {
      expect_error({
        cmaes(fn, x0, lower, upper, ctrl)
      }, regexp = "libcmaesr: objective evaluation failed!")
    })
  })
  expect_string(msg, pattern = "foo")
})


test_that("maximize works", {
  dim = 2
  eval_log <<- data.frame()
  obj = make_logged_sphere_objective(dim, lambda = NA)
  ctrl = cmaes_control(max_fevals = 500, maximize = TRUE)
  res = cmaes(obj$fn, obj$x0, obj$lower, obj$upper, ctrl)
  expect_int(res$status, lower = 0)
  expect_true(all(abs(res$x) > 0.999))
  expect_gt(res$y, 1.9999)
})


test_that("ftarget stops early and meets target (cmaes)", {
  set.seed(1)
  dim = 2
  lambda = 6
  fevals = 2000
  ftarget = 1e-3

  eval_log <<- data.frame()
  obj = make_logged_sphere_objective(dim, lambda)
  ctrl = cmaes_control(algo = "cmaes", max_fevals = fevals, lambda = lambda, seed = 42,
    ftarget = ftarget)
  res = cmaes(obj$fn, obj$x0, obj$lower, obj$upper, ctrl)

  expect_true(res$y <= ftarget + 1e-8)
  expect_gt(nrow(eval_log), 0)
  expect_lte(nrow(eval_log), fevals) # should stop before budget is exhausted
})


test_that("max_iter works", {
  dim = 2
  lambda = 4
  max_iter = 3
  fevals = 10000

  eval_log <<- data.frame()
  obj = make_logged_sphere_objective(dim, lambda)
  ctrl = cmaes_control(algo = "cmaes", max_iter = max_iter, max_fevals = fevals, lambda = lambda, seed = 7)
  res = cmaes(obj$fn, obj$x0, obj$lower, obj$upper, ctrl)

  # number of evaluations should be on the order of iterations * lambda (plus at most one extra batch)
  expect_lte(nrow(eval_log), lambda * (max_iter + 1))
  expect_int(res$status, lower = 0)
})


test_that("seed reproducibility", {
  dim = 3
  lambda = 5
  fevals = 40
  obj = make_logged_sphere_objective(dim, lambda)

  ctrl1 = cmaes_control(algo = "cmaes", max_fevals = fevals, lambda = lambda, seed = 99)
  ctrl2 = cmaes_control(algo = "cmaes", max_fevals = fevals, lambda = lambda, seed = 99)
  ctrl3 = cmaes_control(algo = "cmaes", max_fevals = fevals, lambda = lambda, seed = 100)

  res1 = cmaes(f_sphere, x0, lower, upper, ctrl1)
  res2 = cmaes(f_sphere, x0, lower, upper, ctrl2)
  res3 = cmaes(f_sphere, x0, lower, upper, ctrl3)

  expect_equal(res1$y, res2$y, tolerance = 1e-12)
  expect_equal(res1$x, res2$x, tolerance = 1e-12)
  # different seed should typically differ in either x or y
  expect_false(isTRUE(all.equal(res1$y, res3$y, tolerance = 1e-12)))
  expect_false(isTRUE(all.equal(res1$x, res3$x, tolerance = 1e-12)))
})


test_that("elitism and tpa options run and return valid structure", {
  dim = 2
  x0 = rep(0.5, dim)
  lower = rep(-1, dim)
  upper = rep(1, dim)

  for (elitism in c(0L, 1L, 2L, 3L)) {
    for (tpa in c(0L, 1L, 2L)) {
      ctrl = cmaes_control(algo = "cmaes", max_fevals = 20, lambda = 4,
        elitism = elitism, tpa = tpa, tpa_dsigma = if (tpa == 0L) NA_real_ else 0.1,
        seed = 5)
      res = cmaes(f_sphere, x0, lower, upper, ctrl)
      expect_named(res, c("x", "y", "edm", "time", "status"), ignore.order = TRUE)
      expect_numeric(res$edm, lower = 0, any.missing = FALSE, len = 1)
      expect_numeric(res$time, lower = 0, any.missing = FALSE, len = 1)
      expect_integer(res$status, lower = 0, len = 1)
      expect_numeric(res$x, any.missing = FALSE, len = dim)
      expect_number(res$y)
    }
  }
})
