eval_log = NULL

# helper to create a logged sphere objective with lambda-aware argument checks
make_logged_sphere_single = function(dim, lambda, lower = rep(-1, dim), upper = rep(1, dim), x0 = rep(0.5, dim)) {
  eval_log <<- data.frame(matrix(NA, nrow = 0, ncol = dim + 1))
  names(eval_log) <<- c(sprintf("x%d", seq_len(dim)), "y")

  fn = function(x) {
    assert_numeric(x, len = dim)
    y = sum(x^2)
    row = as.data.frame(as.list(x))
    names(row) = sprintf("x%d", seq_len(dim))
    row$y = y
    eval_log <<- rbind(eval_log, row) # nolint
    y
  }
  list(fn = fn, x0 = x0, lower = lower, upper = upper)
}


test_that("single: finds minimum of sphere function", {
  seed = 123
  # make sure to test the special case dim=1
  for (dim in c(1, 3)) {
    fevals = 100 * dim
    for (lambda in c(3, NA)) {
      for (algo in cmaes_algos) {
        obj = make_logged_sphere_single(dim, lambda)
        ctrl = cmaes_control(algo = algo, max_fevals = fevals, lambda = lambda, seed = seed, max_restarts = 2L)
        res = cmaes(obj$fn, obj$x0, obj$lower, obj$upper, ctrl, batch = FALSE)

        ctx = sprintf("dim=%d, lambda=%s, algo=%s", dim, ifelse(is.na(lambda), "NA", lambda), algo)
        # basic structure checks
        expect_list(res)
        expect_named(res, c("x", "y", "edm", "time", "status_code", "status_msg"), ignore.order = TRUE)
        expect_int(res$status_code, lower = 0)
        expect_string(res$status_msg)
        # solution quality (should get very close to the optimum)
        expect_true(all(abs(res$x) < 1e-2), info = ctx)
        expect_lt(res$y, 1e-3)
        ee = eval_log
        # the log must have at least one evaluation and the right shape
        expect_gt(nrow(ee), 0)
        expect_equal(ncol(ee), dim + 1) # dim cols for x + 1 for y
        expect_true(nrow(ee) <= fevals + 10, info = ctx)

        # all x-evaluations must be within [lower, upper]
        ee$y = NULL
        for (k in seq_len(dim)) {
          expect_true(all(ee[, k] >= obj$lower[k]))
          expect_true(all(ee[, k] <= obj$upper[k]))
        }
      }
    }
  }
})

test_that("single: works with exception in objective", {
  dim = 2
  fn = function(x) stop("foo")

  x0 = rep(0.5, dim)
  lower = rep(-1, dim)
  upper = rep(1, dim)
  ctrl = cmaes_control(max_fevals = 10)
  msg = capture.output(type = "message", {
    x = capture.output(type = "output", {
      expect_error(
        {
          cmaes(fn, x0, lower, upper, ctrl, batch = FALSE)
        },
        regexp = "libcmaesr: objective evaluation failed!"
      )
    })
  })
  expect_string(msg, pattern = "foo")
})

test_that("single: maximize works", {
  dim = 2
  obj = make_logged_sphere_single(dim, lambda = NA)
  ctrl = cmaes_control(max_fevals = 500, maximize = TRUE)
  res = cmaes(obj$fn, obj$x0, obj$lower, obj$upper, ctrl, batch = FALSE)
  expect_int(res$status_code, lower = 0)
  expect_string(res$status_msg)
  expect_true(all(abs(res$x) > 0.999))
  expect_gt(res$y, 1.9999)
})

test_that("single: ftarget stops early and meets target", {
  set.seed(1)
  dim = 2
  lambda = 6
  fevals = 2000
  ftarget = 1

  obj = make_logged_sphere_single(dim, lambda)
  ctrl = cmaes_control(algo = "cmaes", max_fevals = fevals, lambda = lambda, seed = 42, ftarget = ftarget)
  res = cmaes(obj$fn, obj$x0, obj$lower, obj$upper, ctrl, batch = FALSE)

  expect_true(res$y <= ftarget + 1e-8)
  expect_equal(res$status_code, 10) # ftarget status
  expect_gt(nrow(eval_log), 0)
  expect_lte(nrow(eval_log), 100) # should stop before budget is exhausted
})

test_that("single: max_iter works", {
  dim = 2
  lambda = 4
  max_iter = 3
  fevals = 10000

  obj = make_logged_sphere_single(dim, lambda)
  ctrl = cmaes_control(algo = "cmaes", max_iter = max_iter, max_fevals = fevals, lambda = lambda, seed = 7)
  res = cmaes(obj$fn, obj$x0, obj$lower, obj$upper, ctrl, batch = FALSE)

  # number of evaluations should be on the order of iterations * lambda (plus at most one extra batch)
  expect_lte(nrow(eval_log), lambda * (max_iter + 1))
  expect_int(res$status_code, lower = 0)
  expect_string(res$status_msg)
})

test_that("single: seed reproducibility", {
  dim = 3
  lambda = 5
  fevals = 40
  obj = make_logged_sphere_single(dim, lambda)

  ctrl1 = cmaes_control(algo = "cmaes", max_fevals = fevals, lambda = lambda, seed = 99)
  ctrl2 = cmaes_control(algo = "cmaes", max_fevals = fevals, lambda = lambda, seed = 99)
  ctrl3 = cmaes_control(algo = "cmaes", max_fevals = fevals, lambda = lambda, seed = 100)

  res1 = cmaes(obj$fn, obj$x0, obj$lower, obj$upper, ctrl1, batch = FALSE)
  res2 = cmaes(obj$fn, obj$x0, obj$lower, obj$upper, ctrl2, batch = FALSE)
  res3 = cmaes(obj$fn, obj$x0, obj$lower, obj$upper, ctrl3, batch = FALSE)

  expect_equal(res1$y, res2$y, tolerance = 1e-12)
  expect_equal(res1$x, res2$x, tolerance = 1e-12)
  # different seed should typically differ in either x or y
  expect_false(isTRUE(all.equal(res1$y, res3$y, tolerance = 1e-12)))
  expect_false(isTRUE(all.equal(res1$x, res3$x, tolerance = 1e-12)))
})

test_that("single: elitism and tpa options run and return valid structure", {
  dim = 2
  for (elitism in c(0L, 1L, 2L, 3L)) {
    for (tpa in c(0L, 1L, 2L)) {
      obj = make_logged_sphere_single(dim, lambda = 4)
      ctrl = cmaes_control(
        algo = "cmaes",
        max_fevals = 20,
        lambda = 4,
        elitism = elitism,
        tpa = tpa,
        tpa_dsigma = if (tpa == 0L) NA_real_ else 0.1,
        seed = 5
      )
      res = cmaes(obj$fn, obj$x0, obj$lower, obj$upper, ctrl, batch = FALSE)
      expect_named(res, c("x", "y", "edm", "time", "status_code", "status_msg"), ignore.order = TRUE)
      expect_numeric(res$edm, lower = 0, any.missing = FALSE, len = 1)
      expect_numeric(res$time, lower = 0, any.missing = FALSE, len = 1)
      expect_integer(res$status_code, lower = 0, len = 1)
      expect_string(res$status_msg)
      expect_numeric(res$x, any.missing = FALSE, len = dim)
      expect_number(res$y)
    }
  }
})

test_that("ipop/bipop non-batch logs show multiple restarts when budget allows", {
  library(callr)
  for (algo in c("ipop", "bipop")) {
    p = r_bg(
      function(algo) {
        library(libcmaesr)
        dim = 2
        fn = function(x) sum(x^2)
        x0 = rep(0.5, dim)
        lower = rep(-1, dim)
        upper = rep(1, dim)
        ctrl = cmaes_control(
          algo = algo,
          max_fevals = 300,
          max_iter = 10,
          lambda = 2,
          max_restarts = 3,
          seed = 123,
          quiet = FALSE
        )
        cmaes(fn, x0, lower, upper, ctrl, batch = FALSE)
      },
      args = list(algo = algo)
    )
    out = p$read_all_output()
    p$wait()
    res = p$get_result()
    expect_named(res, c("x", "y", "edm", "time", "status_code", "status_msg"), ignore.order = TRUE)
    expect_integer(res$status_code, lower = 0, len = 1)
    expect_string(res$status_msg)
    n_restart = length(gregexpr("restart", out, ignore.case = TRUE)[[1]])
    expect_gte(n_restart, 2) # should indicate multiple restarts in logs
  }
})

test_that("single: x0 randomization accepts vector bounds and runs", {
  dim = 2
  obj = make_logged_sphere_single(dim, lambda = NA)
  # use very small interval for x0 and nearly sigma=0
  # so all initial points should be in this interval
  x0_lower = c(0.2, -0.2)
  x0_upper = c(0.3, -0.1)
  lambda = 10
  ctrl = cmaes_control(
    max_fevals = 10L,
    seed = 999,
    sigma = 1e-6,
    lambda = lambda,
    x0_lower = x0_lower,
    x0_upper = x0_upper
  )
  res = cmaes(obj$fn, obj$x0, obj$lower, obj$upper, ctrl, batch = FALSE)
  ee = eval_log
  expect_true(all(ee$x1 >= x0_lower[1] & ee$x1 <= x0_upper[1]))
})

test_that("single: objective must return numeric scalar (type check)", {
  dim = 2
  fn = function(x) "oops" # wrong type
  x0 = rep(0.5, dim)
  lower = rep(-1, dim)
  upper = rep(1, dim)
  ctrl = cmaes_control(max_fevals = 5)
  expect_error(
    cmaes(fn, x0, lower, upper, ctrl, batch = FALSE),
    regexp = "objective must return.*numeric"
  )
})

test_that("single: objective must return length 1 (length check)", {
  dim = 2
  fn = function(x) c(1, 2) # wrong length
  x0 = rep(0.5, dim)
  lower = rep(-1, dim)
  upper = rep(1, dim)
  ctrl = cmaes_control(max_fevals = 5)
  expect_error(
    cmaes(fn, x0, lower, upper, ctrl, batch = FALSE),
    regexp = "length 1.*got 2"
  )
})
