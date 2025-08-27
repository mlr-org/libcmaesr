eval_log = NULL

# helper to create a logged sphere objective with lambda-aware argument checks
make_logged_sphere_batch = function(dim, algo, lambda,
  lower = rep(-1, dim), upper = rep(1, dim), x0 = rep(0.5, dim)) {

  eval_log <<- data.frame(matrix(NA, nrow = 0, ncol = dim + 1))
  names(eval_log) <<- c(sprintf("x%d", seq_len(dim)), "y")

  fn = function(x) {
    if (!is.na(lambda) && algo %nin% c("ipop", "sepipop", "aipop", "bipop", "sepbipop", "vdbipopcma", "abipop")) {
      assert_matrix(x, nrows = lambda, ncols = dim)
    } else {
      assert_matrix(x, ncols = dim)
    }
    ys = apply(x, 1, function(row) sum(row^2))
    dd = as.data.frame(x)
    names(dd) = sprintf("x%d", seq_len(dim))
    dd$y = ys
    eval_log <<- rbind(eval_log, dd) # nolint
    ys
  }
  list(fn = fn, x0 = x0, lower = lower, upper = upper)
}


test_that("cmaes finds minimum of sphere function", {
  seed = 123
  # make sure to test the special case dim=1
  for (dim in c(1, 3)) {
    fevals = 200 * dim
    for (lambda in c(3, NA)) {
      for (algo in cmaes_algos) {
        obj = make_logged_sphere_batch(dim, algo, lambda)
        ctrl = cmaes_control(algo = algo, max_fevals = fevals, lambda = lambda, seed = seed,
          max_restarts = 2L)
        res = cmaes(obj$fn, obj$x0, obj$lower, obj$upper, ctrl, batch = TRUE)

        ctx = sprintf("dim=%d, lambda=%s, algo=%s", dim, ifelse(is.na(lambda), "NA", lambda), algo)
        ee = eval_log

        # basic structure checks
        expect_list(res)
        expect_named(res, c("x", "y", "edm", "time", "status_code", "status_msg"), ignore.order = TRUE)
        expect_int(res$status_code, lower = 0)
        expect_string(res$status_msg)
        expect_true(all(abs(res$x) < 1e-3))
        # solution quality (should get very close to the optimum)
        expect_lt(res$y, 1e-6)

        # the log must have at least one evaluation and the right shape
        expect_gt(nrow(ee), 0)
        expect_equal(ncol(ee), dim + 1)  # dim cols for x + 1 for y

        # FIXME: bipop does not respect max_fevals, opened an issue
        if (algo %nin% c("bipop", "sepbipop", "vdbipopcma", "abipop", "ipop", "sepipop", "aipop")) {
          expect_true(nrow(ee) <= fevals + 20, info = ctx) # last pop could be over fevals
        }

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
        cmaes(fn, x0, lower, upper, ctrl, batch = TRUE)
      }, regexp = "libcmaesr: objective evaluation failed!")
    })
  })
  expect_string(msg, pattern = "foo")
})


test_that("maximize works", {
  dim = 2
  obj = make_logged_sphere_batch(dim, "cmaes", lambda = NA)
  ctrl = cmaes_control(max_fevals = 500, maximize = TRUE)
  res = cmaes(obj$fn, obj$x0, obj$lower, obj$upper, ctrl, batch = TRUE)
  expect_int(res$status_code, lower = 0)
  expect_string(res$status_msg)
  expect_true(all(abs(res$x) > 0.999))
  expect_gt(res$y, 1.9999)
})


test_that("ftarget stops early and meets target (cmaes)", {
  set.seed(1)
  dim = 2
  lambda = 6
  fevals = 2000
  ftarget = 1

  algo = "cmaes"
  obj = make_logged_sphere_batch(dim, algo, lambda)
  ctrl = cmaes_control(algo = algo, max_fevals = fevals, lambda = lambda, seed = 42,
    ftarget = ftarget)
  res = cmaes(obj$fn, obj$x0, obj$lower, obj$upper, ctrl, batch = TRUE)

  expect_true(res$y <= ftarget + 1e-8)
  expect_equal(res$status_code, 10) # ftarget status
  expect_string(res$status_msg)
  expect_gt(nrow(eval_log), 0)
  expect_lte(nrow(eval_log), 100) # should stop before budget is exhausted
})


test_that("max_iter works", {
  dim = 2
  lambda = 4
  max_iter = 3
  fevals = 10000

  algo = "cmaes"
  obj = make_logged_sphere_batch(dim, algo, lambda)
  ctrl = cmaes_control(algo = algo, max_iter = max_iter, max_fevals = fevals, lambda = lambda, seed = 7)
  res = cmaes(obj$fn, obj$x0, obj$lower, obj$upper, ctrl, batch = TRUE)

  # number of evaluations should be on the order of iterations * lambda (plus at most one extra batch)
  expect_lte(nrow(eval_log), lambda * (max_iter + 1))
  expect_int(res$status_code, lower = 0)
  expect_string(res$status_msg)
})


test_that("seed reproducibility", {
  dim = 3
  lambda = 5
  fevals = 40
  obj = make_logged_sphere_batch(dim = dim, algo = "cmaes", lambda = lambda)

  ctrl1 = cmaes_control(max_fevals = fevals, lambda = lambda, seed = 99)
  ctrl2 = cmaes_control(max_fevals = fevals, lambda = lambda, seed = 99)
  ctrl3 = cmaes_control(max_fevals = fevals, lambda = lambda, seed = 100)

  res1 = cmaes(obj$fn, obj$x0, obj$lower, obj$upper, ctrl1, batch = TRUE)
  res2 = cmaes(obj$fn, obj$x0, obj$lower, obj$upper, ctrl2, batch = TRUE)
  res3 = cmaes(obj$fn, obj$x0, obj$lower, obj$upper, ctrl3, batch = TRUE)

  expect_equal(res1$y, res2$y, tolerance = 1e-12)
  expect_equal(res1$x, res2$x, tolerance = 1e-12)
  # different seed should typically differ in either x or y
  expect_false(isTRUE(all.equal(res1$y, res3$y, tolerance = 1e-12)))
  expect_false(isTRUE(all.equal(res1$x, res3$x, tolerance = 1e-12)))
})


test_that("elitism and tpa options run and return valid structure", {
  f_sphere = function(x) {
    apply(x, 1, function(row) sum(row^2))
  }

  dim = 2
  x0 = rep(0.5, dim)
  lower = rep(-1, dim)
  upper = rep(1, dim)

  for (elitism in c(0L, 1L, 2L, 3L)) {
    for (tpa in c(0L, 1L, 2L)) {
      ctrl = cmaes_control(max_fevals = 20, lambda = 4,
        elitism = elitism, tpa = tpa, tpa_dsigma = if (tpa == 0L) NA_real_ else 0.1,
        seed = 5)
      res = cmaes(f_sphere, x0, lower, upper, ctrl, batch = TRUE)
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


test_that("ipop restarts double lambda across restarts when budget allows", {
  dim = 2
  seed = 42
  init_lambda = 2L # small initial lambda
  max_restarts = 3L
  fevals = 1500
  x0 = rep(0.5, dim)
  lower = rep(-1, dim)
  upper = rep(1, dim)

  # log batch sizes seen by objective
  fn = function(x) {
    batch_sizes <<- c(batch_sizes, nrow(x))
    apply(x, 1, function(row) sum(row^2))
  }

  for (algo in c("ipop", "bipop")) {
    batch_sizes <<- integer(0)
    ctrl = cmaes_control(algo = algo, max_fevals = fevals, lambda = init_lambda,
      max_restarts = max_restarts, seed = seed)
    res = cmaes(fn, x0, lower, upper, ctrl, batch = TRUE)
    bs = unique(batch_sizes)
    expect_contains(bs, c(2, 4))
  }
})
