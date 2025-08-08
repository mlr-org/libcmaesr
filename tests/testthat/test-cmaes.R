test_that("cmaes finds minimum of sphere function", {
  seed = 123
  # make sure to test the special case dim=1
  for (dim in 1:3) {
    fevals = 500 * dim 
    for (lambda in c(3 , 5, NA)) {
      # global log for all objective evaluations in this iteration
      eval_log <<- data.frame()

      fn = function(x) {
        if (!is.na(lambda)) {
          assert_matrix(x, nrows = lambda, ncols = dim)
        } else {
          assert_matrix(x, ncols = dim)
        }
        ys = apply(x, 1, function(row) sum(row^2))
        # append the current evaluation batch to the global log
        df = as.data.frame(x)
        names(df) = sprintf("x%d", seq_len(dim))
        df$y = ys
        eval_log <<- rbind(eval_log, df)
        ys
      }

      x0 = rep(0.5, dim)
      lower = rep(-1, dim)
      upper = rep(1, dim)
      ctrl = cmaes_control(max_fevals = fevals, lambda = lambda, seed = seed)
      res = cmaes(fn, x0, lower, upper, ctrl)

      # basic structure checks
      expect_list(res)
      expect_named(res, c("x", "y", "edm", "time", "status"), ignore.order = TRUE)
      expect_int(res$status, lower = 0)
      print(res$x)
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
        expect_true(all(eval_x[,k] >= lower[k]))
        expect_true(all(eval_x[,k] <= upper[k]))
      }
    }
  }
})

