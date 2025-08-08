library(devtools)
library(roxygen2)
document()
roxygenize()
load_all()

dim = 3
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
lambda = 3
fevals = 5000 * dim
ctrl = cmaes_control(max_fevals = fevals, seed = 123, lambda = lambda)
res = cmaes(fn, x0, lower, upper, ctrl)

print(dim(eval_log))
print(res)
