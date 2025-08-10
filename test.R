library(devtools)
library(roxygen2)
# document()
# roxygenize()
load_all()
seed = 1
set.seed(seed)

dim = 100
eval_log <<- data.frame()
batch_sizes <<- integer(0)

fn = function(x) {
  n = ncol(x)
  ys = 10 * n + rowSums(x^2 - 10 * cos(2 * pi * x))
  df = as.data.frame(x)
  names(df) = sprintf("x%d", seq_len(dim))
  df$y = ys
  eval_log <<- rbind(eval_log, df)
  batch_sizes <<- c(batch_sizes, nrow(x))
  print(sum(batch_sizes))
  ys
}

max_restarts = 3L
fevals = 500


x0 = rep(3, dim)
lower = rep(-5.12, dim)
upper = rep(5.12, dim)
ctrl = cmaes_control(algo = "bipop", max_fevals = fevals,
  max_restarts = max_restarts, seed = seed)
z = cmaes(fn, x0, lower, upper, ctrl)
print(unique(batch_sizes))
print(z)

