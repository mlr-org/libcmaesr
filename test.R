library(devtools)
library(roxygen2)
# document()
# roxygenize()
load_all()
seed = 1
set.seed(seed)

dim = 3
eval_log <<- data.frame()
batch_sizes <<- integer(0)

fn_single = function(x) {
  #print(x)
  #ys = sum(x^2 - 10 * cos(2 * pi * x))
  ys = sum(x^2)
  df = as.data.frame(as.list(x))
  names(df) = sprintf("x%d", seq_len(dim))
  df$y = ys
  eval_log <<- rbind(eval_log, df)
  ys
}

fn_batch = function(x) {
  n = ncol(x)
  ys = 10 * n + rowSums(x^2 - 10 * cos(2 * pi * x))
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
z = cmaes(fn_single, x0, lower, upper, ctrl, batch = FALSE)
#print(unique(batch_sizes))
print(z)

