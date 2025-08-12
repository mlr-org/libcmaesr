library(mlr3learners)
library(mlr3)
library(mlbench)
source("benchmarks/speed_helpers.R")

mytask = tsk("california_housing")
mytask = mytask$select(setdiff(mytask$feature_names, "ocean_proximity"))
mytask = mytask$filter(1:1000)

mylrn = lrn("regr.ranger", num.trees = 100)
mylrn$train(mytask)

DIM = length(mytask$feature_names)


obj_single = function(x) {
  pd = as.data.table(as.list(x))
  colnames(pd) = mytask$feature_names
  mylrn$predict_newdata(pd)$response
}

obj_batch = function(x) {
  pd = as.data.table(x)
  colnames(pd) = mytask$feature_names
  mylrn$predict_newdata(pd)$response
}

REPLS = 3L
BUDGETS = c(10L, 20L)

results = list()
for (budget in BUDGETS) {
  print(paste("Running budget", budget))
  dd = run_speedtest(obj_single = obj_single, obj_batch = obj_batch,
    dim = DIM, budget = budget, reps = REPLS, do_check = FALSE)
  results[[length(results) + 1]] = dd
}
dd = rbindlist(results)
print(dd)

# Aggregate over reps and plot mean time by budget, facetted by dimension
long = melt(dd,
  id.vars = c("rep", "dim", "budget"),
  measure.vars = c("cmaes_single", "cmaes_batch", "libcmaesr_single", "libcmaesr_batch"),
  variable.name = "algorithm",
  value.name = "time"
)

# Count how often NA occurred per dim/budget/algorithm and overall per algorithm
na_counts = long[, .(na_count = sum(is.na(time)), total = .N), by = .(dim, algorithm)]
print(na_counts)

avg = long[, .(time_mean = mean(time, na.rm = TRUE)), by = .(dim, budget, algorithm)]

p = ggplot(avg, aes(x = budget, y = time_mean, color = algorithm)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ dim, scales = "free") +
  labs(x = "Budget (function evaluations)", y = "Mean elapsed time [s]", color = "Algorithm") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("benchmarks/speed_ml.png", p, width = 8, height = 5, dpi = 120)


