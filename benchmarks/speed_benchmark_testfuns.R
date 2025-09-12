REPLS = 5L
DIMS = c(2L, 5L, 10L, 20L)
BUDGETS = c(100L, 200L, 500L)

results = list()
for (dim in DIMS) {
  for (budget in BUDGETS) {
    print(paste("Running", dim, "D,", budget, "B"))
    obj_single = addCountingWrapper(makeRastriginFunction(dim))
    obj_batch = function(x) rowSums(x^2)
    dd = run_speedtest(obj_single, obj_batch, dim, budget, reps = REPLS)
    results[[length(results) + 1]] = dd
  }
}
dd = rbindlist(results)
print(dd)

# Aggregate over reps and plot mean time by budget, facetted by dimension
long = melt(
  dd,
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
  facet_wrap(~dim, scales = "free") +
  labs(x = "Budget (function evaluations)", y = "Mean elapsed time [s]", color = "Algorithm") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("benchmarks/speed_testfuns.png", p, width = 8, height = 5, dpi = 120)
