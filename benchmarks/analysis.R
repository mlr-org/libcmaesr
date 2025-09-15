library(batchtools)
library(ggplot2)
library(data.table)
source("benchmarks/setup.R")

reg = loadRegistry(file.dir = expsetup_reg_path, writeable = FALSE)

res = reduceResultsDataTable(reg = reg)
res = ijoin(getJobPars(reg = reg), res)
jt = getJobTable(reg = reg)
res = unwrap(res)
res$repl = jt$repl
print(res)

# Unnest ys into long format and compute cumulative min per job
long = res[,
  {
    y = as.numeric(ys[[1]])
    data.table(eval = seq_along(y), y = y, algorithm = algorithm, dim = dim, fid = fid, iid = iid, repl = repl)
  },
  by = .(job.id)
]

long[, ybest := cummin(y), by = job.id]

# Start plotting from eval >= 500
plot_long = long[eval >= 500 & iid == 1 & repl == 1]

# Average the curve over repls and iids: mean per algorithm/dim/fid/eval
avg_long = plot_long[, .(ybest_mean = mean(ybest, na.rm = TRUE)), by = .(algorithm, dim, fid, eval)]

# Plot averaged cummin curves per algorithm, facetted by fid
p = ggplot(avg_long, aes(x = eval, y = ybest_mean, color = algorithm)) +
  geom_line() +
  scale_y_log10() +
  facet_wrap(~fid, scales = "free") +
  labs(x = "Evals", y = "y (mean cummin across repls & iids)") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("benchmarks/benchmark.png", p, width = 9, height = 6, dpi = 120)
