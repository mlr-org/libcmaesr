REG_PATH = "benchmarks/registry"

FLOG_SIZE = 10000

DIMS = 5
FIDS = 1:24
IIDS = 1:3
REPLS = 4L
MAX_FEVALS = 5000
MAX_ITERS = 500


GET_RESULT = function(loggedfun, ybest) {
  ys = environment(loggedfun)$obj.vals
  ys = as.numeric(ys)
  max_idx = environment(loggedfun)$curr.idx
  ys = ys[seq_len(max_idx)]
  ys = ys - ybest + 1e-10
  list(ys = ys)
}
