expsetup_reg_path = "benchmarks/registry"

expsetup_flog_size = 10000

expsetup_dims = 5
expsetup_fids = 1:24
expsetup_iids = 1:3
expsetup_reps = 4L
expsetup_max_fevals = 5000
expsetup_max_iters = 500


get_result = function(loggedfun, ybest) {
  ys = environment(loggedfun)$obj.vals
  ys = as.numeric(ys)
  max_idx = environment(loggedfun)$curr.idx
  ys = ys[seq_len(max_idx)]
  ys = ys - ybest + 1e-10
  list(ys = ys)
}
