library(smoof)
library(batchtools)
source("benchmarks/setup.R")

unlink(expsetup_reg_path, recursive = TRUE)
reg = makeExperimentRegistry(file.dir = expsetup_reg_path, packages = "smoof")

addProblem(reg, name = "bbob", fun = function(job, data, dim, fid, iid) {
  bfun = smoof::makeBBOBFunction(dim = dim, fid = fid, iid = iid)
  lower = smoof::getLowerBoxConstraints(bfun)
  upper = smoof::getUpperBoxConstraints(bfun)
  list(bfun = bfun, dim = dim, lower = lower, upper = upper, ybest = smoof::getGlobalOptimum(bfun)$value)
})

addAlgorithm(reg, name = "cma_es", fun = function(job, data, instance) {
  library(cmaes)
  loggedfun = addLoggingWrapper(instance$bfun, size = expsetup_flog_size)
  x0 = runif(instance$dim, min = instance$lower, max = instance$upper)
  ctrl = list(maxit = expsetup_max_iters)
  res = cmaes::cma_es(par = x0, fn = loggedfun, lower = instance$lower, upper = instance$upper, control = ctrl)
  get_result(loggedfun, instance$ybest)
})

addAlgorithm(reg, name = "libcmaesr", fun = function(job, data, instance) {
  library(libcmaesr)
  loggedfun = addLoggingWrapper(instance$bfun, size = expsetup_flog_size)
  x0 = runif(instance$dim, min = instance$lower, max = instance$upper)
  ctrl = cmaes_control(max_fevals = expsetup_max_fevals, algo = "bipop")
  res = libcmaesr::cmaes(
    x0 = x0,
    objective = loggedfun,
    lower = instance$lower,
    upper = instance$upper,
    control = ctrl,
    batch = FALSE
  )
  get_result(loggedfun, instance$ybest)
})


pdes = list(bbob = expand.grid(dim = expsetup_dims, fid = expsetup_fids, iid = expsetup_iids))
ades = list(cma_es = data.frame(), libcmaesr = data.frame())
addExperiments(reg = reg, prob.designs = pdes, algo.designs = ades, repls = expsetup_reps)


submitJobs(reg = reg)
waitForJobs(reg = reg)

# x = testJob(reg = reg, id = 1)
# print()
