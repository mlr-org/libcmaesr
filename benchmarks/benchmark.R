library(smoof)
library(batchtools)
source("benchmarks/setup.R")

unlink(REG_PATH, recursive = TRUE)
reg = makeExperimentRegistry(file.dir = REG_PATH, packages = "smoof")

addProblem(reg, name = "bbob", fun = function(job, data, dim, fid, iid) {
  bfun = smoof::makeBBOBFunction(dim = dim, fid = fid, iid = iid)
  lower = smoof::getLowerBoxConstraints(bfun)
  upper = smoof::getUpperBoxConstraints(bfun)
  list(bfun = bfun, dim = dim, lower = lower, upper = upper, ybest = smoof::getGlobalOptimum(bfun)$value)
})

addAlgorithm(reg, name = "cma_es", fun = function(job, data, instance) {
  library(cmaes)
  loggedfun = addLoggingWrapper(instance$bfun, size = FLOG_SIZE)
  x0 = runif(instance$dim, min = instance$lower, max = instance$upper)
  ctrl = list(maxit = MAX_ITERS)
  res = cmaes::cma_es(par = x0, fn = loggedfun, lower = instance$lower, upper = instance$upper, control = ctrl)
  GET_RESULT(loggedfun, instance$ybest)
})

addAlgorithm(reg, name = "libcmaesr", fun = function(job, data, instance) {
  library(libcmaesr)
  loggedfun = addLoggingWrapper(instance$bfun, size = FLOG_SIZE)
  x0 = runif(instance$dim, min = instance$lower, max = instance$upper)
  ctrl = cmaes_control(max_fevals = MAX_FEVALS, algo = "bipop")
  res = libcmaesr::cmaes(
    x0 = x0,
    objective = loggedfun,
    lower = instance$lower,
    upper = instance$upper,
    control = ctrl,
    batch = FALSE
  )
  GET_RESULT(loggedfun, instance$ybest)
})


pdes = list(bbob = expand.grid(dim = DIMS, fid = FIDS, iid = IIDS))
ades = list(cma_es = data.frame(), libcmaesr = data.frame())
addExperiments(reg = reg, prob.designs = pdes, algo.designs = ades, repls = REPLS)


submitJobs(reg = reg)
waitForJobs(reg = reg)

# x = testJob(reg = reg, id = 1)
# print()
