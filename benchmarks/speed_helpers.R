library(libcmaesr)
library(cmaes)
library(data.table)

check_fevals = function(smoof_fun, budget, lambda) {
  assert_class(smoof_fun, "smoof_counting_function")
  nevals = getNumberOfEvaluations(smoof_fun)
  print(nevals)
  if (nevals < budget - lambda*10 || nevals > budget + lambda*10) {
    stop("maxfevals not respected")
  }
  resetEvaluationCounter(smoof_fun)
  invisible(NULL)
}

run_speedtest = function(obj_single, obj_batch, dim, budget, reps, do_check) {
  budget = budget * dim

  lower = rep(-5, dim)
  upper = rep(5, dim)

  # fix lambda so we can calc cmaes maxit / budget
  lambda = 4L + floor(3 * log(dim))
  results = NULL

  for (r in seq_len(reps)) {
    x0 = runif(dim, min = lower, max = upper)

    print("Running libcmaesr single")
    libcmaesr_ctrl = libcmaesr::cmaes_control(
      algo = "ipop", # run ipop so we cannot stop early due to convergence
      max_fevals = budget,
      max_restarts = 1000,
      lambda = lambda
    )

    libcmaesr_single = system.time({
      libcmaesr::cmaes(
        objective = obj_single,
        x0 = x0,
        lower = lower,
        upper = upper,
        control = libcmaesr_ctrl,
        batch = FALSE
      )
    })["elapsed"]
    if (do_check) {
      check_fevals(obj_single, budget, lambda)
    }

    print("Running libcmaesr batch")
    libcmaesr_batch = system.time({
      libcmaesr::cmaes(
        objective = obj_batch,
        x0 = x0,
        lower = lower,
        upper = upper,
        control = libcmaesr_ctrl,
        batch = TRUE
      )
    })["elapsed"]

    cmaes_maxit = ceiling(budget / lambda)
    cmaes_ctrl = list(
      maxit = cmaes_maxit,
      lambda = lambda,
      sc_tolx = 0
    )

    print("Running cmaes single")
    cmaes_ctrl$vectorized = FALSE
    cmaes_single = tryCatch({
      system.time({
        cmaes::cma_es(
          par = x0,
          fn = obj_single,
          lower = lower,
          upper = upper,
          control = cmaes_ctrl
        )
      })["elapsed"]
    }, error = function(e) NA)
    if (do_check && !is.na(cmaes_single)) {
      check_fevals(obj_single, budget, lambda)
    }

    print("Running cmaes batch")
    cmaes_ctrl$vectorized = TRUE
    cmaes_batch = tryCatch({
    system.time({
      cmaes::cma_es(
        par = x0,
        fn = obj_batch,
        lower = lower,
        upper = upper,
        control = cmaes_ctrl
        )
      })["elapsed"]
    }, error = function(e) NA)

    row = data.table(
      rep = r,
      dim = dim,
      budget = budget,
      cmaes_single = as.numeric(cmaes_single),
      cmaes_batch = as.numeric(cmaes_batch),
      libcmaesr_single = as.numeric(libcmaesr_single),
      libcmaesr_batch = as.numeric(libcmaesr_batch)
    )
    results = if (is.null(results)) row else rbindlist(list(results, row))
  }
  results
}

