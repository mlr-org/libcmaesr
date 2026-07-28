test_that("print.cmaes_control only shows non-default fields", {
  expect_output(print(cmaes_control()), "cmaes_control()", fixed = TRUE)

  out = capture.output(print(cmaes_control(algo = "bipop", max_fevals = 1000, seed = 42)))
  expect_match(out[2L], "algo = \"bipop\"", fixed = TRUE)
  expect_match(out[2L], "max_fevals = 1000L", fixed = TRUE)
  expect_match(out[2L], "seed = 42L", fixed = TRUE)
  expect_no_match(out[2L], "maximize", fixed = TRUE)
  expect_no_match(out[2L], "quiet", fixed = TRUE)
})

test_that("print.cmaes_control shows fields mutated after construction", {
  control = cmaes_control()
  control$quiet = FALSE
  control$x0_lower = c(-1, -1)
  control$x0_upper = c(1, 1)
  out = paste(capture.output(print(control)), collapse = "")
  expect_match(out, "quiet = FALSE", fixed = TRUE)
  expect_match(out, "x0_lower = c(-1, -1)", fixed = TRUE)
  expect_no_match(out, "algo", fixed = TRUE)
})
