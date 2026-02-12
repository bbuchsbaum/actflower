test_that("estimate_fc_multreg matches OLS slopes", {
  set.seed(11)
  n_nodes <- 6
  n_time <- 120
  x <- matrix(rnorm(n_nodes * n_time), nrow = n_nodes, ncol = n_time)

  b <- estimate_fc_multreg(x)
  expect_equal(dim(b), c(n_nodes, n_nodes))
  expect_true(all(diag(b) == 0))

  # Compare against direct regressions with intercept.
  for (target in seq_len(n_nodes)) {
    src <- setdiff(seq_len(n_nodes), target)
    fit <- stats::lm.fit(x = cbind(1, t(x[src, , drop = FALSE])), y = as.numeric(x[target, ]))
    coef_ref <- fit$coefficients[-1]
    expect_equal(as.numeric(b[target, src]), as.numeric(coef_ref), tolerance = 1e-6)
  }
})

test_that("estimate_fc_multreg use_cpp path matches R path", {
  set.seed(12)
  x <- matrix(rnorm(7 * 140), nrow = 7, ncol = 140)
  b_r <- estimate_fc_multreg(x, use_cpp = FALSE)
  b_cpp <- estimate_fc_multreg(x, use_cpp = TRUE)
  expect_equal(b_cpp, b_r, tolerance = 1e-6)
})

test_that("estimate_fc_multreg validates dimensions and ridge", {
  x_bad <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)
  expect_error(estimate_fc_multreg(x_bad), "nodes < timepoints")

  x <- matrix(rnorm(6 * 80), nrow = 6, ncol = 80)
  expect_error(estimate_fc_multreg(x, ridge = -0.1), "non-negative scalar")
})

test_that("estimate_fc_multreg supports orientation aliases", {
  set.seed(13)
  x <- matrix(rnorm(6 * 90), nrow = 6, ncol = 90)
  b_a <- estimate_fc_multreg(x, orientation = "nodes_by_time")
  b_b <- estimate_fc_multreg(x, orientation = "nodes_x_time")
  b_c <- estimate_fc_multreg(t(x), orientation = "time_x_nodes")
  expect_equal(b_b, b_a, tolerance = 1e-8)
  expect_equal(b_c, b_a, tolerance = 1e-8)
})
