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
