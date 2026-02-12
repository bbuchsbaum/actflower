test_that("estimate_fc_combined returns square matrix with zero diagonal", {
  set.seed(31)
  x <- matrix(rnorm(10 * 150), nrow = 10, ncol = 150)
  fc <- estimate_fc_combined(x)
  expect_equal(dim(fc), c(10, 10))
  expect_true(all(diag(fc) == 0))
  expect_true(all(is.finite(fc)))
})

test_that("estimate_fc with method='combined' routes correctly", {
  set.seed(32)
  x <- matrix(rnorm(8 * 120), nrow = 8, ncol = 120)
  fc <- estimate_fc(x, method = "combined")
  expect_equal(dim(fc), c(8, 8))
})

test_that("estimate_fc_combined supports multipleRegression + simpleRegression branches", {
  set.seed(33)
  x <- matrix(rnorm(9 * 140), nrow = 9, ncol = 140)
  fc <- estimate_fc_combined(
    x,
    method_cond_asso = "multipleRegression",
    method_asso = "simpleRegression"
  )
  expect_equal(dim(fc), c(9, 9))
  expect_true(all(diag(fc) == 0))
})

test_that("estimate_fc_combined supports equivalence-test association branch", {
  set.seed(34)
  x <- matrix(rnorm(7 * 160), nrow = 7, ncol = 160)
  fc <- estimate_fc_combined(
    x,
    method_asso = "correlation",
    equivalence_test_asso = TRUE,
    lower_bound = -0.2,
    upper_bound = 0.2
  )
  expect_equal(dim(fc), c(7, 7))
  expect_true(all(diag(fc) == 0))
})

test_that("combinedFC helper guards behave as expected", {
  expect_error(.af_combined_zcutoff(kind = "bad-kind"), "Unknown zcutoff kind")

  z <- .af_combined_fisher_z(r = c(-1, 0, 1), n_datapoints = 20, ho = 0)
  expect_true(all(is.finite(z)))

  sr <- .af_combined_simple_regression_sig(y = rnorm(60), x = rnorm(60), sig_test = FALSE)
  expect_true(is.list(sr))
  expect_true("b" %in% names(sr))
})

test_that("combinedFC helper internals return finite matrices/vectors", {
  set.seed(35)
  d <- matrix(rnorm(140 * 6), nrow = 140, ncol = 6)

  pc_reg <- .af_combined_partial_correlation_sig(d, alpha = 0.05, method = "regression")
  pc_inv <- .af_combined_partial_correlation_sig(d, alpha = 0.05, method = "inverseCovariance")
  expect_equal(dim(pc_reg), c(6, 6))
  expect_equal(dim(pc_inv), c(6, 6))

  mr_sig <- .af_combined_multiple_regression_sig(d, alpha = 0.05, sig_test = TRUE)
  mr_nosig <- .af_combined_multiple_regression_sig(d, alpha = 0.05, sig_test = FALSE)
  expect_equal(dim(mr_sig), c(6, 6))
  expect_equal(dim(mr_nosig), c(6, 6))

  c_sig <- .af_combined_correlation_sig(d, alpha = 0.05, equivalence_test = FALSE)
  c_eq <- .af_combined_correlation_sig(
    d, alpha = 0.05, lower_bound = -0.2, upper_bound = 0.2, equivalence_test = TRUE
  )
  expect_equal(dim(c_sig), c(6, 6))
  expect_equal(dim(c_eq), c(6, 6))

  inv <- .af_combined_pinv_sym(crossprod(d))
  expect_equal(dim(inv), c(6, 6))
  expect_true(all(is.finite(inv)))
})
