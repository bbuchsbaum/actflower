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
