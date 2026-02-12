test_that("estimate_fc routes corr and multreg methods", {
  set.seed(500)
  x <- matrix(rnorm(8 * 120), nrow = 8, ncol = 120)

  fc_corr <- estimate_fc(x, method = "corr")
  fc_multreg <- estimate_fc(x, method = "multreg")

  expect_equal(dim(fc_corr), c(8, 8))
  expect_equal(dim(fc_multreg), c(8, 8))
  expect_true(all(diag(fc_corr) == 0))
  expect_true(all(diag(fc_multreg) == 0))
})
