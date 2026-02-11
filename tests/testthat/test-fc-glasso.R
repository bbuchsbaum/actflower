test_that("estimate_fc_glasso works when glasso is installed", {
  skip_if_not_installed("glasso")
  set.seed(61)
  x <- matrix(rnorm(8 * 160), nrow = 8, ncol = 160)
  fc <- estimate_fc_glasso(x, lambda = 0.05)
  expect_equal(dim(fc), c(8, 8))
  expect_true(all(diag(fc) == 0))
})

test_that("estimate_fc_glasso_cv returns fc and cv metadata", {
  skip_if_not_installed("glasso")
  set.seed(62)
  x <- matrix(rnorm(7 * 150), nrow = 7, ncol = 150)
  out <- estimate_fc_glasso_cv(
    x,
    lambda = c(0.2, 0.1, 0.05),
    k_folds = 4,
    opt_method = "loglikelihood",
    folds_scheme = "blocked"
  )
  expect_true(all(c("fc", "cv") %in% names(out)))
  expect_equal(dim(out$fc), c(7, 7))
  expect_true(length(out$cv$bestParam) == 1)
})
