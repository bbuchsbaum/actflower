test_that("estimate_fc_lasso works when glmnet is installed", {
  skip_if_not_installed("glmnet")
  set.seed(51)
  x <- matrix(rnorm(8 * 150), nrow = 8, ncol = 150)
  fc <- estimate_fc_lasso(x, lambda = 0.05)
  expect_equal(dim(fc), c(8, 8))
  expect_true(all(diag(fc) == 0))
})

test_that("estimate_fc_lasso_cv returns fc and cv metadata", {
  skip_if_not_installed("glmnet")
  set.seed(52)
  x <- matrix(rnorm(7 * 140), nrow = 7, ncol = 140)
  out <- estimate_fc_lasso_cv(
    x,
    lambda = c(0.2, 0.1, 0.05),
    k_folds = 4,
    nodewise_hyperparams = TRUE,
    folds_scheme = "blocked"
  )
  expect_true(all(c("fc", "cv") %in% names(out)))
  expect_equal(dim(out$fc), c(7, 7))
  expect_true(length(out$cv$bestParam) == 7)
})
