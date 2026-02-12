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

test_that("estimate_fc routes partial and combined methods", {
  set.seed(501)
  x <- matrix(rnorm(9 * 140), nrow = 9, ncol = 140)

  fc_partial <- estimate_fc(x, method = "partial", partial_method = "ridge", lambda = 0.1)
  fc_combined <- estimate_fc(x, method = "combined")

  expect_equal(dim(fc_partial), c(9, 9))
  expect_equal(dim(fc_combined), c(9, 9))
  expect_true(all(diag(fc_partial) == 0))
  expect_true(all(diag(fc_combined) == 0))
})

test_that("estimate_fc routes lasso_cv and glasso_cv methods", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("glasso")

  set.seed(502)
  x <- matrix(rnorm(8 * 90), nrow = 8, ncol = 90)

  fc_lasso <- estimate_fc(
    x,
    method = "lasso_cv",
    lambda = c(0.1, 0.03),
    k_folds = 3,
    nodewise_hyperparams = FALSE,
    folds_scheme = "interleaving"
  )
  fc_glasso <- estimate_fc(
    x,
    method = "glasso_cv",
    lambda = c(0.1, 0.03),
    k_folds = 3,
    opt_method = "R2",
    folds_scheme = "interleaving"
  )

  expect_true(is.list(fc_lasso))
  expect_true(is.list(fc_glasso))
  expect_equal(dim(fc_lasso$fc), c(8, 8))
  expect_equal(dim(fc_glasso$fc), c(8, 8))
})
