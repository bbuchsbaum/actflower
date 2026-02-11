test_that("estimate_fc_corr returns square matrix with zero diagonal", {
  set.seed(1)
  x <- matrix(rnorm(300), nrow = 10, ncol = 30)
  fc <- estimate_fc_corr(x)
  expect_equal(dim(fc), c(10, 10))
  expect_true(all(diag(fc) == 0))
})

test_that("estimate_fc_corr use_cpp path matches R path", {
  set.seed(11)
  x <- matrix(rnorm(12 * 80), nrow = 12, ncol = 80)
  fc_r <- estimate_fc_corr(x, use_cpp = FALSE)
  fc_cpp <- estimate_fc_corr(x, use_cpp = TRUE)
  expect_equal(fc_cpp, fc_r, tolerance = 1e-12)
})

test_that("estimate_fc_corr zero-variance rows align between R and cpp paths", {
  x <- matrix(rnorm(6 * 40), nrow = 6, ncol = 40)
  x[2, ] <- 1
  x[5, ] <- -3
  fc_r <- suppressWarnings(estimate_fc_corr(x, use_cpp = FALSE))
  fc_cpp <- estimate_fc_corr(x, use_cpp = TRUE)
  expect_equal(fc_cpp, fc_r, tolerance = 1e-12)
})

test_that("estimate_fc routes partial methods", {
  set.seed(4)
  x <- matrix(rnorm(8 * 120), nrow = 8, ncol = 120)
  fc <- estimate_fc(x, method = "partial", partial_method = "ridge", lambda = 0.1)
  expect_equal(dim(fc), c(8, 8))
  expect_true(all(diag(fc) == 0))
})
