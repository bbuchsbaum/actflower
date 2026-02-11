test_that("estimate_fc_partial ridge returns finite square matrix", {
  set.seed(21)
  x <- matrix(rnorm(10 * 200), nrow = 10, ncol = 200)
  p <- estimate_fc_partial(x, method = "ridge", lambda = 0.1)
  expect_equal(dim(p), c(10, 10))
  expect_true(all(diag(p) == 0))
  expect_true(all(is.finite(p)))
})

test_that("estimate_fc_partial supports lambda grid CV", {
  set.seed(22)
  x <- matrix(rnorm(8 * 180), nrow = 8, ncol = 180)
  p <- estimate_fc_partial(x, method = "ridge", lambda = c(0.01, 0.1, 1), cv_folds = 4, cv_seed = 7)
  expect_equal(dim(p), c(8, 8))
  expect_true(attr(p, "lambda") %in% c(0.01, 0.1, 1))
  expect_true(is.list(attr(p, "cv")))
})

test_that("estimate_fc_partial pc returns directed coefficient matrix", {
  set.seed(23)
  x <- matrix(rnorm(7 * 160), nrow = 7, ncol = 160)
  b <- estimate_fc_partial(x, method = "pc", k = c(2, 3, 4), cv_folds = 3, cv_seed = 9)
  expect_equal(dim(b), c(7, 7))
  expect_true(all(diag(b) == 0))
  expect_true(is.list(attr(b, "cv")))
  expect_true(attr(b, "k") %in% c(2, 3, 4))
})

test_that("estimate_fc_partial glasso works when package installed", {
  skip_if_not_installed("glasso")
  set.seed(24)
  x <- matrix(rnorm(6 * 150), nrow = 6, ncol = 150)
  p <- estimate_fc_partial(x, method = "glasso", lambda = 0.05)
  expect_equal(dim(p), c(6, 6))
  expect_true(all(diag(p) == 0))
})
