test_that("actflow_predict returns length-n vector", {
  set.seed(2)
  n <- 8
  act <- rnorm(n)
  fc <- matrix(rnorm(n * n), n, n)
  pred <- actflow_predict(act, fc)
  expect_length(pred, n)
})

test_that("actflow_predict excludes self-connections for vector path", {
  n <- 6
  act <- rep(1, n)
  fc <- matrix(0, n, n)
  diag(fc) <- 5
  pred <- actflow_predict(act, fc, separate_by_target = FALSE, transfer = "linear")
  expect_equal(pred, rep(0, n))
})

test_that("actflow_predict supports separate_by_target path", {
  n <- 5
  act <- matrix(seq_len(n * n), nrow = n, ncol = n)
  fc <- matrix(1, n, n)
  diag(fc) <- 2

  pred <- actflow_predict(act, fc, separate_by_target = TRUE, transfer = "linear")
  expect_length(pred, n)
  expect_true(all(is.finite(pred)))
})

test_that("actflow_predict validates act shape by mode", {
  n <- 4
  fc <- matrix(1, n, n)
  expect_error(
    actflow_predict(act = matrix(1, 3, 2), fc = fc, separate_by_target = FALSE),
    "must have length"
  )
  expect_error(
    actflow_predict(act = rep(1, n), fc = fc, separate_by_target = TRUE),
    "must be"
  )
})

test_that("transfer_function covers relu, sigmoid, and logit branches", {
  x <- c(-1, 0, 1)
  expect_equal(transfer_function(x, transfer = "relu", threshold = 0), c(0, 0, 1))

  s <- transfer_function(c(-2, 2), transfer = "sigmoid")
  expect_true(all(s > 0 & s < 1))

  y <- c(0.2, 0.8)
  lg <- transfer_function(y, transfer = "logit", a = 2)
  expect_equal(lg, (1 / 2) * log(y / (1 - y)))
})
