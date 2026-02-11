test_that("actflow_predict returns length-n vector", {
  set.seed(2)
  n <- 8
  act <- rnorm(n)
  fc <- matrix(rnorm(n * n), n, n)
  pred <- actflow_predict(act, fc)
  expect_length(pred, n)
})
