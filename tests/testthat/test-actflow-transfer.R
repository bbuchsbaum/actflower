test_that("transfer_function supports linear/relu/sigmoid/logit", {
  x <- c(-1, 0, 0.5, 1)

  expect_equal(transfer_function(x, transfer = "linear"), x)
  expect_equal(transfer_function(x, transfer = "relu", threshold = 0), c(0, 0, 0.5, 1))

  sig <- transfer_function(x, transfer = "sigmoid")
  expect_true(all(sig > 0 & sig < 1))

  p <- c(0.2, 0.5, 0.8)
  expect_equal(
    transfer_function(p, transfer = "logit", a = 1),
    qlogis(p),
    tolerance = 1e-12
  )
})

test_that("transfer_function validates transfer argument", {
  expect_error(transfer_function(1:3, transfer = "unknown"), "should be one of")
})
