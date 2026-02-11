test_that("calcactivity_parcelwise_noncircular builds target-source-condition array", {
  set.seed(71)
  x <- matrix(rnorm(6 * 4), nrow = 6, ncol = 4)
  ex <- lapply(seq_len(6), function(i) if (i %% 2 == 0) c(1, 2) else integer(0))

  arr <- calcactivity_parcelwise_noncircular(x, parcelstoexclude_bytarget = ex)
  expect_equal(dim(arr), c(6, 6, 4))
  # Diagonal target activity preserved.
  expect_equal(arr[3, 3, ], x[3, ])
})

test_that("calcconn_parcelwise_noncircular supports multreg and corr", {
  set.seed(72)
  x <- matrix(rnorm(8 * 180), nrow = 8, ncol = 180)
  ex <- replicate(8, integer(0), simplify = FALSE)
  ex[[1]] <- c(2, 3)

  fc_m <- calcconn_parcelwise_noncircular(x, connmethod = "multreg", parcelstoexclude_bytarget = ex)
  fc_c <- calcconn_parcelwise_noncircular(x, connmethod = "pearsoncorr", parcelstoexclude_bytarget = ex)

  expect_equal(dim(fc_m), c(8, 8))
  expect_equal(dim(fc_c), c(8, 8))
  expect_equal(fc_m[1, c(2, 3)], c(0, 0))
  expect_equal(fc_c[1, c(2, 3)], c(0, 0))
})

test_that("calcconn_parcelwise_noncircular supports combinedFC", {
  set.seed(73)
  x <- matrix(rnorm(7 * 200), nrow = 7, ncol = 200)
  ex <- replicate(7, integer(0), simplify = FALSE)

  fc <- calcconn_parcelwise_noncircular(x, connmethod = "combinedFC", parcelstoexclude_bytarget = ex)
  expect_equal(dim(fc), c(7, 7))
  expect_true(all(diag(fc) == 0))
})
