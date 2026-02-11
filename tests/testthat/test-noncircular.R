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

test_that("calcconn_parcelwise_noncircular sparse path matches dense path", {
  set.seed(74)
  x <- matrix(rnorm(9 * 220), nrow = 9, ncol = 220)
  ex <- lapply(seq_len(9), function(i) sample(setdiff(seq_len(9), i), size = 2))

  methods <- c("multreg", "pearsoncorr", "pc_multregconn", "combinedFC")
  for (m in methods) {
    fc_dense <- calcconn_parcelwise_noncircular(
      x,
      connmethod = m,
      parcelstoexclude_bytarget = ex,
      sparse = FALSE
    )
    fc_sparse_list <- calcconn_parcelwise_noncircular(
      x,
      connmethod = m,
      parcelstoexclude_bytarget = ex,
      sparse = TRUE,
      mask_format = "list"
    )
    fc_sparse_dgc <- calcconn_parcelwise_noncircular(
      x,
      connmethod = m,
      parcelstoexclude_bytarget = ex,
      sparse = TRUE,
      mask_format = "dgCMatrix"
    )

    expect_equal(fc_sparse_list, fc_dense, tolerance = 1e-10, ignore_attr = TRUE)
    expect_equal(fc_sparse_dgc, fc_dense, tolerance = 1e-10, ignore_attr = TRUE)
    expect_equal(attr(fc_sparse_dgc, "source_mask_format"), "dgCMatrix")
  }
})

test_that("calcconn_parcelwise_noncircular accepts sparse exclusion matrices", {
  set.seed(75)
  x <- matrix(rnorm(8 * 200), nrow = 8, ncol = 200)
  ex_dense <- matrix(0, 8, 8)
  ex_dense[1, c(2, 3)] <- 1
  ex_dense[3, 7] <- 1
  ex_sparse <- Matrix::sparseMatrix(
    i = which(ex_dense != 0, arr.ind = TRUE)[, 1],
    j = which(ex_dense != 0, arr.ind = TRUE)[, 2],
    x = 1,
    dims = c(8, 8)
  )

  fc_from_dense <- calcconn_parcelwise_noncircular(
    x,
    connmethod = "multreg",
    parcelstoexclude_bytarget = ex_dense
  )
  fc_from_sparse <- calcconn_parcelwise_noncircular(
    x,
    connmethod = "multreg",
    parcelstoexclude_bytarget = ex_sparse
  )
  expect_equal(fc_from_sparse, fc_from_dense, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("calcactivity_parcelwise_noncircular sparse path matches dense path", {
  set.seed(76)
  x <- matrix(rnorm(7 * 5), nrow = 7, ncol = 5)
  ex <- lapply(seq_len(7), function(i) if (i %% 3 == 0) c(1, 2) else integer(0))

  arr_dense <- calcactivity_parcelwise_noncircular(x, parcelstoexclude_bytarget = ex, sparse = FALSE)
  arr_sparse_list <- calcactivity_parcelwise_noncircular(
    x,
    parcelstoexclude_bytarget = ex,
    sparse = TRUE,
    mask_format = "list"
  )
  arr_sparse_dgc <- calcactivity_parcelwise_noncircular(
    x,
    parcelstoexclude_bytarget = ex,
    sparse = TRUE,
    mask_format = "dgCMatrix"
  )

  expect_equal(arr_sparse_list, arr_dense, tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(arr_sparse_dgc, arr_dense, tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(attr(arr_sparse_dgc, "source_mask_format"), "dgCMatrix")
})
