test_that("hdf5 write/read roundtrip", {
  skip_if_not_installed("hdf5r")
  tf <- tempfile(fileext = ".h5")
  on.exit(unlink(tf), add = TRUE)
  x <- matrix(1:12, nrow = 3)
  write_actflower_h5(tf, "x", x)
  y <- read_actflower_h5(tf, "x")
  expect_equal(y, x)
})

test_that("hdf5 write/read supports nested dataset paths", {
  skip_if_not_installed("hdf5r")
  tf <- tempfile(fileext = ".h5")
  on.exit(unlink(tf), add = TRUE)

  x <- array(seq_len(24), dim = c(2, 3, 4))
  write_actflower_h5(tf, "metrics/fullcompare/corr_vals", x)
  y <- read_actflower_h5(tf, "metrics/fullcompare/corr_vals")
  expect_equal(y, x)
})

test_that("hdf5 read reports missing dataset cleanly for nested paths", {
  skip_if_not_installed("hdf5r")
  tf <- tempfile(fileext = ".h5")
  on.exit(unlink(tf), add = TRUE)
  write_actflower_h5(tf, "x", 1:3)

  expect_error(
    read_actflower_h5(tf, "metrics/fullcompare/corr_vals"),
    "Dataset not found in HDF5 file"
  )
})
