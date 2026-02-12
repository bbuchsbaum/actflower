test_that("estimate_fc_corr returns square matrix with zero diagonal", {
  set.seed(1)
  x <- matrix(rnorm(300), nrow = 10, ncol = 30)
  fc <- estimate_fc_corr(x)
  expect_equal(dim(fc), c(10, 10))
  expect_true(all(diag(fc) == 0))
})

test_that("estimate_fc_corr accepts orientation aliases", {
  set.seed(2)
  x <- matrix(rnorm(10 * 30), nrow = 10, ncol = 30)
  fc_a <- estimate_fc_corr(x, orientation = "nodes_by_time", use_cpp = TRUE)
  fc_b <- estimate_fc_corr(x, orientation = "nodes_x_time", use_cpp = TRUE)
  fc_c <- estimate_fc_corr(t(x), orientation = "time_x_nodes", use_cpp = TRUE)
  expect_equal(fc_b, fc_a, tolerance = 1e-12)
  expect_equal(fc_c, fc_a, tolerance = 1e-12)
})

test_that("estimate_fc_corr use_cpp path matches R path", {
  set.seed(11)
  x <- matrix(rnorm(12 * 80), nrow = 12, ncol = 80)
  fc_r <- estimate_fc_corr(x, use_cpp = FALSE)
  fc_cpp <- estimate_fc_corr(x, use_cpp = TRUE)
  expect_equal(fc_cpp, fc_r, tolerance = 1e-12)
})

test_that("estimate_fc_corr supports 3D arrays and matches per-subject path", {
  set.seed(12)
  x <- array(rnorm(10 * 60 * 4), dim = c(10, 60, 4))
  fc_batch <- estimate_fc_corr(x, orientation = "nodes_by_time", use_cpp = TRUE)

  fc_loop <- array(NA_real_, dim = c(10, 10, 4))
  for (s in seq_len(4)) {
    fc_loop[, , s] <- estimate_fc_corr(x[, , s], orientation = "nodes_by_time", use_cpp = TRUE)
  }

  expect_equal(fc_batch, fc_loop, tolerance = 1e-12)
})

test_that("estimate_fc_corr 3D orientation aliases are consistent", {
  set.seed(13)
  x <- array(rnorm(8 * 50 * 3), dim = c(8, 50, 3))
  x_t <- aperm(x, c(2, 1, 3))

  fc_a <- estimate_fc_corr(x, orientation = "nodes_by_time", use_cpp = TRUE)
  fc_b <- estimate_fc_corr(x_t, orientation = "time_by_nodes", use_cpp = TRUE)

  expect_equal(fc_b, fc_a, tolerance = 1e-12)
})

test_that("estimate_fc_corr zero-variance rows align between R and cpp paths", {
  x <- matrix(rnorm(6 * 40), nrow = 6, ncol = 40)
  x[2, ] <- 1
  x[5, ] <- -3
  fc_r <- suppressWarnings(estimate_fc_corr(x, use_cpp = FALSE))
  fc_cpp <- estimate_fc_corr(x, use_cpp = TRUE)
  expect_equal(fc_cpp, fc_r, tolerance = 1e-12)
})

test_that("estimate_fc_corr 3D uses pairwise-complete fallback with NA", {
  set.seed(14)
  x <- array(rnorm(7 * 45 * 3), dim = c(7, 45, 3))
  x[1, 1, 1] <- NA_real_
  x[4, 10, 2] <- NA_real_

  fc_batch <- estimate_fc_corr(x, orientation = "nodes_by_time", use_cpp = TRUE)
  fc_ref <- array(NA_real_, dim = c(7, 7, 3))
  for (s in seq_len(3)) {
    fc_ref[, , s] <- estimate_fc_corr(x[, , s], orientation = "nodes_by_time", use_cpp = FALSE)
  }

  expect_equal(fc_batch, fc_ref, tolerance = 1e-12)
})

test_that("estimate_fc routes partial methods", {
  set.seed(4)
  x <- matrix(rnorm(8 * 120), nrow = 8, ncol = 120)
  fc <- estimate_fc(x, method = "partial", partial_method = "ridge", lambda = 0.1)
  expect_equal(dim(fc), c(8, 8))
  expect_true(all(diag(fc) == 0))
})
