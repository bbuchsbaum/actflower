test_that("compare_fullcomp_cpp handles NA and dimension mismatch guards", {
  target <- array(rnorm(5 * 4 * 2), dim = c(5, 4, 2))
  pred <- target
  pred[1, 1, 1] <- NA_real_

  out <- compare_fullcomp_cpp(target, pred)
  expect_true(all(c("corr_vals", "R2_vals", "mae_vals") %in% names(out)))
  expect_equal(length(out$corr_vals), 2)

  bad_pred <- array(rnorm(5 * 4 * 3), dim = c(5, 4, 3))
  expect_error(compare_fullcomp_cpp(target, bad_pred), "dimensions must match")
})

test_that("compare_fullcomp_cpp handles low-information edge cases", {
  target <- array(NA_real_, dim = c(3, 2, 2))
  pred <- array(NA_real_, dim = c(3, 2, 2))
  target[1, 1, 1] <- 1
  pred[1, 1, 1] <- 1

  out <- compare_fullcomp_cpp(target, pred)
  expect_true(is.na(out$corr_vals[1]))
  expect_true(is.na(out$R2_vals[1]))
  expect_true(is.na(out$mae_vals[1]))

  # Constant vectors: correlation denominator and SS_tot are zero.
  target2 <- array(2, dim = c(3, 2, 1))
  pred2 <- array(2, dim = c(3, 2, 1))
  out2 <- compare_fullcomp_cpp(target2, pred2)
  expect_true(is.na(out2$corr_vals[1]))
  expect_true(is.na(out2$R2_vals[1]))
  expect_equal(out2$mae_vals[1], 0)
})

test_that("noncircular native activity kernels match R path and validate inputs", {
  set.seed(400)
  n_nodes <- 6
  n_cond <- 4
  mat <- matrix(rnorm(n_nodes * n_cond), nrow = n_nodes, ncol = n_cond)

  exclusions <- lapply(seq_len(n_nodes), function(i) as.integer(setdiff(seq_len(n_nodes), c(i, ((i %% n_nodes) + 1L)))))
  sources <- lapply(seq_len(n_nodes), function(i) setdiff(seq_len(n_nodes), unique(c(i, exclusions[[i]]))))

  out_sparse_cpp <- noncircular_activity_sparse_cpp(mat, sources, fill_value = 0)
  out_dense_cpp <- noncircular_activity_dense_cpp(mat, exclusions, fill_value = 0)

  out_r_sparse <- calcactivity_parcelwise_noncircular(
    mat,
    parcelstoexclude_bytarget = exclusions,
    orientation = "nodes_by_conditions",
    fill_value = 0,
    sparse = TRUE,
    mask_format = "list"
  )

  attr(out_r_sparse, "source_mask_format") <- NULL
  attr(out_r_sparse, "sparsity_stats") <- NULL

  expect_equal(out_sparse_cpp, out_r_sparse, tolerance = 1e-12)
  expect_equal(out_dense_cpp, out_r_sparse, tolerance = 1e-12)

  expect_error(noncircular_activity_sparse_cpp(mat, sources[1:2], fill_value = 0), "one entry per target node")
  expect_error(noncircular_activity_dense_cpp(mat, exclusions[1:2], fill_value = 0), "one entry per target node")
})
