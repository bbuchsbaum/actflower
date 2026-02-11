test_that("HCP parity fixtures match corr, multreg, and combined outputs", {
  skip_if_not_installed("hdf5r")

  fixture <- testthat::test_path("..", "..", "inst", "extdata", "parity_refs", "hcp_parity_refs.h5")
  if (!file.exists(fixture)) {
    skip("Parity fixture not found; run inst/scripts/build_parity_refs.py first.")
  }

  from_py <- function(x) {
    if (is.array(x) && length(dim(x)) >= 2) {
      return(aperm(x, rev(seq_along(dim(x)))))
    }
    x
  }

  rest <- from_py(read_actflower_h5(fixture, "restdata"))
  taskbeta <- from_py(read_actflower_h5(fixture, "taskbeta"))

  fc_corr_ref <- from_py(read_actflower_h5(fixture, "fc_corr"))
  fc_multreg_ref <- from_py(read_actflower_h5(fixture, "fc_multreg"))
  fc_combined_ref <- from_py(read_actflower_h5(fixture, "fc_combined"))
  pred_corr_ref <- from_py(read_actflower_h5(fixture, "pred_corr"))
  pred_multreg_ref <- from_py(read_actflower_h5(fixture, "pred_multreg"))
  pred_combined_ref <- from_py(read_actflower_h5(fixture, "pred_combined"))

  n_subj <- dim(rest)[3]
  n_nodes <- dim(rest)[1]

  fc_corr_est <- array(0, dim = c(n_nodes, n_nodes, n_subj))
  fc_multreg_est <- array(0, dim = c(n_nodes, n_nodes, n_subj))
  fc_combined_est <- array(0, dim = c(n_nodes, n_nodes, n_subj))

  for (s in seq_len(n_subj)) {
    fc_corr_est[, , s] <- estimate_fc_corr(rest[, , s])
    fc_multreg_est[, , s] <- estimate_fc_multreg(rest[, , s])
    fc_combined_est[, , s] <- estimate_fc_combined(rest[, , s])
  }

  expect_lt(max(abs(fc_corr_est - fc_corr_ref)), 1e-10)
  expect_lt(max(abs(fc_multreg_est - fc_multreg_ref)), 1e-5)
  expect_lt(max(abs(fc_combined_est - fc_combined_ref)), 1e-5)

  out_corr <- actflow_test(taskbeta, fc_corr_est, comparison = "fullcompare_compthenavg")
  out_multreg <- actflow_test(taskbeta, fc_multreg_est, comparison = "fullcompare_compthenavg")
  out_combined <- actflow_test(taskbeta, fc_combined_est, comparison = "fullcompare_compthenavg")

  expect_lt(max(abs(out_corr$actPredVector_bytask_bysubj - pred_corr_ref)), 1e-8)
  expect_lt(max(abs(out_multreg$actPredVector_bytask_bysubj - pred_multreg_ref)), 1e-5)
  expect_lt(max(abs(out_combined$actPredVector_bytask_bysubj - pred_combined_ref)), 1e-5)

  corr_ref_metric <- read_actflower_h5(fixture, "metrics_corr/corr_vals")
  multreg_ref_metric <- read_actflower_h5(fixture, "metrics_multreg/corr_vals")
  combined_ref_metric <- read_actflower_h5(fixture, "metrics_combined/corr_vals")

  expect_lt(max(abs(out_corr$model_compare_output$corr_vals - corr_ref_metric)), 1e-8)
  expect_lt(max(abs(out_multreg$model_compare_output$corr_vals - multreg_ref_metric)), 1e-5)
  expect_lt(max(abs(out_combined$model_compare_output$corr_vals - combined_ref_metric)), 1e-5)
})
