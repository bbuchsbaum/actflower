test_that("actflow_test supports nodewise comparison mode", {
  set.seed(41)
  n_nodes <- 6
  n_cond <- 4
  n_subj <- 3

  act <- array(rnorm(n_nodes * n_cond * n_subj), dim = c(n_nodes, n_cond, n_subj))
  fc <- array(0, dim = c(n_nodes, n_nodes, n_subj))
  for (s in seq_len(n_subj)) {
    m <- matrix(rnorm(n_nodes * n_nodes), n_nodes, n_nodes)
    diag(m) <- 0
    fc[, , s] <- m
  }

  out <- actflow_test(act, fc, comparison = "nodewise_compthenavg")
  expect_true(all(c("actPredVector_bytask_bysubj", "actVect_actual_group", "model_compare_output") %in% names(out)))
  expect_equal(dim(out$actPredVector_bytask_bysubj), c(n_nodes, n_cond, n_subj))
  expect_equal(dim(out$model_compare_output$corr_vals), c(n_cond, n_subj))
})

test_that("actflow_test use_cpp path matches R path for fullcompare", {
  set.seed(42)
  n_nodes <- 5
  n_cond <- 3
  n_subj <- 2

  act <- array(rnorm(n_nodes * n_cond * n_subj), dim = c(n_nodes, n_cond, n_subj))
  fc <- array(0, dim = c(n_nodes, n_nodes, n_subj))
  for (s in seq_len(n_subj)) {
    m <- matrix(rnorm(n_nodes * n_nodes), n_nodes, n_nodes)
    diag(m) <- 0
    fc[, , s] <- m
  }

  out_r <- actflow_test(act, fc, comparison = "fullcompare_compthenavg", use_cpp = FALSE)
  out_cpp <- actflow_test(act, fc, comparison = "fullcompare_compthenavg", use_cpp = TRUE)

  expect_equal(out_cpp$actPredVector_bytask_bysubj, out_r$actPredVector_bytask_bysubj, tolerance = 1e-10)
  expect_equal(out_cpp$model_compare_output$corr_vals, out_r$model_compare_output$corr_vals, tolerance = 1e-10)
})

test_that("actflow_test fused fullcompare path respects act_group_test target", {
  set.seed(43)
  n_nodes <- 6
  n_cond <- 4
  n_subj <- 3

  act <- array(rnorm(n_nodes * n_cond * n_subj), dim = c(n_nodes, n_cond, n_subj))
  act_test <- act + 0.05

  fc <- array(0, dim = c(n_nodes, n_nodes, n_subj))
  for (s in seq_len(n_subj)) {
    m <- matrix(rnorm(n_nodes * n_nodes), n_nodes, n_nodes)
    diag(m) <- 0
    fc[, , s] <- m
  }

  out_r <- actflow_test(
    act,
    fc,
    act_group_test = act_test,
    comparison = "fullcompare_compthenavg",
    use_cpp = FALSE
  )
  out_cpp <- actflow_test(
    act,
    fc,
    act_group_test = act_test,
    comparison = "fullcompare_compthenavg",
    use_cpp = TRUE
  )

  expect_equal(out_cpp$actVect_actual_group, act_test, tolerance = 1e-12)
  expect_equal(out_cpp$actPredVector_bytask_bysubj, out_r$actPredVector_bytask_bysubj, tolerance = 1e-10)
  expect_equal(out_cpp$model_compare_output$corr_vals, out_r$model_compare_output$corr_vals, tolerance = 1e-10)
  expect_equal(out_cpp$model_compare_output$R2_vals, out_r$model_compare_output$R2_vals, tolerance = 1e-10)
  expect_equal(out_cpp$model_compare_output$mae_vals, out_r$model_compare_output$mae_vals, tolerance = 1e-10)
})
