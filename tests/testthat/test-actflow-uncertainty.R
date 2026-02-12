test_that("actflow_uncertainty returns expected structure and metrics", {
  set.seed(81)
  n_nodes <- 6
  n_cond <- 5
  n_subj <- 8

  act <- array(rnorm(n_nodes * n_cond * n_subj), dim = c(n_nodes, n_cond, n_subj))
  fc <- array(0, dim = c(n_nodes, n_nodes, n_subj))
  for (s in seq_len(n_subj)) {
    m <- matrix(rnorm(n_nodes * n_nodes), n_nodes, n_nodes)
    diag(m) <- 0
    fc[, , s] <- m
  }

  out <- actflow_uncertainty(
    act,
    fc,
    metric = c("corr", "R2", "mae"),
    n_boot = 120,
    conf_level = 0.95,
    seed = 1001
  )

  expect_true(all(c("point", "interval", "distribution_summary", "bootstrap_config", "bootstrap_draws") %in% names(out)))
  expect_equal(sort(names(out$point)), c("R2", "corr", "mae"))
  expect_equal(nrow(out$interval), 3)
  expect_equal(nrow(out$distribution_summary), 3)
  expect_equal(dim(out$bootstrap_draws), c(120, 3))
  expect_true(all(is.finite(out$point)))

  # Point estimates should align with the direct fullcompare output on the full sample.
  ref <- actflow_test(act, fc, comparison = "fullcompare_compthenavg", use_cpp = TRUE)
  expect_equal(unname(out$point["corr"]), mean(ref$model_compare_output$corr_vals), tolerance = 1e-10)
  expect_equal(unname(out$point["R2"]), mean(ref$model_compare_output$R2_vals), tolerance = 1e-10)
  expect_equal(unname(out$point["mae"]), mean(ref$model_compare_output$mae_vals), tolerance = 1e-10)
})

test_that("actflow_uncertainty is reproducible with seed", {
  set.seed(82)
  n_nodes <- 5
  n_cond <- 4
  n_subj <- 7

  act <- array(rnorm(n_nodes * n_cond * n_subj), dim = c(n_nodes, n_cond, n_subj))
  fc <- array(rnorm(n_nodes * n_nodes * n_subj), dim = c(n_nodes, n_nodes, n_subj))
  for (s in seq_len(n_subj)) diag(fc[, , s]) <- 0

  out_a <- actflow_uncertainty(act, fc, n_boot = 80, seed = 2026)
  out_b <- actflow_uncertainty(act, fc, n_boot = 80, seed = 2026)
  out_c <- actflow_uncertainty(act, fc, n_boot = 80, seed = 2027)

  expect_equal(out_a$bootstrap_draws, out_b$bootstrap_draws, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(out_a$bootstrap_draws, out_c$bootstrap_draws, tolerance = 1e-12)))
})

test_that("actflow_uncertainty validates dimensions and comparison mode", {
  set.seed(83)
  act <- array(rnorm(5 * 4 * 6), dim = c(5, 4, 6))
  fc <- array(rnorm(5 * 5 * 6), dim = c(5, 5, 6))

  bad_target <- array(rnorm(5 * 4 * 5), dim = c(5, 4, 5))
  expect_error(
    actflow_uncertainty(act, fc, act_group_test = bad_target),
    "same dimensions"
  )

  expect_error(
    actflow_uncertainty(act, fc, comparison = "nodewise_compthenavg"),
    "currently supports only comparison='fullcompare_compthenavg'"
  )
})

test_that("actflow_uncertainty MAE intervals show reasonable synthetic coverage", {
  set.seed(84)
  n_rep <- 20L
  n_nodes <- 5L
  n_cond <- 6L
  n_subj <- 16L
  sigma <- 0.7
  true_mae <- sigma * sqrt(2 / pi)

  covered <- logical(n_rep)
  for (r in seq_len(n_rep)) {
    act <- array(rnorm(n_nodes * n_cond * n_subj), dim = c(n_nodes, n_cond, n_subj))
    fc <- array(rnorm(n_nodes * n_nodes * n_subj), dim = c(n_nodes, n_nodes, n_subj))
    for (s in seq_len(n_subj)) diag(fc[, , s]) <- 0

    pred <- array(0, dim = dim(act))
    for (s in seq_len(n_subj)) {
      pred[, , s] <- fc[, , s] %*% act[, , s]
    }
    noise <- array(rnorm(length(pred), sd = sigma), dim = dim(pred))
    target <- pred + noise

    out <- actflow_uncertainty(
      act_group = act,
      fc_group = fc,
      act_group_test = target,
      metric = "mae",
      n_boot = 80,
      conf_level = 0.95,
      seed = 9000 + r,
      use_cpp = TRUE
    )

    lo <- out$interval$lower[out$interval$metric == "mae"]
    hi <- out$interval$upper[out$interval$metric == "mae"]
    covered[r] <- is.finite(lo) && is.finite(hi) && lo <= true_mae && true_mae <= hi
  }

  coverage <- mean(covered)
  expect_gte(coverage, 0.75)
  expect_lte(coverage, 1.00)
})
