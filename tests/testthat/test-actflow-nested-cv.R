test_that("actflow_nested_cv returns expected artifacts and is reproducible", {
  set.seed(91)
  n_nodes <- 5
  n_time <- 80
  n_cond <- 4
  n_subj <- 9

  rest <- array(rnorm(n_nodes * n_time * n_subj), dim = c(n_nodes, n_time, n_subj))
  task <- array(rnorm(n_nodes * n_cond * n_subj), dim = c(n_nodes, n_cond, n_subj))

  grid <- list(
    list(method = "multreg", params = list(ridge = 0)),
    list(method = "multreg", params = list(ridge = 0.2)),
    list(method = "corr", params = list())
  )

  out_a <- actflow_nested_cv(
    data = list(rest_group = rest, task_group = task),
    outer_folds = 3,
    inner_folds = 2,
    estimator_grid = grid,
    seed = 2026,
    save_splits = TRUE,
    selection_metric = "corr"
  )
  out_b <- actflow_nested_cv(
    data = list(rest_group = rest, task_group = task),
    outer_folds = 3,
    inner_folds = 2,
    estimator_grid = grid,
    seed = 2026,
    save_splits = TRUE,
    selection_metric = "corr"
  )

  expect_true(all(c("outer_metrics", "selected_hyperparams", "split_artifacts", "reproducibility_manifest") %in% names(out_a)))
  expect_equal(nrow(out_a$outer_metrics), 3)
  expect_equal(nrow(out_a$selected_hyperparams), 3)
  expect_equal(length(out_a$split_artifacts$outer), 3)
  expect_true(all(is.finite(out_a$outer_metrics$corr)))

  expect_equal(out_a$outer_metrics, out_b$outer_metrics, tolerance = 1e-12)
  expect_equal(out_a$selected_hyperparams$candidate_id, out_b$selected_hyperparams$candidate_id)
  expect_equal(out_a$split_artifacts$outer_fold_ids, out_b$split_artifacts$outer_fold_ids)
})

test_that("actflow_nested_cv outer-test perturbation does not alter selected hyperparams", {
  set.seed(92)
  n_nodes <- 6
  n_time <- 90
  n_cond <- 5
  n_subj <- 12

  rest <- array(rnorm(n_nodes * n_time * n_subj), dim = c(n_nodes, n_time, n_subj))
  task <- array(rnorm(n_nodes * n_cond * n_subj), dim = c(n_nodes, n_cond, n_subj))
  target <- task

  grid <- list(
    list(method = "multreg", params = list(ridge = 0)),
    list(method = "multreg", params = list(ridge = 0.6))
  )

  base <- actflow_nested_cv(
    data = list(rest_group = rest, task_group = task, task_group_test = target),
    outer_folds = 4,
    inner_folds = 2,
    estimator_grid = grid,
    seed = 3030,
    save_splits = TRUE,
    selection_metric = "corr"
  )

  test_idx_fold1 <- base$split_artifacts$outer[[1]]$test_indices
  target_perturbed <- target
  target_perturbed[, , test_idx_fold1] <- array(
    rnorm(length(target_perturbed[, , test_idx_fold1]), sd = 25),
    dim = dim(target_perturbed[, , test_idx_fold1])
  )

  pert <- actflow_nested_cv(
    data = list(rest_group = rest, task_group = task, task_group_test = target_perturbed),
    outer_folds = 4,
    inner_folds = 2,
    estimator_grid = grid,
    seed = 3030,
    save_splits = TRUE,
    selection_metric = "corr"
  )

  # Selection should not change because outer test data are excluded from model selection.
  expect_equal(base$selected_hyperparams$candidate_id, pert$selected_hyperparams$candidate_id)
  expect_equal(base$selected_hyperparams$method, pert$selected_hyperparams$method)

  # Outer-fold performance for perturbed fold should change substantially.
  c1_base <- base$outer_metrics$corr[base$outer_metrics$outer_fold == 1]
  c1_pert <- pert$outer_metrics$corr[pert$outer_metrics$outer_fold == 1]
  expect_gt(abs(c1_base - c1_pert), 0.1)
})

test_that("actflow_nested_cv validates data and estimator grid", {
  set.seed(93)
  rest <- array(rnorm(5 * 40 * 8), dim = c(5, 40, 8))
  task <- array(rnorm(5 * 3 * 8), dim = c(5, 3, 8))

  expect_error(
    actflow_nested_cv(
      data = list(task_group = task),
      outer_folds = 3,
      inner_folds = 2,
      estimator_grid = list(list(method = "multreg"))
    ),
    "must include `rest_group` and `task_group`"
  )

  expect_error(
    actflow_nested_cv(
      data = list(rest_group = rest, task_group = task),
      outer_folds = 3,
      inner_folds = 2,
      estimator_grid = list(list(method = "unknown"))
    ),
    "Unsupported estimator method"
  )
})
