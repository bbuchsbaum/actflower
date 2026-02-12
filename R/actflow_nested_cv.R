#' Leakage-proof nested CV for actflow model selection
#' @param data A list with `rest_group` and `task_group` arrays, each 3D.
#'   Optional `task_group_test` can be supplied for evaluation targets.
#' @param outer_folds Number of outer subject folds.
#' @param inner_folds Number of inner subject folds used for model selection.
#' @param estimator_grid Candidate estimators/hyperparameters as a list or data frame.
#'   For list input, each candidate should contain a `method` field and optional
#'   `params` list. For data frame input, include a `method` column and parameter columns.
#' @param seed Optional RNG seed for reproducible folds and bootstrap-like sampling.
#' @param save_splits If `TRUE`, include detailed fold assignments in output.
#' @param artifact_path Optional path to persist nested-CV artifacts as JSON.
#' @param selection_metric Metric used to select the winning candidate in inner CV.
#'   One of `corr`, `R2`, or `mae`.
#' @param comparison Comparison mode passed to [actflow_test()]. Currently
#'   nested CV supports `fullcompare_compthenavg`.
#' @param transfer Transfer function passed to [actflow_test()].
#' @param use_cpp Use native kernels when available.
#' @param ... Reserved for forward compatibility.
#' @return A list with outer metrics, selected hyperparameters, split artifacts,
#'   and a reproducibility manifest.
#' @export
actflow_nested_cv <- function(
  data,
  outer_folds = 5,
  inner_folds = 5,
  estimator_grid,
  seed = NULL,
  save_splits = TRUE,
  artifact_path = NULL,
  selection_metric = c("corr", "R2", "mae"),
  comparison = "fullcompare_compthenavg",
  transfer = c("linear", "relu", "sigmoid", "logit"),
  use_cpp = TRUE,
  ...
) {
  transfer <- match.arg(transfer)
  selection_metric <- match.arg(selection_metric)

  parsed <- .af_parse_nested_data(data)
  rest_group <- parsed$rest_group
  task_group <- parsed$task_group
  task_group_test <- parsed$task_group_test

  if (!identical(comparison, "fullcompare_compthenavg")) {
    actflower_abort("`actflow_nested_cv()` currently supports only comparison='fullcompare_compthenavg'.")
  }

  n_subj <- dim(task_group)[3]
  if (!is.numeric(outer_folds) || length(outer_folds) != 1L || !is.finite(outer_folds)) {
    actflower_abort("`outer_folds` must be a finite scalar.")
  }
  if (!is.numeric(inner_folds) || length(inner_folds) != 1L || !is.finite(inner_folds)) {
    actflower_abort("`inner_folds` must be a finite scalar.")
  }

  k_outer <- max(2L, min(as.integer(outer_folds), n_subj))
  k_inner_req <- max(2L, as.integer(inner_folds))

  candidates <- .af_normalize_estimator_grid(estimator_grid)
  n_candidates <- length(candidates)

  outer_fold_ids <- .af_subject_fold_ids(n_subj, k = k_outer, seed = seed)

  outer_metrics_rows <- vector("list", k_outer)
  selected_rows <- vector("list", k_outer)
  split_artifacts <- if (isTRUE(save_splits)) vector("list", k_outer) else NULL

  for (fo in seq_len(k_outer)) {
    test_idx <- which(outer_fold_ids == fo)
    train_idx <- setdiff(seq_len(n_subj), test_idx)

    if (length(train_idx) < 2L) {
      actflower_abort("Not enough training subjects for inner CV. Reduce `outer_folds`.")
    }

    k_inner <- max(2L, min(k_inner_req, length(train_idx)))
    inner_seed <- if (is.null(seed)) NULL else as.integer(seed) + fo
    inner_fold_ids_local <- .af_subject_fold_ids(length(train_idx), k = k_inner, seed = inner_seed)
    inner_fold_ids <- setNames(inner_fold_ids_local, train_idx)

    inner_scores <- list(
      corr = matrix(NA_real_, nrow = n_candidates, ncol = k_inner),
      R2 = matrix(NA_real_, nrow = n_candidates, ncol = k_inner),
      mae = matrix(NA_real_, nrow = n_candidates, ncol = k_inner)
    )
    for (nm in names(inner_scores)) {
      rownames(inner_scores[[nm]]) <- paste0("candidate_", seq_len(n_candidates))
      colnames(inner_scores[[nm]]) <- paste0("inner_fold_", seq_len(k_inner))
    }

    for (ci in seq_len(n_candidates)) {
      cand <- candidates[[ci]]

      for (fi in seq_len(k_inner)) {
        val_local <- which(inner_fold_ids_local == fi)
        val_idx <- train_idx[val_local]

        eval_out <- .af_eval_candidate_subjects(
          rest_group = rest_group,
          task_group = task_group,
          task_group_test = task_group_test,
          subject_idx = val_idx,
          candidate = cand,
          comparison = comparison,
          transfer = transfer,
          use_cpp = use_cpp
        )
        inner_scores$corr[ci, fi] <- eval_out$summary$corr
        inner_scores$R2[ci, fi] <- eval_out$summary$R2
        inner_scores$mae[ci, fi] <- eval_out$summary$mae
      }
    }

    mean_scores <- rowMeans(inner_scores[[selection_metric]], na.rm = TRUE)
    best_ci <- if (identical(selection_metric, "mae")) {
      unname(as.integer(which.min(mean_scores)))
    } else {
      unname(as.integer(which.max(mean_scores)))
    }
    best <- candidates[[best_ci]]

    outer_eval <- .af_eval_candidate_subjects(
      rest_group = rest_group,
      task_group = task_group,
      task_group_test = task_group_test,
      subject_idx = test_idx,
      candidate = best,
      comparison = comparison,
      transfer = transfer,
      use_cpp = use_cpp
    )

    outer_metrics_rows[[fo]] <- data.frame(
      outer_fold = fo,
      n_test_subjects = length(test_idx),
      corr = outer_eval$summary$corr,
      R2 = outer_eval$summary$R2,
      mae = outer_eval$summary$mae,
      selected_candidate = best_ci,
      selection_metric = selection_metric,
      selection_score = mean_scores[best_ci],
      stringsAsFactors = FALSE
    )

    selected_rows[[fo]] <- list(
      outer_fold = fo,
      candidate_id = best_ci,
      method = best$method,
      params = best$params,
      selection_metric = selection_metric,
      selection_score = mean_scores[best_ci]
    )

    if (isTRUE(save_splits)) {
      split_artifacts[[fo]] <- list(
        outer_fold = fo,
        train_indices = train_idx,
        test_indices = test_idx,
        inner_fold_ids = as.integer(inner_fold_ids),
        inner_scores = inner_scores,
        inner_mean_scores = list(
          corr = rowMeans(inner_scores$corr, na.rm = TRUE),
          R2 = rowMeans(inner_scores$R2, na.rm = TRUE),
          mae = rowMeans(inner_scores$mae, na.rm = TRUE)
        ),
        selected_candidate = best_ci,
        outer_test_metrics = outer_eval$metrics
      )
      names(split_artifacts[[fo]]$inner_fold_ids) <- names(inner_fold_ids)
    }
  }

  outer_metrics <- do.call(rbind, outer_metrics_rows)
  selected_hyperparams <- .af_selected_rows_to_df(selected_rows)

  manifest <- list(
    seed = if (is.null(seed)) NA_integer_ else as.integer(seed),
    outer_folds = k_outer,
    inner_folds = k_inner_req,
    selection_metric = selection_metric,
    comparison = comparison,
    transfer = transfer,
    use_cpp = isTRUE(use_cpp),
    n_subjects = n_subj,
    n_nodes = dim(task_group)[1],
    n_conditions = dim(task_group)[2],
    estimator_grid_signature = .af_grid_signature(candidates),
    created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )

  out <- list(
    outer_metrics = outer_metrics,
    selected_hyperparams = selected_hyperparams,
    split_artifacts = if (isTRUE(save_splits)) {
      list(
        outer_fold_ids = outer_fold_ids,
        estimator_grid = candidates,
        outer = split_artifacts
      )
    } else {
      NULL
    },
    reproducibility_manifest = manifest
  )

  if (!is.null(artifact_path)) {
    .af_write_nested_cv_artifacts(out, artifact_path)
    out$reproducibility_manifest$artifact_path <- normalizePath(artifact_path, winslash = "/", mustWork = FALSE)
  }

  out
}

.af_parse_nested_data <- function(data) {
  if (!is.list(data)) {
    actflower_abort("`data` must be a list with `rest_group` and `task_group`.")
  }
  if (is.null(data$rest_group) || is.null(data$task_group)) {
    actflower_abort("`data` must include `rest_group` and `task_group`.")
  }

  rest_group <- data$rest_group
  task_group <- data$task_group
  task_group_test <- if (is.null(data$task_group_test)) task_group else data$task_group_test

  if (!is.numeric(rest_group) || length(dim(rest_group)) != 3L) {
    actflower_abort("`data$rest_group` must be a numeric 3D array [nodes x time x subjects].")
  }
  validate_array3(task_group, arg = "data$task_group")
  validate_array3(task_group_test, arg = "data$task_group_test")

  if (dim(rest_group)[1] != dim(task_group)[1]) {
    actflower_abort("`rest_group` and `task_group` must share node dimension.")
  }
  if (dim(rest_group)[3] != dim(task_group)[3]) {
    actflower_abort("`rest_group` and `task_group` must share subject dimension.")
  }
  if (!all(dim(task_group_test) == dim(task_group))) {
    actflower_abort("`task_group_test` must have the same dimensions as `task_group`.")
  }

  list(
    rest_group = rest_group,
    task_group = task_group,
    task_group_test = task_group_test
  )
}

.af_subject_fold_ids <- function(n_subjects, k = 5, seed = NULL) {
  k <- max(2L, min(as.integer(k), n_subjects))
  if (is.null(seed)) {
    perm <- sample.int(n_subjects, n_subjects, replace = FALSE)
  } else {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(seed))
    perm <- sample.int(n_subjects, n_subjects, replace = FALSE)
  }

  out <- integer(n_subjects)
  out[perm] <- rep_len(seq_len(k), n_subjects)
  out
}

.af_normalize_estimator_grid <- function(estimator_grid) {
  if (is.data.frame(estimator_grid)) {
    if (!"method" %in% names(estimator_grid)) {
      actflower_abort("When `estimator_grid` is a data.frame it must include a `method` column.")
    }
    out <- vector("list", nrow(estimator_grid))
    for (i in seq_len(nrow(estimator_grid))) {
      row <- estimator_grid[i, , drop = FALSE]
      method <- as.character(row$method[[1]])
      params <- as.list(row[setdiff(names(row), "method")])
      params <- params[!vapply(params, function(x) is.null(x) || (length(x) == 1L && is.na(x)), logical(1))]
      out[[i]] <- list(method = method, params = params)
    }
    return(.af_validate_candidates(out))
  }

  if (!is.list(estimator_grid) || length(estimator_grid) < 1L) {
    actflower_abort("`estimator_grid` must be a non-empty list or data.frame.")
  }

  if (!is.null(estimator_grid$method) && is.character(estimator_grid$method)) {
    estimator_grid <- list(estimator_grid)
  }

  out <- vector("list", length(estimator_grid))
  for (i in seq_along(estimator_grid)) {
    cand <- estimator_grid[[i]]
    if (!is.list(cand) || is.null(cand$method)) {
      actflower_abort("Each estimator candidate must be a list containing `method`.")
    }
    method <- as.character(cand$method[[1]])
    if (!is.null(cand$params)) {
      params <- cand$params
      if (!is.list(params)) {
        actflower_abort("Candidate `params` must be a list when provided.")
      }
    } else {
      params <- cand[setdiff(names(cand), "method")]
      if (!length(params)) params <- list()
    }
    out[[i]] <- list(method = method, params = params)
  }

  .af_validate_candidates(out)
}

.af_validate_candidates <- function(candidates) {
  valid_methods <- c("corr", "multreg", "partial", "combined", "lasso_cv", "glasso_cv")
  for (i in seq_along(candidates)) {
    m <- candidates[[i]]$method
    if (!is.character(m) || length(m) != 1L || !nzchar(m)) {
      actflower_abort("Each candidate must define a non-empty scalar `method`.")
    }
    if (!(m %in% valid_methods)) {
      actflower_abort(sprintf("Unsupported estimator method in grid: %s", m))
    }
    if (!is.list(candidates[[i]]$params)) {
      actflower_abort("Each candidate `params` must be a list.")
    }
  }
  candidates
}

.af_estimate_fc_for_subjects <- function(rest_group, subject_idx, candidate, use_cpp = TRUE) {
  n_nodes <- dim(rest_group)[1]
  out <- array(NA_real_, dim = c(n_nodes, n_nodes, length(subject_idx)))

  for (i in seq_along(subject_idx)) {
    s <- subject_idx[[i]]
    args <- c(
      list(
        x = rest_group[, , s, drop = TRUE],
        method = candidate$method,
        orientation = "nodes_by_time"
      ),
      candidate$params
    )

    if (candidate$method %in% c("corr", "multreg") && is.null(args$use_cpp)) {
      args$use_cpp <- use_cpp
    }

    est <- do.call(estimate_fc, args)
    fc <- if (is.list(est) && !is.null(est$fc)) est$fc else est
    validate_square_matrix(fc, arg = "candidate FC")
    out[, , i] <- as.matrix(fc)
  }

  out
}

.af_eval_candidate_subjects <- function(
  rest_group,
  task_group,
  task_group_test,
  subject_idx,
  candidate,
  comparison,
  transfer,
  use_cpp
) {
  fc_group <- .af_estimate_fc_for_subjects(rest_group, subject_idx, candidate = candidate, use_cpp = use_cpp)
  task <- task_group[, , subject_idx, drop = FALSE]
  target <- task_group_test[, , subject_idx, drop = FALSE]

  out <- actflow_test(
    act_group = task,
    fc_group = fc_group,
    act_group_test = target,
    comparison = comparison,
    transfer = transfer,
    use_cpp = use_cpp
  )
  cmp <- out$model_compare_output
  summary <- list(
    corr = mean(as.numeric(cmp$corr_vals), na.rm = TRUE),
    R2 = mean(as.numeric(cmp$R2_vals), na.rm = TRUE),
    mae = mean(as.numeric(cmp$mae_vals), na.rm = TRUE)
  )

  list(
    summary = summary,
    metrics = list(
      corr_vals = as.numeric(cmp$corr_vals),
      R2_vals = as.numeric(cmp$R2_vals),
      mae_vals = as.numeric(cmp$mae_vals)
    )
  )
}

.af_selected_rows_to_df <- function(rows) {
  outer_fold <- vapply(rows, `[[`, integer(1), "outer_fold")
  candidate_id <- vapply(rows, `[[`, integer(1), "candidate_id")
  method <- vapply(rows, `[[`, character(1), "method")
  selection_metric <- vapply(rows, `[[`, character(1), "selection_metric")
  selection_score <- vapply(rows, `[[`, numeric(1), "selection_score")
  params <- lapply(rows, `[[`, "params")

  data.frame(
    outer_fold = outer_fold,
    candidate_id = candidate_id,
    method = method,
    selection_metric = selection_metric,
    selection_score = selection_score,
    params = I(params),
    stringsAsFactors = FALSE
  )
}

.af_grid_signature <- function(candidates) {
  raw <- serialize(candidates, NULL, ascii = FALSE)
  paste0("bytes=", length(raw), ";sum=", sum(as.integer(raw)))
}

.af_write_nested_cv_artifacts <- function(x, path) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    actflower_abort("`artifact_path` must be a non-empty character scalar.")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    actflower_abort("Install 'jsonlite' to persist nested-CV artifacts.")
  }

  dirp <- dirname(path)
  if (!dir.exists(dirp)) {
    dir.create(dirp, recursive = TRUE, showWarnings = FALSE)
  }
  jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE)
  invisible(path)
}
