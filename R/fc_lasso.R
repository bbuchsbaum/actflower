.af_zscore_rows <- function(x) {
  m <- as.matrix(x)
  mu <- rowMeans(m, na.rm = TRUE)
  sdv <- apply(m, 1, stats::sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  (m - mu) / sdv
}

.af_make_folds <- function(n_time, k_folds = 10, scheme = c("blocked", "interleaving")) {
  scheme <- match.arg(scheme)
  k_folds <- max(2L, min(as.integer(k_folds), n_time))
  fold_ids <- integer(n_time)

  if (scheme == "blocked") {
    tr_per_fold <- n_time / k_folds
    t1 <- 1L
    for (k in seq_len(k_folds)) {
      t2 <- as.integer(round(k * tr_per_fold))
      t2 <- min(max(t2, t1), n_time)
      fold_ids[t1:t2] <- k
      t1 <- t2 + 1L
      if (t1 > n_time) break
    }
    fold_ids[fold_ids == 0] <- k_folds
  } else {
    fold_ids <- rep_len(seq_len(k_folds), n_time)
  }

  fold_ids
}

.af_lasso_fit_all_targets <- function(data_nodes_by_time, lambda_grid, standardize = TRUE) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    actflower_abort("Install 'glmnet' to use lasso estimators.")
  }

  data <- as.matrix(data_nodes_by_time)
  n_nodes <- nrow(data)
  n_lambda <- length(lambda_grid)

  conn <- array(0, dim = c(n_nodes, n_nodes, n_lambda))
  all_nodes <- seq_len(n_nodes)

  for (target in all_nodes) {
    sources <- all_nodes[all_nodes != target]
    x <- t(data[sources, , drop = FALSE])
    y <- as.numeric(data[target, ])

    src_sd <- apply(x, 2, stats::sd)
    keep <- is.finite(src_sd) & src_sd > 0
    if (!any(keep)) {
      next
    }
    x_fit <- x[, keep, drop = FALSE]
    sources_fit <- sources[keep]

    fit <- glmnet::glmnet(
      x = x_fit,
      y = y,
      alpha = 1,
      lambda = lambda_grid,
      standardize = standardize,
      intercept = TRUE
    )

    beta <- as.matrix(fit$beta)
    for (l in seq_len(n_lambda)) {
      conn[target, sources_fit, l] <- beta[, l]
    }
  }

  conn
}

.af_activity_prediction_scores <- function(activity_nodes_by_time, conn, nodewise = TRUE) {
  activity <- as.matrix(activity_nodes_by_time)
  if (length(dim(conn)) == 2) {
    conn <- array(conn, dim = c(nrow(conn), ncol(conn), 1L))
  }

  n_nodes <- nrow(activity)
  n_models <- dim(conn)[3]
  pred <- array(0, dim = c(n_nodes, ncol(activity), n_models))
  nodes <- seq_len(n_nodes)

  for (n in nodes) {
    src <- nodes[nodes != n]
    x <- activity[src, , drop = FALSE]
    for (m in seq_len(n_models)) {
      beta <- conn[n, src, m]
      pred[n, , m] <- colSums(x * beta)
    }
  }

  if (nodewise) {
    r2 <- matrix(NA_real_, nrow = n_nodes, ncol = n_models)
    pear <- matrix(NA_real_, nrow = n_nodes, ncol = n_models)
    for (m in seq_len(n_models)) {
      for (n in seq_len(n_nodes)) {
        r2[n, m] <- r2_vec(activity[n, ], pred[n, , m])
        pear[n, m] <- pearson_vec(activity[n, ], pred[n, , m])
      }
    }
  } else {
    r2 <- matrix(NA_real_, nrow = 1, ncol = n_models)
    pear <- matrix(NA_real_, nrow = 1, ncol = n_models)
    for (m in seq_len(n_models)) {
      r2[1, m] <- r2_vec(as.vector(activity), as.vector(pred[, , m]))
      pear[1, m] <- pearson_vec(as.vector(activity), as.vector(pred[, , m]))
    }
  }

  list(R2 = r2, pearson = pear)
}

#' Compute lasso FC with fixed lambda
#' @param x Numeric matrix.
#' @param lambda Lambda scalar or vector.
#' @param orientation Matrix orientation.
#' @param zscore Whether to z-score node time series before fitting.
#' @return FC matrix or FC array with dimensions (nodes x nodes x n_lambda).
#' @export
estimate_fc_lasso <- function(
  x,
  lambda,
  orientation = c("nodes_by_time", "time_by_nodes"),
  zscore = TRUE
) {
  data <- as_nodes_by_time(x, orientation = orientation)
  data <- if (zscore) .af_zscore_rows(data) else data

  lambda_grid <- as.numeric(lambda)
  if (length(lambda_grid) < 1 || any(!is.finite(lambda_grid)) || any(lambda_grid <= 0)) {
    actflower_abort("`lambda` must contain positive finite values.")
  }

  conn <- .af_lasso_fit_all_targets(data, lambda_grid = lambda_grid, standardize = FALSE)
  if (length(lambda_grid) == 1L) {
    return(conn[, , 1L, drop = TRUE])
  }
  conn
}

#' Cross-validated lasso FC
#' @param x Numeric matrix.
#' @param lambda Optional lambda grid.
#' @param k_folds Number of CV folds.
#' @param opt_method Optimization method (`"R2"` or `NULL`).
#' @param orientation Matrix orientation.
#' @param nodewise_hyperparams Select lambda per target node.
#' @param folds_scheme Fold assignment strategy.
#' @param zscore Whether to z-score node time series before fitting and scoring.
#' @return List with `fc` matrix and `cv` metadata.
#' @export
estimate_fc_lasso_cv <- function(
  x,
  lambda = NULL,
  k_folds = 10,
  opt_method = "R2",
  orientation = c("nodes_by_time", "time_by_nodes"),
  nodewise_hyperparams = TRUE,
  folds_scheme = c("blocked", "interleaving"),
  zscore = TRUE
) {
  data <- as_nodes_by_time(x, orientation = orientation)
  folds_scheme <- match.arg(folds_scheme)

  if (is.null(lambda)) {
    lambda <- round(10 ^ seq(-0.5, -3.0, by = -0.1), 6)
  }
  lambda_grid <- as.numeric(lambda)
  if (length(lambda_grid) < 1 || any(!is.finite(lambda_grid)) || any(lambda_grid <= 0)) {
    actflower_abort("`lambda` must contain positive finite values.")
  }

  if (!is.null(opt_method) && !identical(opt_method, "R2")) {
    actflower_abort("`opt_method` must be 'R2' or NULL.")
  }

  data <- if (zscore) .af_zscore_rows(data) else data
  n_nodes <- nrow(data)
  n_time <- ncol(data)

  fold_ids <- .af_make_folds(n_time, k_folds = k_folds, scheme = folds_scheme)
  k_folds_eff <- max(fold_ids)

  if (nodewise_hyperparams) {
    scores <- array(NA_real_, dim = c(length(lambda_grid), k_folds_eff, n_nodes))
  } else {
    scores <- array(NA_real_, dim = c(length(lambda_grid), k_folds_eff, 1L))
  }

  for (k in seq_len(k_folds_eff)) {
    te <- fold_ids == k
    tr <- !te

    conn_k <- .af_lasso_fit_all_targets(data[, tr, drop = FALSE], lambda_grid = lambda_grid, standardize = FALSE)

    if (identical(opt_method, "R2")) {
      sc <- .af_activity_prediction_scores(data[, te, drop = FALSE], conn_k, nodewise = nodewise_hyperparams)$R2
      if (nodewise_hyperparams) {
        scores[, k, ] <- t(sc)
      } else {
        scores[, k, 1L] <- as.vector(sc)
      }
    }
  }

  score_mean <- apply(scores, c(1, 3), function(v) mean(v, na.rm = TRUE))

  if (nodewise_hyperparams) {
    best_idx <- apply(score_mean, 2, which.max)
    best_param <- lambda_grid[best_idx]
    conn_full_all <- .af_lasso_fit_all_targets(data, lambda_grid = unique(best_param), standardize = FALSE)

    fc <- matrix(0, n_nodes, n_nodes)
    uniq <- unique(best_param)
    for (u in seq_along(uniq)) {
      lam <- uniq[u]
      nodes_here <- which(best_param == lam)
      lam_idx <- which(unique(best_param) == lam)
      fc[nodes_here, ] <- conn_full_all[nodes_here, , lam_idx]
    }
  } else {
    best_idx <- which.max(score_mean[, 1L])
    best_param <- lambda_grid[best_idx]
    fc <- estimate_fc_lasso(data, lambda = best_param, orientation = "nodes_by_time", zscore = FALSE)
  }

  list(
    fc = fc,
    cv = list(
      bestParam = best_param,
      L1s = lambda_grid,
      R2 = scores,
      fold_ids = fold_ids,
      nodewise = nodewise_hyperparams
    )
  )
}
