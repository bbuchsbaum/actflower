.af_glasso_parcorr <- function(data_nodes_by_time, lambda, zscore = TRUE) {
  if (!requireNamespace("glasso", quietly = TRUE)) {
    actflower_abort("Install 'glasso' to use graphical lasso estimators.")
  }

  data <- as.matrix(data_nodes_by_time)
  data <- if (zscore) .af_zscore_rows(data) else data

  emp_cov <- stats::cov(t(data))
  fit <- glasso::glasso((emp_cov + t(emp_cov)) * 0.5, rho = lambda)
  prec <- (fit$wi + t(fit$wi)) * 0.5

  d <- sqrt(pmax(diag(prec), .Machine$double.eps))
  parcorr <- -prec / (d %o% d)
  diag(parcorr) <- 0

  list(parcorr = parcorr, precision = prec)
}

.af_glasso_neg_loglik <- function(emp_cov, precision) {
  s <- (emp_cov + t(emp_cov)) * 0.5
  p <- (precision + t(precision)) * 0.5
  eig <- eigen(p, symmetric = TRUE, only.values = TRUE)$values
  eig[eig <= .Machine$double.eps] <- .Machine$double.eps
  logdet <- sum(log(eig))
  tr <- sum(s * p)
  # Constant omitted (does not affect optimization over lambda grid)
  0.5 * (tr - logdet)
}

#' Compute graphical-lasso partial-correlation FC
#' @param x Numeric matrix.
#' @param lambda Positive scalar lambda.
#' @param orientation Matrix orientation.
#' @param zscore Whether to z-score node time series before estimation.
#' @return Partial-correlation FC matrix.
#' @export
estimate_fc_glasso <- function(
  x,
  lambda,
  orientation = c("nodes_by_time", "time_by_nodes"),
  zscore = TRUE
) {
  data <- as_nodes_by_time(x, orientation = orientation)
  if (!is.numeric(lambda) || length(lambda) != 1 || !is.finite(lambda) || lambda <= 0) {
    actflower_abort("`lambda` must be a positive scalar.")
  }
  .af_glasso_parcorr(data, lambda = lambda, zscore = zscore)$parcorr
}

#' Cross-validated graphical-lasso FC
#' @param x Numeric matrix.
#' @param lambda Optional lambda grid.
#' @param k_folds Number of CV folds.
#' @param opt_method Optimization method (`"loglikelihood"` or `"R2"`).
#' @param orientation Matrix orientation.
#' @param folds_scheme Fold assignment strategy.
#' @param zscore Whether to z-score node time series before fitting and scoring.
#' @return List with `fc` matrix and `cv` metadata.
#' @export
estimate_fc_glasso_cv <- function(
  x,
  lambda = NULL,
  k_folds = 10,
  opt_method = c("loglikelihood", "R2"),
  orientation = c("nodes_by_time", "time_by_nodes"),
  folds_scheme = c("blocked", "interleaving"),
  zscore = TRUE
) {
  opt_method <- match.arg(opt_method)
  folds_scheme <- match.arg(folds_scheme)

  data <- as_nodes_by_time(x, orientation = orientation)
  data <- if (zscore) .af_zscore_rows(data) else data

  if (is.null(lambda)) {
    lambda <- 10 ^ seq(-0.5, -3.0, by = -0.1)
  }
  lambda_grid <- as.numeric(lambda)
  if (length(lambda_grid) < 1 || any(!is.finite(lambda_grid)) || any(lambda_grid <= 0)) {
    actflower_abort("`lambda` must contain positive finite values.")
  }

  n_time <- ncol(data)
  fold_ids <- .af_make_folds(n_time, k_folds = k_folds, scheme = folds_scheme)
  k_folds_eff <- max(fold_ids)

  scores <- matrix(NA_real_, nrow = length(lambda_grid), ncol = k_folds_eff)

  for (k in seq_len(k_folds_eff)) {
    te <- fold_ids == k
    tr <- !te

    for (l in seq_along(lambda_grid)) {
      est <- .af_glasso_parcorr(data[, tr, drop = FALSE], lambda = lambda_grid[l], zscore = FALSE)

      if (opt_method == "loglikelihood") {
        emp_cov_test <- stats::cov(t(data[, te, drop = FALSE]))
        scores[l, k] <- .af_glasso_neg_loglik(emp_cov_test, est$precision)
      } else {
        sc <- .af_activity_prediction_scores(data[, te, drop = FALSE], est$parcorr, nodewise = FALSE)$R2
        scores[l, k] <- as.numeric(sc[1, 1])
      }
    }
  }

  mean_scores <- rowMeans(scores, na.rm = TRUE)
  best_idx <- if (opt_method == "loglikelihood") which.min(mean_scores) else which.max(mean_scores)
  best_param <- lambda_grid[best_idx]

  fc <- .af_glasso_parcorr(data, lambda = best_param, zscore = FALSE)$parcorr

  list(
    fc = fc,
    cv = list(
      bestParam = best_param,
      L1s = lambda_grid,
      metric = opt_method,
      scores = scores,
      fold_ids = fold_ids
    )
  )
}
