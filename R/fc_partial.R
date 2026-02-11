#' Partial-correlation or PC-regression connectivity
#' @param x Numeric matrix.
#' @param method One of `"ridge"`, `"glasso"`, or `"pc"`.
#' @param orientation Matrix orientation.
#' @param lambda Regularization parameter (scalar or grid for CV).
#' @param k Number of components for `method = "pc"` (scalar or grid for CV).
#' @param cv_folds Number of folds for CV when grid-searching `lambda` or `k`.
#' @param cv_seed Seed for fold assignment.
#' @param penalty_weights Optional penalty matrix for glasso.
#' @return FC matrix (nodes x nodes).
#' @export
estimate_fc_partial <- function(
  x,
  method = c("ridge", "glasso", "pc"),
  orientation = c("nodes_by_time", "time_by_nodes"),
  lambda = 0.1,
  k = 5,
  cv_folds = 5,
  cv_seed = 1,
  penalty_weights = NULL
) {
  method <- match.arg(method)
  x <- as_nodes_by_time(x, orientation = orientation)
  x_time_by_nodes <- t(x)

  if (ncol(x_time_by_nodes) < 2) {
    actflower_abort("Need at least 2 nodes to estimate partial connectivity.")
  }

  if (method %in% c("ridge", "glasso")) {
    lambda_grid <- sort(unique(as.numeric(lambda)))
    if (length(lambda_grid) < 1 || any(!is.finite(lambda_grid)) || any(lambda_grid <= 0)) {
      actflower_abort("`lambda` must contain positive finite values.")
    }

    cv_meta <- NULL
    lambda_use <- lambda_grid
    if (length(lambda_grid) > 1) {
      cv_meta <- .af_pcor_cv_lambda(
        x_time_by_nodes,
        lambdas = lambda_grid,
        method = method,
        folds = cv_folds,
        seed = cv_seed,
        penalty_weights = penalty_weights
      )
      lambda_use <- cv_meta$best
    }

    r <- .af_weighted_cor(x_time_by_nodes)
    kappa <- if (method == "ridge") {
      .af_ridge_precision(r, lambda = lambda_use)
    } else {
      .af_glasso_precision(r, lambda = lambda_use, penalty_weights = penalty_weights)
    }
    p <- .af_precision_to_partial(kappa)
    diag(p) <- 0

    attr(p, "method") <- method
    attr(p, "lambda") <- lambda_use
    if (!is.null(cv_meta)) {
      attr(p, "cv") <- cv_meta
    }
    return(p)
  }

  # method == "pc": returns directed coefficient matrix (target rows x source columns).
  k_grid <- sort(unique(as.integer(k)))
  k_grid <- k_grid[is.finite(k_grid) & k_grid > 0]
  if (length(k_grid) < 1) {
    actflower_abort("`k` must contain positive integers.")
  }

  cv_meta <- NULL
  k_use <- k_grid
  if (length(k_grid) > 1) {
    cv_meta <- .af_pcor_cv_pc(x_time_by_nodes, k_grid = k_grid, folds = cv_folds, seed = cv_seed)
    k_use <- cv_meta$best
  }

  b <- .af_pcor_pc_coef(x_time_by_nodes, k = k_use)
  attr(b, "method") <- method
  attr(b, "k") <- k_use
  if (!is.null(cv_meta)) {
    attr(b, "cv") <- cv_meta
  }
  b
}

.af_weighted_cor <- function(x_time_by_nodes) {
  x <- as.matrix(x_time_by_nodes)
  x <- scale(x, center = TRUE, scale = FALSE)
  s <- stats::cov(x)
  d <- sqrt(pmax(diag(s), .Machine$double.eps))
  r <- s / (d %o% d)
  r[!is.finite(r)] <- 0
  diag(r) <- 1
  r
}

.af_ridge_precision <- function(r, lambda) {
  a <- (r + t(r)) * 0.5 + lambda * diag(nrow(r))
  .af_safe_inverse_sym(a)
}

.af_glasso_precision <- function(r, lambda, penalty_weights = NULL) {
  if (!requireNamespace("glasso", quietly = TRUE)) {
    actflower_abort("Install 'glasso' to use method='glasso'.")
  }
  rho <- if (is.null(penalty_weights)) {
    lambda
  } else {
    pw <- as.matrix(penalty_weights)
    if (!all(dim(pw) == dim(r))) {
      actflower_abort("`penalty_weights` must match FC matrix dimensions.")
    }
    out <- abs(pw) * lambda
    diag(out) <- 0
    out
  }
  fit <- try(glasso::glasso((r + t(r)) * 0.5, rho = rho), silent = TRUE)
  if (inherits(fit, "try-error")) {
    fit <- glasso::glasso((r + t(r)) * 0.5, rho = lambda)
  }
  (fit$wi + t(fit$wi)) * 0.5
}

.af_precision_to_partial <- function(kappa) {
  d <- sqrt(pmax(diag(kappa), .Machine$double.eps))
  p <- -kappa / (d %o% d)
  diag(p) <- 1
  p
}

.af_precision_to_coef <- function(kappa) {
  d <- diag(kappa)
  d[abs(d) < .Machine$double.eps] <- .Machine$double.eps
  b <- -kappa %*% diag(1 / d, nrow(kappa))
  diag(b) <- 0
  b
}

.af_pcor_predictive_r2 <- function(y_true, y_hat) {
  y_true <- as.matrix(y_true)
  y_hat <- as.matrix(y_hat)
  mu <- colMeans(y_true)
  ss_tot <- colSums((y_true - rep(mu, each = nrow(y_true)))^2)
  ss_res <- colSums((y_true - y_hat)^2)
  out <- 1 - ss_res / pmax(ss_tot, .Machine$double.eps)
  out[!is.finite(out)] <- NA_real_
  out
}

.af_fold_assignments <- function(n, folds, seed = 1) {
  k <- max(2L, min(as.integer(folds), n))
  set.seed(seed)
  rep_len(sample(seq_len(k)), n)
}

.af_pcor_cv_lambda <- function(
  x_time_by_nodes,
  lambdas,
  method,
  folds = 5,
  seed = 1,
  penalty_weights = NULL
) {
  x <- as.matrix(x_time_by_nodes)
  fold_id <- .af_fold_assignments(nrow(x), folds, seed)
  scores <- numeric(length(lambdas))

  for (i in seq_along(lambdas)) {
    lam <- lambdas[i]
    fold_scores <- numeric(max(fold_id))

    for (f in seq_len(max(fold_id))) {
      idx_tr <- fold_id != f
      idx_te <- !idx_tr
      xtr <- x[idx_tr, , drop = FALSE]
      xte <- x[idx_te, , drop = FALSE]

      rtr <- .af_weighted_cor(xtr)
      kappa <- if (method == "ridge") {
        .af_ridge_precision(rtr, lam)
      } else {
        .af_glasso_precision(rtr, lam, penalty_weights = penalty_weights)
      }
      b <- .af_precision_to_coef(kappa)
      yhat <- xte %*% b
      fold_scores[f] <- mean(.af_pcor_predictive_r2(xte, yhat), na.rm = TRUE)
    }

    scores[i] <- mean(fold_scores, na.rm = TRUE)
  }

  best_idx <- which.max(scores)
  list(
    best = lambdas[best_idx],
    grid = lambdas,
    scores = scores,
    metric = "predictive_r2"
  )
}

.af_pcor_cv_pc <- function(x_time_by_nodes, k_grid, folds = 5, seed = 1) {
  x <- as.matrix(x_time_by_nodes)
  fold_id <- .af_fold_assignments(nrow(x), folds, seed)
  scores <- numeric(length(k_grid))

  for (i in seq_along(k_grid)) {
    kval <- k_grid[i]
    fold_scores <- numeric(max(fold_id))

    for (f in seq_len(max(fold_id))) {
      idx_tr <- fold_id != f
      idx_te <- !idx_tr
      xtr <- x[idx_tr, , drop = FALSE]
      xte <- x[idx_te, , drop = FALSE]

      b <- .af_pcor_pc_coef(xtr, k = kval)
      yhat <- xte %*% b
      fold_scores[f] <- mean(.af_pcor_predictive_r2(xte, yhat), na.rm = TRUE)
    }

    scores[i] <- mean(fold_scores, na.rm = TRUE)
  }

  best_idx <- which.max(scores)
  list(
    best = k_grid[best_idx],
    grid = k_grid,
    scores = scores,
    metric = "predictive_r2"
  )
}

.af_pcor_pc_coef <- function(x_time_by_nodes, k = 5) {
  x <- as.matrix(x_time_by_nodes)
  x <- scale(x, center = TRUE, scale = FALSE)
  n_nodes <- ncol(x)
  n_time <- nrow(x)

  b <- matrix(0, n_nodes, n_nodes)
  for (target in seq_len(n_nodes)) {
    y <- x[, target]
    src <- x[, -target, drop = FALSE]

    k_use <- max(1L, min(as.integer(k), ncol(src), n_time - 1L))
    sv <- svd(src, nu = k_use, nv = k_use)
    scores <- sv$u[, seq_len(k_use), drop = FALSE] %*% diag(sv$d[seq_len(k_use)], nrow = k_use)

    xtx <- crossprod(scores)
    xty <- crossprod(scores, y)
    beta_pc <- try(.af_safe_inverse_sym(xtx) %*% xty, silent = TRUE)
    if (inherits(beta_pc, "try-error")) {
      beta_pc <- matrix(0, nrow = k_use, ncol = 1)
    }
    beta_src <- sv$v[, seq_len(k_use), drop = FALSE] %*% beta_pc

    b[target, -target] <- as.numeric(beta_src)
  }

  diag(b) <- 0
  b
}
