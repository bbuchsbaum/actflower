.af_native_symbol_cache <- new.env(parent = emptyenv())

.af_has_native_symbol <- function(symbol) {
  if (exists(symbol, envir = .af_native_symbol_cache, inherits = FALSE)) {
    return(get(symbol, envir = .af_native_symbol_cache, inherits = FALSE))
  }
  out <- try(getNativeSymbolInfo(symbol, PACKAGE = "actflower"), silent = TRUE)
  ok <- !inherits(out, "try-error")
  assign(symbol, ok, envir = .af_native_symbol_cache)
  ok
}

.af_can_call_wrapper <- function(fname) {
  exists(fname, mode = "function", inherits = TRUE)
}

.af_native_enabled <- function(fname, symbol) {
  .af_can_call_wrapper(fname) && .af_has_native_symbol(symbol)
}

.af_multreg_fc <- function(x_nodes_by_time, ridge = 0, use_cpp = TRUE) {
  if (use_cpp && .af_native_enabled("multreg_fc_cpp", "_actflower_multreg_fc_cpp")) {
    out <- try(multreg_fc_cpp(x_nodes_by_time, ridge = ridge), silent = TRUE)
    if (!inherits(out, "try-error")) {
      return(out)
    }
  }
  .af_multreg_fc_r(x_nodes_by_time, ridge = ridge)
}

.af_corr_fc <- function(x_nodes_by_time, use_cpp = TRUE) {
  # Keep exact pairwise-complete behavior for missing data via R fallback.
  if (anyNA(x_nodes_by_time) || !all(is.finite(x_nodes_by_time))) {
    return(.af_corr_fc_r(x_nodes_by_time))
  }

  if (use_cpp && .af_native_enabled("corr_fc_cpp", "_actflower_corr_fc_cpp")) {
    out <- try(corr_fc_cpp(x_nodes_by_time), silent = TRUE)
    if (!inherits(out, "try-error")) {
      return(out)
    }
  }
  .af_corr_fc_r(x_nodes_by_time)
}

.af_corr_fc_r <- function(x_nodes_by_time) {
  fc <- stats::cor(t(x_nodes_by_time), use = "pairwise.complete.obs")
  diag(fc) <- 0
  fc
}

.af_multreg_fc_r <- function(x_nodes_by_time, ridge = 0) {
  n_nodes <- nrow(x_nodes_by_time)

  x_centered <- x_nodes_by_time - rowMeans(x_nodes_by_time)
  sigma <- stats::cov(t(x_centered))
  if (ridge > 0) {
    sigma <- sigma + diag(ridge, n_nodes)
  }

  theta <- .af_safe_inverse_sym(sigma)
  theta <- (theta + t(theta)) * 0.5

  denom <- diag(theta)
  denom[abs(denom) < .Machine$double.eps] <- .Machine$double.eps
  beta <- -diag(1 / denom, n_nodes) %*% theta
  diag(beta) <- 0
  beta
}

.af_actflow_predict_batch <- function(act_group, fc_group, use_cpp = TRUE) {
  if (use_cpp && .af_native_enabled("actflow_predict_batch_cpp", "_actflower_actflow_predict_batch_cpp")) {
    out <- try(actflow_predict_batch_cpp(act_group, fc_group, remove_diag = TRUE), silent = TRUE)
    if (!inherits(out, "try-error")) {
      return(out)
    }
  }

  # R fallback: per-subject matrix multiply with explicit self-connection exclusion.
  d <- dim(act_group)
  n_nodes <- d[1]
  n_cond <- d[2]
  n_subj <- d[3]

  pred <- array(NA_real_, dim = c(n_nodes, n_cond, n_subj))
  for (s in seq_len(n_subj)) {
    a <- act_group[, , s, drop = TRUE]
    f <- fc_group[, , s, drop = TRUE]
    p <- f %*% a
    p <- p - diag(f) * a
    pred[, , s] <- p
  }
  pred
}

.af_compare_fullcomp <- function(target, pred, use_cpp = TRUE) {
  if (use_cpp && .af_native_enabled("compare_fullcomp_cpp", "_actflower_compare_fullcomp_cpp")) {
    out <- try(compare_fullcomp_cpp(target, pred), silent = TRUE)
    if (!inherits(out, "try-error")) {
      return(out)
    }
  }

  n_subj <- dim(target)[3]
  corr_vals <- numeric(n_subj)
  r2_vals <- numeric(n_subj)
  mae_vals <- numeric(n_subj)
  for (s in seq_len(n_subj)) {
    yt <- as.vector(target[, , s])
    yp <- as.vector(pred[, , s])
    corr_vals[s] <- pearson_vec(yt, yp)
    r2_vals[s] <- r2_vec(yt, yp)
    mae_vals[s] <- mae_vec(yt, yp)
  }
  list(corr_vals = corr_vals, R2_vals = r2_vals, mae_vals = mae_vals)
}
