#' Combined functional connectivity (combinedFC)
#'
#' Implements the two-step combinedFC procedure from ActflowToolbox defaults:
#' (1) conditional association (partial correlation or multiple regression), then
#' (2) collider check via unconditional association.
#'
#' @param x Numeric matrix.
#' @param orientation Orientation of matrix.
#' @param method_cond_asso Conditional association method.
#' @param method_parcorr Partial-correlation backend when
#'   `method_cond_asso = "partialCorrelation"`.
#' @param alpha_cond_asso Alpha for conditional association significance.
#' @param method_asso Unconditional association method.
#' @param alpha_asso Alpha for unconditional association significance.
#' @param equivalence_test_asso Use equivalence test for unconditional association.
#' @param lower_bound Lower equivalence bound.
#' @param upper_bound Upper equivalence bound.
#' @return FC matrix (nodes x nodes).
#' @export
estimate_fc_combined <- function(
  x,
  orientation = c("nodes_by_time", "time_by_nodes"),
  method_cond_asso = c("partialCorrelation", "multipleRegression"),
  method_parcorr = c("inverseCovariance", "regression"),
  alpha_cond_asso = 0.01,
  method_asso = c("correlation", "simpleRegression"),
  alpha_asso = 0.01,
  equivalence_test_asso = FALSE,
  lower_bound = -0.1,
  upper_bound = 0.1
) {
  orientation <- match.arg(orientation)
  method_cond_asso <- match.arg(method_cond_asso)
  method_parcorr <- match.arg(method_parcorr)
  method_asso <- match.arg(method_asso)

  x_nodes_by_time <- as_nodes_by_time(x, orientation = orientation)
  # combinedFC convention: datapoints x nodes.
  d <- t(x_nodes_by_time)

  if (ncol(d) < 2) {
    actflower_abort("estimate_fc_combined() requires at least 2 nodes.")
  }

  if (method_cond_asso == "partialCorrelation") {
    mca <- .af_combined_partial_correlation_sig(
      d,
      alpha = alpha_cond_asso,
      method = method_parcorr
    )
  } else {
    mca <- .af_combined_multiple_regression_sig(
      d,
      alpha = alpha_cond_asso,
      sig_test = TRUE
    )
  }

  n_nodes <- ncol(d)
  m <- mca

  if (method_asso == "correlation") {
    mcorr <- .af_combined_correlation_sig(
      d,
      alpha = alpha_asso,
      lower_bound = lower_bound,
      upper_bound = upper_bound,
      equivalence_test = equivalence_test_asso
    )

    for (xidx in seq_len(n_nodes - 1L)) {
      for (yidx in seq.int(xidx + 1L, n_nodes)) {
        if (mca[xidx, yidx] != 0 && mcorr[xidx, yidx] == 0) {
          m[xidx, yidx] <- 0
        }
        if (mca[yidx, xidx] != 0 && mcorr[yidx, xidx] == 0) {
          m[yidx, xidx] <- 0
        }
      }
    }
  } else {
    for (xidx in seq_len(n_nodes - 1L)) {
      for (yidx in seq.int(xidx + 1L, n_nodes)) {
        if (mca[xidx, yidx] != 0) {
          bxy <- .af_combined_simple_regression_sig(
            y = d[, xidx],
            x = d[, yidx],
            alpha = alpha_asso,
            sig_test = TRUE
          )$b
          if (bxy == 0) {
            m[xidx, yidx] <- 0
          }
        }
        if (mca[yidx, xidx] != 0) {
          byx <- .af_combined_simple_regression_sig(
            y = d[, yidx],
            x = d[, xidx],
            alpha = alpha_asso,
            sig_test = TRUE
          )$b
          if (byx == 0) {
            m[yidx, xidx] <- 0
          }
        }
      }
    }
  }

  # Remove negative zero artifacts and enforce zero diagonal.
  m <- m + 0
  diag(m) <- 0
  m
}

.af_combined_partial_correlation_sig <- function(dataset, alpha = 0.01, method = "inverseCovariance") {
  d <- as.matrix(dataset)
  n_nodes <- ncol(d)
  n_datapoints <- nrow(d)

  m_parcorr <- if (method == "regression") {
    .af_combined_parcorr_regression(d)
  } else {
    .af_combined_parcorr_invcov(d)
  }
  cond_set_size <- n_nodes - 2

  z_alpha <- .af_combined_zcutoff(alpha = alpha, kind = "two-sided")
  fz <- .af_combined_fisher_z(
    r = m_parcorr,
    n_datapoints = n_datapoints,
    ho = 0,
    cond_set_size = cond_set_size
  )
  m <- m_parcorr * (abs(fz) >= z_alpha)
  m + 0
}

.af_combined_parcorr_regression <- function(dataset) {
  d <- as.matrix(dataset)
  n_nodes <- ncol(d)
  m <- matrix(0, n_nodes, n_nodes)

  for (xidx in seq_len(n_nodes - 1L)) {
    for (yidx in seq.int(xidx + 1L, n_nodes)) {
      idx <- rep(TRUE, n_nodes)
      idx[xidx] <- FALSE
      idx[yidx] <- FALSE

      x_rest <- d[, idx, drop = FALSE]
      fit_x <- stats::lm.fit(x = cbind(1, x_rest), y = d[, xidx])
      fit_y <- stats::lm.fit(x = cbind(1, x_rest), y = d[, yidx])

      # Match Python logic: subtract only the slope component.
      beta_x <- fit_x$coefficients[-1]
      beta_y <- fit_y$coefficients[-1]
      beta_x[is.na(beta_x)] <- 0
      beta_y[is.na(beta_y)] <- 0

      res_x <- d[, xidx] - as.vector(x_rest %*% beta_x)
      res_y <- d[, yidx] - as.vector(x_rest %*% beta_y)
      parcorr <- suppressWarnings(stats::cor(res_x, res_y, use = "pairwise.complete.obs"))
      if (!is.finite(parcorr)) {
        parcorr <- 0
      }
      m[xidx, yidx] <- parcorr
      m[yidx, xidx] <- parcorr
    }
  }

  m
}

.af_combined_parcorr_invcov <- function(dataset) {
  d <- as.matrix(dataset)
  inv_cov <- .af_combined_pinv_sym(stats::cov(d))
  denom <- sqrt(pmax(diag(inv_cov), .Machine$double.eps))
  m <- -inv_cov / (denom %o% denom)
  diag(m) <- 0
  m
}

.af_combined_multiple_regression_sig <- function(dataset, alpha = 0.01, sig_test = FALSE) {
  d <- as.matrix(dataset)
  n_nodes <- ncol(d)
  n_datapoints <- nrow(d)
  m <- matrix(0, n_nodes, n_nodes)

  for (xidx in seq_len(n_nodes)) {
    idx <- rep(TRUE, n_nodes)
    idx[xidx] <- FALSE
    xr <- d[, idx, drop = FALSE]

    fit <- stats::lm.fit(x = cbind(1, xr), y = d[, xidx])
    params <- fit$coefficients
    params[is.na(params)] <- 0

    if (!sig_test) {
      m[xidx, idx] <- params[-1]
      next
    }

    n_params <- length(params)
    y_hat <- as.vector(cbind(1, xr) %*% params)
    mse <- sum((d[, xidx] - y_hat) ^ 2) / max(1, n_datapoints - n_params)
    xtx_inv <- .af_combined_pinv_sym(crossprod(cbind(1, xr)))
    var_params <- mse * diag(xtx_inv)
    std_params <- sqrt(pmax(var_params, .Machine$double.eps))
    t_stats <- params / std_params
    p_vals <- 2 * (1 - stats::pt(abs(t_stats), df = n_datapoints - 1))

    keep <- as.numeric(p_vals[-1] < alpha)
    m[xidx, idx] <- params[-1] * keep
  }

  m
}

.af_combined_correlation_sig <- function(
  dataset,
  alpha = 0.01,
  lower_bound = -0.1,
  upper_bound = 0.1,
  equivalence_test = FALSE
) {
  d <- as.matrix(dataset)
  n_datapoints <- nrow(d)

  mcorr <- stats::cor(d, use = "pairwise.complete.obs")
  diag(mcorr) <- 0

  if (!equivalence_test) {
    z_alpha <- .af_combined_zcutoff(alpha = alpha, kind = "two-sided")
    fz <- .af_combined_fisher_z(r = mcorr, n_datapoints = n_datapoints, ho = 0)
    m <- mcorr * (abs(fz) >= z_alpha)
    return(m + 0)
  }

  z_alpha_u <- .af_combined_zcutoff(alpha = alpha, kind = "one-sided-left")
  z_alpha_l <- .af_combined_zcutoff(alpha = alpha, kind = "one-sided-right")
  fz_u <- .af_combined_fisher_z(r = mcorr, n_datapoints = n_datapoints, ho = upper_bound)
  fz_l <- .af_combined_fisher_z(r = mcorr, n_datapoints = n_datapoints, ho = lower_bound)

  m <- mcorr * ((fz_u <= z_alpha_u) & (fz_l >= z_alpha_l))
  m + 0
}

.af_combined_simple_regression_sig <- function(y, x, alpha = 0.01, sig_test = TRUE) {
  y <- as.numeric(y)
  x <- as.numeric(x)
  n_datapoints <- length(y)

  fit <- stats::lm.fit(x = cbind(1, x), y = y)
  params <- fit$coefficients
  params[is.na(params)] <- 0

  if (!sig_test) {
    return(list(b = params[2]))
  }

  n_params <- length(params)
  y_hat <- as.vector(cbind(1, x) %*% params)
  mse <- sum((y - y_hat) ^ 2) / max(1, n_datapoints - n_params)
  xtx_inv <- .af_combined_pinv_sym(crossprod(cbind(1, x)))
  var_params <- mse * diag(xtx_inv)
  std_params <- sqrt(pmax(var_params, .Machine$double.eps))
  t_stats <- params / std_params
  p_vals <- 2 * (1 - stats::pt(abs(t_stats), df = n_datapoints - 1))
  p_beta <- p_vals[2]
  b <- params[2] * as.numeric(p_beta < alpha)

  list(b = as.numeric(b), p_value = as.numeric(p_beta))
}

.af_combined_zcutoff <- function(alpha = 0.01, kind = "two-sided") {
  if (kind == "two-sided") {
    return(stats::qnorm(1 - alpha / 2))
  }
  if (kind == "one-sided-right") {
    return(stats::qnorm(1 - alpha))
  }
  if (kind == "one-sided-left") {
    return(stats::qnorm(alpha))
  }
  actflower_abort(sprintf("Unknown zcutoff kind: %s", kind))
}

.af_combined_fisher_z <- function(r, n_datapoints, ho = 0, cond_set_size = 0) {
  clamp <- function(x) pmax(pmin(x, 1 - 1e-12), -1 + 1e-12)
  r <- clamp(r)
  ho <- clamp(ho)
  (atanh(r) - atanh(ho)) * sqrt(n_datapoints - cond_set_size - 3)
}

.af_combined_pinv_sym <- function(m) {
  m <- as.matrix(m)
  m <- (m + t(m)) * 0.5
  eig <- eigen(m, symmetric = TRUE)
  vals <- eig$values
  tol <- max(vals) * .Machine$double.eps * length(vals)
  vals_inv <- ifelse(vals > tol, 1 / vals, 0)
  eig$vectors %*% (diag(vals_inv, nrow = length(vals_inv)) %*% t(eig$vectors))
}
