#' Multiple-regression connectivity
#' @param x Numeric matrix.
#' @param orientation Orientation of matrix.
#' @param ridge Optional ridge loading added to covariance diagonal.
#' @param use_cpp Use native kernel when available.
#' @return FC matrix (nodes x nodes) with rows as targets and columns as sources.
#' @export
estimate_fc_multreg <- function(
  x,
  orientation = c("nodes_by_time", "time_by_nodes"),
  ridge = 0,
  use_cpp = TRUE
) {
  x <- as_nodes_by_time(x, orientation = orientation)
  n_nodes <- nrow(x)
  n_time <- ncol(x)
  if (n_nodes >= n_time) {
    actflower_abort(sprintf("multreg requires nodes < timepoints; got nodes=%d timepoints=%d", n_nodes, n_time))
  }
  if (!is.numeric(ridge) || length(ridge) != 1 || ridge < 0) {
    actflower_abort("`ridge` must be a non-negative scalar.")
  }
  .af_multreg_fc(x, ridge = ridge, use_cpp = use_cpp)
}

.af_safe_inverse_sym <- function(m) {
  m <- (m + t(m)) * 0.5
  out <- try(solve(m), silent = TRUE)
  if (!inherits(out, "try-error")) {
    return(out)
  }

  # Fallback for ill-conditioned covariance.
  eig <- eigen(m, symmetric = TRUE)
  vals <- pmax(eig$values, .Machine$double.eps)
  eig$vectors %*% (diag(1 / vals, nrow = length(vals)) %*% t(eig$vectors))
}
