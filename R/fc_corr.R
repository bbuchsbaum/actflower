#' Correlation connectivity
#' @param x Numeric matrix.
#' @param orientation Orientation of matrix.
#' @return FC matrix (nodes x nodes).
#' @export
estimate_fc_corr <- function(x, orientation = c("nodes_by_time", "time_by_nodes")) {
  x <- as_nodes_by_time(x, orientation = orientation)
  fc <- stats::cor(t(x), use = "pairwise.complete.obs")
  diag(fc) <- 0
  fc
}
