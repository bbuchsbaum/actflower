#' Correlation connectivity
#' @param x Numeric matrix.
#' @param orientation Orientation of matrix.
#' @param use_cpp Use native kernel when available.
#' @return FC matrix (nodes x nodes).
#' @export
estimate_fc_corr <- function(
  x,
  orientation = c("nodes_by_time", "time_by_nodes"),
  use_cpp = TRUE
) {
  x <- as_nodes_by_time(x, orientation = orientation)
  .af_corr_fc(x, use_cpp = use_cpp)
}
