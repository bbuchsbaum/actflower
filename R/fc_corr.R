#' Correlation connectivity
#' @param x Numeric matrix or 3D array.
#' @param orientation Orientation of matrix.
#' @param use_cpp Use native kernel when available.
#' @return FC matrix (nodes x nodes) or 3D FC array (nodes x nodes x subjects).
#' @export
estimate_fc_corr <- function(
  x,
  orientation = c("nodes_by_time", "time_by_nodes"),
  use_cpp = TRUE
) {
  orientation <- .af_normalize_orientation(orientation[[1]])

  if (is.array(x) && is.numeric(x) && length(dim(x)) == 3L) {
    x_group <- if (orientation == "nodes_by_time") x else aperm(x, c(2, 1, 3))
    return(.af_corr_fc_group(x_group, use_cpp = use_cpp))
  }

  x <- as_nodes_by_time(x, orientation = orientation)
  .af_corr_fc(x, use_cpp = use_cpp)
}
