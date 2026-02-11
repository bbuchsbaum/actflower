# neuroim2 integration helpers
#
# This module is intentionally lightweight and defensive so actflower can run
# without neuroim2 while still providing a clear adapter path when neuroim2
# objects are supplied.

has_neuroim2 <- function() {
  requireNamespace("neuroim2", quietly = TRUE)
}

is_neuroim2_object <- function(x) {
  inherits(x, c("NeuroVec", "ClusteredNeuroVec", "ROIVec", "NeuroVol"))
}

as_nodes_by_time_neuroim2 <- function(
  x,
  orientation = c("nodes_by_time", "time_by_nodes")
) {
  if (!has_neuroim2()) {
    actflower_abort("Install 'neuroim2' to convert neuroim2 objects.")
  }
  if (!is_neuroim2_object(x)) {
    actflower_abort("`x` must inherit from a neuroim2 class (e.g., NeuroVec, ClusteredNeuroVec).")
  }

  orientation <- match.arg(orientation)

  mat <- try(as.matrix(x), silent = TRUE)
  if (inherits(mat, "try-error") || !is.matrix(mat)) {
    # Fallbacks for list-like objects
    if (!is.null(x$X)) {
      mat <- as.matrix(x$X)
    } else if (!is.null(x$data)) {
      mat <- as.matrix(x$data)
    } else {
      actflower_abort("Could not coerce neuroim2 object to matrix.")
    }
  }

  as_nodes_by_time(mat, orientation = orientation)
}
