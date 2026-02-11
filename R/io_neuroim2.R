# Neuro object integration helpers.
#
# neuroim2 is a required dependency for actflower and is treated as a first-
# class data source for matrix conversion.

is_neuro_object <- function(x) {
  inherits(x, c("NeuroVec", "ClusteredNeuroVec", "ROIVec", "NeuroVol"))
}

as_nodes_by_time_object <- function(
  x,
  orientation = c("nodes_by_time", "time_by_nodes")
) {
  if (!is_neuro_object(x)) {
    actflower_abort("`x` must inherit from a supported neuro object class (e.g., NeuroVec, ClusteredNeuroVec).")
  }

  orientation <- match.arg(orientation)

  mat <- try(as.matrix(x), silent = TRUE)
  if (inherits(mat, "try-error") || !is.matrix(mat) || !is.numeric(mat) || is_neuro_object(mat)) {
    # Fallbacks for list-like objects
    if (!is.null(x$X)) {
      mat <- as.matrix(x$X)
    } else if (!is.null(x$data)) {
      mat <- as.matrix(x$data)
    } else {
      actflower_abort("Could not coerce neuro object to matrix.")
    }
  }

  if (!is.matrix(mat) || !is.numeric(mat)) {
    actflower_abort("Neuro object conversion produced a non-numeric matrix.")
  }

  if (orientation == "nodes_by_time") mat else t(mat)
}
