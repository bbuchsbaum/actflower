as_nodes_by_time <- function(x, orientation = c("nodes_by_time", "time_by_nodes")) {
  orientation <- match.arg(orientation)
  if (!is.matrix(x) || !is.numeric(x)) {
    if (is_neuro_object(x)) {
      return(as_nodes_by_time_object(x, orientation = orientation))
    }
    actflower_abort("`x` must be a numeric matrix or a supported neuro object.")
  }
  if (orientation == "nodes_by_time") {
    return(x)
  }
  t(x)
}
