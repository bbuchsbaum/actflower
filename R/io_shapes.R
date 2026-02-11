as_nodes_by_time <- function(x, orientation = c("nodes_by_time", "time_by_nodes")) {
  orientation <- match.arg(orientation)
  if (!is.matrix(x) || !is.numeric(x)) {
    if (is_neuroim2_object(x)) {
      return(as_nodes_by_time_neuroim2(x, orientation = orientation))
    }
    actflower_abort("`x` must be a numeric matrix or a supported neuroim2 object.")
  }
  if (orientation == "nodes_by_time") {
    return(x)
  }
  t(x)
}
