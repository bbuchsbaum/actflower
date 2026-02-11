.af_normalize_orientation <- function(orientation) {
  if (!is.character(orientation) || length(orientation) != 1L || !nzchar(orientation)) {
    actflower_abort("`orientation` must be a non-empty character scalar.")
  }

  key <- tolower(gsub("[^a-z]", "", orientation))
  if (key %in% c("nodesbytime", "nodesxtime", "nodestime")) return("nodes_by_time")
  if (key %in% c("timebynodes", "timexnodes", "timenodes")) return("time_by_nodes")

  actflower_abort(
    paste0(
      "Unknown `orientation`: ", orientation,
      ". Expected one of: nodes_by_time, time_by_nodes, nodes_x_time, time_x_nodes."
    )
  )
}

as_nodes_by_time <- function(x, orientation = c("nodes_by_time", "time_by_nodes")) {
  orientation <- .af_normalize_orientation(orientation[[1]])
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
