validate_matrix <- function(x, arg = "x") {
  if (!is.matrix(x) || !is.numeric(x)) {
    actflower_abort(sprintf("`%s` must be a numeric matrix.", arg))
  }
  invisible(TRUE)
}

validate_array3 <- function(x, arg = "x") {
  if (!is.array(x) || length(dim(x)) != 3 || !is.numeric(x)) {
    actflower_abort(sprintf("`%s` must be a numeric 3D array [nodes x conditions x subjects].", arg))
  }
  invisible(TRUE)
}

validate_square_matrix <- function(x, arg = "x") {
  validate_matrix(x, arg = arg)
  dx <- dim(x)
  if (dx[1] != dx[2]) {
    actflower_abort(sprintf("`%s` must be square; got %d x %d.", arg, dx[1], dx[2]))
  }
  invisible(TRUE)
}
