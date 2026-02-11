# Internal helper: hdf5r::H5File$exists() can throw when intermediate links
# in a nested path are missing.
.af_h5_exists <- function(h5_file, path) {
  out <- try(h5_file$exists(path), silent = TRUE)
  if (inherits(out, "try-error")) return(FALSE)
  isTRUE(out)
}

.af_h5_ensure_parent_groups <- function(h5_file, path) {
  parts <- strsplit(path, "/", fixed = TRUE)[[1]]
  parts <- parts[nzchar(parts)]
  if (length(parts) <= 1L) return(invisible(NULL))

  cur <- character(0)
  for (p in parts[-length(parts)]) {
    cur <- c(cur, p)
    gpath <- paste(cur, collapse = "/")
    if (!.af_h5_exists(h5_file, gpath)) {
      h5_file$create_group(gpath)
    }
  }
  invisible(NULL)
}

#' Read an HDF5 dataset
#' @param path File path.
#' @param dataset Dataset name/path in the HDF5 file.
#' @return Numeric array/matrix.
#' @export
read_actflower_h5 <- function(path, dataset) {
  if (!file.exists(path)) {
    actflower_abort(sprintf("File does not exist: %s", path))
  }
  f <- hdf5r::H5File$new(path, mode = "r")
  on.exit(f$close_all(), add = TRUE)
  if (!.af_h5_exists(f, dataset)) {
    actflower_abort(sprintf("Dataset not found in HDF5 file: %s", dataset))
  }
  f[[dataset]]$read()
}

#' Write an HDF5 dataset
#' @param path File path.
#' @param name Dataset name/path.
#' @param x Object to write.
#' @export
write_actflower_h5 <- function(path, name, x) {
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    actflower_abort("`name` must be a non-empty character scalar.")
  }
  mode <- if (file.exists(path)) "a" else "w"
  f <- hdf5r::H5File$new(path, mode = mode)
  on.exit(f$close_all(), add = TRUE)
  .af_h5_ensure_parent_groups(f, name)
  if (.af_h5_exists(f, name)) {
    f$link_delete(name)
  }
  f[[name]] <- x
  invisible(path)
}
