.af_normalize_parcel_exclusions <- function(parcelstoexclude_bytarget, n_nodes) {
  out <- vector("list", n_nodes)
  for (i in seq_len(n_nodes)) out[[i]] <- integer(0)

  if (is.null(parcelstoexclude_bytarget)) {
    return(out)
  }

  if (is.list(parcelstoexclude_bytarget)) {
    if (length(parcelstoexclude_bytarget) != n_nodes) {
      actflower_abort("`parcelstoexclude_bytarget` list must have one entry per target node.")
    }
    for (i in seq_len(n_nodes)) {
      v <- as.integer(parcelstoexclude_bytarget[[i]])
      v <- v[is.finite(v)]
      v <- unique(v)
      v <- v[v >= 1 & v <= n_nodes]
      out[[i]] <- sort(v)
    }
    return(out)
  }

  if (is.matrix(parcelstoexclude_bytarget)) {
    if (!all(dim(parcelstoexclude_bytarget) == c(n_nodes, n_nodes))) {
      actflower_abort("`parcelstoexclude_bytarget` matrix must be [nodes x nodes].")
    }
    for (i in seq_len(n_nodes)) {
      out[[i]] <- which(as.logical(parcelstoexclude_bytarget[i, ]))
    }
    return(out)
  }

  actflower_abort("`parcelstoexclude_bytarget` must be NULL, list, or matrix.")
}

.af_read_parcel_exclusions_h5 <- function(path, n_nodes, group = "parcels_to_remove_indices") {
  if (!file.exists(path)) {
    actflower_abort(sprintf("Exclusion file does not exist: %s", path))
  }
  f <- hdf5r::H5File$new(path, mode = "r")
  on.exit(f$close_all(), add = TRUE)

  if (!f$exists(group)) {
    actflower_abort(sprintf("Group '%s' not found in exclusion file.", group))
  }

  out <- vector("list", n_nodes)
  g <- f[[group]]
  for (i in seq_len(n_nodes)) {
    nm <- as.character(i - 1L)
    if (!g$exists(nm)) {
      out[[i]] <- integer(0)
      next
    }
    v <- as.integer(g[[nm]]$read()) + 1L
    v <- v[is.finite(v)]
    v <- unique(v)
    v <- v[v >= 1 & v <= n_nodes]
    out[[i]] <- sort(v)
  }
  out
}

.af_pcr_single_target <- function(source_nodes_by_time, target_ts, n_components = NULL) {
  x <- t(source_nodes_by_time)
  y <- as.numeric(target_ts)

  max_comp <- min(ncol(x), nrow(x) - 1L)
  if (is.null(n_components)) {
    n_components <- max_comp
  }
  n_components <- max(1L, min(as.integer(n_components), max_comp))

  pca <- stats::prcomp(x, center = TRUE, scale. = FALSE)
  scores <- pca$x[, seq_len(n_components), drop = FALSE]
  fit <- stats::lm.fit(x = cbind(1, scores), y = y)
  beta_pc <- fit$coefficients[-1]
  beta_pc[is.na(beta_pc)] <- 0

  loadings <- pca$rotation[, seq_len(n_components), drop = FALSE]
  as.vector(loadings %*% beta_pc)
}

.af_multreg_masked <- function(data_nodes_by_time, mask) {
  data <- as.matrix(data_nodes_by_time)
  n_nodes <- nrow(data)
  fc <- matrix(0, n_nodes, n_nodes)
  all_nodes <- seq_len(n_nodes)

  for (target in all_nodes) {
    src <- which(mask[target, ] != 0)
    src <- setdiff(src, target)
    if (length(src) == 0L) next

    x <- t(data[src, , drop = FALSE])
    y <- as.numeric(data[target, ])
    fit <- stats::lm.fit(x = cbind(1, x), y = y)
    beta <- fit$coefficients[-1]
    beta[is.na(beta)] <- 0
    fc[target, src] <- beta
  }

  fc
}

#' Parcel-level non-circular connectivity estimation
#' @param data Numeric matrix.
#' @param connmethod Connectivity method.
#' @param parcelstoexclude_bytarget Optional per-target exclusions.
#' @param exclusion_h5 Optional HDF5 path containing per-target exclusions.
#' @param orientation Matrix orientation.
#' @param verbose Print progress.
#' @param ... Method-specific args.
#' @return FC matrix (target x source).
#' @export
calcconn_parcelwise_noncircular <- function(
  data,
  connmethod = c("multreg", "pearsoncorr", "pc_multregconn", "combinedFC"),
  parcelstoexclude_bytarget = NULL,
  exclusion_h5 = NULL,
  orientation = c("nodes_by_time", "time_by_nodes"),
  verbose = FALSE,
  ...
) {
  connmethod <- match.arg(connmethod)
  data <- as_nodes_by_time(data, orientation = orientation)
  n_nodes <- nrow(data)

  exclusions <- if (!is.null(exclusion_h5)) {
    .af_read_parcel_exclusions_h5(exclusion_h5, n_nodes = n_nodes)
  } else {
    .af_normalize_parcel_exclusions(parcelstoexclude_bytarget, n_nodes = n_nodes)
  }

  fc <- matrix(0, n_nodes, n_nodes)
  all_nodes <- seq_len(n_nodes)
  dots <- list(...)

  if (connmethod == "combinedFC") {
    net_mask <- matrix(0, n_nodes, n_nodes)

    for (target in all_nodes) {
      if (isTRUE(verbose)) message("combinedFC mask target ", target, "/", n_nodes)
      excluded <- sort(unique(c(target, exclusions[[target]])))
      src <- setdiff(all_nodes, excluded)
      if (length(src) == 0L) next

      sub_data <- rbind(data[target, , drop = FALSE], data[src, , drop = FALSE])
      sub_fc <- estimate_fc_combined(
        sub_data,
        orientation = "nodes_by_time",
        method_cond_asso = dots$method_cond_asso %||% "partialCorrelation",
        method_parcorr = dots$method_parcorr %||% "inverseCovariance",
        alpha_cond_asso = dots$alpha_cond_asso %||% 0.01,
        method_asso = dots$method_asso %||% "correlation",
        alpha_asso = dots$alpha_asso %||% 0.01,
        equivalence_test_asso = dots$equivalence_test_asso %||% FALSE,
        lower_bound = dots$lower_bound %||% -0.1,
        upper_bound = dots$upper_bound %||% 0.1
      )
      net_mask[target, src] <- sub_fc[1, -1]
    }

    fc <- .af_multreg_masked(data, mask = net_mask != 0)
    return(fc)
  }

  for (target in all_nodes) {
    if (isTRUE(verbose)) message("target ", target, "/", n_nodes)

    excluded <- sort(unique(c(target, exclusions[[target]])))
    src <- setdiff(all_nodes, excluded)
    if (length(src) == 0L) next

    y <- as.numeric(data[target, ])
    xsrc <- data[src, , drop = FALSE]

    if (connmethod == "pearsoncorr") {
      vals <- apply(xsrc, 1, function(v) pearson_vec(v, y))
      vals[is.na(vals)] <- 0
      fc[target, src] <- vals
      next
    }

    if (connmethod == "multreg") {
      fit <- stats::lm.fit(x = cbind(1, t(xsrc)), y = y)
      beta <- fit$coefficients[-1]
      beta[is.na(beta)] <- 0
      fc[target, src] <- beta
      next
    }

    # pc_multregconn
    n_comp <- dots$n_components %||% NULL
    fc[target, src] <- .af_pcr_single_target(xsrc, y, n_components = n_comp)
  }

  fc
}

#' Parcel-level non-circular activity tensor
#' @param data Numeric matrix of parcel activity.
#' @param parcelstoexclude_bytarget Optional per-target exclusions.
#' @param exclusion_h5 Optional HDF5 path containing per-target exclusions.
#' @param orientation Matrix orientation.
#' @param fill_value Value used for excluded source entries.
#' @return Array with dimensions (target x source x condition).
#' @export
calcactivity_parcelwise_noncircular <- function(
  data,
  parcelstoexclude_bytarget = NULL,
  exclusion_h5 = NULL,
  orientation = c("nodes_by_conditions", "conditions_by_nodes"),
  fill_value = 0
) {
  orientation <- match.arg(orientation)
  mat <- as.matrix(data)
  if (orientation == "conditions_by_nodes") {
    mat <- t(mat)
  }
  validate_matrix(mat, arg = "data")

  n_nodes <- nrow(mat)
  n_cond <- ncol(mat)

  exclusions <- if (!is.null(exclusion_h5)) {
    .af_read_parcel_exclusions_h5(exclusion_h5, n_nodes = n_nodes)
  } else {
    .af_normalize_parcel_exclusions(parcelstoexclude_bytarget, n_nodes = n_nodes)
  }

  out <- array(fill_value, dim = c(n_nodes, n_nodes, n_cond))
  all_nodes <- seq_len(n_nodes)

  for (target in all_nodes) {
    excluded <- sort(unique(exclusions[[target]]))
    src <- setdiff(all_nodes, excluded)
    if (length(src) > 0) {
      out[target, src, ] <- mat[src, , drop = FALSE]
    }
    out[target, target, ] <- mat[target, ]
  }

  out
}
