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

  if (inherits(parcelstoexclude_bytarget, "dgCMatrix")) {
    if (!all(dim(parcelstoexclude_bytarget) == c(n_nodes, n_nodes))) {
      actflower_abort("`parcelstoexclude_bytarget` sparse matrix must be [nodes x nodes].")
    }

    sm <- Matrix::summary(parcelstoexclude_bytarget)
    if (nrow(sm) > 0) {
      idx <- split(as.integer(sm$j), as.integer(sm$i))
      for (nm in names(idx)) {
        i <- as.integer(nm)
        v <- idx[[nm]]
        v <- v[is.finite(v)]
        v <- unique(v)
        v <- v[v >= 1 & v <= n_nodes]
        out[[i]] <- sort(v)
      }
    }
    return(out)
  }

  actflower_abort("`parcelstoexclude_bytarget` must be NULL, list, matrix, or dgCMatrix.")
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

.af_sources_from_exclusions <- function(exclusions, n_nodes) {
  all_nodes <- seq_len(n_nodes)
  out <- vector("list", n_nodes)
  for (target in all_nodes) {
    excluded <- sort(unique(c(target, exclusions[[target]])))
    out[[target]] <- setdiff(all_nodes, excluded)
  }
  out
}

.af_dense_source_mask_from_exclusions <- function(exclusions, n_nodes) {
  mask <- matrix(TRUE, n_nodes, n_nodes)
  for (target in seq_len(n_nodes)) {
    excluded <- sort(unique(c(target, exclusions[[target]])))
    if (length(excluded) > 0L) {
      mask[target, excluded] <- FALSE
    }
  }
  mask
}

.af_sources_from_dense_mask <- function(mask) {
  n_nodes <- nrow(mask)
  out <- vector("list", n_nodes)
  for (target in seq_len(n_nodes)) {
    out[[target]] <- which(mask[target, ] != 0)
  }
  out
}

.af_sources_to_dgc <- function(sources_by_target, n_nodes) {
  i <- integer(0)
  j <- integer(0)

  for (target in seq_len(n_nodes)) {
    src <- sources_by_target[[target]]
    if (!length(src)) next
    i <- c(i, rep.int(target, length(src)))
    j <- c(j, src)
  }

  Matrix::sparseMatrix(
    i = i,
    j = j,
    x = rep.int(1, length(i)),
    dims = c(n_nodes, n_nodes),
    repr = "C"
  )
}

.af_build_source_index <- function(exclusions, n_nodes, sparse = TRUE, mask_format = c("list", "dgCMatrix")) {
  mask_format <- match.arg(mask_format)

  if (!isTRUE(sparse)) {
    mask <- .af_dense_source_mask_from_exclusions(exclusions, n_nodes = n_nodes)
    sources <- .af_sources_from_dense_mask(mask)
    nnz <- sum(mask)
    return(list(
      sources = sources,
      mask = mask,
      sparsity = list(
        nnz = nnz,
        total = n_nodes * n_nodes,
        density = if ((n_nodes * n_nodes) > 0) nnz / (n_nodes * n_nodes) else 0
      )
    ))
  }

  sources <- .af_sources_from_exclusions(exclusions, n_nodes = n_nodes)
  if (identical(mask_format, "dgCMatrix")) {
    mask <- .af_sources_to_dgc(sources, n_nodes = n_nodes)
    nnz <- length(mask@x)
  } else {
    mask <- NULL
    nnz <- sum(vapply(sources, length, integer(1)))
  }

  list(
    sources = sources,
    mask = mask,
    sparsity = list(
      nnz = nnz,
      total = n_nodes * n_nodes,
      density = if ((n_nodes * n_nodes) > 0) nnz / (n_nodes * n_nodes) else 0
    )
  )
}

.af_multreg_sources <- function(data_nodes_by_time, sources_by_target) {
  data <- as.matrix(data_nodes_by_time)
  n_nodes <- nrow(data)
  fc <- matrix(0, n_nodes, n_nodes)

  for (target in seq_len(n_nodes)) {
    src <- sources_by_target[[target]]
    if (!length(src)) next

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
#' @param sparse If `TRUE`, build and use sparse/mask-aware source indexing.
#' @param mask_format Sparse mask representation when `sparse=TRUE`.
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
  sparse = TRUE,
  mask_format = c("list", "dgCMatrix"),
  verbose = FALSE,
  ...
) {
  connmethod <- match.arg(connmethod)
  data <- as_nodes_by_time(data, orientation = orientation)
  n_nodes <- nrow(data)
  mask_format <- match.arg(mask_format)

  exclusions <- if (!is.null(exclusion_h5)) {
    .af_read_parcel_exclusions_h5(exclusion_h5, n_nodes = n_nodes)
  } else {
    .af_normalize_parcel_exclusions(parcelstoexclude_bytarget, n_nodes = n_nodes)
  }

  source_index <- .af_build_source_index(
    exclusions,
    n_nodes = n_nodes,
    sparse = sparse,
    mask_format = mask_format
  )

  fc <- matrix(0, n_nodes, n_nodes)
  all_nodes <- seq_len(n_nodes)
  dots <- list(...)

  if (connmethod == "combinedFC") {
    active_sources <- vector("list", n_nodes)
    for (i in seq_len(n_nodes)) active_sources[[i]] <- integer(0)

    for (target in all_nodes) {
      if (isTRUE(verbose)) message("combinedFC mask target ", target, "/", n_nodes)
      src <- if (!isTRUE(sparse) && !is.null(source_index$mask)) {
        which(source_index$mask[target, ] != 0)
      } else {
        source_index$sources[[target]]
      }
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
      nz <- which(sub_fc[1, -1] != 0)
      if (length(nz) > 0L) {
        active_sources[[target]] <- src[nz]
      }
    }

    fc <- .af_multreg_sources(data, sources_by_target = active_sources)
    attr(fc, "sparsity_stats") <- source_index$sparsity
    attr(fc, "source_mask_format") <- if (isTRUE(sparse)) mask_format else "dense"
    return(fc)
  }

  for (target in all_nodes) {
    if (isTRUE(verbose)) message("target ", target, "/", n_nodes)
    src <- if (!isTRUE(sparse) && !is.null(source_index$mask)) {
      which(source_index$mask[target, ] != 0)
    } else {
      source_index$sources[[target]]
    }
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

  attr(fc, "sparsity_stats") <- source_index$sparsity
  attr(fc, "source_mask_format") <- if (isTRUE(sparse)) mask_format else "dense"
  fc
}

#' Parcel-level non-circular activity tensor
#' @param data Numeric matrix of parcel activity.
#' @param parcelstoexclude_bytarget Optional per-target exclusions.
#' @param exclusion_h5 Optional HDF5 path containing per-target exclusions.
#' @param orientation Matrix orientation.
#' @param fill_value Value used for excluded source entries.
#' @param sparse If `TRUE`, build and use sparse/mask-aware source indexing.
#' @param mask_format Sparse mask representation when `sparse=TRUE`.
#' @return Array with dimensions (target x source x condition).
#' @export
calcactivity_parcelwise_noncircular <- function(
  data,
  parcelstoexclude_bytarget = NULL,
  exclusion_h5 = NULL,
  orientation = c("nodes_by_conditions", "conditions_by_nodes"),
  fill_value = 0,
  sparse = TRUE,
  mask_format = c("list", "dgCMatrix")
) {
  orientation <- match.arg(orientation)
  mask_format <- match.arg(mask_format)
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

  source_index <- .af_build_source_index(
    exclusions,
    n_nodes = n_nodes,
    sparse = sparse,
    mask_format = mask_format
  )

  if (isTRUE(sparse)) {
    density <- source_index$sparsity$density
    out <- array(fill_value, dim = c(n_nodes, n_nodes, n_cond))

    if (is.finite(density) && density <= 0.35) {
      # Low-density path: write only allowed source-target pairs.
      for (target in seq_len(n_nodes)) {
        src <- source_index$sources[[target]]
        if (length(src) > 0L) {
          out[target, src, ] <- mat[src, , drop = FALSE]
        }
      }
    } else {
      # Dense-ish sparse masks are faster to materialize then blank-out.
      out <- array(rep(mat, each = n_nodes), dim = c(n_nodes, n_nodes, n_cond))
      all_nodes <- seq_len(n_nodes)
      for (target in seq_len(n_nodes)) {
        src <- source_index$sources[[target]]
        excluded <- setdiff(all_nodes, src)
        if (length(excluded) > 0L) {
          out[target, excluded, ] <- fill_value
        }
      }
    }
  } else {
    # Dense baseline path: materialize full target-source tensor, then apply mask.
    out <- array(rep(mat, each = n_nodes), dim = c(n_nodes, n_nodes, n_cond))
    for (target in seq_len(n_nodes)) {
      excluded <- which(source_index$mask[target, ] == 0)
      if (length(excluded) > 0L) {
        out[target, excluded, ] <- fill_value
      }
    }
  }

  for (target in seq_len(n_nodes)) {
    out[target, target, ] <- mat[target, ]
  }

  attr(out, "sparsity_stats") <- source_index$sparsity
  attr(out, "source_mask_format") <- if (isTRUE(sparse)) mask_format else "dense"
  out
}
