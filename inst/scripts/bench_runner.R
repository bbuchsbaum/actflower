#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1L]]
}

as_flag <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(x) %in% c("1", "true", "yes", "y")
}

fixture <- arg_value("--fixture", "inst/extdata/parity_refs/hcp_parity_refs.h5")
iters <- as.integer(arg_value("--iterations", "3"))
out_file <- arg_value("--output", "inst/extdata/benchmarks/latest.json")
threshold_file <- arg_value("--threshold-file", "tools/benchmark_thresholds.json")
enforce <- as_flag(arg_value("--enforce", "false"), default = FALSE)

if (is.na(iters) || iters < 1) {
  stop("--iterations must be a positive integer")
}

load_package <- as_flag(Sys.getenv("ACTFLOWER_LOAD_PACKAGE", "false"), default = FALSE)
if (load_package) {
  suppressPackageStartupMessages(library(actflower))
} else {
  r_files <- list.files("R", pattern = "[.]R$", full.names = TRUE)
  for (f in r_files) source(f)
}

if (!file.exists(fixture)) {
  stop("Benchmark fixture not found: ", fixture)
}

from_py <- function(x) {
  if (is.array(x) && length(dim(x)) >= 2) {
    return(aperm(x, rev(seq_along(dim(x)))))
  }
  x
}

rest <- from_py(read_actflower_h5(fixture, "restdata"))
taskbeta <- from_py(read_actflower_h5(fixture, "taskbeta"))

# Keep benchmark runtime bounded.
n_subj <- min(dim(rest)[3], 3L)
rest <- rest[, , seq_len(n_subj), drop = FALSE]
taskbeta <- taskbeta[, , seq_len(n_subj), drop = FALSE]

time_per_rep <- function(expr_sub, eval_env, min_elapsed = 0.05, max_reps = 512L) {
  reps <- 1L
  repeat {
    gc(FALSE)
    t0 <- Sys.time()
    for (k in seq_len(reps)) {
      eval(expr_sub, envir = eval_env)
    }
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (elapsed >= min_elapsed || reps >= max_reps) {
      if (reps <= 0 || !is.finite(elapsed)) return(NA_real_)
      return(elapsed / reps)
    }
    reps <- reps * 2L
  }
}

median_time <- function(expr, n = 3L, min_elapsed = 0.05) {
  expr_sub <- substitute(expr)
  eval_env <- parent.frame()
  vals <- numeric(n)
  for (i in seq_len(n)) {
    vals[i] <- time_per_rep(expr_sub, eval_env, min_elapsed = min_elapsed)
  }
  as.numeric(stats::median(vals, na.rm = TRUE))
}

make_sparse_exclusions <- function(n_nodes, density = 0.08, seed = 1L) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  out <- vector("list", n_nodes)
  all_nodes <- seq_len(n_nodes)
  n_keep <- max(1L, as.integer(floor((n_nodes - 1L) * density)))

  for (target in all_nodes) {
    pool <- all_nodes[all_nodes != target]
    keep <- sort(sample(pool, size = n_keep, replace = FALSE))
    out[[target]] <- setdiff(pool, keep)
  }
  out
}

alloc_bytes <- function(expr) {
  tf <- tempfile("actflower-mem-", fileext = ".log")
  on.exit(unlink(tf), add = TRUE)
  out <- NA_real_
  ok <- TRUE
  tryCatch(Rprofmem(tf), error = function(e) {
    ok <<- FALSE
  })
  if (!ok) return(NA_real_)
  on.exit(Rprofmem(NULL), add = TRUE)

  force(expr)
  Rprofmem(NULL)

  lines <- readLines(tf, warn = FALSE)
  vals <- suppressWarnings(as.numeric(sub(" .*", "", lines)))
  sum(vals, na.rm = TRUE)
}

dense_noncircular_multreg_reference <- function(data_nodes_by_time, exclusions) {
  data <- as.matrix(data_nodes_by_time)
  n_nodes <- nrow(data)
  all_nodes <- seq_len(n_nodes)
  fc <- matrix(0, n_nodes, n_nodes)

  for (target in all_nodes) {
    other <- all_nodes[all_nodes != target]
    x_full <- t(data[other, , drop = FALSE])
    y <- as.numeric(data[target, ])

    excluded <- exclusions[[target]]
    if (length(excluded) > 0L) {
      ex_idx <- match(excluded, other, nomatch = 0L)
      ex_idx <- ex_idx[ex_idx > 0L]
      if (length(ex_idx) > 0L) {
        x_full[, ex_idx] <- 0
      }
    }

    fit <- stats::lm.fit(x = cbind(1, x_full), y = y)
    beta <- fit$coefficients[-1]
    beta[is.na(beta)] <- 0
    fc[target, other] <- beta
    if (length(excluded) > 0L) {
      fc[target, excluded] <- 0
    }
  }

  fc
}

dense_noncircular_activity_reference <- function(data_nodes_by_conditions, exclusions, fill_value = 0) {
  mat <- as.matrix(data_nodes_by_conditions)
  n_nodes <- nrow(mat)
  n_cond <- ncol(mat)

  scratch <- array(rep(mat, each = n_nodes), dim = c(n_nodes, n_nodes, n_cond))
  mask3 <- array(TRUE, dim = c(n_nodes, n_nodes, n_cond))
  for (target in seq_len(n_nodes)) {
    excluded <- sort(unique(c(target, exclusions[[target]])))
    if (length(excluded) > 0L) {
      mask3[target, excluded, ] <- FALSE
    }
  }
  out <- scratch
  out[!mask3] <- fill_value
  rm(mask3)

  for (target in seq_len(n_nodes)) out[target, target, ] <- mat[target, ]

  out
}

native_status <- list(
  multreg = !inherits(try(getNativeSymbolInfo("_actflower_multreg_fc_cpp", PACKAGE = "actflower"), silent = TRUE), "try-error"),
  actflow = !inherits(try(getNativeSymbolInfo("_actflower_actflow_predict_batch_cpp", PACKAGE = "actflower"), silent = TRUE), "try-error"),
  compare = !inherits(try(getNativeSymbolInfo("_actflower_compare_fullcomp_cpp", PACKAGE = "actflower"), silent = TRUE), "try-error"),
  fullcompare_fused = !inherits(try(getNativeSymbolInfo("_actflower_actflow_fullcomp_batch_cpp", PACKAGE = "actflower"), silent = TRUE), "try-error")
)

# Benchmark multreg FC estimation
multreg_r_sec <- median_time({
  for (s in seq_len(n_subj)) {
    estimate_fc_multreg(rest[, , s], use_cpp = FALSE)
  }
}, n = iters)

multreg_cpp_sec <- median_time({
  for (s in seq_len(n_subj)) {
    estimate_fc_multreg(rest[, , s], use_cpp = TRUE)
  }
}, n = iters)

fc_for_test <- array(NA_real_, dim = c(dim(rest)[1], dim(rest)[1], n_subj))
for (s in seq_len(n_subj)) {
  fc_for_test[, , s] <- estimate_fc_multreg(rest[, , s], use_cpp = TRUE)
}

# Benchmark end-to-end actflow_test (fixed FC, fullcompare mode)
actflow_r_sec <- median_time({
  actflow_test(taskbeta, fc_for_test, comparison = "fullcompare_compthenavg", use_cpp = FALSE)
}, n = iters)

actflow_cpp_sec <- median_time({
  actflow_test(taskbeta, fc_for_test, comparison = "fullcompare_compthenavg", use_cpp = TRUE)
}, n = iters)

## Sparse noncircular benchmark profile (low mask density)
nc_nodes <- 240L
nc_time <- 260L
nc_cond <- 28L
set.seed(2026)
rest_nc <- matrix(rnorm(nc_nodes * nc_time), nrow = nc_nodes, ncol = nc_time)
task_nc <- matrix(rnorm(nc_nodes * nc_cond), nrow = nc_nodes, ncol = nc_cond)
exclusions_nc <- make_sparse_exclusions(nc_nodes, density = 0.08, seed = 3030)

conn_sparse_probe <- calcconn_parcelwise_noncircular(
  rest_nc,
  connmethod = "multreg",
  parcelstoexclude_bytarget = exclusions_nc,
  orientation = "nodes_by_time",
  sparse = TRUE,
  mask_format = "list"
)
nc_density <- as.numeric(attr(conn_sparse_probe, "sparsity_stats")$density)

conn_sparse_sec <- median_time({
  calcconn_parcelwise_noncircular(
    rest_nc,
    connmethod = "multreg",
    parcelstoexclude_bytarget = exclusions_nc,
    orientation = "nodes_by_time",
    sparse = TRUE,
    mask_format = "list"
  )
}, n = iters)

conn_dense_sec <- median_time({
  dense_noncircular_multreg_reference(rest_nc, exclusions_nc)
}, n = iters)

activity_sparse_sec <- median_time({
  calcactivity_parcelwise_noncircular(
    task_nc,
    parcelstoexclude_bytarget = exclusions_nc,
    orientation = "nodes_by_conditions",
    fill_value = 0,
    sparse = TRUE,
    mask_format = "list"
  )
}, n = iters)

activity_dense_sec <- median_time({
  dense_noncircular_activity_reference(task_nc, exclusions_nc, fill_value = 0)
}, n = iters)

activity_sparse_alloc <- alloc_bytes({
  calcactivity_parcelwise_noncircular(
    task_nc,
    parcelstoexclude_bytarget = exclusions_nc,
    orientation = "nodes_by_conditions",
    fill_value = 0,
    sparse = TRUE,
    mask_format = "list"
  )
})

activity_dense_alloc <- alloc_bytes({
  dense_noncircular_activity_reference(task_nc, exclusions_nc, fill_value = 0)
})

speedup <- function(r_sec, cpp_sec) {
  if (!is.finite(r_sec) || !is.finite(cpp_sec) || cpp_sec <= 0) return(NA_real_)
  r_sec / cpp_sec
}

sparse_speedup <- function(dense_sec, sparse_sec) {
  if (!is.finite(dense_sec) || !is.finite(sparse_sec) || sparse_sec <= 0) return(NA_real_)
  dense_sec / sparse_sec
}

ratio <- function(num, den) {
  if (!is.finite(num) || !is.finite(den) || den <= 0) return(NA_real_)
  num / den
}

results <- list(
  timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
  load_mode = if (load_package) "package" else "source",
  fixture = fixture,
  n_subjects = n_subj,
  iterations = iters,
  native_available = native_status,
  benchmarks = list(
    multreg = list(r_sec = multreg_r_sec, cpp_sec = multreg_cpp_sec, speedup = speedup(multreg_r_sec, multreg_cpp_sec)),
    actflow_test = list(r_sec = actflow_r_sec, cpp_sec = actflow_cpp_sec, speedup = speedup(actflow_r_sec, actflow_cpp_sec)),
    noncircular_sparse = list(
      density = nc_density,
      conn = list(
        sparse_sec = conn_sparse_sec,
        dense_sec = conn_dense_sec,
        speedup_dense_over_sparse = sparse_speedup(conn_dense_sec, conn_sparse_sec)
      ),
      activity = list(
        sparse_sec = activity_sparse_sec,
        dense_sec = activity_dense_sec,
        speedup_dense_over_sparse = sparse_speedup(activity_dense_sec, activity_sparse_sec),
        sparse_alloc_bytes = activity_sparse_alloc,
        dense_alloc_bytes = activity_dense_alloc,
        alloc_ratio_sparse_over_dense = ratio(activity_sparse_alloc, activity_dense_alloc)
      )
    )
  )
)

if (!dir.exists(dirname(out_file))) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
}

if (requireNamespace("jsonlite", quietly = TRUE)) {
  jsonlite::write_json(results, out_file, auto_unbox = TRUE, pretty = TRUE)
} else {
  saveRDS(results, sub("[.]json$", ".rds", out_file))
}

cat("Wrote benchmark results to", out_file, "\n")

if (enforce) {
  if (!file.exists(threshold_file)) {
    stop("--enforce requested but threshold file missing: ", threshold_file)
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("--enforce requires jsonlite to parse threshold file.")
  }

  th <- jsonlite::fromJSON(threshold_file)
  req_multreg <- if (!is.null(th$multreg_min_speedup)) as.numeric(th$multreg_min_speedup) else NA_real_
  req_actflow <- if (!is.null(th$actflow_test_min_speedup)) as.numeric(th$actflow_test_min_speedup) else NA_real_
  req_nc_density <- if (!is.null(th$noncircular_sparse_density_max)) as.numeric(th$noncircular_sparse_density_max) else NA_real_
  req_nc_conn <- if (!is.null(th$noncircular_conn_sparse_min_speedup)) as.numeric(th$noncircular_conn_sparse_min_speedup) else NA_real_
  req_nc_act <- if (!is.null(th$noncircular_activity_sparse_min_speedup)) as.numeric(th$noncircular_activity_sparse_min_speedup) else NA_real_
  req_nc_mem <- if (!is.null(th$noncircular_activity_sparse_max_mem_ratio)) as.numeric(th$noncircular_activity_sparse_max_mem_ratio) else NA_real_

  if (isTRUE(native_status$multreg) && is.finite(req_multreg)) {
    sp <- as.numeric(results$benchmarks$multreg$speedup)
    if (!is.finite(sp)) {
      stop("multreg speedup could not be estimated (non-finite value).")
    }
    if (sp < req_multreg) {
      stop(sprintf("multreg speedup %.3f < required %.3f", results$benchmarks$multreg$speedup, req_multreg))
    }
  }

  if (isTRUE(native_status$fullcompare_fused) && is.finite(req_actflow)) {
    sp <- as.numeric(results$benchmarks$actflow_test$speedup)
    if (!is.finite(sp)) {
      stop("actflow_test speedup could not be estimated (non-finite value).")
    }
    if (sp < req_actflow) {
      stop(sprintf("actflow_test speedup %.3f < required %.3f", results$benchmarks$actflow_test$speedup, req_actflow))
    }
  }

  if (!isTRUE(native_status$multreg) || !isTRUE(native_status$fullcompare_fused)) {
    cat("Native symbols not available; threshold checks skipped for missing kernels.\n")
  }

  if (is.finite(req_nc_density)) {
    got <- as.numeric(results$benchmarks$noncircular_sparse$density)
    if (!is.finite(got)) stop("noncircular sparse density could not be estimated.")
    if (got > req_nc_density) {
      stop(sprintf("noncircular sparse density %.3f > required max %.3f", got, req_nc_density))
    }
  }

  if (is.finite(req_nc_conn)) {
    got <- as.numeric(results$benchmarks$noncircular_sparse$conn$speedup_dense_over_sparse)
    if (!is.finite(got)) stop("noncircular sparse conn speedup could not be estimated.")
    if (got < req_nc_conn) {
      stop(sprintf("noncircular conn sparse speedup %.3f < required %.3f", got, req_nc_conn))
    }
  }

  if (is.finite(req_nc_act)) {
    got <- as.numeric(results$benchmarks$noncircular_sparse$activity$speedup_dense_over_sparse)
    if (!is.finite(got)) stop("noncircular sparse activity speedup could not be estimated.")
    if (got < req_nc_act) {
      stop(sprintf("noncircular activity sparse speedup %.3f < required %.3f", got, req_nc_act))
    }
  }

  if (is.finite(req_nc_mem)) {
    got <- as.numeric(results$benchmarks$noncircular_sparse$activity$alloc_ratio_sparse_over_dense)
    if (!is.finite(got)) stop("noncircular sparse activity allocation ratio could not be estimated.")
    if (got > req_nc_mem) {
      stop(sprintf("noncircular activity sparse allocation ratio %.3f > required max %.3f", got, req_nc_mem))
    }
  }
}
