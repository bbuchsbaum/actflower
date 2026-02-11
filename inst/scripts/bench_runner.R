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

native_status <- list(
  multreg = !inherits(try(getNativeSymbolInfo("_actflower_multreg_fc_cpp", PACKAGE = "actflower"), silent = TRUE), "try-error"),
  actflow = !inherits(try(getNativeSymbolInfo("_actflower_actflow_predict_batch_cpp", PACKAGE = "actflower"), silent = TRUE), "try-error"),
  compare = !inherits(try(getNativeSymbolInfo("_actflower_compare_fullcomp_cpp", PACKAGE = "actflower"), silent = TRUE), "try-error")
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

speedup <- function(r_sec, cpp_sec) {
  if (!is.finite(r_sec) || !is.finite(cpp_sec) || cpp_sec <= 0) return(NA_real_)
  r_sec / cpp_sec
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
    actflow_test = list(r_sec = actflow_r_sec, cpp_sec = actflow_cpp_sec, speedup = speedup(actflow_r_sec, actflow_cpp_sec))
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

  if (isTRUE(native_status$multreg) && is.finite(req_multreg)) {
    sp <- as.numeric(results$benchmarks$multreg$speedup)
    if (!is.finite(sp)) {
      stop("multreg speedup could not be estimated (non-finite value).")
    }
    if (sp < req_multreg) {
      stop(sprintf("multreg speedup %.3f < required %.3f", results$benchmarks$multreg$speedup, req_multreg))
    }
  }

  if (isTRUE(native_status$actflow) && is.finite(req_actflow)) {
    sp <- as.numeric(results$benchmarks$actflow_test$speedup)
    if (!is.finite(sp)) {
      stop("actflow_test speedup could not be estimated (non-finite value).")
    }
    if (sp < req_actflow) {
      stop(sprintf("actflow_test speedup %.3f < required %.3f", results$benchmarks$actflow_test$speedup, req_actflow))
    }
  }

  if (!isTRUE(native_status$multreg) || !isTRUE(native_status$actflow)) {
    cat("Native symbols not available; threshold checks skipped for missing kernels.\n")
  }
}
