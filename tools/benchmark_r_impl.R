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

input_h5 <- arg_value("--input-h5", NULL)
out_h5 <- arg_value("--out-h5", NULL)
out_json <- arg_value("--out-json", NULL)
repeats <- as.integer(arg_value("--repeats", "3"))
use_cpp <- as_flag(arg_value("--use-cpp", "true"), default = TRUE)
include_combined <- as_flag(arg_value("--include-combined", "false"), default = FALSE)
include_noncircular <- as_flag(arg_value("--include-noncircular", "false"), default = FALSE)

if (is.null(input_h5) || is.null(out_h5) || is.null(out_json)) {
  stop("Required args: --input-h5, --out-h5, --out-json")
}
if (is.na(repeats) || repeats < 1L) {
  stop("--repeats must be a positive integer")
}

suppressPackageStartupMessages(library(actflower))

timed <- function(fun, repeats = 3L) {
  times <- numeric(repeats)
  result <- NULL
  for (i in seq_len(repeats)) {
    gc(FALSE)
    t0 <- Sys.time()
    result <- fun()
    times[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }
  list(
    result = result,
    times = times,
    median_sec = as.numeric(stats::median(times)),
    mean_sec = as.numeric(mean(times))
  )
}

rest <- read_actflower_h5(input_h5, "restdata")
task <- read_actflower_h5(input_h5, "taskbeta")

from_python <- function(x) {
  if (is.array(x) && length(dim(x)) >= 2L) {
    return(aperm(x, rev(seq_along(dim(x)))))
  }
  x
}

to_python <- from_python

make_ring_exclusions <- function(n_nodes, k = 2L) {
  out <- vector("list", n_nodes)
  for (i in seq_len(n_nodes)) {
    offsets <- setdiff(seq.int(-k, k), 0L)
    idx <- ((i + offsets - 1L) %% n_nodes) + 1L
    out[[i]] <- sort(unique(as.integer(idx)))
  }
  out
}

rest <- from_python(rest)
task <- from_python(task)

if (length(dim(rest)) != 3L || length(dim(task)) != 3L) {
  stop("Expected 3D arrays for restdata and taskbeta.")
}

n_nodes <- dim(rest)[1]
n_cond <- dim(task)[2]
n_subj <- dim(rest)[3]

if (dim(task)[1] != n_nodes || dim(task)[3] != n_subj) {
  stop("Incompatible restdata/taskbeta shapes.")
}

run_fc_corr <- function() {
  out <- array(0, dim = c(n_nodes, n_nodes, n_subj))
  for (s in seq_len(n_subj)) {
    out[, , s] <- estimate_fc_corr(rest[, , s], orientation = "nodes_by_time", use_cpp = use_cpp)
  }
  out
}

run_fc_multreg <- function() {
  out <- array(0, dim = c(n_nodes, n_nodes, n_subj))
  for (s in seq_len(n_subj)) {
    out[, , s] <- estimate_fc_multreg(
      rest[, , s],
      orientation = "nodes_by_time",
      use_cpp = use_cpp
    )
  }
  out
}

pred_from_fc <- function(fc) {
  out <- array(0, dim = c(n_nodes, n_cond, n_subj))
  for (s in seq_len(n_subj)) {
    for (c in seq_len(n_cond)) {
      out[, c, s] <- actflow_predict(
        task[, c, s],
        fc[, , s],
        separate_by_target = FALSE,
        transfer = "linear"
      )
    }
  }
  out
}

corr_run <- timed(run_fc_corr, repeats = repeats)
multreg_run <- timed(run_fc_multreg, repeats = repeats)

fc_corr <- corr_run$result
fc_multreg <- multreg_run$result

pred_corr_run <- timed(function() pred_from_fc(fc_corr), repeats = repeats)
pred_multreg_run <- timed(function() pred_from_fc(fc_multreg), repeats = repeats)
fullcompare_multreg_run <- timed(
  function() {
    actflow_test(
      task,
      fc_multreg,
      comparison = "fullcompare_compthenavg",
      use_cpp = use_cpp
    )
  },
  repeats = repeats
)
fullcompare_multreg_metrics <- fullcompare_multreg_run$result$model_compare_output

combined_payload <- list()
if (include_combined) {
  run_fc_combined <- function() {
    out <- array(0, dim = c(n_nodes, n_nodes, n_subj))
    for (s in seq_len(n_subj)) {
      out[, , s] <- estimate_fc_combined(rest[, , s], orientation = "nodes_by_time")
    }
    out
  }
  combined_run <- timed(run_fc_combined, repeats = repeats)
  fc_combined <- combined_run$result
  pred_combined_run <- timed(function() pred_from_fc(fc_combined), repeats = repeats)

  write_actflower_h5(out_h5, "fc_combined", to_python(fc_combined))
  write_actflower_h5(out_h5, "pred_combined", to_python(pred_combined_run$result))

  combined_payload <- list(
    fc_combined = list(
      median_sec = combined_run$median_sec,
      mean_sec = combined_run$mean_sec,
      times = combined_run$times
    ),
    pred_combined = list(
      median_sec = pred_combined_run$median_sec,
      mean_sec = pred_combined_run$mean_sec,
      times = pred_combined_run$times
    )
  )
}

noncircular_payload <- list()
if (include_noncircular) {
  exclusions <- make_ring_exclusions(n_nodes, k = 2L)

  run_fc_noncircular <- function() {
    out <- array(0, dim = c(n_nodes, n_nodes, n_subj))
    for (s in seq_len(n_subj)) {
      out[, , s] <- calcconn_parcelwise_noncircular(
        rest[, , s],
        connmethod = "multreg",
        parcelstoexclude_bytarget = exclusions,
        orientation = "nodes_by_time"
      )
    }
    out
  }

  run_activity_noncircular <- function() {
    out <- array(0, dim = c(n_nodes, n_nodes, n_cond, n_subj))
    for (s in seq_len(n_subj)) {
      out[, , , s] <- calcactivity_parcelwise_noncircular(
        task[, , s],
        parcelstoexclude_bytarget = exclusions,
        orientation = "nodes_by_conditions",
        fill_value = 0
      )
    }
    out
  }

  noncirc_fc_run <- timed(run_fc_noncircular, repeats = repeats)
  noncirc_act_run <- timed(run_activity_noncircular, repeats = repeats)

  write_actflower_h5(out_h5, "fc_noncircular_multreg", to_python(noncirc_fc_run$result))
  write_actflower_h5(out_h5, "activity_noncircular", to_python(noncirc_act_run$result))

  noncircular_payload <- list(
    fc_noncircular_multreg = list(
      median_sec = noncirc_fc_run$median_sec,
      mean_sec = noncirc_fc_run$mean_sec,
      times = noncirc_fc_run$times
    ),
    activity_noncircular = list(
      median_sec = noncirc_act_run$median_sec,
      mean_sec = noncirc_act_run$mean_sec,
      times = noncirc_act_run$times
    )
  )
}

write_actflower_h5(out_h5, "fc_corr", to_python(fc_corr))
write_actflower_h5(out_h5, "fc_multreg", to_python(fc_multreg))
write_actflower_h5(out_h5, "pred_corr", to_python(pred_corr_run$result))
write_actflower_h5(out_h5, "pred_multreg", to_python(pred_multreg_run$result))
write_actflower_h5(out_h5, "fullcompare_multreg_corr", fullcompare_multreg_metrics$corr_vals)
write_actflower_h5(out_h5, "fullcompare_multreg_R2", fullcompare_multreg_metrics$R2_vals)
write_actflower_h5(out_h5, "fullcompare_multreg_mae", fullcompare_multreg_metrics$mae_vals)

timings <- list(
  engine = "R_actflower",
  use_cpp = use_cpp,
  include_combined = include_combined,
  include_noncircular = include_noncircular,
  repeats = repeats,
  dimensions = list(
    nodes = n_nodes,
    time = dim(rest)[2],
    conditions = n_cond,
    subjects = n_subj
  ),
  operations = list(
    fc_corr = list(median_sec = corr_run$median_sec, mean_sec = corr_run$mean_sec, times = corr_run$times),
    fc_multreg = list(median_sec = multreg_run$median_sec, mean_sec = multreg_run$mean_sec, times = multreg_run$times),
    pred_corr = list(median_sec = pred_corr_run$median_sec, mean_sec = pred_corr_run$mean_sec, times = pred_corr_run$times),
    pred_multreg = list(median_sec = pred_multreg_run$median_sec, mean_sec = pred_multreg_run$mean_sec, times = pred_multreg_run$times),
    fullcompare_multreg = list(
      median_sec = fullcompare_multreg_run$median_sec,
      mean_sec = fullcompare_multreg_run$mean_sec,
      times = fullcompare_multreg_run$times
    )
  )
)

if (length(combined_payload)) {
  timings$operations <- c(timings$operations, combined_payload)
}
if (length(noncircular_payload)) {
  timings$operations <- c(timings$operations, noncircular_payload)
}

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("jsonlite is required to write benchmark JSON.")
}
jsonlite::write_json(timings, out_json, pretty = TRUE, auto_unbox = TRUE)
cat("Wrote R benchmark outputs to ", out_h5, " and ", out_json, "\n", sep = "")
