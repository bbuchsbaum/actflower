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

mode <- arg_value("--mode", "synthetic")
out_path <- arg_value("--out", "inst/extdata/parity_refs/synthetic_refs.h5")
overwrite <- as_flag(arg_value("--overwrite", "false"), default = FALSE)

if (file.exists(out_path) && !overwrite) {
  stop("Output already exists. Re-run with --overwrite true: ", out_path)
}

source_local_package <- function() {
  r_files <- list.files("R", pattern = "[.]R$", full.names = TRUE)
  for (f in r_files) source(f)
}

if (identical(mode, "python_parity")) {
  py <- arg_value("--python", Sys.which("python3"))
  n_subj <- as.integer(arg_value("--n-subj", "5"))
  script <- "inst/scripts/build_parity_refs.py"

  if (!nzchar(py)) stop("Python executable not found. Set --python explicitly.")
  if (!file.exists(script)) stop("Parity builder not found: ", script)
  if (is.na(n_subj) || n_subj < 1) stop("--n-subj must be a positive integer.")

  cmd <- c(
    script,
    "--repo-root", normalizePath(".", mustWork = TRUE),
    "--n-subj", as.character(n_subj),
    "--out", normalizePath(out_path, winslash = "/", mustWork = FALSE)
  )

  status <- system2(py, cmd)
  if (!identical(status, 0L)) stop("Python parity fixture build failed with status ", status)
  cat("Wrote parity fixture with Python builder to ", out_path, "\n", sep = "")
  quit(status = 0L)
}

if (!identical(mode, "synthetic")) {
  stop("Unknown --mode value: ", mode, ". Expected 'synthetic' or 'python_parity'.")
}

n_nodes <- as.integer(arg_value("--nodes", "64"))
n_time <- as.integer(arg_value("--time", "300"))
n_cond <- as.integer(arg_value("--conditions", "24"))
n_subj <- as.integer(arg_value("--subjects", "5"))
seed <- as.integer(arg_value("--seed", "42"))

if (any(is.na(c(n_nodes, n_time, n_cond, n_subj, seed)))) {
  stop("Numeric arguments must be parseable integers.")
}
if (any(c(n_nodes, n_time, n_cond, n_subj) < 1L)) {
  stop("--nodes, --time, --conditions, and --subjects must all be >= 1.")
}

set.seed(seed)

load_package <- as_flag(Sys.getenv("ACTFLOWER_LOAD_PACKAGE", "false"), default = FALSE)
if (load_package) {
  suppressPackageStartupMessages(library(actflower))
} else {
  source_local_package()
}

if (!dir.exists(dirname(out_path))) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
}

restdata <- array(rnorm(n_nodes * n_time * n_subj), dim = c(n_nodes, n_time, n_subj))
taskbeta <- array(rnorm(n_nodes * n_cond * n_subj), dim = c(n_nodes, n_cond, n_subj))

fc_corr <- array(0, dim = c(n_nodes, n_nodes, n_subj))
fc_multreg <- array(0, dim = c(n_nodes, n_nodes, n_subj))
pred_corr <- array(0, dim = c(n_nodes, n_cond, n_subj))
pred_multreg <- array(0, dim = c(n_nodes, n_cond, n_subj))

for (s in seq_len(n_subj)) {
  ts_s <- restdata[, , s, drop = TRUE]
  fc_corr[, , s] <- estimate_fc_corr(ts_s, orientation = "nodes_by_time")
  fc_multreg[, , s] <- estimate_fc_multreg(ts_s, orientation = "nodes_by_time", use_cpp = TRUE)

  for (c in seq_len(n_cond)) {
    a <- taskbeta[, c, s]
    pred_corr[, c, s] <- actflow_predict(a, fc_corr[, , s], transfer = "linear")
    pred_multreg[, c, s] <- actflow_predict(a, fc_multreg[, , s], transfer = "linear")
  }
}

metrics_corr <- model_compare(taskbeta, pred_corr, comparison = "fullcompare_compthenavg", use_cpp = TRUE)
metrics_multreg <- model_compare(taskbeta, pred_multreg, comparison = "fullcompare_compthenavg", use_cpp = TRUE)

write_actflower_h5(out_path, "meta/seed", as.integer(seed))
write_actflower_h5(out_path, "meta/nodes", as.integer(n_nodes))
write_actflower_h5(out_path, "meta/time", as.integer(n_time))
write_actflower_h5(out_path, "meta/conditions", as.integer(n_cond))
write_actflower_h5(out_path, "meta/subjects", as.integer(n_subj))

write_actflower_h5(out_path, "restdata", restdata)
write_actflower_h5(out_path, "taskbeta", taskbeta)
write_actflower_h5(out_path, "fc_corr", fc_corr)
write_actflower_h5(out_path, "fc_multreg", fc_multreg)
write_actflower_h5(out_path, "pred_corr", pred_corr)
write_actflower_h5(out_path, "pred_multreg", pred_multreg)

write_actflower_h5(out_path, "metrics_corr/corr_vals", metrics_corr$corr_vals)
write_actflower_h5(out_path, "metrics_corr/R2_vals", metrics_corr$R2_vals)
write_actflower_h5(out_path, "metrics_corr/mae_vals", metrics_corr$mae_vals)

write_actflower_h5(out_path, "metrics_multreg/corr_vals", metrics_multreg$corr_vals)
write_actflower_h5(out_path, "metrics_multreg/R2_vals", metrics_multreg$R2_vals)
write_actflower_h5(out_path, "metrics_multreg/mae_vals", metrics_multreg$mae_vals)

cat("Wrote synthetic fixture to ", out_path, "\n", sep = "")
