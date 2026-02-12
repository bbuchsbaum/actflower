#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1L]]
}

report_json <- arg_value("--report-json", "inst/extdata/benchmarks/stochastic_parity/report.json")
seed_base <- as.integer(arg_value("--seed", "2026"))

suppressPackageStartupMessages(library(actflower))
suppressPackageStartupMessages(library(jsonlite))

set.seed(seed_base)
n_nodes <- 24L
n_time <- 180L
n_cond <- 10L
n_subj <- 12L

rest_group <- array(rnorm(n_nodes * n_time * n_subj), dim = c(n_nodes, n_time, n_subj))
task_group <- array(rnorm(n_nodes * n_cond * n_subj), dim = c(n_nodes, n_cond, n_subj))

grid <- list(
  list(method = "corr"),
  list(method = "multreg")
)

run_nested <- function(seed) {
  actflow_nested_cv(
    data = list(rest_group = rest_group, task_group = task_group),
    outer_folds = 3,
    inner_folds = 3,
    estimator_grid = grid,
    seed = seed,
    save_splits = TRUE,
    selection_metric = "corr",
    use_cpp = TRUE
  )
}

nested_a <- run_nested(seed_base + 10L)
nested_b <- run_nested(seed_base + 10L)
nested_c <- run_nested(seed_base + 11L)

run_uncertainty <- function(seed) {
  fc <- estimate_fc_corr(rest_group, use_cpp = TRUE)
  actflow_uncertainty(
    act_group = task_group,
    fc_group = fc,
    metric = c("corr", "R2", "mae"),
    n_boot = 120,
    conf_level = 0.95,
    seed = seed,
    use_cpp = TRUE
  )
}

unc_a <- run_uncertainty(seed_base + 20L)
unc_b <- run_uncertainty(seed_base + 20L)
unc_c <- run_uncertainty(seed_base + 21L)

checks <- list(
  nested_cv = list(
    same_seed_split_ids_identical = identical(
      nested_a$split_artifacts$outer_fold_ids,
      nested_b$split_artifacts$outer_fold_ids
    ),
    same_seed_selected_identical = isTRUE(all.equal(
      nested_a$selected_hyperparams,
      nested_b$selected_hyperparams,
      tolerance = 0
    )),
    same_seed_outer_metrics_identical = isTRUE(all.equal(
      nested_a$outer_metrics,
      nested_b$outer_metrics,
      tolerance = 0
    )),
    different_seed_split_ids_differ = !identical(
      nested_a$split_artifacts$outer_fold_ids,
      nested_c$split_artifacts$outer_fold_ids
    )
  ),
  uncertainty = list(
    same_seed_draws_identical = isTRUE(all.equal(
      unc_a$bootstrap_draws,
      unc_b$bootstrap_draws,
      tolerance = 0
    )),
    same_seed_interval_identical = isTRUE(all.equal(
      unc_a$interval,
      unc_b$interval,
      tolerance = 0
    )),
    different_seed_draws_differ = !isTRUE(all.equal(
      unc_a$bootstrap_draws,
      unc_c$bootstrap_draws,
      tolerance = 0
    ))
  )
)

flat <- unlist(checks, recursive = TRUE, use.names = TRUE)
flat_logical <- as.logical(flat)
checks_total <- length(flat_logical)
checks_passed <- sum(flat_logical, na.rm = TRUE)

report <- list(
  schema_version = "1.0.0",
  tool = "stochastic_parity_runner",
  created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  seed_base = seed_base,
  dimensions = list(
    nodes = n_nodes,
    timepoints = n_time,
    conditions = n_cond,
    subjects = n_subj
  ),
  checks = checks,
  summary = list(
    checks_total = checks_total,
    checks_passed = checks_passed,
    pass_ratio = if (checks_total > 0) checks_passed / checks_total else NaN
  )
)

dir.create(dirname(report_json), recursive = TRUE, showWarnings = FALSE)
write_json(report, path = report_json, pretty = TRUE, auto_unbox = TRUE, na = "null")

cat("Wrote stochastic parity report: ", normalizePath(report_json, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat(toJSON(report$summary, pretty = TRUE, auto_unbox = TRUE), "\n")
