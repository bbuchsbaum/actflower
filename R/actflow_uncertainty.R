#' Bootstrap uncertainty for actflow metrics
#' @param act_group Array with dimensions (nodes x conditions x subjects).
#' @param fc_group Array with dimensions (nodes x nodes x subjects) or
#'   (nodes x nodes x conditions x subjects).
#' @param metric Metrics to summarize. Any of `corr`, `R2`, `mae`.
#' @param n_boot Number of bootstrap resamples over subjects.
#' @param conf_level Confidence level for percentile intervals.
#' @param seed Optional RNG seed for reproducible bootstrap draws.
#' @param parallel If `TRUE`, use `future.apply::future_lapply()` when available.
#' @param act_group_test Optional test activations (same dimensions as `act_group`).
#' @param use_cpp Use native kernels when available.
#' @param comparison Comparison type; currently only `fullcompare_compthenavg`.
#' @param transfer Transfer function for actflow prediction.
#' @param ... Reserved for forward compatibility.
#' @return A list with point estimates, interval bounds, bootstrap summaries,
#'   and bootstrap configuration metadata.
#' @export
actflow_uncertainty <- function(
  act_group,
  fc_group,
  metric = c("corr", "R2", "mae"),
  n_boot = 1000,
  conf_level = 0.95,
  seed = NULL,
  parallel = FALSE,
  act_group_test = NULL,
  use_cpp = TRUE,
  comparison = "fullcompare_compthenavg",
  transfer = c("linear", "relu", "sigmoid", "logit"),
  ...
) {
  validate_array3(act_group, arg = "act_group")
  transfer <- match.arg(transfer)
  metric <- unique(match.arg(metric, choices = c("corr", "R2", "mae"), several.ok = TRUE))

  if (!is.numeric(n_boot) || length(n_boot) != 1L || !is.finite(n_boot) || n_boot < 1L) {
    actflower_abort("`n_boot` must be a positive scalar.")
  }
  n_boot <- as.integer(n_boot)

  if (!is.numeric(conf_level) || length(conf_level) != 1L || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1) {
    actflower_abort("`conf_level` must be a scalar in (0, 1).")
  }

  if (!identical(comparison, "fullcompare_compthenavg")) {
    actflower_abort("`actflow_uncertainty()` currently supports only comparison='fullcompare_compthenavg'.")
  }

  fc_dim <- dim(fc_group)
  if (length(fc_dim) < 3L || length(fc_dim) > 4L) {
    actflower_abort("`fc_group` must be a 3D or 4D numeric array.")
  }
  if (fc_dim[1] != dim(act_group)[1] || fc_dim[2] != dim(act_group)[1]) {
    actflower_abort("`fc_group` first two dimensions must match node dimension of `act_group`.")
  }
  if (fc_dim[length(fc_dim)] != dim(act_group)[3]) {
    actflower_abort("`fc_group` subject dimension must match `act_group`.")
  }
  if (length(fc_dim) == 4L && fc_dim[3] != dim(act_group)[2]) {
    actflower_abort("4D `fc_group` must have condition dimension matching `act_group`.")
  }

  target <- if (is.null(act_group_test)) act_group else act_group_test
  if (!is.null(act_group_test)) {
    validate_array3(target, arg = "act_group_test")
    if (!all(dim(target) == dim(act_group))) {
      actflower_abort("`act_group_test` must have the same dimensions as `act_group`.")
    }
  }

  .af_set_seed_if_needed(seed)

  observed <- .af_collect_metrics_for_uncertainty(
    act = act_group,
    fc = fc_group,
    target = target,
    use_cpp = use_cpp,
    comparison = comparison,
    transfer = transfer
  )

  observed_metric_vals <- observed[metric]
  point <- vapply(observed_metric_vals, function(v) mean(v, na.rm = TRUE), numeric(1))

  n_subj <- dim(act_group)[3]
  boot_fun <- function(i) {
    idx <- sample.int(n_subj, size = n_subj, replace = TRUE)
    act_b <- act_group[, , idx, drop = FALSE]
    target_b <- target[, , idx, drop = FALSE]
    fc_b <- if (length(fc_dim) == 3L) {
      fc_group[, , idx, drop = FALSE]
    } else {
      fc_group[, , , idx, drop = FALSE]
    }

    met <- .af_collect_metrics_for_uncertainty(
      act = act_b,
      fc = fc_b,
      target = target_b,
      use_cpp = use_cpp,
      comparison = comparison,
      transfer = transfer
    )
    vapply(met[metric], function(v) mean(v, na.rm = TRUE), numeric(1))
  }

  draws <- .af_bootstrap_draws(boot_fun, n_boot = n_boot, n_metric = length(metric), parallel = parallel)
  colnames(draws) <- metric

  alpha <- (1 - conf_level) / 2
  interval <- data.frame(
    metric = metric,
    lower = apply(draws, 2, function(x) stats::quantile(x, probs = alpha, na.rm = TRUE, names = FALSE)),
    upper = apply(draws, 2, function(x) stats::quantile(x, probs = 1 - alpha, na.rm = TRUE, names = FALSE)),
    conf_level = rep(conf_level, length(metric)),
    method = rep("percentile-bootstrap", length(metric)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  distribution_summary <- data.frame(
    metric = metric,
    mean = apply(draws, 2, function(x) mean(x, na.rm = TRUE)),
    median = apply(draws, 2, function(x) stats::median(x, na.rm = TRUE)),
    sd = apply(draws, 2, function(x) stats::sd(x, na.rm = TRUE)),
    se = apply(draws, 2, function(x) stats::sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))),
    n_valid = apply(draws, 2, function(x) sum(is.finite(x))),
    n_boot = rep(n_boot, length(metric)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  list(
    point = point,
    interval = interval,
    distribution_summary = distribution_summary,
    bootstrap_config = list(
      n_boot = n_boot,
      conf_level = conf_level,
      seed = if (is.null(seed)) NA_integer_ else as.integer(seed),
      resample_unit = "subject",
      parallel = isTRUE(parallel),
      use_cpp = isTRUE(use_cpp),
      comparison = comparison,
      transfer = transfer
    ),
    bootstrap_draws = draws
  )
}

.af_collect_metrics_for_uncertainty <- function(act, fc, target, use_cpp, comparison, transfer) {
  out <- actflow_test(
    act_group = act,
    fc_group = fc,
    act_group_test = target,
    comparison = comparison,
    transfer = transfer,
    use_cpp = use_cpp
  )
  cmp <- out$model_compare_output
  list(
    corr = as.numeric(cmp$corr_vals),
    R2 = as.numeric(cmp$R2_vals),
    mae = as.numeric(cmp$mae_vals)
  )
}

.af_bootstrap_draws <- function(boot_fun, n_boot, n_metric, parallel = FALSE) {
  idx <- seq_len(n_boot)
  if (isTRUE(parallel)) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      warning(
        "parallel=TRUE requested but 'future.apply' is not installed; falling back to sequential bootstrap.",
        call. = FALSE
      )
      parallel <- FALSE
    }
  }

  if (isTRUE(parallel)) {
    out <- future.apply::future_lapply(idx, boot_fun)
  } else {
    out <- lapply(idx, boot_fun)
  }

  draws <- matrix(NA_real_, nrow = n_boot, ncol = n_metric)
  for (i in idx) {
    draws[i, ] <- as.numeric(out[[i]])
  }
  draws
}

.af_set_seed_if_needed <- function(seed) {
  if (is.null(seed)) return(invisible(NULL))
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
    actflower_abort("`seed` must be a finite numeric scalar when provided.")
  }
  set.seed(as.integer(seed))
  invisible(NULL)
}
