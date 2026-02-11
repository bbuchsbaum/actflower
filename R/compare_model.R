#' Compare model predictions to target activations
#' @param target Array with dimensions (nodes x conditions x subjects).
#' @param model1 Array with dimensions (nodes x conditions x subjects).
#' @param model2 Optional model2 array.
#' @param comparison One of `fullcompare_compthenavg`,
#'   `conditionwise_compthenavg`, `conditionwise_avgthencomp`,
#'   `nodewise_compthenavg`, or `nodewise_avgthencomp`.
#' @param use_cpp Use native kernels when available.
#' @export
model_compare <- function(
  target,
  model1,
  model2 = NULL,
  comparison = c(
    "fullcompare_compthenavg",
    "conditionwise_compthenavg",
    "conditionwise_avgthencomp",
    "nodewise_compthenavg",
    "nodewise_avgthencomp"
  ),
  use_cpp = TRUE
) {
  comparison <- match.arg(comparison)

  validate_array3(target, arg = "target")
  validate_array3(model1, arg = "model1")
  if (!all(dim(target) == dim(model1))) {
    actflower_abort("`target` and `model1` must have identical dimensions.")
  }

  if (!is.null(model2)) {
    validate_array3(model2, arg = "model2")
    if (!all(dim(target) == dim(model2))) {
      actflower_abort("`target` and `model2` must have identical dimensions.")
    }
  }

  out1 <- .af_model_compare_predicted_to_actual(
    target,
    model1,
    comparison = comparison,
    use_cpp = use_cpp
  )
  if (is.null(model2)) {
    return(out1)
  }

  out2 <- .af_model_compare_predicted_to_actual(
    target,
    model2,
    comparison = comparison,
    use_cpp = use_cpp
  )
  for (nm in names(out2)) {
    out1[[paste0(nm, "_model2")]] <- out2[[nm]]
  }
  out1
}

.af_model_compare_predicted_to_actual <- function(
  target,
  pred,
  comparison = "conditionwise_compthenavg",
  use_cpp = TRUE
) {
  d <- dim(target)
  n_nodes <- d[1]
  n_conds <- d[2]
  n_subjs <- d[3]

  if (comparison == "fullcompare_compthenavg") {
    return(.af_compare_fullcomp(target, pred, use_cpp = use_cpp))
  }

  if (comparison == "conditionwise_compthenavg") {
    corr_vals <- matrix(NA_real_, nrow = n_nodes, ncol = n_subjs)
    r2_vals <- matrix(NA_real_, nrow = n_nodes, ncol = n_subjs)
    mae_vals <- matrix(NA_real_, nrow = n_nodes, ncol = n_subjs)

    for (n in seq_len(n_nodes)) {
      for (s in seq_len(n_subjs)) {
        yt <- target[n, , s]
        yp <- pred[n, , s]
        corr_vals[n, s] <- pearson_vec(yt, yp)
        r2_vals[n, s] <- r2_vec(yt, yp)
        mae_vals[n, s] <- mae_vec(yt, yp)
      }
    }
    return(list(
      corr_vals = corr_vals,
      R2_vals = r2_vals,
      mae_vals = mae_vals,
      corr_conditionwise_compthenavg_bynode = corr_vals,
      R2_conditionwise_compthenavg_bynode = r2_vals,
      mae_conditionwise_compthenavg_bynode = mae_vals
    ))
  }

  if (comparison == "conditionwise_avgthencomp") {
    corr_vals <- numeric(n_nodes)
    r2_vals <- numeric(n_nodes)
    mae_vals <- numeric(n_nodes)

    for (n in seq_len(n_nodes)) {
      yt <- rowMeans(target[n, , , drop = TRUE], na.rm = TRUE)
      yp <- rowMeans(pred[n, , , drop = TRUE], na.rm = TRUE)
      corr_vals[n] <- pearson_vec(yt, yp)
      r2_vals[n] <- r2_vec(yt, yp)
      mae_vals[n] <- mae_vec(yt, yp)
    }

    return(list(
      corr_conditionwise_avgthencomp_bynode = corr_vals,
      R2_conditionwise_avgthencomp_bynode = r2_vals,
      maeAcc_bynode_avgthencomp = mae_vals
    ))
  }

  if (comparison == "nodewise_compthenavg") {
    corr_vals <- matrix(NA_real_, nrow = n_conds, ncol = n_subjs)
    r2_vals <- matrix(NA_real_, nrow = n_conds, ncol = n_subjs)
    mae_vals <- matrix(NA_real_, nrow = n_conds, ncol = n_subjs)

    for (c in seq_len(n_conds)) {
      for (s in seq_len(n_subjs)) {
        yt <- target[, c, s]
        yp <- pred[, c, s]
        corr_vals[c, s] <- pearson_vec(yt, yp)
        r2_vals[c, s] <- r2_vec(yt, yp)
        mae_vals[c, s] <- mae_vec(yt, yp)
      }
    }
    return(list(
      corr_vals = corr_vals,
      R2_vals = r2_vals,
      mae_vals = mae_vals,
      corr_nodewise_compthenavg_bycond = corr_vals,
      R2_nodewise_compthenavg_bycond = r2_vals,
      mae_nodewise_compthenavg_bycond = mae_vals
    ))
  }

  if (comparison == "nodewise_avgthencomp") {
    corr_vals <- numeric(n_conds)
    r2_vals <- numeric(n_conds)
    mae_vals <- numeric(n_conds)

    for (c in seq_len(n_conds)) {
      yt <- rowMeans(target[, c, , drop = TRUE], na.rm = TRUE)
      yp <- rowMeans(pred[, c, , drop = TRUE], na.rm = TRUE)
      corr_vals[c] <- pearson_vec(yt, yp)
      r2_vals[c] <- r2_vec(yt, yp)
      mae_vals[c] <- mae_vec(yt, yp)
    }

    return(list(
      corr_nodewise_avgthencomp_bycond = corr_vals,
      R2_nodewise_avgthencomp_bycond = r2_vals,
      maeAcc_nodewise_avgthencomp_bycond = mae_vals
    ))
  }

  actflower_abort(sprintf("Unknown comparison type: %s", comparison))
}
