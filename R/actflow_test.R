#' Test actflow predictions across conditions/subjects
#' @param act_group Array with dimensions (nodes x conditions x subjects).
#' @param fc_group Array with dimensions (nodes x nodes x subjects) or
#'   (nodes x nodes x conditions x subjects).
#' @param act_group_test Optional test activations.
#' @param comparison Comparison type.
#' @param transfer Transfer function for actflow prediction.
#' @param use_cpp Use native kernels when available.
#' @export
actflow_test <- function(
  act_group,
  fc_group,
  act_group_test = NULL,
  comparison = "fullcompare_compthenavg",
  transfer = c("linear", "relu", "sigmoid", "logit"),
  use_cpp = TRUE
) {
  validate_array3(act_group, arg = "act_group")
  transfer <- match.arg(transfer)
  d <- dim(act_group)
  n_nodes <- d[1]
  n_cond <- d[2]
  n_subj <- d[3]

  fc_dim <- dim(fc_group)
  if (length(fc_dim) < 3) {
    actflower_abort("`fc_group` must have at least 3 dimensions.")
  }

  if (length(fc_dim) == 3 && transfer == "linear") {
    pred <- .af_actflow_predict_batch(act_group, fc_group, use_cpp = use_cpp)
  } else {
    pred <- array(NA_real_, dim = c(n_nodes, n_cond, n_subj))
    for (s in seq_len(n_subj)) {
      fc_s <- if (length(fc_dim) == 4) NULL else fc_group[, , s, drop = TRUE]
      for (c in seq_len(n_cond)) {
        fc_sc <- if (length(fc_dim) == 4) fc_group[, , c, s, drop = TRUE] else fc_s
        pred[, c, s] <- actflow_predict(act_group[, c, s], fc_sc, transfer = transfer)
      }
    }
  }

  target <- if (is.null(act_group_test)) act_group else act_group_test
  cmp <- model_compare(target = target, model1 = pred, comparison = comparison, use_cpp = use_cpp)

  list(
    actPredVector_bytask_bysubj = pred,
    actVect_actual_group = target,
    model_compare_output = cmp
  )
}
