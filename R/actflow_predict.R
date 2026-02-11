#' Activity flow prediction
#' @param act Activation vector/matrix.
#' @param fc Connectivity matrix (nodes x nodes).
#' @param separate_by_target If TRUE, act is interpreted as (target x source).
#' @param transfer Transfer function name.
#' @export
actflow_predict <- function(act, fc, separate_by_target = FALSE, transfer = c("linear", "relu", "sigmoid", "logit")) {
  validate_square_matrix(fc, arg = "fc")
  transfer <- match.arg(transfer)
  n <- nrow(fc)

  if (!separate_by_target) {
    if (!is.numeric(act) || length(act) != n) {
      actflower_abort(sprintf("When separate_by_target=FALSE, `act` must have length %d.", n))
    }
    a <- transfer_function(as.numeric(act), transfer = transfer)
    pred <- as.vector(fc %*% a)
    # Match ActflowToolbox behavior: always exclude self-connections even if diagonal != 0.
    pred <- pred - diag(fc) * a
    return(pred)
  }

  validate_matrix(act, arg = "act")
  if (!all(dim(act) == c(n, n))) {
    actflower_abort(sprintf("When separate_by_target=TRUE, `act` must be [%d x %d].", n, n))
  }

  a <- transfer_function(act, transfer = transfer)
  pred <- rowSums(a * fc)
  pred <- pred - diag(fc) * diag(a)
  pred
}
