#' Noise ceiling estimate
#' @param act_run1 Array with dimensions (nodes x conditions x subjects).
#' @param act_run2 Array with dimensions (nodes x conditions x subjects).
#' @param ... Passed to model_compare.
#' @export
noise_ceiling <- function(act_run1, act_run2, ...) {
  model_compare(target = act_run2, model1 = act_run1, ...)
}
