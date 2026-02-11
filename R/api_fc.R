#' Estimate connectivity matrix
#' @param x Numeric matrix.
#' @param method FC method.
#' @param partial_method Partial estimator when `method = "partial"`.
#' @param ... Method-specific args.
#' @export
estimate_fc <- function(
  x,
  method = c("corr", "multreg", "partial", "combined", "lasso_cv", "glasso_cv"),
  partial_method = c("ridge", "glasso", "pc"),
  ...
) {
  method <- match.arg(method)
  partial_method <- match.arg(partial_method)
  switch(
    method,
    corr = estimate_fc_corr(x, ...),
    multreg = estimate_fc_multreg(x, ...),
    partial = estimate_fc_partial(x, method = partial_method, ...),
    combined = estimate_fc_combined(x, ...),
    lasso_cv = estimate_fc_lasso_cv(x, ...),
    glasso_cv = estimate_fc_glasso_cv(x, ...)
  )
}
