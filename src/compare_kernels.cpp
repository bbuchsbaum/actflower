#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

inline void fullcomp_metrics_slice(
  const arma::mat& target,
  const arma::mat& pred,
  double& corr,
  double& r2,
  double& mae
) {
  const arma::uword n_rows = target.n_rows;
  const arma::uword n_cols = target.n_cols;

  double n = 0.0;
  double sum_y = 0.0;
  double sum_p = 0.0;
  double sum_yy = 0.0;
  double sum_pp = 0.0;
  double sum_yp = 0.0;
  double ss_res = 0.0;
  double sum_abs_err = 0.0;

  for (arma::uword j = 0; j < n_cols; ++j) {
    for (arma::uword i = 0; i < n_rows; ++i) {
      const double y = target(i, j);
      const double p = pred(i, j);
      if (std::isfinite(y) && std::isfinite(p)) {
        const double e = y - p;
        n += 1.0;
        sum_y += y;
        sum_p += p;
        sum_yy += y * y;
        sum_pp += p * p;
        sum_yp += y * p;
        ss_res += e * e;
        sum_abs_err += std::abs(e);
      }
    }
  }

  if (n < 2.0) {
    corr = NA_REAL;
    r2 = NA_REAL;
    mae = NA_REAL;
    return;
  }

  double ss_tot = sum_yy - (sum_y * sum_y) / n;
  double ss_p = sum_pp - (sum_p * sum_p) / n;
  double cov = sum_yp - (sum_y * sum_p) / n;

  if (ss_tot < 0.0 && ss_tot > -1e-12) ss_tot = 0.0;
  if (ss_p < 0.0 && ss_p > -1e-12) ss_p = 0.0;

  const double den = std::sqrt(ss_tot * ss_p);
  corr = (!std::isfinite(den) || den <= 0.0) ? NA_REAL : (cov / den);
  r2 = (ss_tot > 0.0) ? (1.0 - ss_res / ss_tot) : NA_REAL;
  mae = sum_abs_err / n;
}

// [[Rcpp::export]]
Rcpp::List compare_fullcomp_cpp(const arma::cube& target, const arma::cube& pred) {
  if (target.n_rows != pred.n_rows || target.n_cols != pred.n_cols || target.n_slices != pred.n_slices) {
    Rcpp::stop("target and pred dimensions must match.");
  }

  const arma::uword n_subj = target.n_slices;
  Rcpp::NumericVector corr_vals(n_subj);
  Rcpp::NumericVector r2_vals(n_subj);
  Rcpp::NumericVector mae_vals(n_subj);

  for (arma::uword s = 0; s < n_subj; ++s) {
    double corr = NA_REAL;
    double r2 = NA_REAL;
    double mae = NA_REAL;
    fullcomp_metrics_slice(target.slice(s), pred.slice(s), corr, r2, mae);
    corr_vals[s] = corr;
    r2_vals[s] = r2;
    mae_vals[s] = mae;
  }

  return Rcpp::List::create(
    Rcpp::Named("corr_vals") = corr_vals,
    Rcpp::Named("R2_vals") = r2_vals,
    Rcpp::Named("mae_vals") = mae_vals
  );
}
