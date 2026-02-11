#include <RcppArmadillo.h>
#include <cmath>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]

inline double corr_vec(const arma::vec& x, const arma::vec& y) {
  std::vector<arma::uword> idx;
  idx.reserve(x.n_elem);
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    if (std::isfinite(x[i]) && std::isfinite(y[i])) {
      idx.push_back(i);
    }
  }
  if (idx.size() < 2) return NA_REAL;

  arma::uvec ok(idx.size());
  for (arma::uword i = 0; i < idx.size(); ++i) ok[i] = idx[i];

  arma::vec xx = x.elem(ok);
  arma::vec yy = y.elem(ok);
  xx -= arma::mean(xx);
  yy -= arma::mean(yy);

  double den = std::sqrt(arma::dot(xx, xx) * arma::dot(yy, yy));
  if (!std::isfinite(den) || den <= 0) return NA_REAL;
  return arma::dot(xx, yy) / den;
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
    arma::vec yt = arma::vectorise(target.slice(s));
    arma::vec yp = arma::vectorise(pred.slice(s));

    corr_vals[s] = corr_vec(yt, yp);

    std::vector<arma::uword> idx;
    idx.reserve(yt.n_elem);
    for (arma::uword i = 0; i < yt.n_elem; ++i) {
      if (std::isfinite(yt[i]) && std::isfinite(yp[i])) {
        idx.push_back(i);
      }
    }
    if (idx.size() < 2) {
      r2_vals[s] = NA_REAL;
      mae_vals[s] = NA_REAL;
      continue;
    }

    arma::uvec ok(idx.size());
    for (arma::uword i = 0; i < idx.size(); ++i) ok[i] = idx[i];

    arma::vec y = yt.elem(ok);
    arma::vec p = yp.elem(ok);
    arma::vec e = y - p;

    double ss_res = arma::dot(e, e);
    arma::vec yc = y - arma::mean(y);
    double ss_tot = arma::dot(yc, yc);

    r2_vals[s] = (ss_tot > 0) ? (1.0 - ss_res / ss_tot) : NA_REAL;
    mae_vals[s] = arma::mean(arma::abs(e));
  }

  return Rcpp::List::create(
    Rcpp::Named("corr_vals") = corr_vals,
    Rcpp::Named("R2_vals") = r2_vals,
    Rcpp::Named("mae_vals") = mae_vals
  );
}
