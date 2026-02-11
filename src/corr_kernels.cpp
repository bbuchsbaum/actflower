#include <RcppArmadillo.h>
#include <cmath>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat corr_fc_cpp(const arma::mat& x_nodes_by_time) {
  const arma::uword n_nodes = x_nodes_by_time.n_rows;
  const arma::uword n_time = x_nodes_by_time.n_cols;

  if (n_nodes < 1) {
    return arma::mat();
  }
  if (n_time < 2) {
    Rcpp::stop("corr_fc_cpp requires at least 2 timepoints.");
  }

  arma::mat X = x_nodes_by_time;
  X.each_col() -= arma::mean(X, 1);

  arma::mat cov = (X * X.t()) / static_cast<double>(n_time - 1);
  arma::vec var = cov.diag();
  arma::vec sd = arma::sqrt(arma::clamp(var, 0.0, std::numeric_limits<double>::infinity()));

  arma::mat corr(n_nodes, n_nodes, arma::fill::zeros);
  const double eps = std::numeric_limits<double>::epsilon();

  for (arma::uword i = 0; i < n_nodes; ++i) {
    for (arma::uword j = i; j < n_nodes; ++j) {
      const double denom = sd[i] * sd[j];
      double val = std::numeric_limits<double>::quiet_NaN();
      if (denom > eps) {
        val = cov(i, j) / denom;
      }
      corr(i, j) = val;
      corr(j, i) = val;
    }
  }

  corr.diag().zeros();
  return corr;
}
