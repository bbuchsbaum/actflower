#include <RcppArmadillo.h>
#include <cmath>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat multreg_fc_cpp(const arma::mat& x_nodes_by_time, double ridge = 0.0) {
  const arma::uword n_nodes = x_nodes_by_time.n_rows;
  const arma::uword n_time = x_nodes_by_time.n_cols;

  if (n_nodes < 2) {
    Rcpp::stop("Need at least 2 nodes.");
  }
  if (n_nodes >= n_time) {
    Rcpp::stop("multreg_fc_cpp requires nodes < timepoints.");
  }
  if (ridge < 0) {
    Rcpp::stop("ridge must be non-negative.");
  }

  arma::mat X = x_nodes_by_time;
  X.each_col() -= arma::mean(X, 1);

  arma::mat Sigma = (X * X.t()) / static_cast<double>(n_time - 1);
  if (ridge > 0) {
    Sigma.diag() += ridge;
  }
  Sigma = 0.5 * (Sigma + Sigma.t());

  arma::vec eigval;
  arma::mat eigvec;
  bool ok = arma::eig_sym(eigval, eigvec, Sigma);
  if (!ok) {
    Rcpp::stop("Failed eigendecomposition for covariance matrix.");
  }

  const double tol = std::max(1e-12, eigval.max() * static_cast<double>(n_nodes) * std::numeric_limits<double>::epsilon());
  arma::vec inv_vals = eigval;
  for (arma::uword i = 0; i < inv_vals.n_elem; ++i) {
    inv_vals[i] = (inv_vals[i] > tol) ? 1.0 / inv_vals[i] : 0.0;
  }

  arma::mat Theta = eigvec * arma::diagmat(inv_vals) * eigvec.t();
  Theta = 0.5 * (Theta + Theta.t());

  arma::vec d = Theta.diag();
  for (arma::uword i = 0; i < d.n_elem; ++i) {
    if (std::abs(d[i]) < 1e-12) d[i] = (d[i] >= 0 ? 1e-12 : -1e-12);
  }

  arma::mat B = -arma::diagmat(1.0 / d) * Theta;
  B.diag().zeros();

  return B;
}
