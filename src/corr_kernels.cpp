#include <RcppArmadillo.h>
#include <cmath>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

inline arma::mat corr_fc_from_nodes_by_time(const arma::mat& x_nodes_by_time) {
  const arma::uword n_nodes = x_nodes_by_time.n_rows;
  const arma::uword n_time = x_nodes_by_time.n_cols;

  if (n_nodes < 1) {
    return arma::mat();
  }
  if (n_time < 2) {
    Rcpp::stop("corr_fc_cpp requires at least 2 timepoints.");
  }

  arma::mat X = x_nodes_by_time;
  const double denom = static_cast<double>(n_time - 1);
  X.each_col() -= arma::mean(X, 1);

  // Standardize each row in-place; then correlation is normalized cross-product.
  arma::vec ss = arma::sum(arma::square(X), 1); // sum of squared deviations per row
  const double eps = std::numeric_limits<double>::epsilon() * denom;
  arma::uvec good = arma::find(ss > eps);
  arma::uvec bad = arma::find(ss <= eps);

  arma::vec inv_sd(n_nodes, arma::fill::zeros);
  if (!good.empty()) {
    inv_sd.elem(good) = arma::sqrt(denom / ss.elem(good));
  }
  X.each_col() %= inv_sd;

  arma::mat corr = (X * X.t()) / denom;

  // Match correlation semantics for zero-variance rows.
  for (arma::uword k = 0; k < bad.n_elem; ++k) {
    arma::uword i = bad[k];
    corr.row(i).fill(arma::datum::nan);
    corr.col(i).fill(arma::datum::nan);
  }
  corr.diag().zeros();
  return corr;
}

// [[Rcpp::export]]
arma::mat corr_fc_cpp(const arma::mat& x_nodes_by_time) {
  return corr_fc_from_nodes_by_time(x_nodes_by_time);
}

// [[Rcpp::export]]
arma::cube corr_fc_batch_cpp(const arma::cube& x_nodes_by_time) {
  const arma::uword n_nodes = x_nodes_by_time.n_rows;
  const arma::uword n_time = x_nodes_by_time.n_cols;
  const arma::uword n_subj = x_nodes_by_time.n_slices;

  if (n_nodes < 1) {
    return arma::cube();
  }
  if (n_time < 2) {
    Rcpp::stop("corr_fc_batch_cpp requires at least 2 timepoints.");
  }

  arma::cube out(n_nodes, n_nodes, n_subj);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int s = 0; s < static_cast<int>(n_subj); ++s) {
    out.slice(static_cast<arma::uword>(s)) = corr_fc_from_nodes_by_time(
      x_nodes_by_time.slice(static_cast<arma::uword>(s))
    );
  }
  return out;
}
