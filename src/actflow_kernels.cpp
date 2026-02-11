#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::cube actflow_predict_batch_cpp(const arma::cube& act, const arma::cube& fc, bool remove_diag = true) {
  if (act.n_rows != fc.n_rows || fc.n_rows != fc.n_cols) {
    Rcpp::stop("Dimension mismatch: act must be [nodes x cond x subj] and fc [nodes x nodes x subj].");
  }
  if (act.n_slices != fc.n_slices) {
    Rcpp::stop("Subject dimension mismatch between act and fc.");
  }

  const arma::uword n_nodes = act.n_rows;
  const arma::uword n_cond = act.n_cols;
  const arma::uword n_subj = act.n_slices;

  arma::cube out(n_nodes, n_cond, n_subj, arma::fill::zeros);

  for (arma::uword s = 0; s < n_subj; ++s) {
    const arma::mat& A = act.slice(s);
    const arma::mat& F = fc.slice(s);

    arma::mat P = F * A;
    if (remove_diag) {
      arma::vec d = F.diag();
      P -= arma::diagmat(d) * A;
    }
    out.slice(s) = P;
  }

  return out;
}
