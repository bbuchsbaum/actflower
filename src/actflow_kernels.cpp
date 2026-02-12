#include <RcppArmadillo.h>
#include <cmath>
#include <vector>

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

inline void fullcomp_metrics_slice(
  const arma::mat& target,
  const arma::mat& pred,
  double& corr,
  double& r2,
  double& mae
) {
  arma::vec yt = arma::vectorise(target);
  arma::vec yp = arma::vectorise(pred);

  std::vector<arma::uword> idx;
  idx.reserve(yt.n_elem);
  for (arma::uword i = 0; i < yt.n_elem; ++i) {
    if (std::isfinite(yt[i]) && std::isfinite(yp[i])) {
      idx.push_back(i);
    }
  }

  if (idx.size() < 2) {
    corr = NA_REAL;
    r2 = NA_REAL;
    mae = NA_REAL;
    return;
  }

  arma::uvec ok(idx.size());
  for (arma::uword i = 0; i < idx.size(); ++i) ok[i] = idx[i];

  arma::vec y = yt.elem(ok);
  arma::vec p = yp.elem(ok);
  arma::vec yc = y - arma::mean(y);
  arma::vec pc = p - arma::mean(p);

  double den = std::sqrt(arma::dot(yc, yc) * arma::dot(pc, pc));
  if (!std::isfinite(den) || den <= 0) {
    corr = NA_REAL;
  } else {
    corr = arma::dot(yc, pc) / den;
  }

  arma::vec e = y - p;
  double ss_res = arma::dot(e, e);
  double ss_tot = arma::dot(yc, yc);
  r2 = (ss_tot > 0) ? (1.0 - ss_res / ss_tot) : NA_REAL;
  mae = arma::mean(arma::abs(e));
}

// [[Rcpp::export]]
Rcpp::List actflow_fullcomp_batch_cpp(
  const arma::cube& act,
  const arma::cube& fc,
  const arma::cube& target,
  bool remove_diag = true
) {
  if (act.n_rows != fc.n_rows || fc.n_rows != fc.n_cols) {
    Rcpp::stop("Dimension mismatch: act must be [nodes x cond x subj] and fc [nodes x nodes x subj].");
  }
  if (act.n_rows != target.n_rows || act.n_cols != target.n_cols || act.n_slices != target.n_slices) {
    Rcpp::stop("Dimension mismatch: target must match act dimensions.");
  }
  if (act.n_slices != fc.n_slices) {
    Rcpp::stop("Subject dimension mismatch between act and fc.");
  }

  const arma::uword n_nodes = act.n_rows;
  const arma::uword n_cond = act.n_cols;
  const arma::uword n_subj = act.n_slices;

  arma::cube pred(n_nodes, n_cond, n_subj, arma::fill::zeros);
  Rcpp::NumericVector corr_vals(n_subj);
  Rcpp::NumericVector r2_vals(n_subj);
  Rcpp::NumericVector mae_vals(n_subj);

  for (arma::uword s = 0; s < n_subj; ++s) {
    const arma::mat& A = act.slice(s);
    const arma::mat& F = fc.slice(s);

    arma::mat P = F * A;
    if (remove_diag) {
      arma::vec d = F.diag();
      for (arma::uword c = 0; c < n_cond; ++c) {
        P.col(c) -= d % A.col(c);
      }
    }
    pred.slice(s) = P;

    double corr = NA_REAL;
    double r2 = NA_REAL;
    double mae = NA_REAL;
    fullcomp_metrics_slice(target.slice(s), P, corr, r2, mae);
    corr_vals[s] = corr;
    r2_vals[s] = r2;
    mae_vals[s] = mae;
  }

  return Rcpp::List::create(
    Rcpp::Named("pred") = pred,
    Rcpp::Named("corr_vals") = corr_vals,
    Rcpp::Named("R2_vals") = r2_vals,
    Rcpp::Named("mae_vals") = mae_vals
  );
}
