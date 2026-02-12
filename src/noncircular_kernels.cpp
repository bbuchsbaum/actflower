#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::cube noncircular_activity_sparse_cpp(
  const arma::mat& data_nodes_by_conditions,
  const Rcpp::List& sources_by_target,
  double fill_value = 0.0
) {
  const arma::uword n_nodes = data_nodes_by_conditions.n_rows;
  const arma::uword n_cond = data_nodes_by_conditions.n_cols;

  if (sources_by_target.size() != static_cast<R_xlen_t>(n_nodes)) {
    Rcpp::stop("`sources_by_target` must have one entry per target node.");
  }

  arma::cube out(n_nodes, n_nodes, n_cond);
  out.fill(fill_value);

  for (arma::uword target = 0; target < n_nodes; ++target) {
    Rcpp::IntegerVector src = Rcpp::as<Rcpp::IntegerVector>(sources_by_target[target]);
    for (R_xlen_t k = 0; k < src.size(); ++k) {
      const int src_idx = src[k] - 1;
      if (src_idx < 0 || src_idx >= static_cast<int>(n_nodes)) continue;
      for (arma::uword c = 0; c < n_cond; ++c) {
        out(target, static_cast<arma::uword>(src_idx), c) = data_nodes_by_conditions(static_cast<arma::uword>(src_idx), c);
      }
    }
  }

  // Non-circular workflows preserve each target's own activity on the diagonal.
  for (arma::uword target = 0; target < n_nodes; ++target) {
    for (arma::uword c = 0; c < n_cond; ++c) {
      out(target, target, c) = data_nodes_by_conditions(target, c);
    }
  }

  return out;
}

// [[Rcpp::export]]
arma::cube noncircular_activity_dense_cpp(
  const arma::mat& data_nodes_by_conditions,
  const Rcpp::List& exclusions_by_target,
  double fill_value = 0.0
) {
  const arma::uword n_nodes = data_nodes_by_conditions.n_rows;
  const arma::uword n_cond = data_nodes_by_conditions.n_cols;

  if (exclusions_by_target.size() != static_cast<R_xlen_t>(n_nodes)) {
    Rcpp::stop("`exclusions_by_target` must have one entry per target node.");
  }

  arma::cube out(n_nodes, n_nodes, n_cond);

  // Materialize full target-by-source tensor from source activity rows.
  for (arma::uword target = 0; target < n_nodes; ++target) {
    for (arma::uword source = 0; source < n_nodes; ++source) {
      for (arma::uword c = 0; c < n_cond; ++c) {
        out(target, source, c) = data_nodes_by_conditions(source, c);
      }
    }
  }

  // Apply exclusions per target.
  for (arma::uword target = 0; target < n_nodes; ++target) {
    Rcpp::IntegerVector ex = Rcpp::as<Rcpp::IntegerVector>(exclusions_by_target[target]);
    for (R_xlen_t k = 0; k < ex.size(); ++k) {
      const int source = ex[k] - 1;
      if (source < 0 || source >= static_cast<int>(n_nodes)) continue;
      for (arma::uword c = 0; c < n_cond; ++c) {
        out(target, static_cast<arma::uword>(source), c) = fill_value;
      }
    }
  }

  // Preserve target self-activity on diagonal.
  for (arma::uword target = 0; target < n_nodes; ++target) {
    for (arma::uword c = 0; c < n_cond; ++c) {
      out(target, target, c) = data_nodes_by_conditions(target, c);
    }
  }

  return out;
}
