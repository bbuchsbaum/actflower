# actflower news

## actflower 0.1.1
- Added cross-language benchmark harness (`tools/benchmark_r_vs_python.py`) with report checker thresholds (`tools/check_r_vs_python_report.py`).
- Added CI workflow for cross-language regression gating (`.github/workflows/crosslang-benchmark.yaml`).
- Tightened package benchmark gate to enforce thresholds in CI (`.github/workflows/benchmark-gate.yaml`).
- Expanded cross-language parity coverage to include combinedFC timings and similarity checks.
- Improved FC correlation performance path with native dispatch caching and optimized correlation kernel handling.
- Improved API ergonomics by accepting orientation aliases (`nodes_x_time`, `time_x_nodes`).

## actflower 0.1.0
- First functional release with parity-oriented activity-flow workflows.
- Added FC estimators: correlation, multiregression, partial, combinedFC, lasso CV, and graphical-lasso CV.
- Added model comparison and noise-ceiling helpers.
- Added non-circular parcel connectivity/activity utilities.
- Added native Rcpp kernels for multreg FC, batched actflow prediction, and full-compare metrics.
- Added parity fixtures, benchmark harness, and CI gates.
- Added basic and advanced vignettes.
