# actflower news

## actflower 0.2.0
- Added stochastic parity contract runner/checker (`tools/stochastic_parity_runner.R`, `tools/check_stochastic_parity_report.py`) for seed-locked nested-CV and uncertainty artifacts.
- Added failure/special-case parity contract runner/checker (`tools/failure_parity_runner.py`, `tools/check_failure_parity_report.py`) and divergence mapping (`tools/parity_divergence_map.json`).
- Added parity drift checker and thresholds (`tools/check_parity_drift.py`, `tools/parity_drift_thresholds.json`) to detect regressions against baseline envelopes.
- Added parity contracts vignette (`vignettes/parity-contracts.Rmd`) and updated README parity reproduction commands.
- Added parity-nightly CI workflow (`.github/workflows/parity-nightly.yaml`) and strengthened PR parity gate with fuzz, stochastic, and failure contracts.
- Added OpenMP-enabled correlation batch kernel build flags and subject-parallel loop for `corr_fc_batch_cpp`.
- Added coverage threshold checker (`tools/check_coverage.R`) plus CI enforcement (`tools/coverage_thresholds.json`, `.github/workflows/tests-and-coverage.yaml`).
- Added cross-language benchmark harness (`tools/benchmark_r_vs_python.py`) with report checker thresholds (`tools/check_r_vs_python_report.py`).
- Added CI workflow for cross-language regression gating (`.github/workflows/crosslang-benchmark.yaml`).
- Tightened package benchmark gate to enforce thresholds in CI (`.github/workflows/benchmark-gate.yaml`).
- Expanded cross-language parity coverage to include combinedFC and noncircular synthetic targets.
- Added sparse/mask-aware noncircular parcel paths with `list`/`dgCMatrix` mask formats and parity tests.
- Added fused native fullcompare path for `actflow_test()` (linear, 3D FC) and benchmark coverage for end-to-end fullcompare parity/performance.
- Added `actflow_uncertainty()` for bootstrap-based interval estimates and uncertainty summaries over fullcompare actflow metrics.
- Added `actflow_nested_cv()` with strict outer/inner subject splits, fold-level inner metrics, selected-hyperparameter artifacts, optional JSON artifact persistence, and reproducibility metadata.
- Tightened uncertainty calibration tests to target PRD coverage bounds on synthetic MAE benchmarks.
- Added sparse noncircular efficiency gates (low-density synthetic profile) to benchmark enforcement, including runtime and allocation-ratio checks.
- Added a batched native correlation FC kernel for 3D inputs and wired `estimate_fc_corr()` to dispatch to the batch path, improving cross-language `fc_corr` speedup while preserving parity.
- Optimized fused fullcompare metric kernels to avoid allocation-heavy vectorization and tightened the benchmark gate to require `actflow_test` speedup >= 1.5x (PRD target).
- Added native noncircular activity kernels for sparse and dense mask regimes and reduced exclusion/index overhead in repeated calls.
- Hardened cross-language benchmark timing with adaptive per-repetition timing in both R and Python harnesses to reduce noise-sensitive threshold flakiness.
- Expanded unit and parity test coverage for transfer/predict branches, combinedFC internals, noncircular native kernels, and HCP fixture metrics (`corr`, `R2`, `mae`).
- Added Phase-P1 parity fuzz tooling: `tools/parity_fuzz.py`, aggregate checker `tools/check_parity_fuzz_report.py`, case-level and aggregate threshold profiles, and a smoke baseline report artifact.
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
