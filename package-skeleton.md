# actflower Package Skeleton

## Goal
Bootstrap `actflower` as a production-grade R package with clear module boundaries, parity harnesses, performance gates, and extension points for optimized kernels.

## 1) Repository Layout

```text
actflower/
  DESCRIPTION
  NAMESPACE
  LICENSE
  README.md
  NEWS.md
  .Rbuildignore
  .gitignore
  prd.json
  vision.md
  package-skeleton.md

  R/
    zzz.R
    api_fc.R
    api_actflow.R
    api_compare.R
    api_io.R

    fc_corr.R
    fc_multreg.R
    fc_partial.R
    fc_pcreg.R
    fc_combined.R
    fc_lasso.R
    fc_glasso.R

    actflow_predict.R
    actflow_test.R
    actflow_transfer.R

    compare_metrics.R
    compare_model.R
    noise_ceiling.R

    noncircular_parcel.R
    noncircular_cifti.R

    io_hdf5.R
    io_neuroim2.R
    io_shapes.R
    validators.R
    errors.R

  src/
    Makevars
    Makevars.win
    RcppExports.cpp
    multreg_kernels.cpp
    actflow_kernels.cpp
    compare_kernels.cpp

  man/

  tests/
    testthat.R
    testthat/
      test-validators.R
      test-io-hdf5.R
      test-fc-corr.R
      test-fc-multreg.R
      test-fc-partial.R
      test-fc-pcreg.R
      test-fc-combined.R
      test-fc-lasso.R
      test-fc-glasso.R
      test-actflow-predict.R
      test-actflow-test.R
      test-compare-model.R
      test-noise-ceiling.R
      test-parity-hcp.R
      test-reproducibility.R
      test-errors.R

  inst/
    extdata/
      hcp_example/
      parity_refs/
    scripts/
      build_parity_refs.py
      bench_runner.R

  vignettes/
    actflower-basic-workflow.Rmd
    actflower-advanced-fc-comparison.Rmd

  tools/
    benchmark_profiles.R
    generate-fixtures.R

  .github/
    workflows/
      r-cmd-check.yaml
      tests-and-coverage.yaml
      benchmark-gate.yaml
      parity-gate.yaml
```

## 2) Initial DESCRIPTION Blueprint

- `Package`: `actflower`
- `Type`: `Package`
- `Title`: Activity Flow Modeling in R
- `Version`: `0.0.0.9000`
- `Depends`: `R (>= 4.3.0)`
- `Imports`:
  - `hdf5r`
  - `Matrix`
  - `matrixStats`
  - `stats`
  - `methods`
  - `Rcpp`
- `LinkingTo`:
  - `Rcpp`
  - `RcppArmadillo`
- `Suggests`:
  - `testthat`
  - `bench`
  - `covr`
  - `waldo`
  - `knitr`
  - `rmarkdown`
  - `future.apply`
  - `glmnet`
  - `glasso`
  - `corpcor`
  - `ciftiTools`
  - `neuroim2`
    - GitHub: `https://github.com/bbuchsbaum/neuroim2`
    - Local dev path: `/Users/bbuchsbaum/code/neuroim2`

## 3) Public API (v1)

- `read_actflower_h5(path, dataset)`
- `write_actflower_h5(path, name, x)`
- `estimate_fc(x, method, ...)`
- `estimate_fc_multreg(x, ...)`
- `estimate_fc_partial(x, method = c("ridge", "glasso", "pc"), ...)`
- `estimate_fc_combined(x, ...)`
- `actflow_predict(act, fc, separate_by_target = FALSE, transfer = "linear", ...)`
- `actflow_test(act, fc, act_test = NULL, comparison = "fullcompare_compthenavg", ...)`
- `model_compare(target, model1, model2 = NULL, comparison = "fullcompare_compthenavg", ...)`
- `noise_ceiling(act_run1, act_run2, ...)`
- `calcconn_parcelwise_noncircular(data, ...)`
- `calcactivity_parcelwise_noncircular(data, ...)`

Design rules:
- No wildcard exports.
- All exported functions validate shape and orientation.
- Errors include expected and actual dimensions.

## 4) Internal Contracts

- Canonical activity tensor: `[nodes x conditions x subjects]`.
- Canonical FC tensor: `[nodes x nodes x subjects]` or `[nodes x nodes x conditions x subjects]`.
- Time series input accepted as either orientation with explicit `orientation = c("nodes_by_time", "time_by_nodes")`.
- Deterministic random paths require `seed` argument.

## 5) Reuse Plan from `~/code/ariadne`

Adopt/adapt first:
- `weighted_cor` / `cor_event_weighted`
- `pcor_ridge`, `pcor_glasso`, `pcor_pc`
- CV helpers for lambda/k grid search
- C++ implementation patterns from `src/ewma_corr.cpp`

Rules:
- Preserve attribution in source headers.
- Rewrap under actflower API contracts.
- Add parity + regression tests before claiming done.

## 6) Performance Plan (First Kernels)

- `multreg_kernels.cpp`
  - Fast coefficient extraction from precision matrix.
  - Mask-aware targetwise regression fallback.

- `actflow_kernels.cpp`
  - Batched prediction with diagonal exclusion.
  - Target-specific source masking path.

- `compare_kernels.cpp`
  - Vectorized Pearson/R2/MAE across tensor slices.

Benchmark targets (from PRD):
- `>= 2x` faster multreg FC vs baseline Python implementation.
- `>= 1.5x` faster end-to-end actflow test on HCP example workload.

## 7) Test Strategy

- Unit tests per module.
- Shape and validation tests for every exported function.
- Golden parity tests using HCP example fixtures.
- Reproducibility tests for CV and stochastic branches.
- Snapshot tests for error messages.

Coverage target:
- `>= 90%` line coverage on exported API.

## 8) Parity Harness

- `inst/scripts/build_parity_refs.py`
  - Uses original ActflowToolbox to generate reference outputs.
  - Writes HDF5/NPY fixtures into `inst/extdata/parity_refs/`.

- `tests/testthat/test-parity-hcp.R`
  - Loads fixtures via `hdf5r`.
  - Compares actflower outputs to reference with method-specific tolerances.

## 9) CI Gates

- `r-cmd-check.yaml`
  - R CMD check on macOS and Linux.

- `tests-and-coverage.yaml`
  - Full tests + coverage threshold enforcement.

- `parity-gate.yaml`
  - Runs parity suite on fixture data.

- `benchmark-gate.yaml`
  - Compares benchmark outputs against threshold file.
  - Fails on regression beyond configured tolerance.

## 10) Implementation Order

1. Scaffold package + validators + IO (`hdf5r`).
2. Implement `estimate_fc()` for `corr` and `multreg`.
3. Implement `actflow_predict()` + `actflow_test()`.
4. Implement `model_compare()` + `noise_ceiling()`.
5. Add partial-correlation family (`ridge/glasso/pc`) with CV.
6. Add combinedFC + masked multreg path.
7. Add non-circular parcel workflows.
8. Add native kernels and benchmark gates.
9. Finalize docs and vignettes.

## 11) Definition of Ready for Coding Sprint 1

- `DESCRIPTION`, `NAMESPACE`, `R/`, `tests/`, `src/`, `vignettes/`, `inst/extdata/` created.
- One passing test per module namespace.
- Fixture ingestion via `hdf5r` proven with HCP example files.
- CI workflow stubs committed.

## 12) Definition of Done for v0.1.0

- AC-001 through AC-008 from `prd.json` all passing.
- Public API documented with runnable examples.
- Benchmarks show required speedups.
- No `rhdf5` usage anywhere in package.
- Clear migration notes from Python ActflowToolbox.
