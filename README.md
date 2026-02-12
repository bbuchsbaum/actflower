# actflower

R implementation of activity-flow modeling focused on correctness, speed, and modularity.

Status: active development (v0.1.1).

## Uncertainty API

Bootstrap uncertainty summaries for fullcompare metrics:

```r
out <- actflow_uncertainty(
  act_group = taskbeta,
  fc_group = fc_multreg,
  metric = c("corr", "R2", "mae"),
  n_boot = 1000,
  conf_level = 0.95,
  seed = 2026
)

out$point
out$interval
```

## Cross-language Benchmark

Run R vs Python synthetic benchmarks (including similarity checks):

```bash
python tools/benchmark_r_vs_python.py \
  --nodes 60 \
  --timepoints 200 \
  --conditions 10 \
  --subjects 4 \
  --repeats 2 \
  --include-combined \
  --include-noncircular \
  --r-libs /private/tmp/actflower-lib \
  --report-json /tmp/actflower_report.json
```

Validate a benchmark report against project thresholds:

```bash
python tools/check_r_vs_python_report.py \
  --report /tmp/actflower_report.json \
  --thresholds tools/benchmark_r_vs_python_thresholds.json
```
