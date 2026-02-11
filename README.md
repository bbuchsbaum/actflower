# actflower

R implementation of activity-flow modeling focused on correctness, speed, and modularity.

Status: active development (v0.1.1).

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
  --r-libs /private/tmp/actflower-lib \
  --report-json /tmp/actflower_report.json
```

Validate a benchmark report against project thresholds:

```bash
python tools/check_r_vs_python_report.py \
  --report /tmp/actflower_report.json \
  --thresholds tools/benchmark_r_vs_python_thresholds.json
```
