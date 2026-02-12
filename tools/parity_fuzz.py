#!/usr/bin/env python3
"""Differential parity fuzzing harness for actflower (R) vs ActflowToolbox (Python).

Runs many randomized benchmark cases by delegating each case to
tools/benchmark_r_vs_python.py, then aggregates parity and timing statistics.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import random
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List


def quantile(values: List[float], q: float) -> float:
    if not values:
        return float("nan")
    if q <= 0:
        return min(values)
    if q >= 1:
        return max(values)
    s = sorted(values)
    pos = (len(s) - 1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return s[lo]
    frac = pos - lo
    return s[lo] * (1 - frac) + s[hi] * frac


def op_stats(values: List[float]) -> Dict[str, float]:
    finite = [float(v) for v in values if isinstance(v, (int, float)) and math.isfinite(v)]
    if not finite:
        return {"count": 0, "min": float("nan"), "p50": float("nan"), "p95": float("nan"), "p99": float("nan"), "max": float("nan"), "mean": float("nan")}
    return {
        "count": len(finite),
        "min": min(finite),
        "p50": quantile(finite, 0.50),
        "p95": quantile(finite, 0.95),
        "p99": quantile(finite, 0.99),
        "max": max(finite),
        "mean": statistics.fmean(finite),
    }


def load_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text())


def choose_profile(case_idx: int) -> Dict[str, bool]:
    # Cycle profiles to guarantee operation-family coverage.
    profiles = [
        {"include_combined": False, "include_noncircular": False},
        {"include_combined": True, "include_noncircular": False},
        {"include_combined": False, "include_noncircular": True},
        {"include_combined": True, "include_noncircular": True},
    ]
    return profiles[case_idx % len(profiles)]


def make_case_config(rng: random.Random, case_idx: int) -> Dict[str, Any]:
    profile = choose_profile(case_idx)
    cfg = {
        "nodes": rng.randint(24, 84),
        "timepoints": rng.randint(120, 360),
        "conditions": rng.randint(6, 20),
        "subjects": rng.randint(3, 8),
        "repeats": 1,
    }
    cfg.update(profile)
    return cfg


def evaluate_case_against_thresholds(report: Dict[str, Any], thresholds: Dict[str, Any]) -> Dict[str, Any]:
    speedup_min = thresholds.get("speedup_min", {})
    sim_max = thresholds.get("similarity_max_abs_err", {})
    sim_corr = thresholds.get("similarity_min_corr", {})
    speedups = report.get("r_speedup_over_python", {})
    similarity = report.get("similarity", {})

    ok = True
    reasons: List[str] = []

    for op, min_speed in speedup_min.items():
        if op not in speedups:
            continue
        got = speedups.get(op, float("nan"))
        if not isinstance(got, (int, float)) or not math.isfinite(got) or got < float(min_speed):
            ok = False
            reasons.append(f"speedup:{op} got={got} need>={float(min_speed)}")

    for op, max_abs in sim_max.items():
        if op not in similarity:
            continue
        got = similarity.get(op, {}).get("max_abs_err", float("nan"))
        if not isinstance(got, (int, float)) or not math.isfinite(got) or got > float(max_abs):
            ok = False
            reasons.append(f"max_abs_err:{op} got={got} need<={float(max_abs)}")

    for op, min_corr in sim_corr.items():
        if op not in similarity:
            continue
        got = similarity.get(op, {}).get("corr", float("nan"))
        if not isinstance(got, (int, float)) or not math.isfinite(got) or got < float(min_corr):
            ok = False
            reasons.append(f"corr:{op} got={got} need>={float(min_corr)}")

    return {"pass": ok, "reasons": reasons}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", type=int, default=24, help="Number of fuzz cases.")
    parser.add_argument("--seed", type=int, default=2026, help="Master seed.")
    parser.add_argument("--python-exe", default=sys.executable)
    parser.add_argument("--r-libs", default=None)
    parser.add_argument("--benchmark-script", default="tools/benchmark_r_vs_python.py")
    parser.add_argument("--work-dir", default="inst/extdata/benchmarks/parity_fuzz")
    parser.add_argument("--thresholds", default="tools/parity_fuzz_case_thresholds.json")
    parser.add_argument("--report-json", default="inst/extdata/benchmarks/parity_fuzz/report.json")
    parser.add_argument("--keep-case-reports", action="store_true", default=False)
    parser.add_argument("--retain-case-artifacts", action="store_true", default=False)
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    benchmark_script = (repo_root / args.benchmark_script).resolve()
    work_dir = (repo_root / args.work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)
    thresholds_path = (repo_root / args.thresholds).resolve()
    thresholds = load_json(thresholds_path) if thresholds_path.exists() else {}

    rng = random.Random(args.seed)

    cases: List[Dict[str, Any]] = []
    profile_counts: Dict[str, int] = {"core": 0, "combined": 0, "noncircular": 0, "full": 0}

    for i in range(args.cases):
        cfg = make_case_config(rng, i)
        case_seed = args.seed * 100000 + i
        case_name = f"case_{i:04d}"
        case_work_dir = work_dir / case_name
        case_work_dir.mkdir(parents=True, exist_ok=True)
        case_report = case_work_dir / "report.json"

        if cfg["include_combined"] and cfg["include_noncircular"]:
            profile_counts["full"] += 1
        elif cfg["include_combined"]:
            profile_counts["combined"] += 1
        elif cfg["include_noncircular"]:
            profile_counts["noncircular"] += 1
        else:
            profile_counts["core"] += 1

        cmd = [
            args.python_exe,
            str(benchmark_script),
            "--nodes",
            str(cfg["nodes"]),
            "--timepoints",
            str(cfg["timepoints"]),
            "--conditions",
            str(cfg["conditions"]),
            "--subjects",
            str(cfg["subjects"]),
            "--repeats",
            str(cfg["repeats"]),
            "--seed",
            str(case_seed),
            "--work-dir",
            str(case_work_dir),
            "--report-json",
            str(case_report),
        ]
        if cfg["include_combined"]:
            cmd.append("--include-combined")
        if cfg["include_noncircular"]:
            cmd.append("--include-noncircular")
        if args.r_libs:
            cmd.extend(["--r-libs", args.r_libs])

        started = time.time()
        proc = subprocess.run(cmd, cwd=str(repo_root), capture_output=True, text=True)
        elapsed = time.time() - started

        record: Dict[str, Any] = {
            "case_id": i,
            "case_name": case_name,
            "seed": case_seed,
            "config": cfg,
            "runtime_sec": elapsed,
            "command": cmd,
            "report_path": str(case_report),
        }

        if proc.returncode != 0 or not case_report.exists():
            record["status"] = "exec_failed"
            record["return_code"] = proc.returncode
            record["stderr_tail"] = "\n".join(proc.stderr.splitlines()[-40:])
            record["stdout_tail"] = "\n".join(proc.stdout.splitlines()[-20:])
            if not args.retain_case_artifacts:
                shutil.rmtree(case_work_dir, ignore_errors=True)
                record["report_path"] = None
            cases.append(record)
            continue

        case_data = load_json(case_report)
        record["status"] = "ok"
        record["operations"] = sorted(case_data.get("similarity", {}).keys())
        record["speedups"] = case_data.get("r_speedup_over_python", {})
        record["similarity"] = case_data.get("similarity", {})
        eval_out = evaluate_case_against_thresholds(case_data, thresholds) if thresholds else {"pass": True, "reasons": []}
        record["threshold_pass"] = bool(eval_out["pass"])
        record["threshold_fail_reasons"] = eval_out["reasons"]
        if not args.retain_case_artifacts:
            shutil.rmtree(case_work_dir, ignore_errors=True)
            record["report_path"] = None
        cases.append(record)

    successful = [c for c in cases if c["status"] == "ok"]
    failed_exec = [c for c in cases if c["status"] != "ok"]
    parity_pass = [c for c in successful if c.get("threshold_pass", False)]
    parity_fail = [c for c in successful if not c.get("threshold_pass", False)]

    op_sim: Dict[str, Dict[str, List[float]]] = {}
    op_speed: Dict[str, List[float]] = {}
    for c in successful:
        for op, m in c.get("similarity", {}).items():
            op_sim.setdefault(op, {"max_abs_err": [], "mean_abs_err": [], "corr": []})
            for k in ("max_abs_err", "mean_abs_err", "corr"):
                v = m.get(k, float("nan"))
                if isinstance(v, (int, float)) and math.isfinite(v):
                    op_sim[op][k].append(float(v))
        for op, v in c.get("speedups", {}).items():
            if isinstance(v, (int, float)) and math.isfinite(v):
                op_speed.setdefault(op, []).append(float(v))

    by_operation: Dict[str, Any] = {}
    for op in sorted(set(list(op_sim.keys()) + list(op_speed.keys()))):
        by_operation[op] = {
            "similarity": {
                "max_abs_err": op_stats(op_sim.get(op, {}).get("max_abs_err", [])),
                "mean_abs_err": op_stats(op_sim.get(op, {}).get("mean_abs_err", [])),
                "corr": op_stats(op_sim.get(op, {}).get("corr", [])),
            },
            "speedup": op_stats(op_speed.get(op, [])),
        }

    report = {
        "schema_version": "1.0.0",
        "tool": "parity_fuzz",
        "created_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "repo_root": str(repo_root),
        "benchmark_script": str(benchmark_script),
        "thresholds_path": str(thresholds_path) if thresholds_path.exists() else None,
        "master_seed": args.seed,
        "cases_requested": args.cases,
        "summary": {
            "cases_total": len(cases),
            "cases_successful": len(successful),
            "cases_exec_failed": len(failed_exec),
            "cases_threshold_pass": len(parity_pass),
            "cases_threshold_fail": len(parity_fail),
            "successful_case_ratio": (len(successful) / len(cases)) if cases else float("nan"),
            "threshold_pass_ratio": (len(parity_pass) / len(successful)) if successful else float("nan"),
            "profile_counts": profile_counts,
        },
        "by_operation": by_operation,
        "cases": cases if args.keep_case_reports else [],
    }

    out_report = (repo_root / args.report_json).resolve()
    out_report.parent.mkdir(parents=True, exist_ok=True)
    out_report.write_text(json.dumps(report, indent=2))

    print(f"Wrote parity fuzz report: {out_report}")
    print("Summary:")
    print(json.dumps(report["summary"], indent=2))


if __name__ == "__main__":
    main()
