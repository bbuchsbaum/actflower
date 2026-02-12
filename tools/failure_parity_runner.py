#!/usr/bin/env python3
"""Run cross-language failure/special-case contract checks."""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Callable, Dict

import numpy as np


def load_module(module_name: str, path: Path):
    spec = importlib.util.spec_from_file_location(module_name, str(path))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module from {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run_python_safely(fn: Callable[[], Any]) -> Dict[str, Any]:
    try:
        fn()
        return {"ok": True, "error_class": None, "message": None}
    except Exception as e:  # noqa: BLE001
        return {
            "ok": False,
            "error_class": type(e).__name__,
            "message": str(e),
        }


def run_r_expr(expr: str, rscript_exe: str, r_libs: str | None) -> Dict[str, Any]:
    wrapped = f'''
    suppressPackageStartupMessages(library(actflower))
    suppressPackageStartupMessages(library(jsonlite))
    res <- tryCatch({{
      {expr}
      list(ok = TRUE, error_class = NA_character_, message = NA_character_)
    }}, error = function(e) {{
      list(ok = FALSE, error_class = class(e)[1], message = conditionMessage(e))
    }})
    cat(jsonlite::toJSON(res, auto_unbox = TRUE, null = "null"))
    '''
    env = dict(os.environ)
    if r_libs:
        env["R_LIBS"] = r_libs
    proc = subprocess.run(
        [rscript_exe, "-e", wrapped],
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )
    if proc.returncode != 0:
        return {
            "ok": False,
            "error_class": "RscriptProcessError",
            "message": (proc.stderr or proc.stdout).strip()[-5000:],
        }
    txt = proc.stdout.strip()
    try:
        out = json.loads(txt)
    except json.JSONDecodeError:
        return {
            "ok": False,
            "error_class": "RscriptJSONDecodeError",
            "message": txt[-5000:],
        }
    return out


def status_from_pair(py: Dict[str, Any], r: Dict[str, Any]) -> str:
    py_ok = bool(py.get("ok", False))
    r_ok = bool(r.get("ok", False))
    if py_ok and r_ok:
        return "both_ok"
    if (not py_ok) and (not r_ok):
        return "both_error"
    return "mismatch"


def contains_or_skip(message: str | None, needle: str | None) -> bool:
    if not needle:
        return True
    if not message:
        return False
    return needle.lower() in message.lower()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report-json", default="inst/extdata/benchmarks/failure_parity/report.json")
    parser.add_argument("--divergence-map", default="tools/parity_divergence_map.json")
    parser.add_argument("--rscript-exe", default="Rscript")
    parser.add_argument("--r-libs", default=None)
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    toolbox_root = repo_root / "ActflowToolbox"
    divergence_map_path = (repo_root / args.divergence_map).resolve()
    cfg = json.loads(divergence_map_path.read_text())

    mult_mod = load_module("multregconn", toolbox_root / "connectivity_estimation" / "multregconn.py")
    corr_mod = load_module("corrcoefconn", toolbox_root / "connectivity_estimation" / "corrcoefconn.py")
    compare_mod = load_module(
        "model_compare_predicted_to_actual",
        toolbox_root / "model_compare" / "model_compare_predicted_to_actual.py",
    )

    rng = np.random.default_rng(2026)

    py_scenarios: Dict[str, Callable[[], Any]] = {
        "multreg_nodes_gt_time": lambda: mult_mod.multregconn(rng.standard_normal((10, 5))),
        "model_compare_invalid_comparison": lambda: compare_mod.model_compare_predicted_to_actual(
            rng.standard_normal((6, 4, 2)),
            rng.standard_normal((6, 4, 2)),
            comparison_type="not_a_mode",
        ),
        "model_compare_dim_mismatch": lambda: compare_mod.model_compare_predicted_to_actual(
            rng.standard_normal((6, 4, 2)),
            rng.standard_normal((6, 5, 2)),
            comparison_type="conditionwise_compthenavg",
        ),
        "corr_single_timepoint": lambda: corr_mod.corrcoefconn(rng.standard_normal((6, 1))),
    }

    r_exprs: Dict[str, str] = {
        "multreg_nodes_gt_time": "x <- matrix(rnorm(50), nrow = 10, ncol = 5); estimate_fc_multreg(x)",
        "model_compare_invalid_comparison": "target <- array(rnorm(48), c(6, 4, 2)); pred <- array(rnorm(48), c(6, 4, 2)); model_compare(target, pred, comparison = 'not_a_mode')",
        "model_compare_dim_mismatch": "target <- array(rnorm(48), c(6, 4, 2)); pred <- array(rnorm(60), c(6, 5, 2)); model_compare(target, pred, comparison = 'conditionwise_compthenavg')",
        "corr_single_timepoint": "x <- matrix(rnorm(6), nrow = 6, ncol = 1); estimate_fc_corr(x)",
    }

    scenarios = []
    for sc in cfg.get("scenarios", []):
        sid = sc.get("id")
        if sid not in py_scenarios or sid not in r_exprs:
            scenarios.append(
                {
                    "id": sid,
                    "status": "not_implemented",
                    "pass": False,
                    "message": "Scenario missing implementation in failure_parity_runner.py",
                }
            )
            continue

        py_res = run_python_safely(py_scenarios[sid])
        r_res = run_r_expr(r_exprs[sid], rscript_exe=args.rscript_exe, r_libs=args.r_libs)

        actual = status_from_pair(py_res, r_res)
        expected = sc.get("expected_status")
        passed = actual == expected

        if passed and actual == "both_error":
            passed = contains_or_skip(r_res.get("message"), sc.get("r_error_contains")) and contains_or_skip(
                py_res.get("message"), sc.get("python_error_contains")
            )

        scenarios.append(
            {
                "id": sid,
                "expected_status": expected,
                "actual_status": actual,
                "parity_policy": sc.get("parity_policy"),
                "python": py_res,
                "r": r_res,
                "pass": bool(passed),
            }
        )

    total = len(scenarios)
    passed = sum(1 for s in scenarios if s.get("pass", False))

    report = {
        "schema_version": "1.0.0",
        "tool": "failure_parity_runner",
        "created_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "divergence_map": str(divergence_map_path),
        "summary": {
            "scenarios_total": total,
            "scenarios_passed": passed,
            "pass_ratio": (passed / total) if total > 0 else float("nan"),
        },
        "scenarios": scenarios,
    }

    out_report = (repo_root / args.report_json).resolve()
    out_report.parent.mkdir(parents=True, exist_ok=True)
    out_report.write_text(json.dumps(report, indent=2))

    print(f"Wrote failure parity report: {out_report}")
    print(json.dumps(report["summary"], indent=2))


if __name__ == "__main__":
    main()
