#!/usr/bin/env python3
"""Validate aggregated parity fuzz report against acceptance thresholds."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict


def fail(msg: str) -> None:
    print(f"[FAIL] {msg}")
    sys.exit(1)


def ok(msg: str) -> None:
    print(f"[OK] {msg}")


def get_num(x: Any) -> float:
    if isinstance(x, (int, float)) and math.isfinite(x):
        return float(x)
    return float("nan")


def check_ratio(name: str, got: float, need_min: float) -> None:
    if not math.isfinite(got):
        fail(f"{name} is non-finite")
    if got < need_min:
        fail(f"{name} below threshold: got {got:.4f}, need >= {need_min:.4f}")
    ok(f"{name}: {got:.4f} >= {need_min:.4f}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", required=True)
    parser.add_argument("--thresholds", required=True)
    args = parser.parse_args()

    report = json.loads(Path(args.report).read_text())
    th = json.loads(Path(args.thresholds).read_text())

    summary = report.get("summary", {})
    by_op: Dict[str, Any] = report.get("by_operation", {})

    min_success = float(th.get("minimum_successful_case_ratio", 1.0))
    min_pass = float(th.get("minimum_threshold_pass_ratio", 1.0))
    check_ratio("successful_case_ratio", get_num(summary.get("successful_case_ratio")), min_success)
    check_ratio("threshold_pass_ratio", get_num(summary.get("threshold_pass_ratio")), min_pass)

    req_ops = th.get("operation_requirements", {})
    for op, cfg in req_ops.items():
        if op not in by_op:
            fail(f"required operation missing from report: {op}")
        op_data = by_op[op]
        sim = op_data.get("similarity", {})
        speed = op_data.get("speedup", {})

        min_cases = int(cfg.get("min_cases", 1))
        sim_count = int(sim.get("max_abs_err", {}).get("count", 0))
        if sim_count < min_cases:
            fail(f"{op} has insufficient cases: {sim_count} < {min_cases}")
        ok(f"{op} case_count {sim_count} >= {min_cases}")

        if "max_abs_err_p99_max" in cfg:
            got = get_num(sim.get("max_abs_err", {}).get("p99"))
            need = float(cfg["max_abs_err_p99_max"])
            if not math.isfinite(got) or got > need:
                fail(f"{op} max_abs_err p99 too high: got {got}, need <= {need}")
            ok(f"{op} max_abs_err p99 {got:.3e} <= {need:.3e}")

        if "corr_p01_min" in cfg:
            got = get_num(sim.get("corr", {}).get("min"))
            need = float(cfg["corr_p01_min"])
            if not math.isfinite(got) or got < need:
                fail(f"{op} corr min too low: got {got}, need >= {need}")
            ok(f"{op} corr min {got:.6f} >= {need:.6f}")

        if "speedup_p50_min" in cfg:
            got = get_num(speed.get("p50"))
            need = float(cfg["speedup_p50_min"])
            if not math.isfinite(got) or got < need:
                fail(f"{op} speedup p50 too low: got {got}, need >= {need}")
            ok(f"{op} speedup p50 {got:.3f} >= {need:.3f}")

    print("[PASS] parity fuzz report satisfies thresholds.")


if __name__ == "__main__":
    main()
