#!/usr/bin/env python3
"""Validate cross-language benchmark report against thresholds."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path


def fail(msg: str) -> None:
    print(f"[FAIL] {msg}")
    sys.exit(1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", required=True, help="Path to benchmark report JSON.")
    parser.add_argument("--thresholds", required=True, help="Path to threshold JSON.")
    args = parser.parse_args()

    report = json.loads(Path(args.report).read_text())
    th = json.loads(Path(args.thresholds).read_text())

    speedups = report.get("r_speedup_over_python", {})
    similarity = report.get("similarity", {})

    for op, min_val in th.get("speedup_min", {}).items():
        got = speedups.get(op, float("nan"))
        if not isinstance(got, (int, float)) or not math.isfinite(got):
            fail(f"speedup for {op} is non-finite: {got}")
        if got < float(min_val):
            fail(f"speedup for {op} below threshold: got {got:.6f}, need >= {float(min_val):.6f}")
        print(f"[OK] speedup {op}: {got:.6f} >= {float(min_val):.6f}")

    for op, max_abs in th.get("similarity_max_abs_err", {}).items():
        got = similarity.get(op, {}).get("max_abs_err", float("nan"))
        if not isinstance(got, (int, float)) or not math.isfinite(got):
            fail(f"similarity max_abs_err for {op} is non-finite: {got}")
        if got > float(max_abs):
            fail(f"max_abs_err for {op} above threshold: got {got:.3e}, need <= {float(max_abs):.3e}")
        print(f"[OK] max_abs_err {op}: {got:.3e} <= {float(max_abs):.3e}")

    for op, min_corr in th.get("similarity_min_corr", {}).items():
        got = similarity.get(op, {}).get("corr", float("nan"))
        if not isinstance(got, (int, float)) or not math.isfinite(got):
            fail(f"similarity corr for {op} is non-finite: {got}")
        if got < float(min_corr):
            fail(f"corr for {op} below threshold: got {got:.8f}, need >= {float(min_corr):.8f}")
        print(f"[OK] corr {op}: {got:.8f} >= {float(min_corr):.8f}")

    print("[PASS] cross-language benchmark report satisfies all thresholds.")


if __name__ == "__main__":
    main()
