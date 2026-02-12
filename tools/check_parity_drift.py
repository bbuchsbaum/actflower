#!/usr/bin/env python3
"""Detect parity drift relative to a baseline parity-fuzz report."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any


def fail(msg: str) -> None:
    print(f"[FAIL] {msg}")
    sys.exit(1)


def ok(msg: str) -> None:
    print(f"[OK] {msg}")


def get_num(x: Any) -> float:
    if isinstance(x, (int, float)) and math.isfinite(x):
        return float(x)
    return float("nan")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--current", required=True)
    parser.add_argument("--thresholds", required=True)
    args = parser.parse_args()

    base = json.loads(Path(args.baseline).read_text())
    cur = json.loads(Path(args.current).read_text())
    th = json.loads(Path(args.thresholds).read_text())

    req_ops = th.get("required_operations", [])
    mul = float(th.get("max_abs_err_multiplier", 10.0))
    abs_floor = float(th.get("max_abs_err_abs_floor", 1e-8))
    corr_drop = float(th.get("corr_drop_tolerance", 1e-5))
    corr_floor = float(th.get("corr_floor", 0.9999))

    b_ops = base.get("by_operation", {})
    c_ops = cur.get("by_operation", {})

    for op in req_ops:
      if op not in b_ops:
          fail(f"operation missing from baseline report: {op}")
      if op not in c_ops:
          fail(f"operation missing from current report: {op}")

      b_err = get_num(b_ops[op].get("similarity", {}).get("max_abs_err", {}).get("p99"))
      c_err = get_num(c_ops[op].get("similarity", {}).get("max_abs_err", {}).get("p99"))
      err_cap = max(abs_floor, b_err * mul) if math.isfinite(b_err) else abs_floor
      if not math.isfinite(c_err) or c_err > err_cap:
          fail(f"{op} drift on max_abs_err p99: current={c_err}, cap={err_cap}, baseline={b_err}")
      ok(f"{op} max_abs_err p99 {c_err:.3e} within cap {err_cap:.3e}")

      b_corr = get_num(b_ops[op].get("similarity", {}).get("corr", {}).get("min"))
      c_corr = get_num(c_ops[op].get("similarity", {}).get("corr", {}).get("min"))
      corr_min = max(corr_floor, b_corr - corr_drop) if math.isfinite(b_corr) else corr_floor
      if not math.isfinite(c_corr) or c_corr < corr_min:
          fail(f"{op} drift on corr min: current={c_corr}, floor={corr_min}, baseline={b_corr}")
      ok(f"{op} corr min {c_corr:.6f} >= {corr_min:.6f}")

    print("[PASS] parity drift check satisfied.")


if __name__ == "__main__":
    main()
