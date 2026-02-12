#!/usr/bin/env python3
"""Validate stochastic parity report against threshold criteria."""

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


def get_path(dct: dict[str, Any], path: str) -> Any:
    cur: Any = dct
    for p in path.split("."):
        if not isinstance(cur, dict) or p not in cur:
            return None
        cur = cur[p]
    return cur


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", required=True)
    parser.add_argument("--thresholds", required=True)
    args = parser.parse_args()

    report = json.loads(Path(args.report).read_text())
    th = json.loads(Path(args.thresholds).read_text())

    pass_ratio = report.get("summary", {}).get("pass_ratio", float("nan"))
    min_ratio = float(th.get("minimum_pass_ratio", 1.0))
    if not isinstance(pass_ratio, (int, float)) or not math.isfinite(pass_ratio) or pass_ratio < min_ratio:
        fail(f"pass_ratio below threshold: got {pass_ratio}, need >= {min_ratio}")
    ok(f"pass_ratio: {pass_ratio:.4f} >= {min_ratio:.4f}")

    for path in th.get("required_true", []):
        val = get_path(report, path)
        if val is not True:
            fail(f"required check is not TRUE: {path} (got {val})")
        ok(f"{path} == TRUE")

    print("[PASS] stochastic parity report satisfies thresholds.")


if __name__ == "__main__":
    main()
