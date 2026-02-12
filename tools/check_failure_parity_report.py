#!/usr/bin/env python3
"""Validate failure parity report against thresholds."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path


def fail(msg: str) -> None:
    print(f"[FAIL] {msg}")
    sys.exit(1)


def ok(msg: str) -> None:
    print(f"[OK] {msg}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", required=True)
    parser.add_argument("--thresholds", required=True)
    args = parser.parse_args()

    report = json.loads(Path(args.report).read_text())
    th = json.loads(Path(args.thresholds).read_text())

    summary = report.get("summary", {})
    scenarios = report.get("scenarios", [])

    min_scenarios = int(th.get("minimum_scenarios", 1))
    if len(scenarios) < min_scenarios:
        fail(f"too few scenarios: got {len(scenarios)}, need >= {min_scenarios}")
    ok(f"scenario count {len(scenarios)} >= {min_scenarios}")

    pass_ratio = summary.get("pass_ratio", float("nan"))
    min_ratio = float(th.get("minimum_pass_ratio", 1.0))
    if not isinstance(pass_ratio, (int, float)) or not math.isfinite(pass_ratio) or pass_ratio < min_ratio:
        fail(f"pass_ratio below threshold: got {pass_ratio}, need >= {min_ratio}")
    ok(f"pass_ratio: {pass_ratio:.4f} >= {min_ratio:.4f}")

    for sid in th.get("required_scenarios", []):
        matched = [s for s in scenarios if s.get("id") == sid]
        if not matched:
            fail(f"required scenario missing: {sid}")
        if not matched[0].get("pass", False):
            fail(f"required scenario failed: {sid}")
        ok(f"required scenario passed: {sid}")

    print("[PASS] failure parity report satisfies thresholds.")


if __name__ == "__main__":
    main()
