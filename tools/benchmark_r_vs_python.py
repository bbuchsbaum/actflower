#!/usr/bin/env python3
"""Cross-language benchmark: actflower (R) vs ActflowToolbox (Python).

Runs both implementations on the same synthetic data, compares numerical
similarity, and reports runtime speed ratios.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, Tuple

import h5py
import numpy as np


def load_module(module_name: str, path: Path):
    spec = importlib.util.spec_from_file_location(module_name, str(path))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module from {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def timed(func, repeats: int) -> Tuple[Any, Dict[str, Any]]:
    times = []
    result = None
    for _ in range(repeats):
        t0 = time.perf_counter()
        result = func()
        times.append(time.perf_counter() - t0)
    arr = np.asarray(times, dtype=float)
    return result, {
        "median_sec": float(np.median(arr)),
        "mean_sec": float(np.mean(arr)),
        "times": [float(x) for x in arr.tolist()],
    }


def compare_arrays(a: np.ndarray, b: np.ndarray) -> Dict[str, float]:
    a = np.asarray(a, dtype=float).ravel()
    b = np.asarray(b, dtype=float).ravel()
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() == 0:
        return {"max_abs_err": float("nan"), "mean_abs_err": float("nan"), "corr": float("nan")}
    a = a[mask]
    b = b[mask]
    diff = np.abs(a - b)
    if a.size > 1:
        corr = float(np.corrcoef(a, b)[0, 1])
    else:
        corr = float("nan")
    return {
        "max_abs_err": float(np.max(diff)),
        "mean_abs_err": float(np.mean(diff)),
        "corr": corr,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--nodes", type=int, default=120)
    parser.add_argument("--timepoints", type=int, default=300)
    parser.add_argument("--conditions", type=int, default=16)
    parser.add_argument("--subjects", type=int, default=8)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--use-cpp", action="store_true", default=True)
    parser.add_argument("--include-combined", action="store_true", default=False)
    parser.add_argument("--python-exe", default=sys.executable)
    parser.add_argument("--rscript-exe", default="Rscript")
    parser.add_argument("--r-libs", default=None, help="Optional R library path containing installed actflower.")
    parser.add_argument("--work-dir", default="inst/extdata/benchmarks/r_vs_python")
    parser.add_argument("--report-json", default="inst/extdata/benchmarks/r_vs_python/report.json")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    toolbox_root = repo_root / "ActflowToolbox"
    work_dir = (repo_root / args.work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    if not toolbox_root.exists():
        raise FileNotFoundError(f"ActflowToolbox not found at {toolbox_root}")

    corr_mod = load_module("corrcoefconn", toolbox_root / "connectivity_estimation" / "corrcoefconn.py")
    mult_mod = load_module("multregconn", toolbox_root / "connectivity_estimation" / "multregconn.py")
    act_mod = load_module("actflowcalc", toolbox_root / "actflowcomp" / "actflowcalc.py")
    combined_mod = None
    if args.include_combined:
        combined_mod = load_module("combinedFC", toolbox_root / "connectivity_estimation" / "combinedFC.py")

    rng = np.random.default_rng(args.seed)
    rest = rng.standard_normal((args.nodes, args.timepoints, args.subjects)).astype(np.float64)
    task = rng.standard_normal((args.nodes, args.conditions, args.subjects)).astype(np.float64)

    input_h5 = work_dir / "input.h5"
    py_out_h5 = work_dir / "python_outputs.h5"
    py_json = work_dir / "python_timings.json"
    r_out_h5 = work_dir / "r_outputs.h5"
    r_json = work_dir / "r_timings.json"

    with h5py.File(input_h5, "w") as h5:
        h5.create_dataset("restdata", data=rest)
        h5.create_dataset("taskbeta", data=task)

    def run_fc_corr_py() -> np.ndarray:
        out = np.empty((args.nodes, args.nodes, args.subjects), dtype=np.float64)
        for s in range(args.subjects):
            out[:, :, s] = corr_mod.corrcoefconn(rest[:, :, s])
        return out

    def run_fc_multreg_py() -> np.ndarray:
        out = np.empty((args.nodes, args.nodes, args.subjects), dtype=np.float64)
        for s in range(args.subjects):
            out[:, :, s] = mult_mod.multregconn(rest[:, :, s])
        return out

    def pred_from_fc_py(fc: np.ndarray) -> np.ndarray:
        out = np.empty((args.nodes, args.conditions, args.subjects), dtype=np.float64)
        for s in range(args.subjects):
            for c in range(args.conditions):
                out[:, c, s] = act_mod.actflowcalc(task[:, c, s], fc[:, :, s], separate_activations_bytarget=False, transfer_func=None)
        return out

    def run_fc_combined_py() -> np.ndarray:
        if combined_mod is None:
            raise RuntimeError("combinedFC module not loaded.")
        out = np.empty((args.nodes, args.nodes, args.subjects), dtype=np.float64)
        for s in range(args.subjects):
            out[:, :, s] = combined_mod.combinedFC(rest[:, :, s])
        return out

    fc_corr_py, t_corr_py = timed(run_fc_corr_py, args.repeats)
    fc_multreg_py, t_multreg_py = timed(run_fc_multreg_py, args.repeats)
    pred_corr_py, t_pred_corr_py = timed(lambda: pred_from_fc_py(fc_corr_py), args.repeats)
    pred_multreg_py, t_pred_multreg_py = timed(lambda: pred_from_fc_py(fc_multreg_py), args.repeats)
    combined_payload_py = {}
    if args.include_combined:
        fc_combined_py, t_combined_py = timed(run_fc_combined_py, args.repeats)
        pred_combined_py, t_pred_combined_py = timed(lambda: pred_from_fc_py(fc_combined_py), args.repeats)
        combined_payload_py = {
            "fc_combined": t_combined_py,
            "pred_combined": t_pred_combined_py,
        }

    py_timings = {
        "engine": "Python_ActflowToolbox",
        "include_combined": args.include_combined,
        "repeats": args.repeats,
        "dimensions": {
            "nodes": args.nodes,
            "timepoints": args.timepoints,
            "conditions": args.conditions,
            "subjects": args.subjects,
        },
        "operations": {
            "fc_corr": t_corr_py,
            "fc_multreg": t_multreg_py,
            "pred_corr": t_pred_corr_py,
            "pred_multreg": t_pred_multreg_py,
        },
    }
    py_timings["operations"].update(combined_payload_py)

    with h5py.File(py_out_h5, "w") as h5:
        h5.create_dataset("fc_corr", data=fc_corr_py)
        h5.create_dataset("fc_multreg", data=fc_multreg_py)
        h5.create_dataset("pred_corr", data=pred_corr_py)
        h5.create_dataset("pred_multreg", data=pred_multreg_py)
        if args.include_combined:
            h5.create_dataset("fc_combined", data=fc_combined_py)
            h5.create_dataset("pred_combined", data=pred_combined_py)

    py_json.write_text(json.dumps(py_timings, indent=2))

    r_script = repo_root / "tools" / "benchmark_r_impl.R"
    cmd = [
        args.rscript_exe,
        str(r_script),
        "--input-h5",
        str(input_h5),
        "--out-h5",
        str(r_out_h5),
        "--out-json",
        str(r_json),
        "--repeats",
        str(args.repeats),
        "--use-cpp",
        "true" if args.use_cpp else "false",
        "--include-combined",
        "true" if args.include_combined else "false",
    ]
    env = dict(os.environ)
    if args.r_libs:
        env["R_LIBS"] = args.r_libs
    subprocess.run(cmd, cwd=str(repo_root), check=True, env=env)

    r_timings = json.loads(r_json.read_text())
    with h5py.File(r_out_h5, "r") as h5:
        fc_corr_r = h5["fc_corr"][:]
        fc_multreg_r = h5["fc_multreg"][:]
        pred_corr_r = h5["pred_corr"][:]
        pred_multreg_r = h5["pred_multreg"][:]
        combined_payload_r = {}
        if args.include_combined:
            combined_payload_r["fc_combined"] = h5["fc_combined"][:]
            combined_payload_r["pred_combined"] = h5["pred_combined"][:]

    similarity = {
        "fc_corr": compare_arrays(fc_corr_r, fc_corr_py),
        "fc_multreg": compare_arrays(fc_multreg_r, fc_multreg_py),
        "pred_corr": compare_arrays(pred_corr_r, pred_corr_py),
        "pred_multreg": compare_arrays(pred_multreg_r, pred_multreg_py),
    }
    if args.include_combined:
        similarity["fc_combined"] = compare_arrays(combined_payload_r["fc_combined"], fc_combined_py)
        similarity["pred_combined"] = compare_arrays(combined_payload_r["pred_combined"], pred_combined_py)

    def speed_ratio(op: str) -> float:
        py_sec = py_timings["operations"][op]["median_sec"]
        r_sec = r_timings["operations"][op]["median_sec"]
        if not np.isfinite(py_sec) or not np.isfinite(r_sec) or r_sec <= 0:
            return float("nan")
        return float(py_sec / r_sec)

    ops = ["fc_corr", "fc_multreg", "pred_corr", "pred_multreg"]
    if args.include_combined:
        ops.extend(["fc_combined", "pred_combined"])
    speedups = {op: speed_ratio(op) for op in ops}

    report = {
        "seed": args.seed,
        "repeats": args.repeats,
        "dimensions": py_timings["dimensions"],
        "python_timings": py_timings["operations"],
        "r_timings": r_timings["operations"],
        "r_speedup_over_python": speedups,
        "similarity": similarity,
        "artifacts": {
            "input_h5": str(input_h5),
            "python_outputs_h5": str(py_out_h5),
            "r_outputs_h5": str(r_out_h5),
            "python_timings_json": str(py_json),
            "r_timings_json": str(r_json),
        },
    }

    out_report = (repo_root / args.report_json).resolve()
    out_report.parent.mkdir(parents=True, exist_ok=True)
    out_report.write_text(json.dumps(report, indent=2))

    print(f"Wrote report: {out_report}")
    print("R speedup over Python (median time ratio; >1 means R faster):")
    for op, val in speedups.items():
        print(f"  {op:12s}: {val:.3f}x")
    print("Similarity (R vs Python):")
    for op, met in similarity.items():
        print(
            f"  {op:12s}: "
            f"corr={met['corr']:.6f}, "
            f"max_abs_err={met['max_abs_err']:.3e}, "
            f"mean_abs_err={met['mean_abs_err']:.3e}"
        )


if __name__ == "__main__":
    main()
