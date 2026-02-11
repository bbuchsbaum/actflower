#!/usr/bin/env python3
"""Build parity fixtures from Python ActflowToolbox outputs.

This script intentionally avoids importing the full ActflowToolbox package
(whose optional dependencies may be missing) and instead loads specific modules
by file path.
"""

from __future__ import annotations

import argparse
import glob
import importlib.util
from pathlib import Path
from typing import Dict

import h5py
import numpy as np


def load_module(module_name: str, path: Path):
    spec = importlib.util.spec_from_file_location(module_name, str(path))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module from {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def r2_score(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    if ss_tot == 0:
        return np.nan
    return float(1.0 - ss_res / ss_tot)


def fullcompare_metrics(target: np.ndarray, pred: np.ndarray) -> Dict[str, np.ndarray]:
    """Match fullcompare_compthenavg summary behavior per subject."""
    n_subj = target.shape[2]
    corr = np.zeros(n_subj, dtype=np.float64)
    r2 = np.zeros(n_subj, dtype=np.float64)
    mae = np.zeros(n_subj, dtype=np.float64)

    for s in range(n_subj):
        yt = target[:, :, s].reshape(-1)
        yp = pred[:, :, s].reshape(-1)
        corr[s] = np.corrcoef(yt, yp)[0, 1]
        r2[s] = r2_score(yt, yp)
        mae[s] = np.mean(np.abs(yt - yp))

    return {
        "corr_vals": corr,
        "R2_vals": r2,
        "mae_vals": mae,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=None, help="Path to actflower repo root")
    parser.add_argument("--n-subj", type=int, default=5, help="Number of subjects to include")
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output HDF5 path (default: inst/extdata/parity_refs/hcp_parity_refs.h5)",
    )
    args = parser.parse_args()

    repo_root = args.repo_root.resolve() if args.repo_root else Path(__file__).resolve().parents[2]
    toolbox_root = repo_root / "ActflowToolbox"
    data_root = toolbox_root / "examples" / "HCP_example_data"

    if not toolbox_root.exists():
        raise FileNotFoundError(f"ActflowToolbox not found at: {toolbox_root}")
    if not data_root.exists():
        raise FileNotFoundError(f"HCP example data not found at: {data_root}")

    out_path = args.out.resolve() if args.out else (repo_root / "inst" / "extdata" / "parity_refs" / "hcp_parity_refs.h5")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    corr_mod = load_module("corrcoefconn", toolbox_root / "connectivity_estimation" / "corrcoefconn.py")
    mult_mod = load_module("multregconn", toolbox_root / "connectivity_estimation" / "multregconn.py")
    comb_mod = load_module("combinedFC", toolbox_root / "connectivity_estimation" / "combinedFC.py")
    calc_mod = load_module("actflowcalc", toolbox_root / "actflowcomp" / "actflowcalc.py")

    subj_files = sorted(glob.glob(str(data_root / "HCP_example_restrun1_subj*_data.h5")))
    if len(subj_files) == 0:
        raise RuntimeError("No subject rest files found.")

    n_subj = min(args.n_subj, len(subj_files))
    subj_files = subj_files[:n_subj]
    subj_ids = [Path(f).stem.split("subj")[-1].split("_data")[0] for f in subj_files]

    with h5py.File(data_root / "HCP_example_taskactivations_data.h5", "r") as h:
        taskbeta_all = h["taskbeta"][:]
    taskbeta = taskbeta_all[:, :, :n_subj]

    n_nodes = taskbeta.shape[0]
    n_cond = taskbeta.shape[1]

    rest = None
    fc_corr = np.zeros((n_nodes, n_nodes, n_subj), dtype=np.float64)
    fc_multreg = np.zeros((n_nodes, n_nodes, n_subj), dtype=np.float64)
    fc_combined = np.zeros((n_nodes, n_nodes, n_subj), dtype=np.float64)

    print(f"Building parity refs for {n_subj} subjects, {n_nodes} nodes, {n_cond} conditions")

    for s, fpath in enumerate(subj_files):
        with h5py.File(fpath, "r") as h:
            rs = h["restdata"][:]

        if rest is None:
            rest = np.zeros((rs.shape[0], rs.shape[1], n_subj), dtype=np.float64)
        rest[:, :, s] = rs

        print(f"  [{s + 1}/{n_subj}] FC corr")
        fc_corr[:, :, s] = corr_mod.corrcoefconn(rs)

        print(f"  [{s + 1}/{n_subj}] FC multreg")
        fc_multreg[:, :, s] = mult_mod.multregconn(rs)

        print(f"  [{s + 1}/{n_subj}] FC combined")
        fc_combined[:, :, s] = comb_mod.combinedFC(rs)

    pred_corr = np.zeros((n_nodes, n_cond, n_subj), dtype=np.float64)
    pred_multreg = np.zeros((n_nodes, n_cond, n_subj), dtype=np.float64)
    pred_combined = np.zeros((n_nodes, n_cond, n_subj), dtype=np.float64)

    for s in range(n_subj):
        for c in range(n_cond):
            pred_corr[:, c, s] = calc_mod.actflowcalc(taskbeta[:, c, s], fc_corr[:, :, s])
            pred_multreg[:, c, s] = calc_mod.actflowcalc(taskbeta[:, c, s], fc_multreg[:, :, s])
            pred_combined[:, c, s] = calc_mod.actflowcalc(taskbeta[:, c, s], fc_combined[:, :, s])

    metrics_corr = fullcompare_metrics(taskbeta, pred_corr)
    metrics_multreg = fullcompare_metrics(taskbeta, pred_multreg)
    metrics_combined = fullcompare_metrics(taskbeta, pred_combined)

    with h5py.File(out_path, "w") as h:
        h.create_dataset("subj_ids", data=np.array(subj_ids, dtype="S"))
        h.create_dataset("restdata", data=rest)
        h.create_dataset("taskbeta", data=taskbeta)

        h.create_dataset("fc_corr", data=fc_corr)
        h.create_dataset("fc_multreg", data=fc_multreg)
        h.create_dataset("fc_combined", data=fc_combined)
        h.create_dataset("pred_corr", data=pred_corr)
        h.create_dataset("pred_multreg", data=pred_multreg)
        h.create_dataset("pred_combined", data=pred_combined)

        h.create_dataset("metrics_corr/corr_vals", data=metrics_corr["corr_vals"])
        h.create_dataset("metrics_corr/R2_vals", data=metrics_corr["R2_vals"])
        h.create_dataset("metrics_corr/mae_vals", data=metrics_corr["mae_vals"])

        h.create_dataset("metrics_multreg/corr_vals", data=metrics_multreg["corr_vals"])
        h.create_dataset("metrics_multreg/R2_vals", data=metrics_multreg["R2_vals"])
        h.create_dataset("metrics_multreg/mae_vals", data=metrics_multreg["mae_vals"])

        h.create_dataset("metrics_combined/corr_vals", data=metrics_combined["corr_vals"])
        h.create_dataset("metrics_combined/R2_vals", data=metrics_combined["R2_vals"])
        h.create_dataset("metrics_combined/mae_vals", data=metrics_combined["mae_vals"])

    print(f"Wrote parity fixtures: {out_path}")


if __name__ == "__main__":
    main()
