#!/usr/bin/env python
"""Evaluate the breakpoint hot list by TRANSPORT into the restricted subset.

Why this script exists
----------------------
A first attempt refit the model *inside* the restricted subset and produced
near-chance results (SCNA AUC 0.4999 at 3 years). That is not a finding about
the restricted subset — it is a finding about fitting 47 features to 22 events.

The locked model does something different, and the handoff methods are explicit
about it: broad-cohort models are applied to the restricted subset **without
subset-specific refitting**. It tests transport into a cleaner clinical setting,
not the performance of a separately trained restricted model.

This script reproduces that design:

    1. Fit on the BROAD cohort with patient-grouped out-of-fold prediction, so
       every tumour receives a score from a model that never saw it or any other
       tumour from the same patient.
    2. Subset those honest out-of-fold scores to the restricted tumours.
    3. Compute AUC there.

Window selection still happens inside the broad training folds, so discovery
cost is still paid.

Locked-model reference (restricted subset):

    horizon   n    events   SCNA     SCNA+breakpoint
    3 y       387  22       0.7792   0.7866
    5 y       248  37       0.7601   0.8167
    7 y       155  43       0.7205   0.7950

Usage:
    python scripts/hotlist_transport.py --k 8 --repeats 10
"""
from __future__ import annotations
import os

import argparse
import logging
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

from src.evaluation.splits import assert_no_patient_leakage, patient_grouped_folds
from src.features.architecture import build_feature_table, feature_groups
from src.features.clinical import load_restricted_flag
from src.features.hotlist import build_window_matrix, consensus_hotlist, select_in_fold
from src.harmonization.breakpoints import derive_breakpoints, tile_arms
from src.harmonization.institutional import (
    fixed_horizon_labels, load_master_seg, load_outcomes,
)
from src.harmonization.references import ArmReference
from src.harmonization.schema import BoundaryClass

H = os.environ.get("MENIBREAK_DATA_ROOT", "E:/Breakpoint_Architecture_Paper_Analysis_Handoff_2026-08-02/")
REF = H + "02_REFERENCE_AND_PROCESSED_DATA/"
log = logging.getLogger("transport")

LOCKED_RESTRICTED = {
    3: {"n": 387, "events": 22, "scna": 0.7792, "full": 0.7866},
    5: {"n": 248, "events": 37, "scna": 0.7601, "full": 0.8167},
    7: {"n": 155, "events": 43, "scna": 0.7205, "full": 0.7950},
}


def oof_scores(
    base: pd.DataFrame,
    windows: pd.DataFrame | None,
    y: np.ndarray,
    groups: pd.DataFrame,
    *,
    k: int,
    n_repeats: int,
    seed: int,
    n_splits: int = 5,
) -> tuple[np.ndarray, list[list[str]]]:
    """Broad-cohort out-of-fold scores, averaged over repeats."""
    acc = np.full((len(y), n_repeats), np.nan)
    selections: list[list[str]] = []

    for rep in range(n_repeats):
        for tr, te in patient_grouped_folds(groups, n_splits=n_splits, seed=seed + rep):
            assert_no_patient_leakage(groups, tr, te)
            if len(np.unique(y[tr])) < 2:
                continue

            Xtr = base.iloc[tr].reset_index(drop=True)
            Xte = base.iloc[te].reset_index(drop=True)

            if windows is not None and k > 0:
                chosen = select_in_fold(windows.iloc[tr], y[tr], k=k)
                selections.append(chosen)
                if chosen:
                    Xtr = pd.concat([Xtr, windows.iloc[tr][chosen].reset_index(drop=True)],
                                    axis=1)
                    Xte = pd.concat([Xte, windows.iloc[te][chosen].reset_index(drop=True)],
                                    axis=1)

            a, b = Xtr.to_numpy(float), Xte.to_numpy(float)
            mu, sd = a.mean(0), a.std(0)
            sd = np.where(sd > 0, sd, 1.0)
            m = LogisticRegression(penalty="l2", C=1.0, max_iter=3000).fit((a - mu) / sd, y[tr])
            acc[te, rep] = m.predict_proba((b - mu) / sd)[:, 1]

    return np.nanmean(acc, axis=1), selections


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--horizons", nargs="+", type=int, default=[3, 5, 7])
    ap.add_argument("--k", type=int, default=8)
    ap.add_argument("--n-bins", type=int, default=20)
    ap.add_argument("--repeats", type=int, default=10)
    ap.add_argument("--min-prevalence", type=int, default=20)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--out", default="reports/hotlist_transport.csv")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        datefmt="%H:%M:%S")

    seg = load_master_seg(
        H + "01_RAW_DATA/institutional/meningioma_masterseg_ALL_2026.txt",
        forbidden_subpaths=["RESTRICTED"])
    ar = ArmReference.from_csv(REF + "arm_reference_GRCh37_software.csv", "GRCh37")
    tiles, _ = tile_arms(seg, ar)
    bps = derive_breakpoints(tiles, ar, boundary_class=BoundaryClass.SNP_ARRAY_CHANGEPOINT)
    m = seg[["sample_id", "tumor_id"]].drop_duplicates()
    bps = bps.merge(m, on="sample_id", how="left")
    tiles["unit"] = tiles["tumor_id"]
    bps["unit"] = bps["tumor_id"]

    base_feats = build_feature_table(tiles, bps, ar, sample_col="unit")
    scna = base_feats[feature_groups(base_feats.columns)["scna"]]
    windows = build_window_matrix(bps, sample_col="unit", n_bins=args.n_bins,
                                  binary=True, min_prevalence=args.min_prevalence)
    restricted = load_restricted_flag(
        H + "03_ANALYSIS_READY_DATA/ml_breakpoint_model_dataset_2026.csv")
    outcomes = load_outcomes(H + "03_ANALYSIS_READY_DATA/outcomes.csv")

    rows = []
    hot_by_h = {}
    for horizon in args.horizons:
        labels = fixed_horizon_labels(outcomes, horizon)
        d = labels.merge(scna, left_on="tumor_id", right_index=True, how="inner")
        d = d.merge(windows, left_on="tumor_id", right_index=True, how="inner")
        d = d.merge(restricted, left_on="tumor_id", right_index=True, how="left")
        d["is_restricted"] = d["is_restricted"].fillna(False)

        y = d["label"].to_numpy(int)
        groups = d[["patient_id"]].reset_index(drop=True)
        S = d[scna.columns].reset_index(drop=True)
        W = d[[c for c in windows.columns if c in d.columns]].reset_index(drop=True)
        rmask = d["is_restricted"].to_numpy(bool)

        s_scna, _ = oof_scores(S, None, y, groups, k=0,
                               n_repeats=args.repeats, seed=args.seed)
        s_full, sels = oof_scores(S, W, y, groups, k=args.k,
                                  n_repeats=args.repeats, seed=args.seed)
        # Windows ALONE: an empty base, so the model sees only the k selected
        # windows. This isolates how much the panel carries by itself, which is
        # a different question from how much it adds.
        empty = pd.DataFrame(index=range(len(y)))
        empty["_intercept"] = 1.0
        s_win, _ = oof_scores(empty, W, y, groups, k=args.k,
                              n_repeats=args.repeats, seed=args.seed)

        hot_by_h[horizon] = consensus_hotlist(sels, n_bins=args.n_bins)

        for scope, mask in (("broad", np.ones(len(y), bool)), ("restricted", rmask)):
            yy = y[mask]
            if len(np.unique(yy)) < 2:
                continue
            a_s = roc_auc_score(yy, s_scna[mask])
            a_f = roc_auc_score(yy, s_full[mask])
            a_w = roc_auc_score(yy, s_win[mask])
            ref = LOCKED_RESTRICTED[horizon] if scope == "restricted" else None
            rows.append({
                "horizon": horizon, "scope": scope, "n": int(mask.sum()),
                "events": int(yy.sum()),
                "auc_windows_only": a_w,
                "auc_scna": a_s,
                f"auc_scna_plus_{args.k}windows": a_f, "delta": a_f - a_s,
                "locked_scna": ref["scna"] if ref else np.nan,
                "locked_full": ref["full"] if ref else np.nan,
            })
            log.info("%dy %-11s n=%-4d ev=%-3d win_only=%.4f SCNA=%.4f  +%dwin=%.4f  (%+.4f)",
                     horizon, scope, int(mask.sum()), int(yy.sum()),
                     a_w, a_s, args.k, a_f, a_f - a_s)

    res = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(args.out, index=False)

    print("\n" + "=" * 86)
    print(f"TRANSPORT EVALUATION — broad-cohort model, {args.k} in-fold-selected windows")
    print("Locked design: broad model applied to the restricted subset, no refitting")
    print("=" * 86)
    print(res.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print(f"\nsaved {args.out}")

    for h, hl in hot_by_h.items():
        if not hl.empty:
            hl.insert(0, "horizon", h)
    pd.concat([v for v in hot_by_h.values() if not v.empty], ignore_index=True) \
        .to_csv(args.out.replace(".csv", "_windows.csv"), index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
