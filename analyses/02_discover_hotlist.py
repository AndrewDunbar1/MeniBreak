#!/usr/bin/env python
"""Discover a parsimonious breakpoint hot list and test it against WHO / IG / SCNA.

Question: can fewer than ~10 specific breakpoint windows add prognostic
information over WHO grade, integrated grade, and arm-level SCNA?

Selection happens **inside training folds only**. Each fold ranks candidate
windows using its own training rows, fits on those, and scores rows it has never
seen. The reported AUC therefore includes the cost of discovery, unlike any
number computed with the frozen windows.

The consensus hot list is reported separately, with the fraction of folds that
selected each window, so a reader can see which windows are stable and which are
fold-specific noise.

Usage:
    python scripts/discover_hotlist.py --horizons 3 5 7 --k 3 5 10 20
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
from src.features.hotlist import (
    build_window_matrix, consensus_hotlist, select_in_fold,
    windows_to_genomic_coordinates,
)
from src.harmonization.breakpoints import derive_breakpoints, tile_arms
from src.harmonization.institutional import (
    fixed_horizon_labels, load_master_seg, load_outcomes,
)
from src.harmonization.references import ArmReference
from src.harmonization.schema import BoundaryClass

H = os.environ.get("MENIBREAK_DATA_ROOT", "E:/Breakpoint_Architecture_Paper_Analysis_Handoff_2026-08-02/")
REF = H + "02_REFERENCE_AND_PROCESSED_DATA/"
log = logging.getLogger("hotlist")

#: Comparator AUCs from the locked project, on common comparator tumours.
LOCKED_COMPARATORS = {
    3: {"WHO": 0.690, "integrated_grade": 0.766},
    5: {"WHO": 0.656, "integrated_grade": 0.735},
    7: {"WHO": 0.635, "integrated_grade": 0.733},
}


def oof_auc_with_in_fold_selection(
    windows: pd.DataFrame,
    base: pd.DataFrame | None,
    y: np.ndarray,
    groups: pd.DataFrame,
    *,
    k: int,
    n_repeats: int,
    seed: int,
    n_splits: int = 5,
) -> tuple[float, list[list[str]]]:
    """Repeated patient-grouped CV where each fold picks its OWN top-k windows."""
    aucs: list[float] = []
    selections: list[list[str]] = []

    for rep in range(n_repeats):
        scores = np.full(len(y), np.nan)
        for tr, te in patient_grouped_folds(groups, n_splits=n_splits, seed=seed + rep):
            assert_no_patient_leakage(groups, tr, te)
            if len(np.unique(y[tr])) < 2:
                continue

            chosen = select_in_fold(windows.iloc[tr], y[tr], k=k)
            selections.append(chosen)
            if not chosen:
                continue

            Xtr = windows.iloc[tr][chosen]
            Xte = windows.iloc[te][chosen]
            if base is not None:
                Xtr = pd.concat([base.iloc[tr].reset_index(drop=True),
                                 Xtr.reset_index(drop=True)], axis=1)
                Xte = pd.concat([base.iloc[te].reset_index(drop=True),
                                 Xte.reset_index(drop=True)], axis=1)

            a, b = Xtr.to_numpy(float), Xte.to_numpy(float)
            mu, sd = a.mean(0), a.std(0)
            sd = np.where(sd > 0, sd, 1.0)
            m = LogisticRegression(penalty="l2", C=1.0, max_iter=3000).fit((a - mu) / sd, y[tr])
            scores[te] = m.predict_proba((b - mu) / sd)[:, 1]

        ok = ~np.isnan(scores)
        if len(np.unique(y[ok])) > 1:
            aucs.append(roc_auc_score(y[ok], scores[ok]))

    return (float(np.mean(aucs)) if aucs else float("nan")), selections


def plain_auc(X: pd.DataFrame, y: np.ndarray, groups: pd.DataFrame,
              *, n_repeats: int, seed: int, n_splits: int = 5) -> float:
    """Baseline AUC with no feature selection."""
    Xv = X.to_numpy(float)
    aucs = []
    for rep in range(n_repeats):
        scores = np.full(len(y), np.nan)
        for tr, te in patient_grouped_folds(groups, n_splits=n_splits, seed=seed + rep):
            if len(np.unique(y[tr])) < 2:
                continue
            mu, sd = Xv[tr].mean(0), Xv[tr].std(0)
            sd = np.where(sd > 0, sd, 1.0)
            m = LogisticRegression(penalty="l2", C=1.0, max_iter=3000).fit((Xv[tr] - mu) / sd, y[tr])
            scores[te] = m.predict_proba((Xv[te] - mu) / sd)[:, 1]
        ok = ~np.isnan(scores)
        if len(np.unique(y[ok])) > 1:
            aucs.append(roc_auc_score(y[ok], scores[ok]))
    return float(np.mean(aucs)) if aucs else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--horizons", nargs="+", type=int, default=[3, 5, 7])
    ap.add_argument("--k", nargs="+", type=int, default=[3, 5, 10, 20])
    ap.add_argument("--n-bins", type=int, default=20)
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--min-prevalence", type=int, default=20)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--restricted", action="store_true",
                    help="limit to the locked model's restricted clinical subset "
                         "(primary, GTR, no prior/adjuvant RT, not radiation-induced, "
                         "non-angiomatous)")
    ap.add_argument("--out-prefix", default="reports/hotlist")
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
    log.info("candidate windows: %d (%d bins/arm, prevalence >= %d)",
             windows.shape[1], args.n_bins, args.min_prevalence)

    outcomes = load_outcomes(H + "03_ANALYSIS_READY_DATA/outcomes.csv")

    restricted = None
    if args.restricted:
        from src.features.clinical import load_restricted_flag
        restricted = load_restricted_flag(
            H + "03_ANALYSIS_READY_DATA/ml_breakpoint_model_dataset_2026.csv")

    rows = []
    all_lists = []

    for horizon in args.horizons:
        labels = fixed_horizon_labels(outcomes, horizon)
        if restricted is not None:
            labels = labels.merge(restricted, left_on="tumor_id",
                                  right_index=True, how="inner")
            labels = labels[labels["is_restricted"]].drop(columns=["is_restricted"])
            if labels["label"].sum() < 15:
                log.warning("%dy restricted: only %d events -- selection of %s "
                            "windows is severely underpowered; treat as exploratory",
                            horizon, int(labels["label"].sum()), args.k)
        d = labels.merge(scna, left_on="tumor_id", right_index=True, how="inner")
        d = d.merge(windows, left_on="tumor_id", right_index=True, how="inner")
        y = d["label"].to_numpy(int)
        groups = d[["patient_id"]].reset_index(drop=True)
        S = d[scna.columns].reset_index(drop=True)
        W = d[[c for c in windows.columns if c in d.columns]].reset_index(drop=True)

        scna_auc = plain_auc(S, y, groups, n_repeats=args.repeats, seed=args.seed)
        comp = LOCKED_COMPARATORS[horizon]
        log.info("%dy n=%d events=%d | SCNA=%.4f WHO=%.3f IG=%.3f",
                 horizon, len(d), int(y.sum()), scna_auc, comp["WHO"],
                 comp["integrated_grade"])
        rows.append({"horizon": horizon, "model": "SCNA_alone", "k": 0,
                     "n": len(d), "events": int(y.sum()), "auc": scna_auc})
        for nm, v in comp.items():
            rows.append({"horizon": horizon, "model": nm, "k": 0,
                         "n": len(d), "events": int(y.sum()), "auc": v})

        for k in args.k:
            auc_w, sels = oof_auc_with_in_fold_selection(
                W, None, y, groups, k=k, n_repeats=args.repeats, seed=args.seed)
            auc_sw, sels2 = oof_auc_with_in_fold_selection(
                W, S, y, groups, k=k, n_repeats=args.repeats, seed=args.seed)
            rows.append({"horizon": horizon, "model": f"windows_only_k{k}", "k": k,
                         "n": len(d), "events": int(y.sum()), "auc": auc_w})
            rows.append({"horizon": horizon, "model": f"SCNA+windows_k{k}", "k": k,
                         "n": len(d), "events": int(y.sum()), "auc": auc_sw})
            log.info("  k=%-3d windows_only=%.4f  SCNA+windows=%.4f  (vs SCNA %.4f, %+0.4f)",
                     k, auc_w, auc_sw, scna_auc, auc_sw - scna_auc)

            hl = consensus_hotlist(sels2, n_bins=args.n_bins)
            hl.insert(0, "horizon", horizon)
            hl.insert(1, "k", k)
            all_lists.append(hl)

    res = pd.DataFrame(rows)
    Path(args.out_prefix).parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(f"{args.out_prefix}_performance.csv", index=False)

    print("\n" + "=" * 80)
    print("PERFORMANCE vs COMPARATORS  (selection done inside training folds)")
    print("=" * 80)
    piv = res.pivot_table(index="horizon", columns="model", values="auc")
    print(piv.to_string(float_format=lambda x: f"{x:.4f}"))

    if all_lists:
        HL = pd.concat(all_lists, ignore_index=True)
        HL.to_csv(f"{args.out_prefix}_windows.csv", index=False)
        print("\n" + "=" * 80)
        print("CONSENSUS HOT LIST  (k=5, stability across folds)")
        print("=" * 80)
        for horizon in args.horizons:
            sub = HL[(HL.horizon == horizon) & (HL.k == 5)].head(8)
            if sub.empty:
                continue
            coords = windows_to_genomic_coordinates(sub, ar, n_bins=args.n_bins)
            print(f"\n-- {horizon} years --")
            cols = ["window", "chrom", "genomic_start_mb", "genomic_end_mb",
                    "selection_frequency"]
            cols = [c for c in cols if c in coords.columns]
            print(coords[cols].to_string(index=False))

    print(f"\nsaved {args.out_prefix}_performance.csv / _windows.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
