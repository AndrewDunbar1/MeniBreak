#!/usr/bin/env python
"""Reclassification of grade-discordant tumours by the breakpoint panel.

The clinically interesting question is not overall sensitivity or specificity —
it is what happens to the tumours **grading gets wrong**. Two failure modes:

* **Grade over-called** — WHO or integrated grade 2/3, but the tumour never
  recurred through the horizon. These patients are candidates for over-treatment.
  A useful model should place them BELOW its high-risk threshold.

* **Grade under-called** — WHO or integrated grade 1, but the tumour recurred
  within the horizon. These are the misses that matter most. A useful model
  should place them ABOVE its high-risk threshold.

For each model the high-risk threshold is set on the whole cohort to reach
~80% sensitivity, exactly as in the locked project's elevated-grade analysis, so
every model is compared at the same operating sensitivity. The two discordant
groups are then scored separately: fraction correctly reclassified in each.

Model scores are honest out-of-fold predictions with patient-grouped folds, and
window selection happens inside the training folds.

Writes reports/grade_discordance.csv.
"""
from __future__ import annotations

import argparse
import os
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
from sklearn.metrics import roc_curve

from src.evaluation.splits import assert_no_patient_leakage, patient_grouped_folds
from src.features.architecture import build_feature_table, feature_groups
from src.features.clinical import load_restricted_flag
from src.features.hotlist import build_window_matrix, select_in_fold
from src.harmonization.breakpoints import derive_breakpoints, tile_arms
from src.harmonization.institutional import (
    fixed_horizon_labels, load_master_seg, load_outcomes,
)
from src.harmonization.references import ArmReference
from src.harmonization.schema import BoundaryClass

H = os.environ.get("MENIBREAK_DATA_ROOT", "E:/Breakpoint_Architecture_Paper_Analysis_Handoff_2026-08-02/")
REF = H + "02_REFERENCE_AND_PROCESSED_DATA/"
SRC = H + "03_ANALYSIS_READY_DATA/ml_breakpoint_model_dataset_2026.csv"

TARGET_SENSITIVITY = 0.80
log = logging.getLogger("discord")


def oof(base, windows, y, groups, *, k, n_repeats, seed, n_splits=5):
    acc = np.full((len(y), n_repeats), np.nan)
    for rep in range(n_repeats):
        for tr, te in patient_grouped_folds(groups, n_splits=n_splits, seed=seed + rep):
            assert_no_patient_leakage(groups, tr, te)
            if len(np.unique(y[tr])) < 2:
                continue
            Xtr = base.iloc[tr].reset_index(drop=True)
            Xte = base.iloc[te].reset_index(drop=True)
            if windows is not None and k > 0:
                sel = select_in_fold(windows.iloc[tr], y[tr], k=k)
                if sel:
                    Xtr = pd.concat([Xtr, windows.iloc[tr][sel].reset_index(drop=True)], axis=1)
                    Xte = pd.concat([Xte, windows.iloc[te][sel].reset_index(drop=True)], axis=1)
            a, b = Xtr.to_numpy(float), Xte.to_numpy(float)
            mu = np.nanmean(a, axis=0)
            mu = np.where(np.isfinite(mu), mu, 0.0)
            a = np.where(np.isnan(a), mu, a)
            b = np.where(np.isnan(b), mu, b)
            sd = a.std(0)
            sd = np.where(sd > 0, sd, 1.0)
            m = LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.25,
                                   C=0.1, max_iter=3000).fit((a - mu) / sd, y[tr])
            acc[te, rep] = m.predict_proba((b - mu) / sd)[:, 1]
    return np.nanmean(acc, axis=1)


def high_risk_threshold(y, score, target=TARGET_SENSITIVITY) -> float:
    fpr, tpr, thr = roc_curve(y, score)
    ok = np.where(tpr >= target)[0]
    return float(thr[ok[0]]) if len(ok) else float(thr[int(np.argmax(tpr))])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=8)
    ap.add_argument("--n-bins", type=int, default=20)
    ap.add_argument("--repeats", type=int, default=10)
    ap.add_argument("--horizons", nargs="+", type=int, default=[3, 5, 7])
    ap.add_argument("--out", default="reports/grade_discordance.csv")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        datefmt="%H:%M:%S")

    ar = ArmReference.from_csv(REF + "arm_reference_GRCh37_software.csv", "GRCh37")
    seg = load_master_seg(H + "01_RAW_DATA/institutional/meningioma_masterseg_ALL_2026.txt",
                          forbidden_subpaths=["RESTRICTED"])
    tiles, _ = tile_arms(seg, ar)
    bps = derive_breakpoints(tiles, ar, boundary_class=BoundaryClass.SNP_ARRAY_CHANGEPOINT)
    m = seg[["sample_id", "tumor_id"]].drop_duplicates()
    bps = bps.merge(m, on="sample_id", how="left")
    tiles["unit"] = tiles["tumor_id"]
    bps["unit"] = bps["tumor_id"]

    base = build_feature_table(tiles, bps, ar, sample_col="unit")
    scna_cols = feature_groups(base.columns)["scna"]
    windows = build_window_matrix(bps, sample_col="unit", n_bins=args.n_bins,
                                  binary=True, min_prevalence=20)
    grades = pd.read_csv(SRC, usecols=lambda c: c in ("cg", "who_grade", "integrated_grade"),
                         low_memory=False).rename(columns={"cg": "tumor_id"})
    grades["tumor_id"] = grades["tumor_id"].astype(str)
    restricted = load_restricted_flag(SRC)
    outcomes = load_outcomes(H + "03_ANALYSIS_READY_DATA/outcomes.csv")

    rows = []
    for horizon in args.horizons:
        lab = fixed_horizon_labels(outcomes, horizon)
        d = (lab.merge(base[scna_cols], left_on="tumor_id", right_index=True, how="inner")
                .merge(windows, left_on="tumor_id", right_index=True, how="inner")
                .merge(grades, on="tumor_id", how="left")
                .merge(restricted, left_on="tumor_id", right_index=True, how="left"))
        d["is_restricted"] = d["is_restricted"].fillna(False)

        y = d["label"].to_numpy(int)
        groups = d[["patient_id"]].reset_index(drop=True)
        S = d[scna_cols].reset_index(drop=True)
        W = d[[c for c in windows.columns if c in d.columns]].reset_index(drop=True)
        blank = pd.DataFrame({"_intercept": np.ones(len(d))})

        models = {
            "Breakpoint panel (8 windows)": oof(blank, W, y, groups, k=args.k,
                                                n_repeats=args.repeats, seed=11),
            "Arm-level SCNA": oof(S, None, y, groups, k=0,
                                  n_repeats=args.repeats, seed=11),
            "SCNA + breakpoint panel": oof(S, W, y, groups, k=args.k,
                                           n_repeats=args.repeats, seed=11),
        }

        for scope, mask in (("Broad", np.ones(len(d), bool)),
                            ("Restricted", d["is_restricted"].to_numpy(bool))):
            yy = y[mask]
            if len(np.unique(yy)) < 2:
                continue
            for model_name, sc in models.items():
                s = sc[mask]
                ok = np.isfinite(s)
                if ok.sum() < 20:
                    continue
                thr = high_risk_threshold(yy[ok], s[ok])
                high_risk = s[ok] >= thr

                for gsys, col in (("WHO grade", "who_grade"),
                                  ("Integrated grade", "integrated_grade")):
                    g = pd.to_numeric(d[col], errors="coerce").to_numpy(float)[mask][ok]
                    ev = yy[ok]

                    # Grade said high risk, tumour never recurred.
                    over = (g >= 2) & (ev == 0)
                    # Grade said low risk, tumour recurred.
                    under = (g == 1) & (ev == 1)

                    for direction, sel, correct in (
                        ("Grade over-called\n(grade 2/3, no recurrence)", over, ~high_risk),
                        ("Grade under-called\n(grade 1, recurred)", under, high_risk),
                    ):
                        if sel.sum() < 5:
                            continue
                        n_ok = int((sel & correct).sum())
                        rows.append({
                            "horizon": horizon, "scope": scope, "model": model_name,
                            "grade_system": gsys, "direction": direction.split("\n")[0],
                            "direction_full": direction.replace("\n", " "),
                            "n_discordant": int(sel.sum()),
                            "n_correct": n_ok,
                            "fraction_correct": n_ok / int(sel.sum()),
                            "threshold": thr,
                        })

    res = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(args.out, index=False)

    print("\nRECLASSIFICATION OF GRADE-DISCORDANT TUMOURS")
    print("(high-risk threshold set at ~80% sensitivity for every model)\n")
    for gsys in ("WHO grade", "Integrated grade"):
        print(f"--- vs {gsys} ---")
        sub = res[res.grade_system == gsys]
        if sub.empty:
            continue
        piv = sub.pivot_table(index=["scope", "horizon", "direction"],
                              columns="model", values="fraction_correct")
        cnt = sub.pivot_table(index=["scope", "horizon", "direction"],
                              values="n_discordant", aggfunc="max")
        print(piv.join(cnt).to_string(float_format=lambda x: f"{x:.2f}"))
        print()
    print(f"saved {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
