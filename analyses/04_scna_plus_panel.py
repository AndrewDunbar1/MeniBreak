#!/usr/bin/env python
"""Clinical-style evaluation of the 8-window breakpoint panel.

Reproduces the analyses used in the locked project so results are directly
comparable:

* **Discrimination** — AUC for WHO grade, integrated grade and the panel.
* **Matched-sensitivity operating point** — each model gets its own high-risk
  threshold chosen to reach ~80% sensitivity among observed events, then
  specificity is compared at that common sensitivity. Comparing specificity at a
  shared sensitivity is the only way the numbers mean the same thing.
* **Elevated-grade non-recurrence specificity** — among tumours that stayed
  recurrence-free through the horizon *but* were WHO 2/3 or integrated grade
  2/3, the fraction the model places below its high-risk threshold. These are the
  patients conventional grading would over-treat.
* **Event capture** — fraction of observed recurrences falling in the
  highest-risk 10% / 20% / 33% of the cohort.

All model scores are honest out-of-fold predictions with patient-grouped folds,
and window selection happens inside the training folds.

Writes reports/clinical_metrics_*.csv.
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
from sklearn.metrics import roc_auc_score, roc_curve

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
log = logging.getLogger("clinical")


def oof(base, windows, y, groups, *, k, n_repeats, seed, n_splits=5):
    """Patient-grouped out-of-fold scores, averaged over repeats."""
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


def operating_point(y, score, target=TARGET_SENSITIVITY):
    """Threshold reaching >= target sensitivity, with specificity there."""
    fpr, tpr, thr = roc_curve(y, score)
    ok = np.where(tpr >= target)[0]
    i = ok[0] if len(ok) else int(np.argmax(tpr))
    return {"threshold": float(thr[i]), "sensitivity": float(tpr[i]),
            "specificity": float(1 - fpr[i])}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=8)
    ap.add_argument("--n-bins", type=int, default=20)
    ap.add_argument("--repeats", type=int, default=10)
    ap.add_argument("--horizons", nargs="+", type=int, default=[3, 5, 7])
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

    disc, oper, elev, capture = [], [], [], []

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

        scores = {
            "WHO grade": pd.to_numeric(d["who_grade"], errors="coerce").to_numpy(float),
            "Integrated grade": pd.to_numeric(d["integrated_grade"], errors="coerce").to_numpy(float),
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
            for name, sc in scores.items():
                s = sc[mask]
                ok = np.isfinite(s)
                if ok.sum() < 20 or len(np.unique(yy[ok])) < 2:
                    continue
                auc = roc_auc_score(yy[ok], s[ok])
                disc.append({"horizon": horizon, "scope": scope, "model": name,
                             "n": int(ok.sum()), "events": int(yy[ok].sum()), "auc": auc})

                op = operating_point(yy[ok], s[ok])
                oper.append({"horizon": horizon, "scope": scope, "model": name, **op})

                # Event capture in the highest-risk tails.
                order = np.argsort(-s[ok])
                ys = yy[ok][order]
                for frac, tail in ((0.10, "Top 10%"), (0.20, "Top 20%"), (1/3, "Top 33%")):
                    n_tail = max(int(round(frac * len(ys))), 1)
                    capture.append({"horizon": horizon, "scope": scope, "model": name,
                                    "tail": tail,
                                    "capture_fraction": float(ys[:n_tail].sum() / max(ys.sum(), 1))})

                # Elevated-grade non-recurrence specificity, at the matched point.
                for gsys, col in (("WHO grade 2/3", "who_grade"),
                                  ("Integrated grade 2/3", "integrated_grade")):
                    g = pd.to_numeric(d[col], errors="coerce").to_numpy(float)[mask][ok]
                    elevated_nonrec = (yy[ok] == 0) & (g >= 2)
                    if elevated_nonrec.sum() < 5:
                        continue
                    below = (s[ok][elevated_nonrec] < op["threshold"]).sum()
                    elev.append({"horizon": horizon, "scope": scope, "model": name,
                                 "grade_system": gsys,
                                 "low_risk_n": int(below),
                                 "total": int(elevated_nonrec.sum()),
                                 "fraction_below_threshold": float(below / elevated_nonrec.sum()),
                                 "achieved_sensitivity": op["sensitivity"]})

    out = ROOT / "reports"
    pd.DataFrame(disc).to_csv(out / "clinical_metrics_discrimination.csv", index=False)
    pd.DataFrame(oper).to_csv(out / "clinical_metrics_operating.csv", index=False)
    pd.DataFrame(elev).to_csv(out / "clinical_metrics_elevated_grade.csv", index=False)
    pd.DataFrame(capture).to_csv(out / "clinical_metrics_capture.csv", index=False)

    D = pd.DataFrame(disc)
    print("\nDISCRIMINATION (AUC)")
    print(D.pivot_table(index=["scope", "horizon"], columns="model", values="auc")
          .to_string(float_format=lambda x: f"{x:.3f}"))
    O = pd.DataFrame(oper)
    print(f"\nSPECIFICITY at matched sensitivity ~{TARGET_SENSITIVITY:.0%}")
    print(O.pivot_table(index=["scope", "horizon"], columns="model", values="specificity")
          .to_string(float_format=lambda x: f"{x:.3f}"))
    print("\nsaved reports/clinical_metrics_*.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
