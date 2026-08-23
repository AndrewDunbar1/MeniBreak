#!/usr/bin/env python
"""Full model ladder in BOTH cohorts, from WHO grade to the enriched model.

Existing outputs cover only part of this. The maximization search
(reports/maximize_search.csv) ran the enriched feature set but broad cohort only;
the hot-list transport ran both cohorts but only the 8-window panel. This script
computes every rung in both, on identical folds, so the numbers are comparable.

Ladder:
    1  WHO grade
    2  Integrated grade
    3  Arm-level SCNA
    4  SCNA + 8-window breakpoint panel
    5  SCNA + enriched breakpoints   (amplitude, allelic, density, architecture)
    6  Clinical + SCNA + enriched breakpoints

Restricted-cohort figures use the locked project's transport design: the model is
fitted on the broad cohort with patient-grouped out-of-fold prediction, and those
honest scores are then subset to the restricted tumours. Nothing is refitted
inside the restricted subset, which has only 22-43 events.

Writes reports/full_ladder.csv.
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
from src.features.clinical import CLINICAL_FEATURES, load_clinical_covariates, load_restricted_flag
from src.features.enriched import build_enriched
from src.features.hotlist import build_window_matrix, select_in_fold
from src.harmonization.breakpoints import derive_breakpoints, tile_arms
from src.harmonization.institutional import (
    fixed_horizon_labels, load_master_seg, load_outcomes,
)
from src.harmonization.references import ArmReference, load_frozen_hotspot_windows
from src.harmonization.schema import BoundaryClass

H = os.environ.get("MENIBREAK_DATA_ROOT", "E:/Breakpoint_Architecture_Paper_Analysis_Handoff_2026-08-02/")
REF = H + "02_REFERENCE_AND_PROCESSED_DATA/"
SRC = H + "03_ANALYSIS_READY_DATA/ml_breakpoint_model_dataset_2026.csv"
log = logging.getLogger("ladder")


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
            # Clinical covariates have genuine missingness (age, mitotic count,
            # resection and radiotherapy status). Impute from TRAINING rows only,
            # as the locked model does -- imputing from the full cohort would
            # leak test-fold information into the training statistics.
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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=8)
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--horizons", nargs="+", type=int, default=[3, 5, 7])
    ap.add_argument("--out", default="reports/full_ladder.csv")
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

    win_frozen = load_frozen_hotspot_windows(REF + "figure4_hotspot_windows_2026.csv")
    base = build_feature_table(tiles, bps, ar, sample_col="unit",
                               hotspot_windows=win_frozen)
    g = feature_groups(base.columns)
    scna_cols = g["scna"]
    arch_hot = [c for c in base.columns if c not in scna_cols]
    enriched = build_enriched(tiles, bps, sample_col="unit", include_bins=False)
    windows = build_window_matrix(bps, sample_col="unit", n_bins=20,
                                  binary=True, min_prevalence=20)
    clin = load_clinical_covariates(SRC)
    grades = pd.read_csv(SRC, usecols=lambda c: c in ("cg", "who_grade", "integrated_grade"),
                         low_memory=False).rename(columns={"cg": "tumor_id"})
    grades["tumor_id"] = grades["tumor_id"].astype(str)
    restricted = load_restricted_flag(SRC)
    outcomes = load_outcomes(H + "03_ANALYSIS_READY_DATA/outcomes.csv")

    log.info("scna=%d arch+hot=%d enriched=%d windows=%d",
             len(scna_cols), len(arch_hot), enriched.shape[1], windows.shape[1])

    rows = []
    for horizon in args.horizons:
        lab = fixed_horizon_labels(outcomes, horizon)
        d = (lab.merge(base, left_on="tumor_id", right_index=True, how="inner")
                .merge(enriched, left_on="tumor_id", right_index=True, how="left")
                .merge(windows, left_on="tumor_id", right_index=True, how="left")
                .merge(clin, left_on="tumor_id", right_index=True, how="left")
                .merge(grades, on="tumor_id", how="left")
                .merge(restricted, left_on="tumor_id", right_index=True, how="left"))
        d["is_restricted"] = d["is_restricted"].fillna(False)
        d = d.fillna({c: 0.0 for c in list(enriched.columns) + list(windows.columns)})

        y = d["label"].to_numpy(int)
        groups = d[["patient_id"]].reset_index(drop=True)
        S = d[scna_cols].reset_index(drop=True)
        AH = d[[c for c in arch_hot if c in d.columns]].reset_index(drop=True)
        E = d[[c for c in enriched.columns if c in d.columns]].reset_index(drop=True)
        W = d[[c for c in windows.columns if c in d.columns]].reset_index(drop=True)
        C = d[CLINICAL_FEATURES].reset_index(drop=True)

        scores = {
            "WHO grade": pd.to_numeric(d["who_grade"], errors="coerce").to_numpy(float),
            "Integrated grade": pd.to_numeric(d["integrated_grade"], errors="coerce").to_numpy(float),
            "Arm-level SCNA": oof(S, None, y, groups, k=0, n_repeats=args.repeats, seed=11),
            # The panel on its own: 8 binary indicators, no arm-level calls and
            # no clinical context. The floor of the parsimony curve.
            "8-window panel alone": oof(
                pd.DataFrame({"_intercept": np.ones(len(d))}), W, y, groups,
                k=args.k, n_repeats=args.repeats, seed=11),
            "SCNA + 8-window panel": oof(S, W, y, groups, k=args.k,
                                         n_repeats=args.repeats, seed=11),
            "SCNA + enriched breakpoints": oof(pd.concat([S, AH, E], axis=1), None, y,
                                               groups, k=0, n_repeats=args.repeats, seed=11),
            # Parsimonious panel plus clinical covariates: the deployable
            # feature count with the clinical context the locked model uses.
            "Clinical + 8-window panel": oof(C, W, y, groups, k=args.k,
                                             n_repeats=args.repeats, seed=11),
            "Clinical + SCNA + 8-window panel": oof(pd.concat([S, C], axis=1), W, y,
                                                    groups, k=args.k,
                                                    n_repeats=args.repeats, seed=11),
            "Clinical + SCNA + enriched breakpoints":
                oof(pd.concat([S, AH, E, C], axis=1), None, y, groups, k=0,
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
                rows.append({"horizon": horizon, "scope": scope, "model": name,
                             "n": int(ok.sum()), "events": int(yy[ok].sum()),
                             "auc": roc_auc_score(yy[ok], s[ok])})
                log.info("%dy %-11s %-40s %.4f", horizon, scope, name, rows[-1]["auc"])

    res = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(args.out, index=False)
    print("\nFULL LADDER")
    print(res.pivot_table(index=["scope", "horizon"], columns="model", values="auc")
          .to_string(float_format=lambda x: f"{x:.3f}"))
    print(f"\nsaved {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
