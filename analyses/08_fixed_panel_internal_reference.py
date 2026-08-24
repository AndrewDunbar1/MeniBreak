#!/usr/bin/env python
"""Internal reference AUC for the FIXED panel — the number an external cohort compares against.

WHY THIS EXISTS
---------------
The reported panel AUCs (0.702 / 0.735 / 0.752 broad, at 3/5/7 years) come from a
*procedure*: eight windows are re-selected by chi-square inside every training
fold and weighted by elastic net. `apply/validate_fixed_panel.py`, by contrast,
scores a *fixed* window list with an unweighted count and zero fitted
parameters. Those are different models, so the published numbers are NOT the
right yardstick for a fixed-panel validation result.

This script computes the missing yardstick: what the frozen panel in
`reference/menibreak_panel_v1.csv` achieves on the institutional cohort, scored
exactly the way an external cohort will score it. Run it once on the
institutional data and the resulting CSV becomes the internal reference that
`docs/THE_PANEL.md` currently reports as "not yet computed".

It deliberately mirrors `apply/validate_fixed_panel.py`: same fixed windows, same
`count` and `refit` scoring, same strict fixed-horizon endpoint, so the internal
and external numbers are directly comparable.

Requires the institutional dataset:
    export MENIBREAK_DATA_ROOT=/path/to/Breakpoint_Architecture_Paper_Analysis_Handoff/
    python analyses/07_fixed_panel_internal_reference.py

Writes results/reference/fixed_panel_internal_reference.csv
"""
from __future__ import annotations

import argparse
import logging
import os
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from apply.menibreak_panel import oof
from src.features.architecture import build_feature_table, feature_groups
from src.features.hotlist import build_window_matrix
from src.harmonization.breakpoints import derive_breakpoints, tile_arms
from src.harmonization.institutional import (
    fixed_horizon_labels, load_master_seg, load_outcomes,
)
from src.harmonization.references import ArmReference, load_frozen_hotspot_windows
from src.harmonization.schema import BoundaryClass

H = os.environ.get("MENIBREAK_DATA_ROOT",
                   "E:/Breakpoint_Architecture_Paper_Analysis_Handoff_2026-08-02/")
REF = H + "02_REFERENCE_AND_PROCESSED_DATA/"
SRC = H + "03_ANALYSIS_READY_DATA/ml_breakpoint_model_dataset_2026.csv"
log = logging.getLogger("fixed-ref")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--horizons", nargs="+", type=int, default=[3, 5, 7])
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--panel-set", choices=["published", "rule", "both"], default="both")
    ap.add_argument("--out", default="results/reference/fixed_panel_internal_reference.csv")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")

    ar = ArmReference.from_csv(REF + "arm_reference_GRCh37_software.csv", "GRCh37")
    panel_all = pd.read_csv(ROOT / "reference" / "menibreak_panel_v1.csv")

    seg = load_master_seg(H + "01_RAW_DATA/institutional/meningioma_masterseg_ALL_2026.txt",
                          forbidden_subpaths=["RESTRICTED"])
    tiles, _ = tile_arms(seg, ar)
    bps = derive_breakpoints(tiles, ar, boundary_class=BoundaryClass.SNP_ARRAY_CHANGEPOINT)
    m = seg[["sample_id", "tumor_id"]].drop_duplicates()
    bps = bps.merge(m, on="sample_id", how="left")
    tiles["unit"] = tiles["tumor_id"]
    bps["unit"] = bps["tumor_id"]

    # No prevalence filter: a fixed panel window must be scored even when rare,
    # exactly as apply/validate_fixed_panel.py does for an external cohort.
    all_w = build_window_matrix(bps, sample_col="unit", n_bins=20,
                                binary=True, min_prevalence=0)

    base = build_feature_table(tiles, bps, ar, sample_col="unit",
                               hotspot_windows=load_frozen_hotspot_windows(
                                   REF + "figure4_hotspot_windows_2026.csv"))
    scna_cols = feature_groups(base.columns)["scna"]
    outcomes = load_outcomes(H + "03_ANALYSIS_READY_DATA/outcomes.csv")

    sets = ["published", "rule"] if args.panel_set == "both" else [args.panel_set]
    rows = []
    for pset in sets:
        panel = panel_all[panel_all["in_published_table"]] if pset == "published" else panel_all
        cols = list(panel["window"])
        W_all = pd.DataFrame(0.0, index=all_w.index, columns=cols)
        for c in cols:
            if c in all_w.columns:
                W_all[c] = all_w[c].astype(float)
        log.info("%s panel: %d windows, %d present in institutional matrix",
                 pset, len(cols), sum(c in all_w.columns for c in cols))

        for horizon in args.horizons:
            lab = fixed_horizon_labels(outcomes, horizon)
            d = (lab.merge(W_all, left_on="tumor_id", right_index=True, how="inner")
                    .merge(base[scna_cols], left_on="tumor_id", right_index=True, how="left"))
            d = d.fillna({c: 0.0 for c in cols + list(scna_cols)})
            if d.empty or d["label"].nunique() < 2:
                continue
            y = d["label"].to_numpy(int)
            groups = d[["patient_id"]].reset_index(drop=True)
            P = d[cols].reset_index(drop=True)
            S = d[scna_cols].reset_index(drop=True)

            scores = {
                "Fixed panel (count, 0 fitted params)": P.sum(axis=1).to_numpy(float),
                "Fixed panel (refit coefficients)": oof(
                    pd.DataFrame({"_intercept": np.ones(len(d))}), P, y, groups,
                    k=P.shape[1], n_repeats=args.repeats, seed=11),
                "Arm-level SCNA (reference)": oof(
                    S, None, y, groups, k=0, n_repeats=args.repeats, seed=11),
            }
            for name, sc in scores.items():
                ok = np.isfinite(sc)
                if ok.sum() < 20 or len(np.unique(y[ok])) < 2:
                    continue
                auc = roc_auc_score(y[ok], sc[ok])
                rows.append({"horizon": horizon, "panel_set": pset, "model": name,
                             "n": int(ok.sum()), "events": int(y[ok].sum()),
                             "auc": round(float(auc), 4)})
                log.info("%dy %-10s %-38s AUC=%.4f", horizon, pset, name, auc)

    res = pd.DataFrame(rows)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out, index=False)
    print("\nFIXED-PANEL INTERNAL REFERENCE")
    if not res.empty:
        print(res.pivot_table(index=["panel_set", "horizon"], columns="model", values="auc")
              .to_string(float_format=lambda x: f"{x:.3f}"))
    print(f"\nsaved {out}")
    print("\nThis is the number an external cohort's fixed-panel result should be "
          "compared against — NOT the 0.702/0.735/0.752 dynamic-selection figures.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
