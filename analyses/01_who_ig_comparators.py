#!/usr/bin/env python
"""Compute WHO-grade and integrated-grade AUCs separately for each cohort.

Necessary because the locked project's published comparator AUCs
(WHO 0.690 / 0.656 / 0.635, integrated grade 0.766 / 0.735 / 0.733) are
**broad-cohort** figures. Re-using them on the restricted subset would imply
they were measured there, which they were not — and they are materially
different there: integrated grade *improves* in the cleaner cohort
(0.827 vs 0.766 at 3 years) because the confounding it competes against is gone.

Writes reports/comparators.csv.
"""
from __future__ import annotations
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

from src.features.clinical import load_restricted_flag
from src.harmonization.institutional import fixed_horizon_labels, load_outcomes

H = os.environ.get("MENIBREAK_DATA_ROOT", "E:/Breakpoint_Architecture_Paper_Analysis_Handoff_2026-08-02/")
SRC = H + "03_ANALYSIS_READY_DATA/ml_breakpoint_model_dataset_2026.csv"

#: Published broad-cohort values, retained for the reproduction check.
LOCKED_BROAD = {
    3: {"who_grade": 0.690, "integrated_grade": 0.766},
    5: {"who_grade": 0.656, "integrated_grade": 0.735},
    7: {"who_grade": 0.635, "integrated_grade": 0.733},
}


def main() -> int:
    grades = pd.read_csv(
        SRC, usecols=lambda c: c in ("cg", "who_grade", "integrated_grade"),
        low_memory=False,
    ).rename(columns={"cg": "tumor_id"})
    grades["tumor_id"] = grades["tumor_id"].astype(str)

    restricted = load_restricted_flag(SRC)
    outcomes = load_outcomes(H + "03_ANALYSIS_READY_DATA/outcomes.csv")

    rows = []
    for horizon in (3, 5, 7):
        lab = (fixed_horizon_labels(outcomes, horizon)
               .merge(grades, on="tumor_id", how="left")
               .merge(restricted, left_on="tumor_id", right_index=True, how="left"))
        lab["is_restricted"] = lab["is_restricted"].fillna(False)

        for scope, mask in (("broad", np.ones(len(lab), bool)),
                            ("restricted", lab["is_restricted"].to_numpy(bool))):
            sub = lab[mask]
            for marker in ("who_grade", "integrated_grade"):
                v = pd.to_numeric(sub[marker], errors="coerce")
                ok = v.notna()
                if ok.sum() < 20 or sub["label"][ok].nunique() < 2:
                    continue
                rows.append({
                    "horizon": horizon, "scope": scope, "marker": marker,
                    "n": int(ok.sum()), "events": int(sub["label"][ok].sum()),
                    "auc": roc_auc_score(sub["label"][ok], v[ok]),
                    "locked_published": (LOCKED_BROAD[horizon][marker]
                                         if scope == "broad" else np.nan),
                })

    res = pd.DataFrame(rows)
    out = ROOT / "reports" / "comparators.csv"
    res.to_csv(out, index=False)

    print(res.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    broad = res[res.scope == "broad"]
    print("\nreproduction check (broad vs published):")
    for r in broad.itertuples():
        print(f"  {r.horizon}y {r.marker:18} mine={r.auc:.4f} "
              f"published={r.locked_published:.3f} diff={r.auc - r.locked_published:+.4f}")
    print(f"\nsaved {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
