#!/usr/bin/env python
"""External validation of the breakpoint hot list on GSE189521.

This is a genuine external test, and the design is deliberately strict:

* the **8 windows are fixed coordinates** chosen from institutional data only;
* the **model is fitted on the full institutional cohort** and applied frozen;
* **nothing is refitted, rescaled or re-thresholded** on GSE189521;
* GSE189521 has been quarantined for the whole project and never used for
  pretraining, feature discovery, hyperparameter choice or thresholding.

It is also a **cross-platform** test, which makes it harder than an internal
split: institutional data is event-called aCGH, GSE189521 is EPIC
methylation-derived CNA via conumee. Methylation change points are the weakest
evidence class in this project's schema, and the segment count per sample is far
lower, so the window-occupancy features are systematically sparser. A drop
relative to internal performance is expected and is not by itself evidence
against the panel.

Reference: the locked SEG-only model scored 0.763 / 0.704 / 0.722 on this cohort
at 3 / 5 / 7 years.

Usage:
    python scripts/external_validate_hotlist.py --k 8
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

from src.features.architecture import build_feature_table, feature_groups
from src.features.hotlist import build_window_matrix, select_in_fold
from src.harmonization.breakpoints import derive_breakpoints, tile_arms
from src.harmonization.institutional import (
    fixed_horizon_labels, load_master_seg, load_outcomes,
)
from src.harmonization.references import ArmReference
from src.harmonization.schema import SEGMENT_SPEC, BoundaryClass, EventClass, normalise_chromosome

H = os.environ.get("MENIBREAK_DATA_ROOT", "E:/Breakpoint_Architecture_Paper_Analysis_Handoff_2026-08-02/")
REF = H + "02_REFERENCE_AND_PROCESSED_DATA/"
EXT_SEG = REF + "GSE189521_conumee_segments_all.tsv.gz"
EXT_CLIN = H + "01_RAW_DATA/GSE189521/GSE189521_Clinical_data.xlsx"

#: The locked model's prespecified event threshold for this cohort.
EVENT_THRESHOLD = 0.15

#: Locked SEG-only model performance on GSE189521, for reference.
LOCKED_EXTERNAL = {3: 0.763, 5: 0.704, 7: 0.722}

log = logging.getLogger("external")


def load_gse189521_segments(path: str = EXT_SEG) -> pd.DataFrame:
    """conumee segments -> canonical schema, using the locked 0.15 threshold."""
    raw = pd.read_csv(path, sep="\t")
    raw = raw[raw["chrom"].astype(str).str.replace("chr", "", case=False)
              .isin([str(i) for i in range(1, 23)])]

    seg = pd.DataFrame({
        "dataset_id": "gse189521",
        "sample_id": raw["sample_id"].astype("string"),
        "patient_id": raw["sample_id"].astype("string"),
        "tumor_id": raw["sample_id"].astype("string"),
        "aliquot_id": pd.Series(pd.NA, index=raw.index, dtype="string"),
        "chromosome": raw["chrom"].map(normalise_chromosome).astype("string"),
        "start": raw["loc.start"].astype("int64"),
        "end": raw["loc.end"].astype("int64"),
        "genome_build": "GRCh37",
        "segment_mean": raw["seg.mean"].astype("float64"),
        "total_copy_number": np.nan,
        "minor_copy_number": np.nan,
        "major_copy_number": np.nan,
        "probe_count": raw["num.mark"].astype("float64"),
        "platform": "methylation_epic",
        "caller": "conumee",
        "source_file": Path(path).name,
        "source_coordinate_build": "GRCh37",
        "qc_status": "pass",
    })

    m = seg["segment_mean"]
    event = pd.Series(EventClass.NEUTRAL.value, index=seg.index, dtype="object")
    event[m >= EVENT_THRESHOLD] = EventClass.GAIN.value
    event[m <= -EVENT_THRESHOLD] = EventClass.LOSS.value
    seg["event"] = event.astype("string")

    seg["segment_id"] = (
        seg["sample_id"].astype(str) + ":" + seg["chromosome"].astype(str) + ":"
        + seg["start"].astype(str) + "-" + seg["end"].astype(str) + ":cn"
    ).astype("string")

    SEGMENT_SPEC.validate(seg, strict=True)
    return seg.assign(track="cn")


def load_gse189521_outcomes(path: str = EXT_CLIN) -> pd.DataFrame:
    """Recurrence + follow-up, converted to the project's outcome schema."""
    d = pd.read_excel(path, sheet_name=0, header=2)
    d = d[d["Sample"].notna()]

    rec = pd.to_numeric(d["Recurrence"], errors="coerce")
    if rec.isna().all():
        rec = d["Recurrence"].astype(str).str.strip().str.lower().map(
            {"yes": 1, "no": 0, "y": 1, "n": 0, "1": 1, "0": 0})
    months = pd.to_numeric(d["MonthsFollowUp"], errors="coerce")

    out = pd.DataFrame({
        "dataset_id": "gse189521",
        "tumor_id": d["Sample"].astype(str).str.strip(),
        "patient_id": d["Sample"].astype(str).str.strip(),
        "recurrence_event": rec.astype("float64"),
        "recurrence_time": (months / 12.0).astype("float64"),
        "who_grade": pd.to_numeric(d["WHO grade"], errors="coerce"),
        "provenance": Path(path).name,
    })
    return out.dropna(subset=["recurrence_time", "recurrence_event"])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=8)
    ap.add_argument("--n-bins", type=int, default=20)
    ap.add_argument("--horizons", nargs="+", type=int, default=[3, 5, 7])
    ap.add_argument("--min-prevalence", type=int, default=20)
    ap.add_argument("--out", default="reports/external_validation.csv")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        datefmt="%H:%M:%S")

    ar = ArmReference.from_csv(REF + "arm_reference_GRCh37_software.csv", "GRCh37")

    # ---------------- institutional (training) --------------------------
    seg = load_master_seg(
        H + "01_RAW_DATA/institutional/meningioma_masterseg_ALL_2026.txt",
        forbidden_subpaths=["RESTRICTED"])
    tiles, _ = tile_arms(seg, ar)
    bps = derive_breakpoints(tiles, ar, boundary_class=BoundaryClass.SNP_ARRAY_CHANGEPOINT)
    m = seg[["sample_id", "tumor_id"]].drop_duplicates()
    bps = bps.merge(m, on="sample_id", how="left")
    tiles["unit"] = tiles["tumor_id"]
    bps["unit"] = bps["tumor_id"]

    base = build_feature_table(tiles, bps, ar, sample_col="unit")
    scna_cols = feature_groups(base.columns)["scna"]
    win_int = build_window_matrix(bps, sample_col="unit", n_bins=args.n_bins,
                                  binary=True, min_prevalence=args.min_prevalence)
    outcomes = load_outcomes(H + "03_ANALYSIS_READY_DATA/outcomes.csv")

    # ---------------- external (held out all project) -------------------
    ext_seg = load_gse189521_segments()
    ext_tiles, _ = tile_arms(ext_seg, ar)
    ext_bps = derive_breakpoints(
        ext_tiles, ar, boundary_class=BoundaryClass.METHYLATION_ARRAY_CHANGEPOINT)
    ext_tiles["unit"] = ext_tiles["sample_id"]
    ext_bps["unit"] = ext_bps["sample_id"]

    ext_base = build_feature_table(ext_tiles, ext_bps, ar, sample_col="unit")
    ext_win = build_window_matrix(ext_bps, sample_col="unit", n_bins=args.n_bins,
                                  binary=True, min_prevalence=1)
    ext_out = load_gse189521_outcomes()

    log.info("institutional: %d tumours, %d windows | external: %d tumours, %d windows",
             len(base), win_int.shape[1], len(ext_base), ext_win.shape[1])
    log.info("median breakpoints/sample  institutional=%d  external=%d",
             int(bps.groupby("unit").size().median()),
             int(ext_bps.groupby("unit").size().median()))

    rows = []
    for horizon in args.horizons:
        # --- choose the panel from INSTITUTIONAL data only ---------------
        lab = fixed_horizon_labels(outcomes, horizon)
        d = lab.merge(base[scna_cols], left_on="tumor_id", right_index=True, how="inner")
        d = d.merge(win_int, left_on="tumor_id", right_index=True, how="inner")
        y_tr = d["label"].to_numpy(int)
        W_tr = d[[c for c in win_int.columns if c in d.columns]]
        panel = select_in_fold(W_tr, y_tr, k=args.k)
        log.info("%dy panel: %s", horizon, ", ".join(panel))

        S_tr = d[scna_cols]

        # --- external labels, same strict endpoint -----------------------
        e_lab = fixed_horizon_labels(ext_out, horizon)
        e = e_lab.merge(ext_base[scna_cols], left_on="tumor_id",
                        right_index=True, how="inner")
        if e.empty:
            log.warning("%dy: no external overlap", horizon)
            continue
        # Windows absent in the external matrix become 0 (not observed).
        ext_w = ext_win.reindex(columns=panel, fill_value=0.0)
        e = e.merge(ext_w, left_on="tumor_id", right_index=True, how="left")
        e[panel] = e[panel].fillna(0.0)
        y_te = e["label"].to_numpy(int)
        if len(np.unique(y_te)) < 2:
            log.warning("%dy: external has one class only (n=%d)", horizon, len(e))
            continue

        for name, cols_tr, cols_te in (
            ("SCNA", S_tr, e[scna_cols]),
            ("windows_only", W_tr[panel], e[panel]),
            (f"SCNA+{args.k}windows", pd.concat([S_tr, W_tr[panel]], axis=1),
             pd.concat([e[scna_cols], e[panel]], axis=1)),
        ):
            a = cols_tr.to_numpy(float)
            b = cols_te.to_numpy(float)
            mu, sd = a.mean(0), a.std(0)
            sd = np.where(sd > 0, sd, 1.0)
            model = LogisticRegression(penalty="l2", C=1.0, max_iter=3000)
            model.fit((a - mu) / sd, y_tr)          # institutional only
            score = model.predict_proba((b - mu) / sd)[:, 1]   # frozen transfer
            rows.append({"horizon": horizon, "model": name, "n_external": len(e),
                         "events": int(y_te.sum()),
                         "auc": roc_auc_score(y_te, score)})

        who = pd.to_numeric(e["who_grade"], errors="coerce")
        ok = who.notna()
        if ok.sum() > 20 and len(np.unique(y_te[ok])) > 1:
            rows.append({"horizon": horizon, "model": "WHO_grade",
                         "n_external": int(ok.sum()),
                         "events": int(y_te[ok].sum()),
                         "auc": roc_auc_score(y_te[ok], who[ok])})
        rows.append({"horizon": horizon, "model": "locked_SEG_only_reference",
                     "n_external": len(e), "events": int(y_te.sum()),
                     "auc": LOCKED_EXTERNAL[horizon]})

    res = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(args.out, index=False)

    print("\n" + "=" * 74)
    print("EXTERNAL VALIDATION — GSE189521 (quarantined; frozen transfer)")
    print("Trained on institutional data only; nothing refitted externally")
    print("=" * 74)
    print(res.pivot_table(index="horizon", columns="model", values="auc")
          .to_string(float_format=lambda x: f"{x:.4f}"))
    print("\nsample sizes:")
    print(res.groupby("horizon")[["n_external", "events"]].max().to_string())
    print(f"\nsaved {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
