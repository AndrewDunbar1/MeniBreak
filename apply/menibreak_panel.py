#!/usr/bin/env python
"""Run the MeniBreak breakpoint panel on YOUR OWN segmented copy-number data.

This is the entry point for an external cohort. It takes a generic SEG file plus
recurrence outcomes and reproduces the core analysis of the MeniBreak report on
your data:

  1. adapt your SEG to the canonical schema,
  2. tile chromosome arms and derive breakpoints (change points),
  3. build the 20-bin arm-relative window matrix,
  4. discover the top-k prognostic windows INSIDE each training fold,
  5. score them with patient-grouped, repeated out-of-fold cross-validation,
  6. report the panel-alone / SCNA / SCNA+panel AUC ladder and the consensus
     window list translated to genomic coordinates.

Every modelling choice (elastic-net l1_ratio=0.25, C=0.1; in-fold selection;
train-only imputation and standardisation; patient-grouped folds) is identical
to analyses/03_full_ladder.py, which produced the institutional numbers.

INPUT (see docs/INPUT_SPEC.md)
  --segments   CSV: sample_id, chromosome, start, end, event   (+ optional patient_id)
  --outcomes   CSV: sample_id, patient_id, recurrence_time, recurrence_event
  --build      GRCh37 (default) or GRCh38

Example (bundled synthetic cohort):
  python apply/menibreak_panel.py \
      --segments data/example/synthetic_segments.csv \
      --outcomes data/example/synthetic_outcomes.csv \
      --out results/example
"""
from __future__ import annotations

import argparse
import logging
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")  # elastic-net convergence / API-deprecation noise

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
from src.harmonization.references import ArmReference, load_frozen_hotspot_windows
from src.harmonization.schema import (
    COPY_NUMBER_EVENTS, BoundaryClass, EventClass, normalise_chromosome,
)

log = logging.getLogger("menibreak")

#: Permissive mapping of common caller vocabularies to the canonical EventClass.
EVENT_SYNONYMS = {
    "loss": EventClass.LOSS, "del": EventClass.LOSS, "deletion": EventClass.LOSS,
    "cn loss": EventClass.LOSS, "-1": EventClass.LOSS, "loh_loss": EventClass.LOSS,
    "gain": EventClass.GAIN, "amp": EventClass.GAIN, "amplification": EventClass.GAIN,
    "cn gain": EventClass.GAIN, "1": EventClass.GAIN, "dup": EventClass.GAIN,
    "homozygous_loss": EventClass.HOMOZYGOUS_LOSS, "homdel": EventClass.HOMOZYGOUS_LOSS,
    "homozygous copy loss": EventClass.HOMOZYGOUS_LOSS, "-2": EventClass.HOMOZYGOUS_LOSS,
    "high_copy_gain": EventClass.HIGH_COPY_GAIN, "high copy gain": EventClass.HIGH_COPY_GAIN,
    "2": EventClass.HIGH_COPY_GAIN,
    "loh": EventClass.LOH, "cnn-loh": EventClass.LOH, "cnloh": EventClass.LOH,
    "allelic_imbalance": EventClass.ALLELIC_IMBALANCE, "ai": EventClass.ALLELIC_IMBALANCE,
    "neutral": EventClass.NEUTRAL, "0": EventClass.NEUTRAL,
}


def _map_event(v: object) -> EventClass:
    s = str(v).strip().lower()
    if s not in EVENT_SYNONYMS:
        raise ValueError(
            f"unrecognised event label {v!r}. Map it in EVENT_SYNONYMS or "
            f"recode your SEG to one of: loss, gain, homozygous_loss, "
            f"high_copy_gain, loh, allelic_imbalance, neutral."
        )
    return EVENT_SYNONYMS[s]


def adapt_segments(path: str | Path, build: str = "GRCh37") -> pd.DataFrame:
    """Generic SEG -> canonical segment frame (adds `track`, `segment_id`)."""
    raw = pd.read_csv(path)
    need = {"sample_id", "chromosome", "start", "end", "event"}
    missing = need - set(raw.columns)
    if missing:
        raise ValueError(f"{path}: segments missing required columns {sorted(missing)}")
    ev = raw["event"].map(_map_event)
    track = ev.map(lambda e: "cn" if e in COPY_NUMBER_EVENTS else "allelic").astype("string")
    sid = raw["sample_id"].astype("string")
    out = pd.DataFrame({
        "dataset_id": "external_cohort",
        "sample_id": sid,
        # tumour-level key; defaults to sample_id when no separate tumour id given
        "tumor_id": (raw["tumor_id"].astype("string") if "tumor_id" in raw.columns else sid),
        "patient_id": (raw["patient_id"].astype("string")
                       if "patient_id" in raw.columns else sid),
        "chromosome": raw["chromosome"].map(normalise_chromosome).astype("string"),
        "start": raw["start"].astype("int64"),
        "end": raw["end"].astype("int64"),
        "event": ev.map(lambda e: e.value).astype("string"),
        # optional continuous columns: used if present, otherwise unknown
        "segment_mean": (raw["segment_mean"].astype("float64")
                         if "segment_mean" in raw.columns else np.nan),
        "probe_count": (raw["probe_count"].astype("float64")
                        if "probe_count" in raw.columns else np.nan),
        "platform": (raw["platform"].astype("string") if "platform" in raw.columns else "unknown"),
        "caller": (raw["caller"].astype("string") if "caller" in raw.columns else "unknown"),
        "genome_build": raw["genome_build"] if "genome_build" in raw.columns else build,
        "source_file": Path(path).name,
        "track": track,
    })
    out = out[out["event"] != EventClass.NEUTRAL.value].reset_index(drop=True)
    out["segment_id"] = [f"S{i}" for i in range(len(out))]
    return out


# --- verbatim modelling recipe from analyses/03_full_ladder.py --------------
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


def per_fold_consensus(windows, y, groups, *, k, n_repeats, seed, n_splits=5):
    """Which windows recur across training folds (stability = a finding)."""
    selections = []
    for rep in range(n_repeats):
        for tr, _ in patient_grouped_folds(groups, n_splits=n_splits, seed=seed + rep):
            if len(np.unique(y[tr])) < 2:
                continue
            selections.append(select_in_fold(windows.iloc[tr], y[tr], k=k))
    return consensus_hotlist(selections, n_bins=20)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--segments", required=True)
    ap.add_argument("--outcomes", required=True)
    ap.add_argument("--build", default="GRCh37", choices=["GRCh37", "GRCh38"])
    ap.add_argument("--horizons", nargs="+", type=int, default=[3, 5, 7])
    ap.add_argument("--k", type=int, default=8, help="windows selected per fold")
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--min-prevalence", type=int, default=20,
                    help="drop windows present in fewer than this many tumours")
    ap.add_argument("--out", default="results/external_run")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
    outdir = Path(args.out); outdir.mkdir(parents=True, exist_ok=True)

    ref_dir = ROOT / "reference"
    arm_csv = ref_dir / f"arm_reference_{args.build}.csv"
    ar = ArmReference.from_csv(arm_csv, args.build)
    hotspots = load_frozen_hotspot_windows(ref_dir / "hotspot_windows_v1.csv")

    seg = adapt_segments(args.segments, build=args.build)
    log.info("segments: %d rows, %d tumours", len(seg), seg["sample_id"].nunique())
    tiles, _ = tile_arms(seg, ar)
    bps = derive_breakpoints(tiles, ar, boundary_class=BoundaryClass.SNP_ARRAY_CHANGEPOINT)
    log.info("derived %d breakpoints", len(bps))

    # arm-level SCNA feature table + the 20-bin window matrix
    base_feat = build_feature_table(tiles, bps, ar, sample_col="sample_id",
                                    hotspot_windows=hotspots)
    scna_cols = feature_groups(base_feat.columns)["scna"]
    windows = build_window_matrix(bps, sample_col="sample_id", n_bins=20,
                                  binary=True, min_prevalence=args.min_prevalence)
    windows.to_csv(outdir / "window_features.csv")
    log.info("window matrix: %d tumours x %d windows", *windows.shape)

    # Resolution / comparability QC (see docs/INPUT_SPEC.md). The institutional
    # platform detects a median of ~97 breakpoints/sample; a much lower number
    # means a coarser assay on which the ~5 Mb windows are not comparable.
    bp_per_sample = bps.groupby("sample_id").size()
    occupancy = float(windows.to_numpy().mean()) if windows.size else float("nan")
    log.info("QC: median breakpoints/sample = %.0f  |  mean window occupancy = %.1f%%",
             float(bp_per_sample.median()), 100 * occupancy)
    if bp_per_sample.median() < 40:
        log.warning("QC: low breakpoint density (<40/sample). The panel was defined on "
                    "~97/sample data; AUCs here may not be comparable (see INPUT_SPEC.md).")

    # outcomes -> per-horizon labels (strict fixed-horizon endpoint)
    oc = pd.read_csv(args.outcomes)
    need = {"sample_id", "patient_id", "recurrence_time", "recurrence_event"}
    if need - set(oc.columns):
        raise ValueError(f"outcomes missing {sorted(need - set(oc.columns))}")
    oc = oc.rename(columns={"sample_id": "tumor_id"})
    oc["tumor_id"] = oc["tumor_id"].astype(str)

    rows, all_windows = [], []
    for horizon in args.horizons:
        t, e = oc["recurrence_time"], oc["recurrence_event"]
        is_event = (e == 1) & (t <= horizon)
        eligible = is_event | (t >= horizon)
        lab = oc.loc[eligible].copy()
        lab["label"] = is_event.loc[eligible].astype(int)

        d = (lab.merge(windows, left_on="tumor_id", right_index=True, how="inner")
                .merge(base_feat[scna_cols], left_on="tumor_id", right_index=True, how="left"))
        d = d.fillna({c: 0.0 for c in list(windows.columns) + list(scna_cols)})
        if d.empty or d["label"].nunique() < 2:
            log.warning("%dy: too few eligible tumours or one class only; skipped", horizon)
            continue

        y = d["label"].to_numpy(int)
        groups = d[["patient_id"]].reset_index(drop=True)
        S = d[scna_cols].reset_index(drop=True)
        W = d[[c for c in windows.columns if c in d.columns]].reset_index(drop=True)
        intercept = pd.DataFrame({"_intercept": np.ones(len(d))})

        ladder = {
            "Arm-level SCNA": oof(S, None, y, groups, k=0, n_repeats=args.repeats, seed=11),
            "8-window panel alone": oof(intercept, W, y, groups, k=args.k,
                                        n_repeats=args.repeats, seed=11),
            "SCNA + 8-window panel": oof(S, W, y, groups, k=args.k,
                                         n_repeats=args.repeats, seed=11),
        }
        for name, sc in ladder.items():
            ok = np.isfinite(sc)
            if ok.sum() < 20 or len(np.unique(y[ok])) < 2:
                continue
            auc = roc_auc_score(y[ok], sc[ok])
            rows.append({"horizon": horizon, "n": int(ok.sum()),
                         "events": int(y[ok].sum()), "model": name, "auc": round(auc, 4)})
            log.info("%dy  %-24s AUC=%.3f  (n=%d, events=%d)",
                     horizon, name, auc, int(ok.sum()), int(y[ok].sum()))

        cons = per_fold_consensus(W, y, groups, k=args.k, n_repeats=args.repeats, seed=11)
        if not cons.empty:
            cons = windows_to_genomic_coordinates(cons, ar, n_bins=20)
            cons.insert(0, "horizon", horizon)
            all_windows.append(cons)

    perf = pd.DataFrame(rows)
    perf.to_csv(outdir / "performance.csv", index=False)
    if all_windows:
        pd.concat(all_windows, ignore_index=True).to_csv(outdir / "hotlist_windows.csv", index=False)

    print("\n=== PERFORMANCE (patient-grouped out-of-fold AUC) ===")
    if not perf.empty:
        print(perf.pivot_table(index="horizon", columns="model", values="auc")
              .to_string(float_format=lambda x: f"{x:.3f}"))
    print(f"\nwrote: {outdir/'performance.csv'}, {outdir/'window_features.csv'}, "
          f"{outdir/'hotlist_windows.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
