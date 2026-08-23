#!/usr/bin/env python
"""Validate the PUBLISHED MeniBreak panel on your cohort — fixed windows, no discovery.

This is the external-validation entry point. It answers:

    "Does the published MeniBreak panel discriminate recurrence in MY cohort?"

which is a different — and stricter — question than the one answered by
`apply/menibreak_panel.py`, which re-runs the *discovery procedure* and selects
windows from your own data. Here the windows are FIXED, read from
`reference/menibreak_panel_v1.csv`, and nothing about them is learned from your
cohort.

Two scoring modes, testing two different claims:

  --score count  (default)
      Risk score = number of panel windows containing a breakpoint.
      ZERO parameters are fitted on your data. This is a fully frozen marker and
      the cleanest external validation: the panel either transfers or it does not.

  --score refit
      The fixed windows are used as features and coefficients are refitted under
      patient-grouped, repeated cross-validation. This validates the *window set*
      while allowing effect sizes to differ between cohorts. It is a weaker claim
      than `count` but a fairer one if your platform or population differs.

Note there are no shipped coefficients for the panel, so a fully frozen
coefficient-level validation is not offered rather than being faked.

INPUT (see docs/INPUT_SPEC.md) — identical to apply/menibreak_panel.py.

Example:
  python apply/validate_fixed_panel.py \
      --segments data/example/synthetic_segments.csv \
      --outcomes data/example/synthetic_outcomes.csv \
      --out results/example_fixed
"""
from __future__ import annotations

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
from sklearn.metrics import roc_auc_score

from apply.menibreak_panel import adapt_segments, load_outcomes_table, oof
from src.features.architecture import build_feature_table, feature_groups
from src.features.hotlist import build_window_matrix
from src.harmonization.breakpoints import derive_breakpoints, tile_arms
from src.harmonization.references import ArmReference, load_frozen_hotspot_windows
from src.harmonization.schema import BoundaryClass

log = logging.getLogger("validate")


def load_panel(path: Path, panel_set: str) -> pd.DataFrame:
    """Read the frozen panel definition; `published` (9) or `rule` (11)."""
    p = pd.read_csv(path)
    if panel_set == "published":
        p = p[p["in_published_table"]]
    return p.reset_index(drop=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--segments", required=True)
    ap.add_argument("--outcomes", required=True)
    ap.add_argument("--build", default="GRCh37", choices=["GRCh37", "GRCh38"])
    ap.add_argument("--horizons", nargs="+", type=int, default=[3, 5, 7])
    ap.add_argument("--score", choices=["count", "refit"], default="count")
    ap.add_argument("--panel-set", choices=["published", "rule"], default="published",
                    help="published = the 9 windows in the report; rule = the 11 the "
                         "reproducible >=40%% rule admits (see reference/build_panel_definition.py)")
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--out", default="results/fixed_panel")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
    outdir = Path(args.out); outdir.mkdir(parents=True, exist_ok=True)

    ref_dir = ROOT / "reference"
    ar = ArmReference.from_csv(ref_dir / f"arm_reference_{args.build}.csv", args.build)
    panel = load_panel(ref_dir / "menibreak_panel_v1.csv", args.panel_set)
    log.info("panel: %s set, %d fixed windows (version %s)",
             args.panel_set, len(panel), panel["panel_version"].iloc[0])

    seg = adapt_segments(args.segments, build=args.build)
    tiles, _ = tile_arms(seg, ar)
    bps = derive_breakpoints(tiles, ar, boundary_class=BoundaryClass.SNP_ARRAY_CHANGEPOINT)
    log.info("segments: %d rows, %d tumours, %d breakpoints",
             len(seg), seg["sample_id"].nunique(), len(bps))

    # Window matrix with NO prevalence filter: a fixed panel window must be
    # scored even when it is rare in this cohort, otherwise the panel silently
    # changes size between cohorts.
    all_windows = build_window_matrix(bps, sample_col="sample_id", n_bins=20,
                                      binary=True, min_prevalence=0)
    tumours = sorted(seg["sample_id"].unique())
    W = pd.DataFrame(0.0, index=tumours, columns=list(panel["window"]))
    present = [w for w in panel["window"] if w in all_windows.columns]
    if present:
        W.loc[all_windows.index.intersection(tumours), present] = \
            all_windows.loc[all_windows.index.intersection(tumours), present].astype(float)

    # Transferability report: which fixed windows actually fire here?
    occ = W.mean(axis=0).rename("occupancy").reset_index().rename(columns={"index": "window"})
    occ = panel.merge(occ, on="window", how="left")
    occ["fires_in_cohort"] = occ["occupancy"] > 0
    occ.to_csv(outdir / "panel_window_occupancy.csv", index=False)
    dead = occ.loc[~occ["fires_in_cohort"], "window"].tolist()
    log.info("panel occupancy: mean %.1f%% across %d windows; %d fire in zero tumours%s",
             100 * float(occ["occupancy"].mean()), len(occ), len(dead),
             (" -> " + ", ".join(dead)) if dead else "")
    if len(dead) >= max(2, len(occ) // 3):
        log.warning("many panel windows are silent in this cohort; this usually means a "
                    "resolution mismatch rather than absence of signal (docs/INPUT_SPEC.md)")

    # SCNA reference model, so the panel is judged against something
    base_feat = build_feature_table(tiles, bps, ar, sample_col="sample_id",
                                    hotspot_windows=load_frozen_hotspot_windows(
                                        ref_dir / "hotspot_windows_v1.csv"))
    scna_cols = feature_groups(base_feat.columns)["scna"]

    oc = load_outcomes_table(args.outcomes)

    rows = []
    for horizon in args.horizons:
        t, e = oc["recurrence_time"], oc["recurrence_event"]
        is_event = (e == 1) & (t <= horizon)
        eligible = is_event | (t >= horizon)
        lab = oc.loc[eligible].copy()
        lab["label"] = is_event.loc[eligible].astype(int)

        d = (lab.merge(W, left_on="tumor_id", right_index=True, how="inner")
                .merge(base_feat[scna_cols], left_on="tumor_id", right_index=True, how="left"))
        d = d.fillna({c: 0.0 for c in list(W.columns) + list(scna_cols)})
        if d.empty or d["label"].nunique() < 2:
            log.warning("%dy: too few eligible tumours or a single class; skipped", horizon)
            continue

        y = d["label"].to_numpy(int)
        groups = d[["patient_id"]].reset_index(drop=True)
        P = d[list(W.columns)].reset_index(drop=True)
        S = d[scna_cols].reset_index(drop=True)

        scores: dict[str, np.ndarray] = {}
        if args.score == "count":
            # Fully frozen: unweighted count of panel windows present.
            scores["Fixed panel (count, 0 fitted params)"] = P.sum(axis=1).to_numpy(float)
        else:
            scores["Fixed panel (refit coefficients)"] = oof(
                pd.DataFrame({"_intercept": np.ones(len(d))}), P, y, groups,
                k=P.shape[1], n_repeats=args.repeats, seed=11)
        scores["Arm-level SCNA (reference)"] = oof(
            S, None, y, groups, k=0, n_repeats=args.repeats, seed=11)

        for name, sc in scores.items():
            ok = np.isfinite(sc)
            if ok.sum() < 20 or len(np.unique(y[ok])) < 2:
                continue
            auc = roc_auc_score(y[ok], sc[ok])
            rows.append({"horizon": horizon, "panel_set": args.panel_set,
                         "scoring": args.score, "model": name,
                         "n": int(ok.sum()), "events": int(y[ok].sum()),
                         "auc": round(float(auc), 4)})
            log.info("%dy  %-38s AUC=%.3f  (n=%d, events=%d)",
                     horizon, name, auc, int(ok.sum()), int(y[ok].sum()))

    perf = pd.DataFrame(rows)
    perf.to_csv(outdir / "fixed_panel_performance.csv", index=False)
    print("\n=== FIXED-PANEL VALIDATION ===")
    print(f"panel set: {args.panel_set} ({len(panel)} windows) | scoring: {args.score}")
    if not perf.empty:
        print(perf.pivot_table(index="horizon", columns="model", values="auc")
              .to_string(float_format=lambda x: f"{x:.3f}"))
    print(f"\nwrote: {outdir/'fixed_panel_performance.csv'}, "
          f"{outdir/'panel_window_occupancy.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
