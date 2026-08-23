#!/usr/bin/env python
"""Preflight: is this SEG file fine-grained enough for the 8-window panel?

Any SEG file will *run*. That is not the question. The question is whether the
assay behind it localized breakpoints precisely enough for window-level features
to mean anything.

The panel scores breakpoints in bins of 5% of a chromosome arm. With a median
autosomal arm of ~55 Mb that is a ~2.8 Mb window. If the caller's typical event
is smaller than that window, a breakpoint lands in a well-defined bin and the
panel works. If typical events are much larger than the window, breakpoint
positions are smeared across neighbouring bins and the window features degrade
to noise -- while arm-level features, which are 20 bins wide, are unaffected.

Two reference points, both measured (see reports/platform_resolution.csv in the
development tree):

    Agilent aCGH, institutional cohort   median event 1.26 Mb  = 0.45x window
    conumee on EPIC methylation array    median event 5.10 Mb  = 1.80x window

The panel was trained on the first and lost its advantage over arm-level SCNA on
the second, so the practical cutoff sits between them.

CALIBRATION CAVEAT. Those are two platforms, not a validation study. The bands
below are a decision aid, not a validated instrument. AMBER especially means
"run both models and compare", not "expect degraded performance of size X".

Usage:
    python apply/check_resolution.py --segments mycohort_segments.csv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from apply.menibreak_panel import adapt_segments  # noqa: E402
from src.harmonization.references import ArmReference  # noqa: E402

WINDOW_FRAC = 0.05

# Ratio of median called-event size to the panel window width.
GREEN_MAX = 1.0    # events finer than the window
AMBER_MAX = 1.5    # marginal; between the two measured platforms


def median_arm_length(ref_dir: Path, build: str = "GRCh37") -> float:
    ref = ArmReference.from_csv(ref_dir / f"arm_reference_{build}.csv", build)
    return float(np.median([a.length for a in ref._arms if a.chromosome.isdigit()]))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--segments", required=True)
    ap.add_argument("--reference-dir", default=str(ROOT / "reference"))
    ap.add_argument("--build", default="GRCh37")
    args = ap.parse_args()

    seg = adapt_segments(args.segments, build=args.build)
    # adapt_segments already drops NEUTRAL rows, so this is called events only --
    # the like-for-like object across platforms that tile vs. those that do not.
    cn = seg[seg["track"] == "cn"].copy()
    if cn.empty:
        print("No copy-number events found. Nothing to assess.")
        return 1

    cn["length"] = cn["end"] - cn["start"]
    L = cn["length"]
    per_sample = cn.groupby("sample_id").size()

    edges = pd.concat([
        cn[["sample_id", "chromosome", "start"]].rename(columns={"start": "pos"}),
        cn[["sample_id", "chromosome", "end"]].rename(columns={"end": "pos"}),
    ]).drop_duplicates()
    cps = edges.groupby("sample_id").size()

    med_arm = median_arm_length(Path(args.reference_dir), args.build)
    window_mb = WINDOW_FRAC * med_arm / 1e6
    med_mb = float(L.median()) / 1e6
    ratio = med_mb / window_mb

    print(f"file                      {Path(args.segments).name}")
    print(f"samples                   {cn['sample_id'].nunique()}")
    print(f"copy-number events        {len(cn):,}")
    print(f"events per sample         {per_sample.median():.0f} (median)")
    print(f"change-points per sample  {cps.median():.0f} (median)")
    print(f"median event size         {med_mb:.2f} Mb")
    print(f"  events under 1 Mb       {(L < 1e6).mean()*100:.1f}%")
    print(f"  events under 5 Mb       {(L < 5e6).mean()*100:.1f}%")
    print(f"panel window              {window_mb:.2f} Mb (5% of a {med_arm/1e6:.1f} Mb arm)")
    print(f"median event / window     {ratio:.2f}x")
    print()

    if ratio <= GREEN_MAX:
        print("VERDICT: SUITABLE for the 8-window panel.")
        print("  Typical events are finer than the scoring window, so breakpoint")
        print("  positions resolve to individual bins. Run validate_fixed_panel.py.")
        code = 0
    elif ratio <= AMBER_MAX:
        print("VERDICT: MARGINAL.")
        print("  Events are somewhat coarser than the scoring window. The panel may")
        print("  still carry signal, but do not assume it. Run the panel AND the")
        print("  arm-level SCNA reference and compare them on your own outcomes;")
        print("  report whichever wins, and report that you checked.")
        code = 0
    else:
        print("VERDICT: TOO COARSE for the 8-window panel.")
        print(f"  Typical events are {ratio:.1f}x the scoring window, so a breakpoint")
        print("  cannot be localized to a single bin. This is the regime where the")
        print("  panel lost to arm-level SCNA on methylation-array data.")
        print("  Use the arm-level SCNA model, which is unaffected by this and did")
        print("  transport. Treat any window-level result as exploratory.")
        code = 2

    print()
    print("Reference: aCGH 0.45x (panel trained here) | EPIC/conumee 1.80x (panel failed)")
    print("Calibrated on two platforms. A decision aid, not a validated instrument.")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
