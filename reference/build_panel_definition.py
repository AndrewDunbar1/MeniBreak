#!/usr/bin/env python
"""Derive the frozen MeniBreak panel definition from the discovery output.

The panel is the set of windows the in-fold discovery procedure re-selected
often enough to be a finding rather than a coincidence:

    union over horizons of {windows selected in >= 40% of training folds at k=5}

Freezing these turns a *procedure* into a *marker*, so an external cohort can
validate the published panel rather than re-run discovery on its own data.

NOTE — a documented discrepancy. That rule yields ELEVEN windows, while the
original report published NINE. The rule additionally admits `bp_14q_b18`
(0.40 at 7y) and `bp_2q_b14` (0.44 at 5y), yet the published table retains
`bp_3p_b10` and `bp_6q_b19` at the same 0.40 frequency. So the published
nine-window list is not recoverable from any single threshold on the retained
output. Rather than quietly pick one, both sets are emitted and labelled:

    in_published_table = True   -> the 9 windows printed in the report
    in_published_table = False  -> the 2 the reproducible rule also admits

`apply/validate_fixed_panel.py` defaults to the published nine (that is the
claim being validated) and exposes `--panel-set rule` for the eleven.

Reads   results/reference/hotlist_windows.csv
Writes  reference/menibreak_panel_v1.csv
"""
from __future__ import annotations
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import pandas as pd
from src.features.hotlist import windows_to_genomic_coordinates
from src.harmonization.references import ArmReference

K = 5              # the k at which the published stability table was reported
MIN_FREQ = 0.40    # "a finding, not a coincidence"
N_BINS = 20


def main() -> int:
    src = pd.read_csv(ROOT / "results" / "reference" / "hotlist_windows.csv")
    sel = src[(src["k"] == K) & (src["selection_frequency"] >= MIN_FREQ)].copy()

    # one row per window; keep the per-horizon frequencies as evidence
    freq = (sel.pivot_table(index="window", columns="horizon",
                            values="selection_frequency", aggfunc="max")
               .rename(columns={3: "sel_freq_3y", 5: "sel_freq_5y", 7: "sel_freq_7y"}))

    base = (sel.sort_values("selection_frequency", ascending=False)
               .drop_duplicates("window")[["window", "arm", "bin", "start_norm", "end_norm"]]
               .set_index("window"))
    panel = base.join(freq).reset_index()

    ar = ArmReference.from_csv(ROOT / "reference" / "arm_reference_GRCh37.csv", "GRCh37")
    panel = windows_to_genomic_coordinates(panel, ar, n_bins=N_BINS)

    # The nine windows printed in the report (docs/THE_PANEL.md). Held explicitly
    # because they are not recoverable from a single threshold - see module docstring.
    PUBLISHED_NINE = {
        "bp_1p_b18", "bp_1p_b19", "bp_6q_b18", "bp_6q_b19",
        "bp_19q_b9", "bp_19q_b10", "bp_2p_b15", "bp_4p_b19", "bp_3p_b10",
    }
    panel["in_published_table"] = panel["window"].isin(PUBLISHED_NINE)
    missing = PUBLISHED_NINE - set(panel["window"])
    if missing:
        raise SystemExit(f"published window(s) not found in discovery output: {sorted(missing)}")

    panel["n_bins"] = N_BINS
    panel["panel_version"] = "v1"
    cols = ["window", "arm", "bin", "n_bins", "start_norm", "end_norm",
            "chrom", "genomic_start_mb", "genomic_end_mb",
            "sel_freq_3y", "sel_freq_5y", "sel_freq_7y",
            "in_published_table", "panel_version"]
    panel = panel[[c for c in cols if c in panel.columns]]
    panel = panel.sort_values(["arm", "bin"]).reset_index(drop=True)

    out = ROOT / "reference" / "menibreak_panel_v1.csv"
    panel.to_csv(out, index=False)
    print(f"wrote {out}  ({len(panel)} windows, k={K}, min selection frequency {MIN_FREQ})")
    print(panel.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
