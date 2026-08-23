#!/usr/bin/env python
"""Does the external entry point reproduce the published numbers?

`analyses/` reads the institutional master SEG through the internal loader.
`apply/` reads a generic SEG through `adapt_segments`. An external team only
ever touches the second path, so the question that matters for deployment is
whether the two agree -- not whether `apply/` runs without crashing.

This routes the institutional cohort, a DNA-based assay (Agilent aCGH), through
the *external* door: master SEG rewritten into the documented generic format,
handed to `apply/menibreak_panel.py`, and the resulting ladder compared against
`results/reference/full_ladder.csv`.

If a column mapping, an event vocabulary or a coordinate convention differs
between the two paths, the AUCs separate here.

Only the Broad scope is checked. The restricted cohort is defined by clinical
criteria (GTR, primary, no prior RT, not radiation-induced, not angiomatous)
that the generic SEG format does not carry, so an external cohort cannot
reproduce it without supplying those flags itself.

Requires MENIBREAK_DATA_ROOT. Not runnable by external users -- it exists to
prove the path they *do* use is faithful.

Usage:
    python analyses/07_verify_external_path.py
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

DATA_ROOT = os.environ.get("MENIBREAK_DATA_ROOT")
if not DATA_ROOT:
    raise SystemExit(
        "MENIBREAK_DATA_ROOT is not set.\n"
        "This script reproduces the institutional analysis and needs the\n"
        "IRB-restricted dataset, which is not shipped with this repository.\n"
        "External users: see apply/ instead -- it needs no institutional data."
    )

H = Path(DATA_ROOT)
MASTER_SEG = H / "01_RAW_DATA/institutional/meningioma_masterseg_ALL_2026.txt"
OUTCOMES = H / "03_ANALYSIS_READY_DATA/outcomes.csv"

# Institutional Event vocabulary -> the documented generic vocabulary.
EVENT_MAP = {
    "CN Gain": "gain",
    "CN Loss": "loss",
    "Homozygous Copy Loss": "homozygous_loss",
    "High Copy Gain": "high_copy_gain",
    "LOH": "loh",
    "Allelic Imbalance": "allelic_imbalance",
}

# Models the external path can produce. Clinical models need age/mitoses/GTR,
# which the generic SEG format does not carry by design.
COMPARABLE = ["Arm-level SCNA", "8-window panel alone", "SCNA + 8-window panel"]

TOL = 0.01  # AUC agreement required to call the path faithful


def export_external_format(outdir: Path) -> tuple[Path, Path]:
    """Rewrite the institutional data as a generic external submission."""
    raw = pd.read_csv(MASTER_SEG, sep="\t", low_memory=False)

    # `cg` is the deidentified tumour key and the one outcomes.csv joins on;
    # `sample` is the assay run and does not link.
    seg = pd.DataFrame({
        "sample_id": raw["cg"].astype(str),
        "chromosome": raw["Chr"].astype(str),
        "start": raw["Start"].astype("int64"),
        "end": raw["End"].astype("int64"),
        "event": raw["Event"].map(EVENT_MAP),
    })
    unmapped = raw.loc[seg["event"].isna(), "Event"].unique()
    if len(unmapped):
        raise SystemExit(f"unmapped Event values: {list(unmapped)}")

    seg = seg[seg["sample_id"].ne("nan")].reset_index(drop=True)

    out = pd.read_csv(OUTCOMES)
    keep = set(seg["sample_id"]) & set(out["sample_id"].astype(str))
    seg = seg[seg["sample_id"].isin(keep)]
    out = out[out["sample_id"].astype(str).isin(keep)]

    sp, op = outdir / "external_segments.csv", outdir / "external_outcomes.csv"
    seg.to_csv(sp, index=False)
    out.to_csv(op, index=False)
    print(f"exported {len(seg):,} segments / {seg['sample_id'].nunique():,} tumours "
          f"/ {len(out):,} outcome rows")
    return sp, op


def main() -> int:
    work = ROOT / "results" / "_verify_external"
    work.mkdir(parents=True, exist_ok=True)
    sp, op = export_external_format(work)

    print("\nrunning apply/menibreak_panel.py (the external entry point)...")
    r = subprocess.run(
        [sys.executable, str(ROOT / "apply" / "menibreak_panel.py"),
         "--segments", str(sp), "--outcomes", str(op),
         "--k", "8", "--repeats", "5", "--out", str(work / "run")],
        cwd=ROOT, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-3000:]); print(r.stderr[-3000:])
        raise SystemExit("external entry point failed")

    got = pd.read_csv(work / "run" / "performance.csv")
    ref = pd.read_csv(ROOT / "results" / "reference" / "full_ladder.csv")
    ref = ref[ref["scope"] == "Broad"]

    rows = []
    for _, g in got.iterrows():
        m = str(g["model"])
        if m not in COMPARABLE:
            continue
        hit = ref[(ref["model"] == m) & (ref["horizon"] == g["horizon"])]
        if hit.empty:
            continue
        pub = float(hit["auc"].iloc[0])
        rows.append({
            "horizon": int(g["horizon"]), "model": m,
            "published": round(pub, 4), "via_apply": round(float(g["auc"]), 4),
            "delta": round(float(g["auc"]) - pub, 4),
            "n_pub": int(hit["n"].iloc[0]), "n_apply": int(g["n"]),
            "events_pub": int(hit["events"].iloc[0]), "events_apply": int(g["events"]),
        })

    res = pd.DataFrame(rows).sort_values(["horizon", "model"])
    print("\n=== PUBLISHED vs REPRODUCED THROUGH apply/ ===")
    print(res.to_string(index=False))

    res.to_csv(work / "external_path_agreement.csv", index=False)

    worst = res["delta"].abs().max()
    n_ok = (res["n_pub"] == res["n_apply"]).all()
    e_ok = (res["events_pub"] == res["events_apply"]).all()
    print(f"\nlargest |delta| AUC : {worst:.4f}   (tolerance {TOL})")
    print(f"cohort sizes match  : {n_ok}")
    print(f"event counts match  : {e_ok}")

    ok = worst <= TOL and n_ok and e_ok
    print("\nVERDICT:", "FAITHFUL" if ok else "DISCREPANCY -- investigate")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
