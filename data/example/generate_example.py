#!/usr/bin/env python
"""Generate a small, fully synthetic meningioma-like cohort for the tutorial.

NOT real patient data. Purpose: let a new user run the whole panel pipeline
end-to-end offline and see non-trivial output. Recurrence risk is deliberately
made to depend on a breakpoint in the chr1p pericentromeric region (arm-relative
position ~0.95, i.e. the panel's headline window), so the discovery step should
recover that window and the panel should beat chance.

Writes:
  synthetic_segments.csv  (sample_id, patient_id, chromosome, start, end, event)
  synthetic_outcomes.csv  (sample_id, patient_id, recurrence_time, recurrence_event)
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

# GRCh37 autosome lengths (bp), used only to place plausible synthetic segments.
CHR_LEN = {
    1: 249250621, 2: 243199373, 3: 198022430, 4: 191154276, 5: 180915260,
    6: 171115067, 7: 159138663, 8: 146364022, 9: 141213431, 10: 135534747,
    11: 135006516, 12: 133851895, 13: 115169878, 14: 107349540, 15: 102531392,
    16: 90354753, 17: 81195210, 18: 78077248, 19: 59128983, 20: 63025520,
    21: 48129895, 22: 51304566,
}
CHR1_CENT = 121_500_000  # ~ chr1 p/q boundary region (GRCh37)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=300, help="number of tumours")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    out = Path(args.outdir)

    seg_rows, out_rows = [], []
    for i in range(args.n):
        sid = f"TUM{i:04d}"
        pid = f"PAT{i:04d}"  # one tumour per patient in the toy set

        # Driver event: a chr1p pericentromeric loss (breakpoint near the
        # centromere) in ~45% of tumours; this is the planted prognostic window.
        has_1p_peri = rng.random() < 0.45
        if has_1p_peri:
            # loss ending just below the centromere -> a breakpoint at ~0.95
            # arm-relative position on 1p.
            start = int(rng.integers(1_000_000, 60_000_000))
            end = int(rng.integers(112_000_000, 119_000_000))
            seg_rows.append([sid, pid, "1", start, end, "loss",
                             float(np.round(rng.normal(-0.45, 0.12), 3))])

        # Background genome-wide events. Counts and sizes are matched to the
        # institutional aCGH cohort the panel was trained on, because the panel
        # scores breakpoints in ~2.8 Mb windows and a toy set built from uniform
        # 3-40 Mb blocks is far coarser than real array data -- it would read as
        # "too coarse" under apply/check_resolution.py and give a new user a
        # misleading first impression.
        #
        # Event LENGTH is matched to the institutional cohort: median 1.26 Mb,
        # Q1 0.38, Q3 10.5, 47% under 1 Mb. Lognormal in log10 units, sigma
        # fitted from that IQR (sigma = IQR / 1.349).
        #
        # Event COUNT is deliberately below the institutional median of 60 rows
        # per tumour. This is a documented compromise. Background events here
        # are placed uniformly at random, whereas real breakpoints recur at
        # shared loci; at a realistic burden that uniform noise floods the one
        # planted chr1p signal and the fixed panel's AUC *falls* with horizon,
        # which is the opposite of what the example exists to demonstrate.
        #
        # 25 background events gives ~48 breakpoints/sample: above the panel's
        # own low-density QC threshold of 40, so the two QC signals agree, while
        # leaving the planted window recoverable at selection frequency 1.0.
        # Real cohorts should sit higher (institutional: ~97/sample).
        #
        # Events also carry `segment_mean` (log2 ratio) and include an allelic
        # track. Both are optional in the input spec but the enriched model's
        # amplitude and allelic blocks are built from them, so an example
        # lacking them cannot exercise the strongest SEG-only model. In the
        # institutional cohort 63.5% of segment rows are LOH / allelic
        # imbalance; the copy-number:allelic split below matches that.
        n_bg = int(np.clip(rng.lognormal(np.log(25), 1.15), 2, 300))
        for _ in range(n_bg):
            c = int(rng.integers(1, 23))
            L = CHR_LEN[c]
            length = int(np.clip(10 ** rng.normal(6.10, 1.07), 20_000, int(0.9 * L)))
            a = int(rng.integers(1_000_000, max(2_000_000, L - length - 1)))
            b = min(L, a + length)
            if rng.random() < 0.635:
                ev_lab = str(rng.choice(["loh", "allelic_imbalance"], p=[0.82, 0.18]))
                amp = float(np.round(rng.normal(0.0, 0.05), 3))  # copy-neutral
            elif rng.random() < 0.5:
                ev_lab = "loss"
                amp = float(np.round(rng.normal(-0.40, 0.18), 3))
            else:
                ev_lab = "gain"
                amp = float(np.round(rng.normal(0.38, 0.18), 3))
            seg_rows.append([sid, pid, str(c), a, b, ev_lab, amp])

        # Outcome: baseline hazard + boost if the planted window is present.
        base = 0.10
        risk = base + (0.60 if has_1p_peri else 0.0) + rng.normal(0, 0.04)
        recurs = rng.random() < np.clip(risk, 0.02, 0.95)
        if recurs:
            t = float(np.round(rng.uniform(0.3, 6.5), 2))
            ev = 1
        else:
            t = float(np.round(rng.uniform(3.0, 10.0), 2))  # event-free follow-up
            ev = 0
        out_rows.append([sid, pid, t, ev])

    seg = pd.DataFrame(seg_rows, columns=["sample_id", "patient_id", "chromosome",
                                          "start", "end", "event", "segment_mean"])
    oc = pd.DataFrame(out_rows, columns=["sample_id", "patient_id", "recurrence_time", "recurrence_event"])
    out.mkdir(parents=True, exist_ok=True)
    seg.to_csv(out / "synthetic_segments.csv", index=False)
    oc.to_csv(out / "synthetic_outcomes.csv", index=False)
    print(f"wrote {len(seg)} segment rows for {args.n} tumours -> {out/'synthetic_segments.csv'}")
    print(f"wrote {len(oc)} outcome rows -> {out/'synthetic_outcomes.csv'}")
    print(f"5-year event rate: {(oc.recurrence_event.eq(1) & oc.recurrence_time.le(5)).mean():.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
