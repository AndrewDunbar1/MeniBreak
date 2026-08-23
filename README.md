# MeniBreak — a parsimonious breakpoint panel for meningioma recurrence

Eight binary **breakpoint** indicators at fixed genomic coordinates — *where*
copy-number change points fall, not just *which* arms are gained or lost — that
carry recurrence-risk information beyond WHO grade, integrated grade, and
arm-level copy-number alteration (SCNA), and that hold up at the long (5–7 year)
horizons where the conventional markers fade.

This repository contains **the actual code that produced every number in the
MeniBreak report**, the verified result tables, and a tested entry point so an
external cohort can **run the same method on its own SEG file**.

---

## The finding (institutional cohort, patient-grouped out-of-fold AUC)

| Model | Features | Broad 5 yr | Restricted 5 yr |
|---|---:|---:|---:|
| WHO grade | 1 | 0.656 | 0.647 |
| Integrated grade | 1 | 0.735 | 0.782 |
| Arm-level SCNA | 47 | 0.766 | 0.766 |
| **8-window panel alone** | **8** | **0.735** | **0.753** |
| SCNA + panel | 55 | 0.797 | 0.847 |
| Clinical + SCNA + panel | 60 | 0.843 | 0.843 |
| Clinical + SCNA + enriched breakpoints | 102 | 0.857 | 0.836 |

The panel **alone** matches integrated grade at 5 years and beats it at 7; it is
the only marker that *improves* with horizon. The chr1p pericentromeric window
(chr1:118.75–125.00 Mb) is selected in 88–96% of folds across horizons. See
[`docs/THE_PANEL.md`](docs/THE_PANEL.md).

## Three caveats that travel with every number

1. **Internal validation only.** Every AUC here is patient-grouped out-of-fold
   within a single institutional cohort. Validation in an independent cohort of
   comparable copy-number resolution is the necessary next step, and is exactly
   what [`apply/menibreak_panel.py`](apply/menibreak_panel.py) exists to enable.
   See [`docs/LIMITATIONS.md`](docs/LIMITATIONS.md).
2. **The restricted cohort has 22/37/43 events** at 3/5/7 years — its largest
   effects are its least precise.
3. **Research use only.** No claim of clinical utility or prospective validation.

---

## Run it on YOUR data (external cohort)

You need a segmented copy-number file and recurrence outcomes — nothing else.
Full walk-through: [`docs/QUICKSTART.md`](docs/QUICKSTART.md); input format:
[`docs/INPUT_SPEC.md`](docs/INPUT_SPEC.md).

### 0. First, check your data's resolution (30 seconds)

```bash
python apply/check_resolution.py --segments mycohort_segments.csv
```

Any SEG file will run. The question is whether the assay behind it localized
breakpoints finely enough for the panel to mean anything. The panel scores
breakpoints in bins of 5% of a chromosome arm — about **2.8 Mb**. If your typical
event is smaller than that, breakpoints resolve to individual bins and the panel
applies. If events are much larger, positions smear across neighbouring bins and
the window features degrade to noise, while arm-level features are unaffected.

Two measured reference points:

| Platform | Median event | vs 2.8 Mb window | Panel |
|---|---|---|---|
| Agilent aCGH (training cohort) | 1.26 Mb | 0.45× | works |
| conumee on EPIC methylation array | 5.10 Mb | 1.80× | lost to arm-level SCNA |

**Methylation-derived copy number is the common failure case.** If your SEG came
from 450k/EPIC via conumee, expect `TOO COARSE` and use the arm-level SCNA
reference model, which did transport. SEG from aCGH, SNP array, OncoScan,
panel-seq or WGS should read `SUITABLE`.

This is calibrated on two platforms, so treat it as a decision aid rather than a
validated instrument — `MARGINAL` means "run both models and compare on your own
outcomes", not "expect a specific amount of degradation".

There are then **two entry points, answering two different questions.**

### 1. Validate the published panel (fixed windows)

> *Does the published MeniBreak panel discriminate recurrence in my cohort?*

The nine windows are **fixed**, read from
[`reference/menibreak_panel_v1.csv`](reference/menibreak_panel_v1.csv); nothing
about them is learned from your data. Default scoring is a plain count of panel
windows containing a breakpoint — **zero parameters fitted on your cohort**, so
the panel either transfers or it does not.

```bash
pip install -r requirements.txt

python apply/validate_fixed_panel.py \
    --segments YOUR_segments.csv --outcomes YOUR_outcomes.csv \
    --build GRCh37 --out results/your_cohort
```

Writes `fixed_panel_performance.csv` and `panel_window_occupancy.csv` (which of
the fixed windows actually fire in your cohort — check this first). Use
`--score refit` to keep the windows fixed but refit coefficients under
patient-grouped CV, and `--panel-set rule` for the 11-window variant.

### 2. Re-run the discovery method

> *Does the discovery procedure find prognostic windows in my cohort?*

Windows are re-selected inside each training fold of **your** data, so this
replicates the method rather than testing the published marker.

```bash
python apply/menibreak_panel.py \
    --segments YOUR_segments.csv --outcomes YOUR_outcomes.csv \
    --out results/your_cohort_discovery
```

Writes `performance.csv`, `hotlist_windows.csv` (the windows *your* cohort
selects, in genomic coordinates) and `window_features.csv`. Every modelling
choice matches the institutional analysis in `analyses/03_full_ladder.py`.

Both scripts run on the bundled synthetic cohort out of the box — start there to
confirm your install:

```bash
python data/example/generate_example.py --n 300 --seed 7 --outdir data/example
python apply/validate_fixed_panel.py \
    --segments data/example/synthetic_segments.csv \
    --outcomes data/example/synthetic_outcomes.csv --out results/example_fixed
```

## Reproduce the institutional numbers

The scripts in [`analyses/`](analyses/) regenerate each table of the report on the
institutional dataset (IRB-restricted, **not** shipped). Point them at the data:

```bash
export MENIBREAK_DATA_ROOT=/path/to/Breakpoint_Architecture_Paper_Analysis_Handoff/
python analyses/03_full_ladder.py --k 8 --repeats 5
```

WHO/integrated-grade comparators reproduce exactly; the panel/ladder/windows
regenerate from the same code that produced the report.

---

## Repository map

| Path | What it is |
|---|---|
| `apply/check_resolution.py` | **Preflight: is your SEG fine-grained enough for the panel?** Run this first |
| `apply/validate_fixed_panel.py` | **Validate the published panel on your SEG** — fixed windows, no discovery |
| `apply/menibreak_panel.py` | **Re-run the discovery method on your SEG** |
| `src/` | The pipeline package: `harmonization/` (SEG → arms → breakpoints), `features/` (windows, SCNA, enriched), `evaluation/` (patient-grouped folds) |
| `analyses/` | The six scripts that produced the report tables |
| `results/reference/` | The 13 verified result CSVs (match the report to the printed digit) |
| `reference/` | Arm references (GRCh37/GRCh38), frozen hotspot windows, and the **frozen panel definition** + the script that derives it |
| `tests/` | 79 pytest tests — patient-leakage guarantees, tiling, feature isolation, external-input handling, resolution preflight ([details](tests/README.md)) |
| `data/example/` | Synthetic cohort + generator, so everything runs offline |
| `figures/` | Publication figures, [with captions](figures/README.md) |
| `RESULTS_MANIFEST.yaml` / `.json` | **Machine-readable index**: 725 result records → script → output → value → status ([regenerate](build_manifest.py)) |
| `docs/` | Input spec, quickstart, the panel, provenance, limitations |

## Provenance

The report numbers originate in the local repository
`BreakpointsFoundationModel`, packaged to a drive export and reassembled here —
**not** a ChatGPT thread (a prior claim that the code lived only in an online
workspace was checked against every file's metadata and rejected). Full audit:
[`docs/PROVENANCE.md`](docs/PROVENANCE.md).

Research use only. See [`LICENSE.md`](LICENSE.md) and [`docs/LIMITATIONS.md`](docs/LIMITATIONS.md).
