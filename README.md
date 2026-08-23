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

You need a segmented copy-number file and recurrence outcomes. Full walk-through:
[`docs/QUICKSTART.md`](docs/QUICKSTART.md); input format:
[`docs/INPUT_SPEC.md`](docs/INPUT_SPEC.md).

```bash
pip install -r requirements.txt

# 1) sanity-check on the bundled synthetic cohort (runs in ~1 min)
python apply/menibreak_panel.py \
    --segments data/example/synthetic_segments.csv \
    --outcomes data/example/synthetic_outcomes.csv \
    --out results/example

# 2) then point it at your own files (same two-column contract)
python apply/menibreak_panel.py \
    --segments YOUR_segments.csv --outcomes YOUR_outcomes.csv \
    --build GRCh37 --out results/your_cohort
```

It writes `performance.csv` (the AUC ladder on your cohort), `hotlist_windows.csv`
(the windows it discovers, in genomic coordinates), and `window_features.csv`
(per-tumour panel occupancy). Every modelling choice matches the institutional
analysis in `analyses/03_full_ladder.py`.

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
| `apply/menibreak_panel.py` | **Run the panel on your own SEG** (tested entry point) |
| `src/` | The pipeline package: `harmonization/` (SEG → arms → breakpoints), `features/` (windows, SCNA, enriched), `evaluation/` (patient-grouped folds) |
| `analyses/` | The six scripts that produced the report tables |
| `results/reference/` | The verified result CSVs (match the report to the printed digit) |
| `reference/` | Arm reference (GRCh37/GRCh38) and frozen hotspot windows |
| `data/example/` | Synthetic cohort + generator, so everything runs offline |
| `figures/` | Publication figures |
| `RESULTS_MANIFEST.yaml` / `.json` | **Machine-readable index**: every result → script → output → value → status ([regenerate](build_manifest.py)) |
| `docs/` | Input spec, quickstart, the panel, methods, provenance, limitations |

## Provenance

The report numbers originate in the local repository
`BreakpointsFoundationModel`, packaged to a drive export and reassembled here —
**not** a ChatGPT thread (a prior claim that the code lived only in an online
workspace was checked against every file's metadata and rejected). Full audit:
[`docs/PROVENANCE.md`](docs/PROVENANCE.md).

Research use only. See [`LICENSE.md`](LICENSE.md) and [`docs/LIMITATIONS.md`](docs/LIMITATIONS.md).
