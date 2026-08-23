# Quickstart — run the panel on your own cohort

This is the whole workflow for an external group. It needs only a segmented
copy-number file and recurrence outcomes; no institutional data is required.

## 0. Install

```bash
git clone https://github.com/AndrewDunbar1/MeniBreak.git
cd MeniBreak
pip install -r requirements.txt     # pandas, numpy, scipy, scikit-learn, pyyaml
```

## 1. Confirm it runs (bundled synthetic cohort)

```bash
python data/example/generate_example.py --n 300 --seed 7 --outdir data/example

# the validation run: fixed published panel, zero fitted parameters
python apply/validate_fixed_panel.py \
    --segments data/example/synthetic_segments.csv \
    --outcomes data/example/synthetic_outcomes.csv \
    --out results/example_fixed

# the method run: rediscover windows inside your own folds
python apply/menibreak_panel.py \
    --segments data/example/synthetic_segments.csv \
    --outcomes data/example/synthetic_outcomes.csv \
    --out results/example
```

Expected on this synthetic cohort: breakpoints are derived, a resolution-QC line
is printed, and both runs report an AUC ladder that **rises with horizon** — the
fixed panel reaches ≈0.66 → 0.69 → 0.72 and the discovery run ≈0.54 → 0.63 →
0.66, with the planted chr1p pericentromeric window recovered at selection
frequency 1.0. A QC warning fires because the toy data is deliberately sparse;
real data will have far more breakpoints per sample.

You can also run the test suite (53 tests, ~1 s):

```bash
pip install -r requirements-dev.txt && python -m pytest tests/ -q
```

## 2. Format your data

Follow [`INPUT_SPEC.md`](INPUT_SPEC.md). The minimum is:

- `segments.csv`: `sample_id, chromosome, start, end, event` (+ `patient_id` if a
  patient can have more than one tumour).
- `outcomes.csv`: `sample_id, patient_id, recurrence_time, recurrence_event`.

If your caller uses different event words, either recode them to
`loss`/`gain`/… or add them to `EVENT_SYNONYMS` in `apply/menibreak_panel.py`.

## 3. Run on your cohort

### 3a. The validation run — fixed published panel (start here)

This is the one that tests **our** panel. The windows come from
`reference/menibreak_panel_v1.csv` and nothing is learned from your data; with
the default `--score count` no parameters are fitted at all.

```bash
python apply/validate_fixed_panel.py \
    --segments YOUR_segments.csv \
    --outcomes YOUR_outcomes.csv \
    --build GRCh37 \
    --horizons 3 5 7 \
    --out results/your_cohort
```

| File | Contents |
|---|---|
| `fixed_panel_performance.csv` | AUC of the fixed panel vs an arm-level SCNA reference, per horizon. |
| `panel_window_occupancy.csv` | How often each fixed window fires in your cohort. **Read this first** — a panel window that fires in zero tumours cannot contribute, and several silent windows usually means a resolution mismatch rather than absent signal. |

Useful variants: `--score refit` (fixed windows, coefficients refit under
patient-grouped CV — a fairer test if your platform or population differs) and
`--panel-set rule` (11 windows instead of the published 9; see
[`THE_PANEL.md`](THE_PANEL.md)).

### 3b. The method run — rediscover windows in your data

```bash
python apply/menibreak_panel.py \
    --segments YOUR_segments.csv \
    --outcomes YOUR_outcomes.csv \
    --build GRCh37 \
    --horizons 3 5 7 \
    --k 8 --repeats 5 \
    --out results/your_cohort_discovery
```

This selects windows inside your own training folds, so it answers "does the
*method* work here?", not "does the *published panel* work here?". It is the
right run for asking whether your cohort has its own recurrent breakpoint
structure — and whether it agrees with ours.

Outputs in `results/your_cohort_discovery/`:

| File | Contents |
|---|---|
| `performance.csv` | Patient-grouped out-of-fold AUC for `Arm-level SCNA`, `8-window panel alone`, `SCNA + 8-window panel`, at each horizon. |
| `hotlist_windows.csv` | Windows your cohort selects across folds, with `selection_frequency` and **genomic coordinates** in your build. |
| `window_features.csv` | Per-tumour occupancy of every arm-relative window (the model inputs). |

## 4. Read the QC line first

Before trusting the AUCs, check the printed `median breakpoints/sample` and
`mean window occupancy`. The panel was defined on data with a **median ~97
breakpoints/sample and 27.8% window occupancy**. Window-based features assume a
change-point detector of comparable sensitivity: if your platform is much coarser
(the tool warns below 40 breakpoints/sample), the ~5 Mb windows fire too rarely
to be comparable, and a weak result should not be read as a test of the panel.
See [`LIMITATIONS.md`](LIMITATIONS.md).

## 5. Interpret

- Does `8-window panel alone` approach or beat integrated grade in your cohort,
  and does its advantage grow at 5–7 years? That is the institutional finding.
- Do your `hotlist_windows.csv` top windows include the chr1p pericentromeric
  region (chr1:118.75–125.00 Mb, GRCh37)? That is the panel's most stable feature
  ([`THE_PANEL.md`](THE_PANEL.md)).
- `SCNA + panel` vs `SCNA` alone is the incremental-value test.

## Full comparison against grade / clinical models

`apply/menibreak_panel.py` covers the seg-derivable rungs (SCNA, panel,
SCNA+panel). To reproduce the complete ladder including WHO/integrated grade and
the clinical + 13/60/102-feature models, use `analyses/03_full_ladder.py` with
your grade and clinical covariates wired in the same shape as
`src/features/clinical.py` expects (see comments there).
