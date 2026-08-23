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
python apply/menibreak_panel.py \
    --segments data/example/synthetic_segments.csv \
    --outcomes data/example/synthetic_outcomes.csv \
    --out results/example
```

Expected: it derives breakpoints, builds the 20-bin window matrix, prints a
resolution-QC line, and reports a small AUC ladder. On the synthetic data the
panel-alone AUC rises with horizon (≈0.54 → 0.63 → 0.66), the discovery step
recovers a chr1p pericentromeric window at selection frequency 1.0, and a QC
warning fires because the toy data is deliberately sparse. Your real data will
have far more breakpoints per sample.

## 2. Format your data

Follow [`INPUT_SPEC.md`](INPUT_SPEC.md). The minimum is:

- `segments.csv`: `sample_id, chromosome, start, end, event` (+ `patient_id` if a
  patient can have more than one tumour).
- `outcomes.csv`: `sample_id, patient_id, recurrence_time, recurrence_event`.

If your caller uses different event words, either recode them to
`loss`/`gain`/… or add them to `EVENT_SYNONYMS` in `apply/menibreak_panel.py`.

## 3. Run on your cohort

```bash
python apply/menibreak_panel.py \
    --segments YOUR_segments.csv \
    --outcomes YOUR_outcomes.csv \
    --build GRCh37 \
    --horizons 3 5 7 \
    --k 8 --repeats 5 \
    --out results/your_cohort
```

Outputs in `results/your_cohort/`:

| File | Contents |
|---|---|
| `performance.csv` | Patient-grouped out-of-fold AUC for `Arm-level SCNA`, `8-window panel alone`, `SCNA + 8-window panel`, at each horizon. |
| `hotlist_windows.csv` | Windows your cohort selects across folds, with `selection_frequency` and **genomic coordinates** in your build. |
| `window_features.csv` | Per-tumour occupancy of every arm-relative window (the model inputs). |

## 4. Read the QC line first

Before trusting the AUCs, check the printed
`median breakpoints/sample` and `mean window occupancy`. The panel was defined on
data with a **median ~97 breakpoints/sample and 27.8% window occupancy**. If your
platform is much coarser (the tool warns below 40/sample), the ~5 Mb windows are
not comparable and a same-resolution cohort is needed — this is the measured
reason validation failed on an external methylation cohort (see
[`LIMITATIONS.md`](LIMITATIONS.md)).

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
