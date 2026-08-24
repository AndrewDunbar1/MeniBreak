# AGENTS.md — instructions for a coding agent working in this repository

You are most likely here because an external lab wants to **validate the MeniBreak
breakpoint hot list on their own cohort**. This file tells you exactly what to do.
A machine-parseable version of the same protocol is in
[`VALIDATION_PROTOCOL.yaml`](VALIDATION_PROTOCOL.yaml).

## Your job, in one sentence

Get the lab's segmented copy-number data and recurrence outcomes into two CSVs,
run the fixed-panel validation, and return two small result files. **Nothing here
requires modelling expertise or any data from us.**

---

## Setup

```bash
pip install -r requirements.txt
```

Python ≥3.9. Everything runs locally on a laptop in minutes. No network access is
needed after install, and no data leaves the machine.

Sanity-check the install on the bundled synthetic cohort before touching real
data:

```bash
python data/example/generate_example.py --n 300 --seed 7 --outdir data/example
python apply/validate_fixed_panel.py \
    --segments data/example/synthetic_segments.csv \
    --outcomes data/example/synthetic_outcomes.csv \
    --out results/example_fixed
```

Expected: a `SUITABLE` resolution line and an AUC ladder rising with horizon
(~0.61 → 0.64 → 0.72). If that works, the install is fine.

---

## Step 1 — Build the two input files

### `segments.csv` (required)

One row per **altered** copy-number segment. Do not list neutral regions.

| Column | Type | Notes |
|---|---|---|
| `sample_id` | string | Tumour identifier |
| `chromosome` | string/int | `1`–`22` or `chr1`–`chr22`; X/Y ignored |
| `start` | int | bp |
| `end` | int | bp |
| `event` | string | see vocabulary below |
| `patient_id` | string | **optional but important** — supply it if any patient contributed more than one tumour; folds are grouped by patient to prevent leakage. Defaults to `sample_id`. |

Optional extras that unlock more: `segment_mean` (log2 ratio → amplitude
features), `loh` / `allelic_imbalance` rows (→ allelic features), `tumor_id`,
`probe_count`, `platform`, `caller`, `genome_build`.

**Event vocabulary.** Canonical: `loss`, `gain`, `homozygous_loss`,
`high_copy_gain`, `loh`, `allelic_imbalance`, `neutral`. Common synonyms are
mapped automatically (`del`, `deletion`, `-1`, `amp`, `dup`, `homdel`, `-2`,
`cnloh`, …). An **unrecognised label raises an error** rather than being silently
treated as neutral — if that happens, either recode the column or add the synonym
to `EVENT_SYNONYMS` in `apply/menibreak_panel.py`.

### `outcomes.csv` (required)

| Column | Type | Notes |
|---|---|---|
| `sample_id` | string | must match `segments.csv` |
| `recurrence_time` | float | **YEARS**, from treatment to recurrence or last follow-up |
| `recurrence_event` | 0/1 | 1 = recurrence/progression observed |
| `patient_id` | string | optional; defaults to `sample_id` |

> ⚠️ **`recurrence_time` must be in years.** Months or days will not crash — they
> silently push nearly every tumour past the horizon and return a meaningless AUC.
> The loader rejects a median > 25 for this reason. Convert: months ÷ 12,
> days ÷ 365.25.

Full schema, including the optional clinical file: [`docs/INPUT_SPEC.md`](docs/INPUT_SPEC.md).

---

## Step 2 — Resolution gate (run this before anything else)

```bash
python apply/check_resolution.py --segments YOUR_segments.csv
```

Three possible verdicts:

| Verdict | What to do |
|---|---|
| `SUITABLE` | Proceed to step 3. |
| `MARGINAL` | Proceed, but report the verdict alongside the result. |
| `TOO COARSE` | **STOP.** Report this back and do not proceed. |

**Why this gate matters, and why you must not work around it.** The panel scores
breakpoints in windows of ~2.8 Mb. A platform that localises breakpoints coarsely
— methylation-array–derived copy number is the usual case, at roughly one-ninth
the change-point density — cannot test this marker. A weak result on such data
measures the assay, not the biology, and would be a misleading contribution.
If you get `TOO COARSE`, say so; do not tune parameters to get a number.

---

## Step 3 — The primary test (this is the deliverable)

```bash
python apply/validate_fixed_panel.py \
    --segments YOUR_segments.csv \
    --outcomes YOUR_outcomes.csv \
    --build GRCh37 \
    --out results/validation
```

Use `--build GRCh38` if the coordinates are GRCh38.

This scores a **fixed nine-window panel** read from
`reference/menibreak_panel_v1.csv`. Nothing about the windows is learned from the
lab's data, and with the default `--score count` **no parameters are fitted at
all** — the risk score is simply how many panel windows contain a breakpoint.

**Return these two files:**

- `results/validation/fixed_panel_performance.csv`
- `results/validation/panel_window_occupancy.csv`

Neither contains patient-level data. Check `panel_window_occupancy.csv` first: a
window firing in zero tumours cannot contribute, and several silent windows
usually indicates a resolution mismatch rather than absent signal.

---

## Step 4 — How to interpret (do not skip)

The tool prints an arm-level SCNA AUC beside the panel AUC. **That is the
comparison that matters**, because it is computed on the lab's own data and folds.

For context, on the derivation cohort the fixed panel scores **0.690 / 0.772 /
0.785** at 3/5/7 years vs arm-level SCNA at 0.777 / 0.766 / 0.767.

> ⚠️ **Those derivation numbers are an optimistic anchor, not a target.** The nine
> windows were selected using that cohort's outcomes, so scoring them there is
> in-sample with respect to window selection. **A lower number in the lab's cohort
> is the expected outcome, not a failure.** Do not treat falling short of 0.772 as
> a negative result, and do not tune anything to approach it.
> See [`docs/WHICH_MODEL.md`](docs/WHICH_MODEL.md).

---

## Step 5 — Only if the primary test is flat

A null result is a legitimate, publishable outcome and should be reported as
such. Do **not** run these first and report the best one. In order:

1. `--score refit` — same fixed windows, weights re-fitted on the lab's cohort
   under patient-grouped CV. Tests whether the *window set* carries signal even if
   our unweighted scoring is not calibrated to this population.
2. `--panel-set rule` — an 11-window variant of the list.
3. `python apply/menibreak_panel.py --segments … --outcomes … --out results/discovery`
   — ignores our list and rediscovers windows inside the lab's own folds. Answers
   whether their cohort has recurrent prognostic breakpoint structure at all, and
   whether it lands near ours (look for chr1:118.75–125.00 Mb in
   `hotlist_windows.csv`).

Report every run you performed, not only the favourable ones.

---

## Optional — the full model ladder (secondary)

If the lab has clinical fields, this runs all nine models from the report
(grade comparators, arm-level SCNA, hot list, combinations, and the larger
clinical + SCNA + breakpoint models) across broad/restricted scopes and three
horizons:

```bash
python apply/menibreak_panel.py \
    --segments YOUR_segments.csv --outcomes YOUR_outcomes.csv \
    --clinical YOUR_clinical.csv --out results/full_ladder
```

`clinical.csv` needs only `sample_id`; every other column is optional and the tool
scores what it can, logging what it skipped. Useful columns:
`age, mitoses, gtr, prior_rt, adj_rt` (clinical models);
`who_grade, integrated_grade` (comparators);
`primary, rt_induced, angiomatous` (restricted subset). See
[`docs/INPUT_SPEC.md`](docs/INPUT_SPEC.md) §3.

Grade is used as a **comparator only** and never as a predictor inside a genomic
model. This is enforced in code — do not attempt to add it as a feature.

---

## Guardrails

**Do not:**

- Tune parameters, seeds, horizons, or `--k` to improve an AUC. Report what the
  default run produces. Sweeping and reporting the maximum invalidates the test.
- Modify anything in `reference/` — those files define the frozen panel. Changing
  them silently changes what is being validated.
- Commit patient data, SEG/IDAT files, sample-level predictions, or identifiers.
  `results/` run outputs are gitignored for this reason.
- Present the derivation-cohort figures as a target the lab should reach.
- Change `src/` or `analyses/` to make a result come out differently.

**Do:**

- Report the resolution verdict with every result.
- Report the arm-level SCNA comparator alongside the panel figure.
- Say clearly which scoring mode and panel set produced each number.
- Report a null result plainly if that is what you get.

## Where things are

| Path | What |
|---|---|
| `apply/validate_fixed_panel.py` | **The primary test** — fixed panel, no discovery |
| `apply/check_resolution.py` | The resolution gate |
| `apply/menibreak_panel.py` | Discovery mode + full clinical ladder |
| `reference/menibreak_panel_v1.csv` | The frozen nine/eleven windows — **do not edit** |
| `docs/INPUT_SPEC.md` | Full input schema |
| `docs/WHICH_MODEL.md` | Which of the three "panels" is which — read before quoting AUCs |
| `docs/QUICKSTART.md` | Long-form human walkthrough |
| `docs/LIMITATIONS.md` | What this does not support |
| `results/reference/` | Our published result tables |
| `RESULTS_MANIFEST.yaml` | Machine-readable index of every reported number |
| `tests/` | 88 tests; `python -m pytest tests/ -q` |

## Context you should have

These are **internally validated, single-institution** results. The lab's run is
the first genuine external test, which is the entire point of the exercise.
Research use only; no clinical claims.
