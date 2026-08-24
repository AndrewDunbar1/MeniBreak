# Quickstart — run the panel on your own cohort

The whole workflow for an external group. It needs a segmented copy-number file
and recurrence outcomes; clinical data is optional but unlocks the full ladder.
No institutional data is required and nothing leaves your machine.

## 0. Install

```bash
git clone https://github.com/AndrewDunbar1/MeniBreak.git
cd MeniBreak
pip install -r requirements.txt     # pandas, numpy, scipy, scikit-learn, pyyaml
```

## 1. Confirm it runs (bundled synthetic cohort)

```bash
python data/example/generate_example.py --n 300 --seed 7 --outdir data/example

python apply/check_resolution.py \
    --segments data/example/synthetic_segments.csv

python apply/validate_fixed_panel.py \
    --segments data/example/synthetic_segments.csv \
    --outcomes data/example/synthetic_outcomes.csv \
    --out results/example_fixed

python apply/menibreak_panel.py \
    --segments data/example/synthetic_segments.csv \
    --outcomes data/example/synthetic_outcomes.csv \
    --out results/example
```

Expected on this synthetic cohort (seed 7): the preflight reports a median event
of ~1.2 Mb and `SUITABLE`; the fixed panel reports an AUC ladder that **rises
with horizon**, roughly 0.61 → 0.64 → 0.72; the discovery run recovers the
planted chr1p window at **selection frequency 1.0** at all three horizons.

Exact figures shift with `--seed` and `--n`. The pattern to look for is the rise
with horizon and the recovery of 1p, not the third decimal.

The example is synthetic and its background breakpoints are placed uniformly at
random, so its event burden (~48 breakpoints/sample) sits below a real cohort's
(institutional: ~97). It is a smoke test, not a performance benchmark.

Run the test suite (88 tests, ~1 min):

```bash
pip install -r requirements-dev.txt && python -m pytest tests/ -q
```

## 2. Check your resolution before anything else

```bash
python apply/check_resolution.py --segments YOUR_segments.csv
```

Any SEG file will *run*. The question is whether the assay behind it localized
breakpoints finely enough for window-level features to mean anything. The panel
scores breakpoints in bins of 5% of a chromosome arm — about **2.8 Mb**.

| Platform | Median event | vs window | Panel |
|---|---|---|---|
| Agilent aCGH (training cohort) | 1.26 Mb | 0.45× | works |
| conumee on EPIC methylation | 5.10 Mb | 1.80× | lost to arm-level SCNA |

`SUITABLE` → proceed. `MARGINAL` → run the panel *and* the arm-level reference
and compare on your own outcomes. `TOO COARSE` → use arm-level SCNA and treat
window results as exploratory. **Methylation-derived copy number is the common
failure case.** Calibrated on two platforms, so it is a decision aid, not a
validated instrument.

## 3. Format your data

Full schema: [`INPUT_SPEC.md`](INPUT_SPEC.md). The minimum is:

- `segments.csv`: `sample_id, chromosome, start, end, event` (+ `patient_id` if a
  patient can have more than one tumour).
- `outcomes.csv`: `sample_id, recurrence_time, recurrence_event`
  (+ `patient_id`). **`recurrence_time` must be in YEARS** — months are rejected.

Two *optional* segment columns materially change what you get:

| Include | Why |
|---|---|
| `segment_mean` (log2 ratio) | amplitude features in the enriched model |
| `loh` / `allelic_imbalance` rows | allelic features — 63.5% of institutional segment rows |

Without them the enriched model still runs but is hollow; the tool logs which
blocks your data supports and warns on each missing one.

If your caller uses different event words, either recode them to `loss`/`gain`/…
or add them to `EVENT_SYNONYMS` in `apply/menibreak_panel.py`.

## 4. Run on your cohort

### 4a. The validation run — fixed published panel (start here)

Tests **our** panel. Windows come from `reference/menibreak_panel_v1.csv`;
nothing is learned from your data and with the default `--score count` no
parameters are fitted at all.

```bash
python apply/validate_fixed_panel.py \
    --segments YOUR_segments.csv --outcomes YOUR_outcomes.csv \
    --build GRCh37 --horizons 3 5 7 --out results/your_cohort
```

| File | Contents |
|---|---|
| `fixed_panel_performance.csv` | AUC of the fixed panel vs an arm-level SCNA reference, per horizon. |
| `panel_window_occupancy.csv` | How often each fixed window fires. **Read this first** — a window firing in zero tumours cannot contribute, and several silent windows usually means a resolution mismatch rather than absent signal. |

Variants: `--score refit` (fixed windows, coefficients refit under
patient-grouped CV — fairer if your platform or population differs) and
`--panel-set rule` (11 windows instead of 9; see [`THE_PANEL.md`](THE_PANEL.md)).

### 4b. The method run — rediscover windows in your data

```bash
python apply/menibreak_panel.py \
    --segments YOUR_segments.csv --outcomes YOUR_outcomes.csv \
    --build GRCh37 --horizons 3 5 7 --k 8 --repeats 5 \
    --out results/your_cohort_discovery
```

Selects windows inside your own training folds, so it answers "does the *method*
work here?", not "does the *published panel* work here?".

| File | Contents |
|---|---|
| `performance.csv` | Patient-grouped out-of-fold AUC per horizon for `Arm-level SCNA`, `8-window panel alone`, `SCNA + 8-window panel` and `SCNA + enriched breakpoints`. |
| `hotlist_windows.csv` | Windows your cohort selects across folds, with `selection_frequency` and **genomic coordinates** in your build. |
| `window_features.csv` | Per-tumour occupancy of every arm-relative window. |

### 4c. Add clinical data for the full ladder

```bash
python apply/menibreak_panel.py \
    --segments YOUR_segments.csv --outcomes YOUR_outcomes.csv \
    --clinical YOUR_clinical.csv --out results/your_cohort_full
```

This adds the grade comparators, the three `Clinical +` models and the
**Restricted** scope — the complete nine-model report on your cohort. Only
`sample_id` is required in that file; partial input scores what it can and logs
what it skipped. Grade is read as a comparator only and never enters a genomic
model as a predictor.

## 5. Read the QC line first

Check the printed `median breakpoints/sample` and `mean window occupancy`. The
panel was defined on data with a **median ~97 breakpoints/sample and 27.8%
window occupancy**. Below 40 the tool warns: the ~2.8 Mb windows fire too rarely
to be comparable, and a weak result should not be read as a test of the panel.
See [`LIMITATIONS.md`](LIMITATIONS.md).

## 6. Interpret

- Does `8-window panel alone` approach or beat integrated grade in your cohort,
  and does its advantage grow at 5–7 years? That is the institutional finding.
- Do your top `hotlist_windows.csv` entries include the chr1p pericentromeric
  region (chr1:118.75–125.00 Mb, GRCh37)? That is the panel's most stable
  feature ([`THE_PANEL.md`](THE_PANEL.md)).
- `SCNA + panel` vs `SCNA` alone is the incremental-value test.
- For the **fixed panel**, compare against the derivation-cohort anchor in
  `results/reference/fixed_panel_internal_reference.csv` (0.690 / 0.772 / 0.785 at
  3/5/7 y) — but read [`WHICH_MODEL.md`](WHICH_MODEL.md) first: that anchor is
  optimistic because the windows were selected on that cohort, so scoring below
  it is expected. The unbiased comparison is against the arm-level SCNA figure
  computed on *your* folds.
- `SCNA + enriched breakpoints` is the strongest SEG-only model in the report
  (0.819 at 5 y) — worth reporting alongside the parsimonious panel.

A null result is a real and publishable outcome. The institutional numbers are
internally validated only; the single external attempt to date failed on
resolution grounds, and what you are running *is* the external validation.

## Is this path faithful to the report?

Yes, and it is checked rather than asserted.
`analyses/07_verify_external_path.py` rewrites the institutional cohort into
these generic formats and runs it through `apply/`: all **9 models × 2 scopes ×
3 horizons** reproduce the published AUCs to within **0.0031**, with exact
cohort sizes and event counts. That script needs the institutional data and is
not runnable outside — it exists to prove the path you *do* use is faithful.
