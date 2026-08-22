# MeniBreak

MeniBreak is a research toolkit for testing whether the location and geometry
of copy-number breakpoints improve meningioma recurrence prediction. The repo
is intentionally ordered like the MeniBreak scientific narrative:

1. reproduce the outcome-linked dynamic top-eight procedure;
2. evaluate the stronger scale-free breakpoint-geometry definition;
3. test the nine-window stability shortlist in an external cohort;
4. compare breakpoint models with WHO grade and prespecified integrated grade;
5. measure what breakpoints add beyond SCNA;
6. examine grade-discordant reclassification; and
7. validate the full locked SEG-to-score model.

**Start at [`START_HERE.md`](START_HERE.md).** It identifies the correct command,
inputs, output, model status, and interpretation for each stage. Machine readers
should use [`MODEL_INDEX.yml`](MODEL_INDEX.yml) as the authoritative registry.

The recovered historical six-AUC table is not execution-reproduced and is
internally inconsistent with its stated candidate count, top-k description,
and stability table. The provenance audit is in
[`docs/09_INTERNAL_REPRODUCTION_AUDIT.md`](docs/09_INTERNAL_REPRODUCTION_AUDIT.md).

## Outcome-linked dynamic top-eight workflow

```sh
./mbpr top8-candidates --segments cohort.seg --format auto \
  --genome-build GRCh37 --output-dir top8_candidates

python3 workflows/01_dynamic_top8_validation/run_validation.py \
  --features top8_candidates/breakpoint_candidates.csv \
  --outcomes outcomes.csv --output-dir top8_validation
```

The top eight are selected anew inside every training fold. This is a faithful
implementation of the narrative recipe, but it does not reproduce the old Word
report's complete AUC table and is not a fixed score for outcome-free samples.

## Strongest freshly reproduced breakpoint-only analysis

The scale-free geometry model uses 56 fixed, normalized summaries of breakpoint
dispersion, concentration, arm-relative position, and loss/gain route shape.
Under ten repeats of nested patient-grouped validation, broad-cohort AUC was
0.764, 0.769, and 0.755 at 3, 5, and 7 years. Restricted-transport AUC was
0.715, 0.747, and 0.722. This is the strongest freshly reproduced
breakpoint-only internal definition, but remains exploratory rather than a
locked external predictor. See
[`docs/10_BREAKPOINT_ONLY_MODEL_COMPARISON.md`](docs/10_BREAKPOINT_ONLY_MODEL_COMPARISON.md).

## Fastest Heidelberg workflow: stability shortlist

```sh
./mbpr hotlist-extract --segments cohort.seg --format auto \
  --genome-build GRCh37 --output-dir hotlist_results

Rscript workflows/01_hotlist_validation/run_validation.R \
  hotlist_results/hotlist_features.csv outcomes.csv hotlist_validation
```

This tests the nine regions that were selected most consistently by the dynamic
procedure. They should initially be reported individually. The dynamic model's
AUCs must not be assigned to their union.

## Full locked model

```sh
./mbpr inspect --segments cohort.seg --format auto --genome-build GRCh37
./mbpr score --segments cohort.seg --format auto \
  --genome-build GRCh37 --output-dir locked_scores
./mbpr validate --predictions locked_scores/predictions.csv \
  --outcomes outcomes.csv --output-dir locked_validation
```

The locked model is `mbpr-derivation-v1`: 44 SCNA features plus 18 breakpoint
architecture features. Five-year recurrence is primary; 3- and 7-year outputs
are secondary. Its probabilities are internally calibrated and provisional.

## Repository map

- [`docs/`](docs/): the scientific story and exact analysis contracts, in order.
- [`workflows/01_hotlist_validation/`](workflows/01_hotlist_validation/): the
  self-contained first external-validation workflow.
- [`R/`](R/), [`inst/`](inst/), [`mbpr`](mbpr): supported package and CLI.
- [`archive/method-development/`](archive/method-development/): BISCUT,
  TRISCUT, geometry, and other method-development analyses. These are retained
  for auditability but are not part of the shortest collaborator path.

No institutional raw data, sample-level outcomes, linkage keys, or PHI are
included. The software is research-use-only and is not for diagnosis, grading,
treatment selection, or patient management.

Before publishing changes, run:

```sh
R_PROFILE_USER=/dev/null Rscript tests/smoke.R
./tools/audit_no_phi.sh
```

Code is released under Apache-2.0. See [`CITATION.cff`](CITATION.cff).
