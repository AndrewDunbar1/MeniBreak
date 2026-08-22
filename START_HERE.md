# Start here

The historical Word report describes a dynamic procedure that selected windows
inside training folds, but its AUC table, candidate count, top-k statements, and
stability table cannot be reconciled into one retained execution. Read
[`docs/09_INTERNAL_REPRODUCTION_AUDIT.md`](docs/09_INTERNAL_REPRODUCTION_AUDIT.md)
before quoting or attempting to reproduce that table.

## Choose the question before choosing a script

| Scientific question | Correct entry point | What is fixed? | What may be fitted in the external cohort? |
|---|---|---|---|
| Does the dynamic top-eight procedure generalize in an outcome-linked cohort? | `mbpr top8-candidates`, then `workflows/01_dynamic_top8_validation/run_validation.py` | Candidate construction, top-k, association test, estimator, endpoints and resampling | Eight windows and coefficients, inside each training fold only |
| Does normalized breakpoint geometry outperform a short location panel? | `archive/method-development/analysis/03_geometry_route_bakeoff.R`, then `04_scale_free_geometry_multihorizon.R` | 56 scale-free shape features and nested-validation contract | Coefficients and penalty, inside training folds only |
| Do the nine stable report-readable regions reproduce individually? | `mbpr hotlist-extract`, then `workflows/01_hotlist_validation/run_validation.R` | Window locations and boundary rules | Secondary composite sensitivity only; do not assign it the dynamic model's AUCs |
| Does a fully specified SEG model transport without refitting? | `mbpr score`, then `mbpr validate` | Preprocessing, 62 features, coefficients, scaling, calibration, percentiles | Nothing |
| What does breakpoint architecture add beyond SCNA internally? | `docs/04_INCREMENTAL_VALUE_BEYOND_SCNA.md` | Feature families and comparison contract | Models within patient-grouped nested validation |
| Which elevated-grade nonrecurrent or low-grade recurrent tumors are reclassified? | `docs/05_GRADE_DISCORDANT_RECLASSIFICATION.md` | Eligibility and risk-threshold contract | Nothing outside prespecified cross-validation |
| How were alternative breakpoint representations explored? | `archive/method-development/` | Varies by experiment | Method-development only |

Do not call the dynamic top-eight procedure a “locked predictor”: its selected
windows can vary by fold and cohort. Do not assign the unverified historical
AUCs to either the dynamic implementation or the nine-window stability union.
Do not call the full locked model the “hotlist.” It uses a larger 62-feature
SEG representation.

## Required data

For feature extraction:

`sample_id`, `chromosome`, `start`, `end`, `event`

Genome build is mandatory. Event-called data are the primary supported format.
Continuous log2 SEG data require an explicit caller profile.

For fixed-horizon validation:

`sample_id`, `recurrence_time`, `recurrence_event`

Add `patient_id` when one patient can contribute multiple tumors. Optional
comparison columns are `who_grade` and `integrated_grade`; integrated grade must
be the prespecified recorded label, not a newly constructed proxy.

## Endpoint rule

At horizon `h`, recurrence on or before `h` years is an event; event-free
follow-up through `h` is a non-event; earlier censoring is excluded. Every model
comparison must use the identical eligible tumors and patient-grouped folds.

Continue through the numbered files in [`docs/`](docs/).

Read [`docs/08_PROVENANCE_AND_OPEN_GAPS.md`](docs/08_PROVENANCE_AND_OPEN_GAPS.md)
before quoting the document's historical “panel” AUCs.
