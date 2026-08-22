# Stability-shortlist replication sensitivity

This is the recommended first analysis for an external collaborator. It answers:
do the nine stable 5%-of-arm breakpoint windows show prognostic signal in a new
cohort?

## Inputs

1. `hotlist_features.csv` from `mbpr hotlist-extract`.
2. An outcomes CSV with `sample_id`, `recurrence_time`, and
   `recurrence_event`; add `patient_id` for repeated tumors.
3. Optional comparators: `who_grade`, `integrated_grade`.
4. Optional clinical augmentation: `age`, `mitoses`, `gtr`, `prior_rt`,
   `adj_rt`.

## Run

```sh
Rscript workflows/01_hotlist_validation/run_validation.R \
  hotlist_results/hotlist_features.csv outcomes.csv output_directory
```

The script uses the nine-window stability union, creates strict 3-/5-/7-year
endpoints, and performs repeated patient-grouped nested validation. It emits
aggregate results only. It does not rediscover windows, use the test fold for
tuning, or recalculate integrated grade.

## Outputs

- `panel_auc_summary.csv`: mean nested AUC by model and horizon.
- `panel_auc_by_repetition.csv`: resampling stability.
- `individual_window_associations.csv`: descriptive per-window direction,
  standardized odds ratio, and prevalence.
- `validation_manifest.json`: input checksums, fixed settings, and QC exclusions.

Samples with incompatible overlapping loss/gain segments are marked
`FAIL_INCOMPATIBLE_OVERLAP` by extraction and excluded from outcome validation;
one problematic sample does not block clean samples in the same cohort.

Individual-window results are the intended first Heidelberg analysis. The
composite refit is a secondary sensitivity analysis: it is neither the dynamic
top-eight procedure nor a validation of frozen coefficients.
