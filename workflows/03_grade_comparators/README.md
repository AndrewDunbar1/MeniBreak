# Recorded-grade comparator validation

This workflow reproduces the two Word-document comparators that were recovered
exactly to rounding: recorded WHO grade and the **prespecified recorded**
integrated-grade label. Integrated grade must never be recalculated from a proxy.

Required outcome columns are:

`sample_id`, `recurrence_time`, `recurrence_event`, `who_grade`, `integrated_grade`

Add `patient_id` when patients contribute multiple tumors. An optional binary
`restricted_eligible` column produces the restricted-cohort table using a cohort
flag constructed before running this script. Its definition must be reported.

```sh
Rscript workflows/03_grade_comparators/run_validation.R \
  outcomes_with_recorded_grades.csv grade_validation
```

These ordinal labels are used directly as discrimination scores. They are not
fit, tuned, or transformed. The strict 3-/5-/7-year endpoint rule is identical
to the genomic workflows.
