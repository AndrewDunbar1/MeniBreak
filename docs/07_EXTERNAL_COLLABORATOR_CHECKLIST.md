# 7. External collaborator checklist

## First pass: hotlist replication

1. Confirm build and coordinate convention.
2. Run `mbpr inspect`; resolve unknown events and incompatible overlaps.
3. Run `mbpr hotlist-extract`.
4. Review per-window prevalence and zero-count windows before outcomes.
5. Join outcomes by `sample_id`; supply `patient_id` for repeated tumors.
6. Run the fixed-panel validation script unchanged.
7. Return aggregate window associations, panel AUCs, denominators, and manifests.

## Second pass: full model

1. Run locked `mbpr score` on the same normalized input.
2. Inspect QC and out-of-distribution flags.
3. Run `mbpr validate` on the same fixed-horizon outcomes.
4. Compare frozen full and frozen SCNA-only predictions on identical tumors.

## Required reporting

State cohort, genome build, assay/caller, event mapping, coordinate convention,
preprocessing lineage, horizon denominators/events, patient grouping, model ID,
and whether coefficients were refit. Do not share PHI, raw institutional data,
or sample-level results through this repository.
