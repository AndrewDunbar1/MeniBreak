# Script index and lineage

## Supported scoring software

`software/mbpr/` is the portable implementation. It normalizes external SEG files, extracts frozen SCNA and breakpoint features, applies the checksum-locked model, emits QC, and validates outcomes. This is the interface another lab should use first.

## Research analysis scripts

The scripts below are ordered by purpose. They preserve the analysis code but have had local workstation paths replaced by `MENIBREAK_PROJECT_ROOT` and, where needed, `MENIBREAK_SEGMENTS`.

| Script | Purpose | Status |
|---|---|---|
| `01_breakpoint_report_hotlist.R` | Nested selection of report-readable breakpoint regions, architecture thresholds, and paired patterns; survival attribution | Exploratory/internal validation |
| `02_compact_six_rule_hotlist.R` | Fixed-size six-rule approximation and comparison with WHO, prespecified integrated grade, and the full SEG model | Exploratory/internal validation |
| `03_geometry_route_bakeoff.R` | Tests arm-relative distribution, entropy, fragmentation, and route features without thresholded arm calls | Exploratory/internal validation |
| `04_scale_free_geometry_multihorizon.R` | Repeats the scale-free geometry comparison at 3, 5, and 7 years | Exploratory/internal validation |
| `05_triscut_inspired_hotlist.R` | TRISCUT-inspired segmentation geometry, iterative rule selection, ablations, and robustness checks | Method-development sensitivity analysis |
| `06_biscut_breakpoint_bakeoff.R` | BISCUT-inspired fold-local peak discovery versus the current breakpoint architecture model | Method-development sensitivity analysis |
| `07_confirmatory_comparators.R` | Harmonized WHO and prespecified integrated-grade comparisons on identical eligible tumors/folds | Comparator audit |
| `08_clinical_integrated_models.R` | Clinical + SCNA and clinical + SCNA + breakpoint models incorporating age, mitoses, GTR, prior RT, and adjuvant RT | Internal prognostic analysis; not the SEG-only deployable model |

## Required reconstructed analysis layout

Several scripts expect the original project-relative layout under `MENIBREAK_PROJECT_ROOT`:

- `analysis_folders/ml_breakpoint_survival_stratification_2026/ml_breakpoint_model_dataset_2026.csv`: derived feature table, one row per tumor.
- `software/mbpr/analysis/bakeoff_2x2/outcomes.csv`: `sample_id`, `patient_id`, recurrence time, and recurrence event.
- `software/mbpr/inst/extdata/arm_reference_GRCh37.csv`: chromosome-arm coordinates.
- `analysis_folders/Breakpoint_Report_Hotlist_2026/`: intermediate hotlist results required by downstream scripts.
- `software/vendor/BISCUT-py3/`: BISCUT reference assets for script 06.

These private derivation tables are deliberately absent. External teams should create equivalent derived tables from their own data or use the MBPR CLI for frozen scoring.

## Recommended execution order

1. Run MBPR `inspect` and resolve all QC failures.
2. Run locked MBPR `score` for the principal external validation.
3. If outcomes are available, run MBPR `validate` with patient IDs for clustered bootstrap.
4. Only then run the exploratory analysis families. Perform all feature selection inside outer-training folds; never select hotspots on the full external cohort and report performance on that same cohort.
