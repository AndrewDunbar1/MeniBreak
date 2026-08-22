# 10. Breakpoint-only model comparison

## Purpose

This comparison asks whether another legitimate breakpoint representation
outperformed the reconstructed location shortlist. Candidate definitions were
required to keep patient grouping, strict fixed-horizon endpoints, and all
model fitting inside training folds.

## Results

| Definition | Cohort | 3 years | 5 years | 7 years | Interpretation |
|---|---|---:|---:|---:|---|
| Literal dynamic top-eight | Broad | 0.724 | 0.707 | 0.741 | Fresh executable reconstruction |
| Literal dynamic top-eight | Restricted transport | 0.686 | 0.675 | 0.733 | Same broad OOF scores subset to restricted tumors |
| Native-R dynamic top-eight | Broad | 0.724 | 0.716 | 0.745 | Closest method-faithful reconstruction |
| Native-R dynamic top-eight | Restricted transport | 0.703 | 0.711 | 0.757 | Same broad OOF scores subset to restricted tumors |
| Scale-free geometry | Broad | **0.764** | **0.769** | **0.755** | Strongest freshly reproduced breakpoint-only definition |
| Scale-free geometry | Restricted transport | **0.715** | **0.747** | 0.722 | Sparse-event transport estimate |

The geometry model is higher than both reproducible shortlist implementations
at every broad-cohort horizon. It rises slightly from 3 to 5 years and remains
near 0.76 at 7 years, rather than showing the perfectly increasing trajectory
in the unverified Word table.

## Why geometry is a legitimate alternative

The 56 inputs are fixed transformations of breakpoint positions relative to
their chromosome arms. They describe dispersion, concentration, radial
position, entropy, and loss/gain route asymmetry. They exclude clinical
covariates, treatment variables, mitotic count, and thresholded arm-level SCNA
calls. Model coefficients and `lambda.1se` are estimated within nested
patient-grouped training folds.

The direct lineage is:

- `archive/method-development/analysis/03_geometry_route_bakeoff.R`
- `archive/method-development/analysis/04_scale_free_geometry_multihorizon.R`

## Trend sensitivity

A count-aware top-eight sensitivity analysis produced a more visibly increasing
mean-fold pattern, but it was not promoted because it depended on count-valued
segment boundaries, an F-test selector, and fold-level AUC averaging. Those
choices differ from the report-readable binary shortlist and are more sensitive
to assay segmentation density. Selecting it because its curve rises would add
post hoc model-selection bias.

## Status

Scale-free geometry is an exploratory internal model recipe, not the locked
`mbpr-derivation-v1` software score. It should be externally tested without
changing its 56-feature definition or validation contract.
