# Analysis methods and result lineage

## Endpoint and validation framework

Strict fixed-horizon recurrence was used at 3, 5, and 7 years. A recurrence at or before the horizon is an event; recurrence-free follow-up through the horizon is a non-event; earlier censoring is excluded. Resampling is patient-grouped so tumors from the same patient cannot occur in both training and test folds. Model tuning and feature selection occur within training data only. Reported internal AUCs come from held-out out-of-fold predictions, with patient-level bootstrap confidence intervals where generated.

## Locked SEG software model

The production model uses portable event-called SEG inputs, 30% arm calling selected by the prespecified bakeoff, continuous SCNA loss/gain/altered-arm information, and frozen breakpoint-architecture/hotspot features. The model artifact freezes preprocessing lineage, feature order, imputation, centering/scaling, coefficients, calibration, seeds, and checksums. Five-year recurrence is primary; 3- and 7-year results are secondary. The frozen model must not be refit on an external scoring cohort.

## Report-readable hotlist

The hotlist analysis asked whether a small set of easily described breakpoint patterns could add value to WHO or the prespecified integrated-grade label when the full SEG model was unavailable. Candidate features included arm-relative hotspot bins, loss/gain-specific breakpoint locations, architecture counts, and selected pairs. Screening and elastic-net selection were nested within patient-grouped outer folds.

The compact six-rule analysis at 5 years produced AUC 0.729 with WHO plus hotlist versus 0.656 for WHO alone, and 0.763 with integrated grade plus hotlist versus 0.735 for integrated grade alone. The full clinical + SCNA + breakpoint reference remained higher at 0.844. These values are internal and do not establish transportability.

The reduced reportable list in `aggregate_results/breakpoint_hotlist_reduced_reportable.csv` contains evidence tiers. Only Tier 1 findings should be described as high-signal internally; Tier 2 is supported; Tier 3 remains exploratory. P-values and q-values are descriptive after a complex selection procedure and are not substitutes for external validation.

## Geometry/route features

Scale-free features summarized where breakpoints occurred relative to telomere and centromere, how widely they were distributed across arms/chromosomes, entropy of their distribution, clustering, fragmentation, focal-segment fraction, and loss/gain balance. These were designed to capture genome geometry without requiring an exact platform-specific base-pair hotspot.

Internally, geometry alone yielded AUCs 0.764, 0.769, and 0.755 at 3, 5, and 7 years. Adding geometry to WHO yielded 0.773, 0.772, and 0.757. Adding it to integrated grade yielded 0.779, 0.782, and 0.768. These are exploratory nested-validation estimates.

## TRISCUT- and BISCUT-inspired analyses

The TRISCUT-inspired analysis tested iterative geometry/routing rules, background normalization, CNA-only ablations, and stability guardrails. It was an engineering analogue, not a faithful execution of the complete original TRISCUT pipeline. At 5 years, its best explored eight-rule combinations reached 0.771 with WHO and 0.793 with integrated grade, but model search across variants makes these hypothesis-generating estimates.

The BISCUT-inspired analysis learned breakpoint peak features inside each outer-training fold and compared them on the same held-out patients used for the current model. It did not materially improve the current breakpoint model. This negative result is important: more elaborate peak discovery did not provide evidence that the current breakpoint benefit had been underestimated.

## Clinical-integrated internal model

The separate clinical-integrated models added age, log-transformed mitotic count, gross-total resection, prior radiotherapy, and adjuvant radiotherapy to SCNA, then evaluated the incremental addition of breakpoint architecture. These models can be useful for internal prognosis but are not the deployable SEG-only product because external labs may not encode these covariates identically and some variables reflect treatment decisions.

## Grade comparators

WHO grade and the prespecified integrated-grade variable were treated as recorded cohort labels in the final comparator analysis; integrated grade was not re-created from a proxy formula. Comparator models used the same eligible tumors and resampling splits as the molecular models. In the exact restricted cohort (primary, GTR, no prior RT, no adjuvant RT, non-angiomatous), WHO AUCs were 0.686, 0.647, and 0.622 at 3, 5, and 7 years, respectively.

## Interpretation hierarchy

1. External use: locked MBPR score is the primary reproducible product.
2. Internal paper claim: breakpoint augmentation is evaluated against the matched SCNA-containing baseline under nested patient-grouped validation.
3. Hotlist and geometry results: report-readable and mechanistic research layers requiring external confirmation.
4. Clinical-integrated models: higher-information internal analyses, not SEG-only deployment claims.
