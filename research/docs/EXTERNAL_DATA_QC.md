# External dataset QC and validation checklist

## Before scoring

- Confirm the genome build from laboratory metadata. Never infer GRCh37 versus GRCh38 from coordinate ranges.
- Record assay/platform, segmentation software and version, normalization method, caller thresholds, coordinate convention, and whether segments are event-called or continuous.
- Preserve the original file read-only and calculate its SHA-256 checksum.
- Use stable, deidentified sample IDs. Maintain any re-identification key outside this project and outside version control.
- Confirm one genome build and one coordinate convention per run.
- Check required fields: sample ID, chromosome, start, end, and event; continuous formats additionally require a caller profile.
- Normalize chromosome labels and event vocabulary deterministically. Unknown labels must fail rather than being silently discarded.

## Segment-level checks

- Start and end must be finite integers, nonnegative under the declared convention, and within chromosome bounds; end must exceed start.
- Flag exact duplicates, nested segments, conflicting overlaps, gaps, arm-spanning segments, and samples with no usable autosomal segments.
- Report total segments; loss, gain, LOH, and allelic-imbalance counts; chromosome coverage; median segment length; focal-segment fraction; and breakpoint count per sample.
- Treat the first/last segment boundary cautiously: array coverage edges are not necessarily biological breakpoints. Apply the same boundary convention as the frozen model.
- Do not count adjacent segments with identical states as separate biological events; merge them only under a prespecified deterministic rule.
- Confirm the handling of centromeres, acrocentric short arms, sex chromosomes, LOH, and allelic imbalance against `software/mbpr/INPUT_SPEC.md`.

## Cross-platform comparability

- Event-called SEG is the validated v0.1.0 pathway. Continuous log2 data require a named caller profile with frozen thresholds.
- For methylation-derived segments, document the preprocessing pipeline, probe filtering, reference set, segmentation parameters, minimum probes/segment, and absolute log2 call threshold.
- Quantify out-of-distribution behavior versus the derivation reference: segment count, breakpoint count, altered-arm burden, missing features, feature range violations, and prediction percentile distribution.
- Run sensitivity analyses across scientifically plausible caller thresholds without choosing the best threshold based on outcome performance.
- Hotspot rules depend on breakpoint localization resolution. If median segment-boundary uncertainty or array spacing is broader than a hotspot window, classify the hotspot result as unresolved rather than absent.
- GRCh38 uses its own arm reference and arm-relative routing. Do not silently liftover individual breakpoints.

## Outcome validation

- Required outcome fields are deidentified `sample_id`, recurrence time in years, and recurrence event; add `patient_id` when more than one tumor may belong to a patient.
- Construct 3-, 5-, and 7-year endpoints before looking at scores, using the strict rule stated in `docs/ANALYSIS_METHODS_AND_RESULTS.md`.
- Report eligible tumors and events separately at every horizon.
- Use patient-clustered bootstrap confidence intervals when repeated tumors exist.
- Report AUC, average precision, Brier score, calibration intercept/slope, event-capture curves, and decision curves. Do not report AUC alone.
- Compare models on exactly the same tumors. Use paired bootstrap draws for AUC differences.
- Do not recalibrate or refit the frozen score before first reporting its unmodified external performance. Any subsequent recalibration must be labeled and cross-validated or assessed in a separate cohort.

## Stop conditions

Block scoring or interpretation for mixed/unknown genome builds, unresolved caller thresholds, incompatible overlapping segments, systematic chromosome loss, extreme feature-distribution shifts, insufficient outcome events, or sample-ID linkage discrepancies. Document failures in the run manifest rather than manually overriding them.
