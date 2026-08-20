# Model card

## Intended use

Research prediction and external validation of meningioma recurrence from
segment-level copy-number data. Five-year recurrence is primary; three- and
seven-year results are secondary.

## Status

`LOCKED`: `mbpr-derivation-v1`.

Exact legacy reproduction passed for 57 shared legacy features across 1,273
tumors within a tolerance of `1e-12`. The prespecified repeated
patient-grouped nested five-fold 2×2 bakeoff selected legacy preprocessing with
a 30% arm threshold. This is a documented protocol deviation: the best V6
candidate was more than 0.01 AUC below the best candidate.

## Inputs and features

The portable feature set contains arm-level loss/gain calls, loss/gain/altered
arm burden, prespecified breakpoint burden and spatial-distribution features,
and four frozen hotspot-routing features. Assay-derived NF2 status is excluded;
22q loss is retained as a SEG-derived feature.

## Outputs

Once locked, each horizon emits a raw score, derivation-cohort percentile, and
internally calibrated provisional probability. No clinical high/low treatment
label is emitted.

## Internal validation

| Horizon | Samples | Events | Breakpoint AUC | SCNA-only AUC | AUC gain |
|---|---:|---:|---:|---:|---:|
| 3 years | 904 | 204 | 0.779 | 0.772 | 0.007 |
| 5 years | 674 | 264 | 0.788 | 0.753 | 0.035 |
| 7 years | 514 | 296 | 0.780 | 0.747 | 0.033 |

Values are means from 10 repeats of patient-grouped nested five-fold
validation. At five years, breakpoint augmentation improved AUC in all repeats
and improved mean Brier score from 0.191 to 0.176. Three-year discrimination
gain was small and its uncalibrated Brier score was not better, reinforcing its
secondary status.

Fold-derived hotspot windows had median overlap/Jaccard 0.60 with the frozen
windows; 38% were exact matches. The frozen derivation windows are therefore
versioned model assets, not claimed as universal biological boundaries.

## Limitations

- Absolute risk is provisional until independently calibrated on an external,
  outcome-linked segment-level cohort.
- Caller and platform shifts may materially affect features and scores.
- Five-year breakpoint gain remained positive in both derivation data-source
  strata, but the smaller legacy-remapped subset had less precise performance.
- GRCh38 uses its own arm coordinates and normalized arm-relative hotspots; no
  silent breakpoint liftover is performed.
- GRCh38 scoring is technically supported but was not the derivation build and
  requires independent validation.
- The software is not for clinical diagnosis or treatment decisions.
