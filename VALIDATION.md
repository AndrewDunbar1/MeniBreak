# Validation record

## Frozen gates

- Legacy reproduction: passed, 1,273 overlapping tumors, 57 features,
  tolerance `1e-12`.
- Primary candidate bakeoff: 10 repeats of patient-grouped nested five-fold
  validation.
- Selected candidate: legacy preprocessing, 30% arm threshold.
- Primary five-year AUC: 0.788 versus 0.753 for SCNA-only.
- Five-year AUC gain: positive in 10/10 release-audit repeats.
- Five-year Brier score: 0.176 versus 0.191 for SCNA-only.
- Frozen feature count: 62.
- Artifact SHA-256:
  `87bf64eded79f44b0ac805f1cb490f760b7bec3699460183744fdce21ee3093f`.

## Calibration

Calibration uses averaged repeated held-out predictions and is labeled internal
and provisional. On averaged out-of-fold predictions, calibration intercepts
were 0.021, 0.009, and 0.003 and slopes were 1.023, 1.018, and 1.022 for the
three-, five-, and seven-year endpoints, respectively.

## External validation requirement

External laboratories should run `mbpr score` locally and then use
`mbpr validate` with outcome-linked data. Absolute risk should not be treated as
externally calibrated until independent results are available.
