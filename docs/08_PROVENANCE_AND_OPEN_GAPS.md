# 8. Provenance and open gaps

## Recovered source

The historical AUC values were not inferred from the later MeniBreak document.
They are recorded in an older report created August 4–5, 2026:

`breakpoint_panel_word_tables_fixed.docx`

Its reproducibility table names original execution artifacts:

- `discover_hotlist.py` and `hotlist_windows.csv`
- `full_ladder.py` and `full_ladder.csv`
- `hotlist_clinical_metrics.py`
- `grade_discordance.py`
- `external_validate_hotlist.py`
- `compute_comparators.py` and `comparators.csv`

Only the generated Word report was found locally or in the August 10 backup.
The named scripts and CSV outputs were apparently retained in the generating
workspace rather than downloaded.

## What the source claims

The source describes a dynamic top-eight procedure: binary 20-bin arm-relative
breakpoint windows were ranked by chi-square association inside each training
fold and an elastic-net model with mixing parameter 0.25 scored held-out
patients. It does not specify a fixed `C`; that value was introduced by the
later scikit-learn reconstruction. Validation reportedly used five repeats
of patient-grouped 5-fold splitting. Restricted results transported the broad
out-of-fold predictions rather than refitting in the small restricted cohort.

The nine windows displayed later are the union of regions selected in at least
40% of folds at one or more horizons. They are not the model inputs for every
fold. The earlier fixed-nine interpretation was therefore incorrect.

## What is exact and executable

- Extraction of the nine stability windows matches the legacy source columns
  exactly across 1,273 matched tumors.
- `extract_breakpoint_candidates()` creates outcome-free binary arm-bin features.
- `workflows/01_dynamic_top8_validation/run_validation.py` keeps selection,
  scaling, and fitting inside outer training folds.
- WHO and integrated-grade comparator AUCs and strict cohort sizes reproduce the
  historical report to rounding.
- The locked 62-feature `mbpr-derivation-v1` software artifact remains a separate
  checksum-verifiable model.

## Provenance contradiction

The recovered report says the support filter left 745 candidates, while the
retained matrix yields 499 untyped candidates, 750 patient-supported combined
candidates, or 773 tumor-supported combined candidates. It states both top
`k=8` and top `k=5`; its “eight-window” stability union contains nine windows.
Several windows reportedly selected in most folds rank only 49th–102nd under
the stated chi-square procedure, and the named sets do not reproduce the AUCs.
Raw-SEG boundary variants, fold schemes, penalty conventions, and aggregation
rules did not resolve the contradiction.

Historical values are preserved in
[`historical_dynamic_top8_benchmark.csv`](../inst/extdata/historical_dynamic_top8_benchmark.csv)
with an explicit historical/unverified status. They document what the Word
report said, but must not be presented as reproduced results.
