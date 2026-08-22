# 11. External validation ladder

The repository intentionally includes complementary models rather than forcing
one model to answer every question.

| Stage | Model | Minimum input | External action | Status |
|---|---|---|---|---|
| 1 | Nine-window stability shortlist | SEG + build | Extract nine report-readable regions; test each region and a clearly labeled composite | Simplest biological replication |
| 2 | Scale-free breakpoint geometry | Complete SEG + build + outcomes | Extract 56 fixed features; refit only inside nested validation folds | Strongest reproduced breakpoint-only recipe |
| 3 | WHO and prespecified integrated grade | Outcomes + recorded labels | Use labels directly on the identical eligible tumors | Exactly reproduced comparators |
| 4 | Clinical + SCNA and clinical + SCNA + breakpoint architecture | SEG + five clinical variables + outcomes | Repeat the documented patient-grouped nested validation | Reproduced internal higher-information analyses |
| 5 | Locked `mbpr-derivation-v1` | SEG + build | Apply frozen coefficients without refitting, then validate if outcomes exist | Deployable SEG-only software model |

## Why both shortlist and geometry are retained

The shortlist is interpretable enough to check manually and is the easiest
first collaboration analysis. Geometry is operationally simple when a complete
SEG file is available and retains more breakpoint-only discrimination, but its
56 summaries are not a paper-report checklist. The two analyses therefore test
the same biological story at different levels of information.

The SEG-native geometry extractor was verified sample-by-sample against the
retained institutional feature table: all 56 feature columns were numerically
identical for all 1,273 overlapping derivation samples (tolerance `1e-12`).

## Which Word-document values are executable benchmarks

The recorded WHO and prespecified integrated-grade AUCs reproduced to rounding.
The clinical + SCNA and clinical + SCNA + breakpoint-architecture analyses also
have retained source code and fresh reproducible nested-validation outputs, but
their fresh values replace—not retroactively validate—the older Word-table
values. The historical nine-window AUC table remains provenance-only because no
retained execution reproduces all six values under one coherent definition.

Machine-readable fresh benchmarks are in
`inst/extdata/internal_reproduced_model_benchmark_20260822.csv`. Historical
values are preserved separately in `inst/extdata/historical_dynamic_top8_benchmark.csv`.
