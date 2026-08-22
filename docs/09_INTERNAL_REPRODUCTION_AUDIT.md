# Internal reproduction audit: historical shortlist

## Bottom line

The six historical “eight-window” AUCs are recoverable from the August Word
report, but not from a retained execution. They should be treated as
**historical/unverified**, not as results reproduced by this repository.

The report contains mutually incompatible specifications:

- top `k=8` in the methods but `k=5` in the stability-table caption;
- nine unique displayed windows despite the “eight-window” label;
- 745 candidates, whereas the retained matrix yields 499 untyped, 750
  patient-supported combined, or 773 tumor-supported combined candidates; and
- high reported stability for windows that rank approximately 49th–102nd under
  the stated chi-square selector.

## Fresh literal implementation

The executable Python workflow uses binary untyped arm bins, minimum support
20, foldwise chi-square top-eight selection, five repeated patient-grouped
folds, standardization, and elastic-net logistic regression
(`l1_ratio=0.25`, `C=0.1`). It reproduced the strict cohorts and direct grade
comparators, but not the historical panel AUCs.

| Cohort | 3 years | 5 years | 7 years |
|---|---:|---:|---:|
| Broad | 0.724 | 0.707 | 0.741 |
| Restricted transport | 0.686 | 0.675 | 0.733 |

## Closest faithful native-R interpretation

Because the source specifies elastic-net mixing but not fixed penalty strength,
we also tested the contemporaneous R convention: combined untyped/loss/gain
candidates, patient-level support, R-generated patient folds, foldwise
chi-square top 8, and nested `cv.glmnet` using `lambda.1se`. AUC was computed
from predictions averaged across five repeats.

| Cohort | 3 years | 5 years | 7 years |
|---|---:|---:|---:|
| Broad | 0.724 | 0.716 | 0.745 |
| Restricted transport | 0.703 | 0.711 | 0.757 |

This is the closest method-faithful reconstruction, not an exact reproduction.

## Stronger legitimate breakpoint-only definition

The separately prespecified scale-free geometry recipe performed better than
either shortlist reconstruction in the broad cohort: 0.764, 0.769, and 0.755.
It captures normalized breakpoint distribution and route shape rather than only
eight locations. Restricted-transport AUCs were 0.715, 0.747, and 0.722. The
7-year restricted decrease is a reason not to select a definition solely for a
preferred trend.

See `docs/10_BREAKPOINT_ONLY_MODEL_COMPARISON.md` and
`inst/extdata/internal_breakpoint_only_benchmark_20260822.csv`.

## Reporting rule

- Describe the six Word-report values only as historical/unverified provenance.
- Use fresh execution results when claiming reproducibility.
- Call dynamic top-eight a fold-local outcome-linked procedure, not a locked
  patient-level predictor.
- Call the nine-window object a descriptive stability shortlist.
- Call scale-free geometry the strongest freshly reproduced breakpoint-only
  internal analysis, while labeling it exploratory and requiring external
  validation.
