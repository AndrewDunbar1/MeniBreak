# tests/

```bash
pip install -r requirements-dev.txt
python -m pytest tests/ -q          # 88 tests
```

These test the shipped pipeline directly. They need no data — every fixture is
constructed in the test file — and they run in under a second.

| File | Guards |
|---|---|
| `test_leakage.py` | The central methodological claim: patient-grouped cross-validation. A patient's tumours never split across train and test; every tumour is tested exactly once; splits are deterministic under a seed; the leakage detector **catches a deliberately bad split** (negative control) and passes a good one; the segment schema declares no outcome columns. |
| `test_harmonization.py` | SEG → arms → breakpoints. Chromosome normalisation rejects junk rather than dropping it; unmapped event labels raise; arm tiling is gapless and complete; contiguous segments produce no filler; centromere-spanning segments are split and produce no false change point there; true overlaps are truncated *and reported*; the two tracks do not interfere; breakpoint coordinate uncertainty is never zero. |
| `test_features.py` | Feature construction. Arms are called lost only above threshold; gains and losses stay distinct; filler segments never count as alterations; route entropy behaves correctly; positional zones use telomere→centromere orientation; **breakpoint information never leaks into the arm-level SCNA feature group**; the feature table is deterministic. |
| `test_clinical.py` | Exactly five clinical covariates; sex is not among them; **grade is a comparator, never a predictor** (a guard rejects any feature matrix containing it); mitoses are log1p-transformed; duplicate tumour ids are collapsed. |

## Why these particular tests

Two guarantees carry the scientific weight of this repository, and both are
tested rather than asserted in prose:

1. **No patient leakage.** Reported AUCs are only meaningful if a patient's
   tumours never appear in both training and test folds. `test_leakage.py`
   verifies the splitter *and* verifies that the detector would fail loudly on a
   bad split.
2. **No target leakage into the baseline.** The panel's claim is incremental
   value *over* arm-level SCNA, which requires that breakpoint information stay
   out of the SCNA feature group, and that grade never enter a model as a
   predictor. `test_features.py` and `test_clinical.py` enforce both.

## Scope

Tests covering components not shipped here (TCGA/PCAWG harmonization, the
dataset registry, pretraining-corpus separation) were removed when the panel was
split out of the wider project. Each test file notes this at the end.
