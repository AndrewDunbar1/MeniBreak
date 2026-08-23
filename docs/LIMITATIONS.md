# Limitations

These travel with every number in this repository. None of them are footnotes.

## 1. No external validation

The panel was tested on an external methylation-derived cohort (GSE189521) and
**degraded** at every horizon. The reason was measured, not guessed:

- median **97** breakpoints/sample institutional vs **11** external;
- panel window occupancy **27.8%** vs **3.7%**;
- two windows fire in **zero** external tumours.

Arm-level SCNA *does* transfer (AUC 0.755 vs 0.763 published), so the failure is
specific to the fine-resolution windows: a marker defined at ~5 Mb cannot be
evaluated on a platform that detects one-ninth the change points. **A
same-resolution external cohort is required** before any cross-platform claim.
`apply/menibreak_panel.py` prints a resolution-QC line so you can detect this
mismatch in your own data.

## 2. The restricted cohort is small

The restricted (primary, gross-total resection, no radiotherapy, non-angiomatous)
cohort has **22 / 37 / 43 events** at 3 / 5 / 7 years. Its largest reported
effects (e.g. +0.081 AUC from adding the panel to SCNA at 5 years) are also its
least precise. Restricted-cohort numbers use transported broad-cohort out-of-fold
scores (no refitting in the small subset), which is honest but not a substitute
for a properly powered cohort.

## 3. Internal validation only

All institutional AUCs are patient-grouped out-of-fold on a single institution's
cohort. Window selection happens inside training folds, which removes the most
obvious optimism, but this is internal cross-validation, not prospective or
multi-site validation.

## 4. Research use only

Nothing here is a clinical test, a grading criterion, or a treatment rule. No
prospective validation, no regulatory assessment, no claim of clinical utility.

## 5. Known tensions in the historical description

The original report contained internal inconsistencies in its *method
description* (top-k stated as both 8 and 5; "eight windows" vs nine displayed;
a 745-candidate count not exactly recoverable from the retained matrix). The
shipped code (`analyses/`, `src/`) is the authoritative description of what was
actually run; where prose and code disagree, trust the code.
