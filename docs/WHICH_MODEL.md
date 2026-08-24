# Which model am I validating?

Read this before quoting any AUC. The word "panel" refers to **three different
objects** in this project, and they are not interchangeable.

| # | Object | Windows | Weighting | Where it appears |
|---|---|---|---|---|
| **A** | **Dynamic top-8 procedure** | 8, **re-selected inside every training fold** by chi-square | elastic net (l1_ratio 0.25, C 0.1) | The reported AUCs: **0.702 / 0.735 / 0.752** (broad, 3/5/7 y). `analyses/03_full_ladder.py`, `apply/menibreak_panel.py` |
| **B** | **Frozen published panel** | **9, fixed** (`reference/menibreak_panel_v1.csv`, `in_published_table`) | unweighted count (default) or refit coefficients | `apply/validate_fixed_panel.py` — the external-validation default |
| **C** | **Frozen rule panel** | **11, fixed** (the reproducible ≥40% rule) | as B | `apply/validate_fixed_panel.py --panel-set rule` |

## The consequence, stated plainly

**A is what was measured. B is what is portable.** They are different models, so:

> **Do not compare a fixed-panel result (B or C) against 0.702 / 0.735 / 0.752.**
> Those numbers came from a procedure that re-fits window selection *and*
> coefficients inside each fold of the derivation cohort. A fixed, unweighted
> nine-window count is a strictly more constrained model and would not be
> expected to match them.

The honest position is that **A is a demonstration that breakpoint position
carries recurrence signal**, while **B is the deployable marker that demonstration
motivates.** A is internally validated; B is what an external cohort can actually
test.

## The internal reference for model B

Computed by `analyses/08_fixed_panel_internal_reference.py` on the derivation
cohort, scored exactly the way `apply/validate_fixed_panel.py` scores an external
cohort. Full table: `results/reference/fixed_panel_internal_reference.csv`.

**Published nine-window panel, broad cohort:**

| Scoring | 3 y | 5 y | 7 y |
|---|---:|---:|---:|
| Fixed panel, unweighted count (0 fitted params) | 0.690 | **0.772** | **0.785** |
| Fixed panel, refit coefficients | 0.713 | 0.756 | 0.756 |
| Arm-level SCNA (same folds, reference) | 0.777 | 0.766 | 0.767 |

Cohort sizes reproduce the published denominators exactly (904/204, 674/264,
514/296). The frozen count rises with horizon and overtakes arm-level SCNA at
5 and 7 years — the same pattern as model A, from a marker with **no fitted
parameters at all**.

### Read this number as an optimistic anchor, not an expectation

⚠️ **It is biased upward, and the bias cannot be removed retrospectively.** The
count score fits nothing at scoring time, but the *nine windows themselves were
chosen because they recurred across folds of this cohort's outcomes*. Scoring
them on that same cohort is therefore in-sample **with respect to window
selection**.

Two consequences:

1. **Do not read 0.785 as "the frozen panel beats the dynamic 0.752."** Model A's
   figures are honest estimates (selection happens inside each training fold);
   model B's derivation figures are not. They are not comparable in that
   direction either.
2. **An external cohort scoring below these values is the expected outcome, not a
   failure.** Treat 0.690 / 0.772 / 0.785 as an upper reference. The genuinely
   unbiased estimate of what the frozen panel does on data it did not help
   choose is precisely what an external validation produces — which is why it is
   worth running.

The honest comparator inside your own cohort is the **arm-level SCNA AUC printed
alongside your result**, computed on your data and your folds. That comparison is
internal to your cohort and carries none of this bias.

## How to report a validation result

State all four of these, or the number is not interpretable:

1. **Which object** — B (published 9) or C (rule 11).
2. **Which scoring** — `count` (zero fitted parameters) or `refit` (coefficients
   refit under patient-grouped CV).
3. **Panel occupancy** — from `panel_window_occupancy.csv`. Windows that fire in
   zero tumours cannot contribute, and several silent windows usually indicates a
   resolution mismatch rather than absent signal (see
   [`INPUT_SPEC.md`](INPUT_SPEC.md)).
4. **The comparator** — the arm-level SCNA AUC printed alongside, computed in the
   same cohort on the same folds.

## Why the panel was not simply frozen from the start

Freezing the windows on the derivation cohort and then reporting performance on
that same cohort would be circular — the windows were chosen because they
separated those outcomes. Selecting inside each training fold is the honest way
to estimate the *method's* performance. The cost is that the reported number
describes the method rather than a fixed marker, which is precisely why the two
must be labelled separately rather than merged.
