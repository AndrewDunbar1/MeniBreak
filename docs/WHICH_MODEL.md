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

## The gap this repository has not yet closed

**The internal reference AUC for model B has not been computed.** The derivation
cohort has never been scored with the frozen nine-window count, so at present an
external cohort obtains a number for B with no matching internal number to
compare it against.

`analyses/07_fixed_panel_internal_reference.py` computes exactly that, scored the
same way `apply/validate_fixed_panel.py` scores an external cohort. It requires
the institutional dataset (`MENIBREAK_DATA_ROOT`). Until it has been run, treat a
fixed-panel result as a **standalone measurement**, not as agreement or
disagreement with a published figure.

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
