# 1. Concept and study frame

Arm loss records **how much** of the genome is altered. Breakpoint position and
distribution record **how the tumor reached that state**. MeniBreak tests whether
that route information adds prognostic signal beyond grade and SCNA burden.

The derivation analyses use strict 3-, 5-, and 7-year binary recurrence
endpoints. A tumor is eligible only when its status at the horizon is known.
All resampling is patient-grouped so tumors from the same patient cannot appear
on both sides of a train/test split.

Two populations appear in the scientific narrative:

- **Broad cohort:** all tumors meeting the horizon-specific outcome rule.
- **Restricted cohort:** primary, gross-total resection, no prior radiotherapy,
  no adjuvant radiotherapy, and non-angiomatous. The recorded fields—not text
  labels or manual guesses—must implement these restrictions.

The first external collaboration should not repeat feature discovery. It should
apply the nine predefined windows and ask whether their association and panel
performance reproduce.
