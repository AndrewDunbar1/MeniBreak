# A Parsimonious Breakpoint Hot List for Meningioma Recurrence

**Date:** 2026-08-04 · institutional cohort, GRCh37
**Figures:** `reports/figures/fig1–fig4.png`

---

## Summary

Eight breakpoint windows — binary "is there a breakpoint here" indicators,
specified by genomic coordinates — add prognostic information over arm-level
SCNA **within the institutional cohort**. The gain is **+0.025 to +0.036 AUC in
the broad cohort** and rises to **+0.079 in the restricted surgical cohort** at
5 and 7 years.

Every window is rediscovered inside training folds, so the numbers include the
cost of discovery. The panel independently recovers the project's known
**chr1p 112–125 Mb hotspot** at 88–100% selection frequency, and nominates
**chr6q** and **chr19q 41–44.5 Mb** as additional loci.

Two results qualify this substantially, and both are reported in full below:

* **At 7 years the eight windows beat arm-level SCNA on their own** (0.747 vs
  0.738 broad; 0.745 vs 0.692 restricted) and exceed integrated grade in the
  broad cohort. Breakpoint architecture is the only marker examined here that
  *improves* as the horizon lengthens, while WHO grade, integrated grade and
  arm-level SCNA all decay.
* **The panel does not externally validate on GSE189521.** It makes performance
  worse at every horizon. The mechanism is measured and specific: panel windows
  occupy 27.8% of institutional tumours but only 3.7% of the
  methylation-derived external cohort, which detects one-ninth as many
  breakpoints per sample. **No cross-platform claim can be made.**

---

## 1. Performance (Figure 1)

Comparators computed **separately for each cohort** — the published WHO/IG
figures are broad-cohort values and differ materially in the restricted subset.

### Broad cohort

| Horizon | n / events | WHO | Integrated grade | SCNA | **SCNA + 8 windows** | Δ |
|---|---|---:|---:|---:|---:|---:|
| 3 y | 904 / 204 | 0.690 | 0.766 | 0.768 | 0.779 | +0.011 |
| 5 y | 674 / 264 | 0.656 | 0.735 | 0.757 | **0.782** | +0.025 |
| 7 y | 514 / 296 | 0.635 | 0.733 | 0.738 | **0.774** | +0.036 |

### Restricted cohort — primary · GTR · no prior/adjuvant RT · not RT-induced · non-angiomatous

| Horizon | n / events | WHO | Integrated grade | SCNA | **SCNA + 8 windows** | Δ |
|---|---|---:|---:|---:|---:|---:|
| 3 y | 387 / 22 | 0.686 | **0.827** | 0.780 | 0.804 | +0.024 |
| 5 y | 248 / 37 | 0.647 | 0.782 | 0.756 | **0.835** | **+0.079** |
| 7 y | 155 / 43 | 0.622 | 0.752 | 0.692 | **0.770** | **+0.079** |

**Validation:** the SCNA baselines reproduce the locked model's restricted
figures closely (0.780 vs 0.779 at 3 y; 0.756 vs 0.760 at 5 y), and WHO /
integrated grade reproduce the published broad values to within **0.0005**.

---

## 2. The headline finding (Figure 2)

**Breakpoint architecture is roughly three times more informative in a clean
surgical cohort**: +0.079 versus +0.025 at 5 years, +0.079 versus +0.036 at
7 years.

This is biologically coherent. In the broad cohort, recurrence is driven
substantially by treatment and presentation — subtotal resection, prior
radiotherapy, radiation-induced tumours, angiomatous histology. Remove those and
what remains is closer to the tumour's intrinsic biology, which is exactly where
a structural marker should matter most.

At 5 years the 8-window model reaches **0.835**, exceeding the locked model's
restricted result of 0.8167 while using a far simpler genomic representation.

**Where it does not win:** at 3 years in the restricted cohort, integrated grade
(0.827) beats SCNA + windows (0.804). The 3-year horizon is null for breakpoints
in every analysis in this project, consistent with the locked model's
non-significant result (p = 0.051).

---

## 3. The hot list (Figures 3 and 4)

`selection_frequency` = fraction of training folds that independently chose the
window. This is what separates a finding from a coincidence.

### Known hotspot recovered — positive control passes

| Window | Coordinates (GRCh37) | 3 y | 5 y | 7 y |
|---|---|---:|---:|---:|
| chr1p 0.90–0.95 | **112.50–118.75 Mb** | **1.00** | — | 0.44 |
| chr1p 0.95–1.00 | **118.75–125.00 Mb** | **0.96** | **0.88** | **0.92** |

The project's known **chr1p 112–123 Mb hotspot**, recovered at 88–100% across
all three horizons **without the frozen windows ever being an input**. This is
the strongest positive control the project has produced.

### Additional loci

| Window | Coordinates (GRCh37) | 3 y | 5 y | 7 y |
|---|---|---:|---:|---:|
| chr6q 0.90–0.95 | 66.51–72.01 Mb | **1.00** | 0.64 | 0.36 |
| chr19q 0.50–0.55 | 41.18–42.81 Mb | — | **0.80** | **0.68** |
| chr19q 0.45–0.50 | 42.81–44.45 Mb | — | 0.72 | 0.36 |
| chr2p 0.75–0.80 | 69.97–74.64 Mb | — | 0.76 | — |
| chr6q 0.95–1.00 | 61.00–66.51 Mb | 0.40 | 0.24 | **0.68** |

**chr6q 61–72 Mb** and **chr19q 41–44.5 Mb** are the most consistent non-1p
signals. 6q loss is established meningioma biology; the **19q13 window is the
most interesting novel candidate**, selected in 68–80% of folds at the longer
horizons.

Anything below ~0.4 is fold-specific noise and is shown greyed in Figure 3. The
locked model's own audit found its frozen windows exactly reproducible in only
38.3% of folds — the same phenomenon measured differently.

---

## 3b. What the 8 windows carry *by themselves* (Figure 5)

A separate question from "what do they add": how much do eight binary
indicators carry with no arm-level information at all?

| Horizon | scope | WHO | Integrated grade | **Arm SCNA (47 feat.)** | **8 windows ALONE** |
|---|---|---:|---:|---:|---:|
| 3 y | broad | 0.690 | 0.766 | **0.768** | 0.700 |
| 5 y | broad | 0.656 | 0.735 | **0.757** | 0.730 |
| **7 y** | broad | 0.635 | 0.733 | 0.738 | **0.747** ★ |
| 3 y | restricted | 0.686 | **0.827** | 0.780 | 0.716 |
| 5 y | restricted | 0.647 | **0.782** | 0.756 | 0.733 |
| **7 y** | restricted | 0.622 | 0.752 | 0.692 | **0.745** ★ |

★ = beats arm-level SCNA.

Three things follow, and the third is the interesting one:

1. **At 3 and 5 years the windows do not beat SCNA.** They beat WHO grade
   comfortably at every horizon, but 47 arm-level features carry more.
2. **At 7 years they overtake SCNA in both cohorts** — 0.747 vs 0.738 broad,
   and 0.745 vs 0.692 restricted (+0.053). In the broad cohort they also exceed
   **integrated grade** (0.747 vs 0.733).
3. **The direction of travel is opposite for the two marker types.** Across
   3 → 7 years, every conventional marker decays while the breakpoint panel
   strengthens:

   | Marker | 3 y → 7 y (broad) |
   |---|---|
   | WHO grade | 0.690 → 0.635 (**−0.055**) |
   | Integrated grade | 0.766 → 0.733 (**−0.033**) |
   | Arm-level SCNA | 0.768 → 0.738 (**−0.030**) |
   | **8 breakpoint windows** | 0.700 → **0.747** (**+0.047**) |

**Grade predicts early recurrence; breakpoint architecture predicts late
recurrence.** That is a genuinely different prognostic axis, and it is the most
biologically interesting result in this report — eight binary indicators are the
only marker here that improves as the horizon lengthens.

It also reframes the 3-year null. Breakpoints are not "weak at 3 years"; they
are measuring something that expresses itself later.

---

## 3c. EXTERNAL VALIDATION FAILS (Figure 6)

The panel was applied to **GSE189521** (110 tumours, EPIC methylation-derived
CNA), quarantined and untouched for the entire project. Windows fixed from
institutional data, model fitted on institutional data, **nothing refitted**.

| Horizon | n / events | SCNA | SCNA + 8 windows | 8 windows alone | locked SEG-only ref |
|---|---|---:|---:|---:|---:|
| 3 y | 50 / 13 | 0.755 | 0.705 | 0.600 | 0.763 |
| 5 y | 30 / 17 | 0.643 | 0.620 | 0.543 | 0.704 |
| 7 y | 22 / 18 | 0.681 | 0.653 | 0.604 | 0.722 |

**The windows do not transfer. They make performance worse at every horizon**,
and alone they are close to chance.

### Mechanism — this is a resolution failure, not (necessarily) a biology failure

| | institutional (aCGH) | external (EPIC methylation) |
|---|---:|---:|
| median breakpoints per sample | **97** | **11** |
| mean panel-window occupancy | **27.8%** | **3.7%** |

Two of the eight windows fire in **zero** external tumours. The panel is
essentially silent on methylation-derived CNA, so it can only contribute noise.

Note what *does* transfer: **arm-level SCNA holds up** (0.755 externally at
3 years, against the locked model's 0.763). Coarse arm-level calls are
resolution-robust; fine positional breakpoint windows are not. This is the
Phase 3 negative control — platform identifiable at AUC 0.985 from burden
features alone — arriving exactly where it was predicted to matter.

### What this does and does not license

* It does **not** show the windows are biologically spurious. A marker defined
  at ~5 Mb resolution cannot be evaluated on a platform detecting one-ninth as
  many change points.
* It **does** show the panel is **not portable to methylation-derived CNA**, and
  no clinical or cross-platform claim can rest on it as it stands.
* External sample sizes are very small (50 / 30 / 22, with 13 / 17 / 18 events).
  Even a working panel would be hard to confirm here.

A fair external test needs a **same-resolution** cohort — SNP-array or WGS
meningioma with outcomes. GSE147673 (383 tumours, CoreExome genotyping) is
already downloaded and is the obvious candidate, though its outcome fields still
need verification.

---

## 4. Calibration — what this does *not* say

* **The 8 windows do not beat arm-level SCNA on their own.** Alone they reach
  0.726 (5 y) and 0.740 (7 y), below SCNA's 0.757 / 0.738. Their value is
  **additive**. They do beat WHO grade outright, and at 7 years they approach
  integrated grade.
* **The restricted-cohort result rests on 37 and 43 events.** The +0.079
  increment is the largest effect here and the least precisely estimated.
  Confidence intervals will be wide; this needs external replication before it
  can be relied on.
* **Windows are ~5–6 Mb arm-relative deciles**, not breakpoint positions. They
  are regions of interest, not mapped junctions.
* **Selection frequency measures internal stability**, across folds of one
  cohort. It is not external reproducibility.
* Single institution, single platform, retrospective ascertainment.

---

## 5. Next step

The panel is **eight coordinate intervals**, so it transfers to any segmented
copy-number dataset without refitting. **GSE189521** (110 tumours, EPIC
methylation-derived CNA) is downloaded, quarantined, and has never been touched
in this project — making it a genuine external test of exactly this claim.

That is the single highest-value experiment remaining.
