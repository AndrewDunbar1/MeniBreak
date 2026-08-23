# LIMITATIONS

Maintained alongside `DECISIONS.md`. Negative results and known weaknesses are
recorded here permanently, per go/no-go acceptance test #10.

---

## 1. Platform and caller differences

* The institutional cohort is **event-called**, not continuously segmented. Only
  altered regions are listed; 54.1% of the tiled genome is materialised filler.
  Its amplitude field (`Probe_Median`) is a per-segment summary, not a
  probe-level signal.
* Two provenance strata coexist in one file: `master_cohort_CNA_calls`
  (167,839 rows) and `old_remapped_seg` (4,252 rows). The second was remapped to
  GRCh37 from an older build. Remapped coordinates may carry systematic
  offsets; `source_file` is retained as a covariate and a stratified sensitivity
  analysis is required before any routing claim.
* Planned public cohorts span radically different resolutions:
  Affymetrix 50K+50K (GSE16583, ~115K markers) → Illumina CoreExome-24
  (GSE147673, non-uniform, exon-biased) → EPIC methylation (GSE270371,
  GSE183647) → SNP6 (TCGA) → WGS (PCAWG). **Breakpoint counts are not comparable
  across these without explicit resolution conditioning**, and any model
  trained naively across them will learn platform, not biology.

## 2. Boundary uncertainty

* No array breakpoint is a point. Uncertainty is estimated from segment length ÷
  probe count, which assumes **uniform probe spacing within a segment**. This is
  a poor assumption for exon-biased arrays such as CoreExome-24, where true local
  spacing varies by more than an order of magnitude.
* Where probe counts are missing the fallback is a deliberately coarse
  50 kb. This is a placeholder, not a measurement.
* Methylation-derived change points are **inferred**, driven by probe design,
  and are recorded with `boundary_class = METHYLATION_ARRAY_CHANGEPOINT`. They
  must never be pooled with WGS SV junctions.

## 3. Outcome limitations

* The strict fixed-horizon endpoint **excludes** tumours censored before the
  horizon. The denominator therefore shrinks with horizon
  (904 → 674 → 514 tumours at 3 → 5 → 7 years), and the event *fraction* rises
  (22.6% → 39.2% → 57.6%). Longer-horizon cohorts are enriched for
  longer-followed, and plausibly higher-risk, patients. **AUCs at different
  horizons are not directly comparable.**
* 718 of 1,991 tumours with segment data have **no linked outcome record**. They
  are usable for self-supervised pretraining but contribute nothing to outcome
  evaluation, and their exclusion is not random.
* A further **50 assay samples (493 segment rows) carry a null `cg`** and cannot
  be linked to a tumour at all: 17 named `BL-*` (apparently blood/germline) and
  33 named `S##-####` (surgical accessions never assigned a deidentified id).
  They are excluded from every labelled analysis with a logged count and remain
  available for self-supervised pretraining. Left unhandled they would collapse
  into a single NA pseudo-sample pooling 50 unrelated specimens — a silent and
  fairly serious error, now guarded by
  `test_null_units_are_dropped_with_a_count_not_pooled`. Whether any of the 33
  are genuine tumours that merely failed linkage is **unresolved**, so this is a
  small, non-random exclusion rather than a clean filter.
* Recurrence is ascertained retrospectively from a single institution's imaging
  follow-up. Surveillance intensity is not uniform and is plausibly correlated
  with grade and resection extent.

## 4. Cohort overlap

* **GSE16584 is a SuperSeries of GSE16583 (170) + GSE16581 (68) = 238.** These
  are one cohort in three GEO records, not three cohorts. Verified arithmetically
  and by identical leading GSM ids.
* **GSE16583's 170 GEO samples are ~85 tumours.** The Affymetrix 100K product is
  two physical arrays per DNA sample (Hind240 = GPL2004, Xba240 = GPL2005), each
  registered as a separate GSM. Both numbers are retained in the registry;
  neither is discarded. Per-sample pairing still needs confirmation from sample
  titles at Stage A.
* **GSE252291's "~1,298 profiles" is a harmonised compendium**, not 1,298 new
  samples; only 279 are newly profiled in that series. Constituent studies
  overlap GSE136661 / GSE183647 / GSE270375 and must be deduplicated before any
  split.
* TCGA donors appear inside PCAWG, and the Xena PCAWG hub carries **only ICGC
  identifiers** — no TCGA barcode or project code in any of its histology,
  survival, specimen-type, purity/ploidy or whitelist files. The donor overlap is
  therefore **not resolved by identity**. It is instead closed structurally:
  PCAWG contains **0 meningiomas** of 3,205 specimens and TCGA contains none, so
  neither can share a patient with a meningioma evaluation cohort; both are
  constrained to be pretraining-only and to sit on the same side of every split
  (D-011). Residual exposure is limited to held-out *pretraining diagnostics*,
  not the headline outcome estimate. This becomes a live issue only if either
  corpus is promoted to an evaluation role, or if objective N (SV-support
  prediction) is pursued — which needs the mapping to find matched donors.
* **Unresolved:** the institutional cohort's overlap with public meningioma
  series is unknown and, without shared identifiers, may be unresolvable. Any
  external validation claim must state this explicitly.

## 5. Selection bias

* Single-institution retrospective derivation with referral bias toward complex
  and recurrent disease. The locked model's own documentation flags referral,
  treatment and assay-selection bias.
* 29 tumours were assayed more than once. Repeat assays are not independent
  samples and are flagged so they cannot cross a split.
* 57 samples are segment-count outliers within their own platform/caller
  stratum. These are flagged, not dropped; whether they represent genuine
  chromothripsis or assay failure is unresolved and matters, because extreme
  fragmentation is both a biological signal of interest and a classic artefact.

## 6. Caller sensitivity

* Breakpoint counts depend on the event caller's amplitude threshold and
  smoothing. The institutional caller's parameters are not documented in the
  handoff. A different caller on the same DNA would yield a different breakpoint
  architecture.
* **TCGA allele-specific data mixes two callers.** The GDC "Allele-specific Copy
  Number Segment" type is 11,104 ASCAT2 files and 10,632 ASCAT3 files, and
  **10,502 aliquots carry both** — the same physical assay called twice by
  different pipeline versions. Consequences:
  * The 21,736 files are **~11,234 unique aliquots**, not 21,736 samples.
    Summing them would inflate the pretraining corpus roughly two-fold.
  * Two callings of one aliquot must never land in different splits.
  * `caller` must be a modelled covariate, not an ignored detail.
  * **Silver lining:** these 10,502 pairs are a natural experiment. Same DNA,
    two callers, so caller sensitivity of breakpoint architecture can be
    measured directly rather than assumed — and they are an ideal negative
    control for acceptance test #3 (the route embedding must not encode caller).
* **TCGA ASCAT files carry no `Segment_Mean` and no `Num_Probes`**, only absolute
  `Copy_Number` / `Major` / `Minor`. Any log-ratio must be *derived* using a
  ploidy value that does not ship with the file. Ploidy source must be declared,
  because the choice changes every amplitude and therefore every
  amplitude-thresholded breakpoint.
* The locked model's 30% arm-coverage threshold was selected in a prespecified
  bakeoff, but its selection was a **documented protocol deviation** (V6
  preprocessing was the pre-registered preference and lost by more than the
  0.01 AUC tolerance). Any comparison inherits that choice.
* The frozen hotspot windows were only exactly reproducible in **38.3%** of
  training folds (median Jaccard 0.60, derivable at all in 89.8%). They are
  positive controls, **not ground truth**.

## 7. External validity

* GSE189521 (n=110) is already spent as the locked model's external cohort and
  is quarantined. At the 7-year horizon it retained only **four eligible
  non-events**, making that estimate essentially uninformative. It cannot be
  reused for model selection without destroying its evaluation role.
* Most large public meningioma cohorts are **methylation**, whose inferred CNA
  boundaries are the weakest evidence class. The abundance of public meningioma
  data does not translate into abundant high-resolution breakpoint data.
* OGM support rests on **15 curated breakpoints from 6 exports**. This is
  mechanistic corroboration, not a training signal, and cannot support a
  quantitative claim about what fraction of array change points are physical
  rearrangements.

## 8. Tooling and reproducibility gaps

* ~~R is not installed.~~ **Resolved 2026-08-02**: R 4.6.1 + all 16 BISCUT
  dependencies + the reticulate↔Python bridge are installed and verified.
  Locked-model comparison still uses frozen artefacts **by choice**, not
  necessity (D-008).
* BISCUT has **no test suite, no CI, no dependency pin**, and hard-coded
  parameters; it silently drops segments with `Num_Probes < 4`.
* Two historical upstream SEG files named by `phase2_seg_cg_merge.R` were not
  found at handoff, so the consolidated master SEG cannot be rebuilt from source.
* No complete instrument-level OGM BNX/molecule archive was located.

## 9. Claims explicitly NOT made

This is a **research system**. No claim is made of clinical utility, treatment
guidance, diagnostic validity, or prospective validation. "Locked" means
modelling choices were frozen before external scoring; it does **not** mean
clinically validated.
