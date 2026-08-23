# DECISIONS

Running log. Newest first. Every entry: decision, rationale, alternatives,
evidence, date, code affected.

---

## D-015 — Adversarial platform removal abandoned; invariance measured by a FRESH probe
**Date:** 2026-08-04
**Decision:** Stop using adversarial platform removal as the invariance
mechanism. Measure platform identifiability only with an **independently
trained** classifier, never the model's own probe head. Prefer conditional
encoding (platform as an explicit input) plus **within-platform evaluation**.

**Rationale — the adversary fools its own probe without removing information.**
An adversarial-weight sweep at a fixed 12 epochs (0.0 / 0.05 / 0.2 / 0.5):

| weight | fresh-probe platform acc | route→burden R² | route→hotspot R² |
|---:|---:|---:|---:|
| 0.00 | 0.972 | 0.099 | 0.210 |
| 0.05 | 0.980 | 0.198 | 0.192 |
| 0.20 | 0.994 | 0.163 | 0.296 |
| 0.50 | 0.998 | 0.192 | **0.310** |

At 12 epochs the model's **internal** probe reported CE 1.232 against a 1.059
chance level — apparent success — while a freshly trained classifier still
recovered platform at **99.8%**. The encoder had learned to defeat its own
probe, not to discard the information. Only by 25 epochs was platform genuinely
removed (fresh probe 51.8%), and by then route→burden and route→hotspot had both
gone to zero: platform and biology were destroyed together.

**A wrong intermediate conclusion, recorded:** an earlier reading called this a
"cliff between weight 0.2 and 0.5". That compared 12-epoch sweep points against a
25-epoch run and attributed a duration effect to weight. The driver is
accumulated adversarial pressure (weight x epochs), not weight.

**Why removal cannot work here:** the Phase 3 negative control showed breakpoint
burden largely *is* assay resolution (platform AUC 0.985 from ordinary burden
features). Platform and biology are entangled at the feature level, so
subtracting one subtracts the other.

**Positive finding retained:** pan-cancer pretraining genuinely improves the
representation — best route→hotspot R² 0.310 versus 0.269 for the single-cohort
prototype, recovered with the frozen windows never used as an input.

**Affects:** `src/models/encoder.py`, `src/objectives/pretraining.py`,
`src/evaluation/representation.py`; Phase 5 evaluation design.

---

## D-014 — Platform invariance uses a BOUNDED confusion loss, not gradient reversal
**Date:** 2026-08-03
**Decision:** Replace gradient-reversal-on-cross-entropy with
``KL(p || uniform) = log K - H(p)``, bounded in ``[0, log K]``. Train the
platform probe separately on a **detached** embedding. Add a divergence guard
that aborts training when any objective leaves its plausible range.
**Rationale:** The first pan-cancer run diverged. Cross-entropy is unbounded
above, so the encoder maximised it by inflating logits rather than removing
platform information:

| epoch | adversarial | coordinate |
|---:|---:|---:|
| 5 | 5.7 | 0.0231 (healthy) |
| 6 | 351.2 | 0.0601 (route damaged) |
| 24 | 14,804.4 | 0.0654 (diverged) |

Useful learning stopped at epoch 5; the following 19 epochs degraded the route
embedding's positional accuracy — the signal D-012 identified as most valuable.
**A wrong first fix, recorded:** cross-entropy against a uniform target was also
unbounded (1355 at logit scale 1,000). Only the entropy formulation saturates
correctly — measured 1.9428 at scale 1,000 and 1.9459 at 100,000 against the
log(7) = 1.9459 ceiling, and exactly 0.0 for a uniform posterior.
**Result after the fix:** coordinate loss halved (0.0248 → 0.0129) instead of
degrading; every objective improved monotonically over 25 epochs; the detached
probe sat at CE 1.23 against a prior-entropy chance level of 1.06, i.e. it could
not beat the class prior — platform information was genuinely removed.
**Alternatives:** lowering the adversarial weight — rejected, it delays
divergence rather than preventing it, since the objective is unbounded at any
weight.
**Evidence:** `tests/test_models.py::test_adversarial_confusion_loss_is_bounded`,
`::test_divergence_guard_stops_training`.
**Affects:** `src/objectives/pretraining.py`, `src/models/encoder.py`,
`src/models/train.py`.

---

## D-013 — Resolution/credibility fields are NUISANCE covariates, never free predictors
**Date:** 2026-08-02
**Decision:** `boundary_credibility` and `coordinate_uncertainty_bp` are removed
from the architecture feature set and from the route branch of the encoder. They
remain available as nuisance covariates (for conditioning and adversarial
training) and as a TARGET for the SV-support objective — never as free inputs to
a representation intended to be platform-invariant.
**Rationale:** `boundary_credibility` is a deterministic function of
`boundary_class`, which is a deterministic function of platform. The Phase 3
negative control ranked it the single strongest platform discriminator
(|coef| 2.11): it was handing the platform label to the model as a feature.
**What the fix did and did not achieve:** array-vs-array platform AUC moved only
0.9864 → 0.9854. The leak was real but was not the main driver — platform is
also encoded in ordinary burden features (`total_breakpoint_boundaries`,
`arm_gain_count`, `loss_segment_count`), because assay resolution determines how
many change points are detected at all. **Acceptance test #3 therefore cannot be
passed by feature hygiene alone**, which is what justifies the gradient-reversal
adversary rather than treating it as optional.
**Caveat recorded:** the 0.985 figure conflates platform with disease (each
cohort is a different disease mix; `cnv_22q`/`cnv_1p` among the top
discriminators are meningioma biology). It is an upper bound. The clean test is
same-disease/different-platform — institutional aCGH vs methylation-derived
meningioma CNA — and is queued.
**Affects:** `src/features/architecture.py`, `src/models/tokenizer.py`
(`RESOLUTION_FEATURE_INDICES`), `src/models/encoder.py`.
**Evidence:** `tests/test_models.py::test_route_branch_ignores_resolution_features`.

---

## D-012 — The route embedding must preserve POSITION, not just geometry summaries
**Date:** 2026-08-02
**Decision:** The Phase 4 route embedding will retain per-arm positional
structure. It will not compress an arm to burden/entropy/concentration
statistics, and architecture benchmarking will always include a
positional-resolution ablation.
**Rationale — measured on the institutional cohort:** splitting the breakpoint
increment into a *global architecture* component and a *location-specific
routing* component gave, at 5 years:

| Component | AUC | increment |
|---|---:|---:|
| SCNA only | 0.7543 | — |
| + global architecture (burden, arms touched, entropy, concentration, zonal fractions, fragmentation) | 0.7696 | +0.0153 |
| + hotspot routing (position-specific) | 0.7945 | **+0.0249** |

**Where a breakpoint lands carried more prognostic information than every global
architecture summary combined.** Any representation that reduces an arm to
counts or dispersion statistics forfeits the larger half of the known breakpoint
signal before training starts. At 7 years the split was +0.0138 / +0.0174, same
direction.
**Alternatives:** a compact per-arm summary encoder — rejected on this evidence.
**Caveat carried forward:** the frozen windows used to measure the routing
component were derived on this same cohort, so +0.0249 is an optimistic upper
bound. The foundation model must rediscover windows inside training folds only;
the honest number will be smaller. See `docs/PHASE3_BASELINES.md` §3.
**Affects:** Phase 4 encoder design, `src/features/architecture.py`.

---

## D-011 — TCGA↔PCAWG donor overlap is neutralised structurally, not resolved by identity
**Date:** 2026-08-02
**Decision:** Do **not** pursue ICGC↔TCGA donor identity mapping now. Instead
constrain TCGA and PCAWG to be **pretraining-only** corpora that always sit on
the same side of every split and are never used as held-out evaluation sets.
Overlap then cannot cross a split by construction.

**Rationale — measured, not assumed:**
* PCAWG contains **0 meningiomas** of 3,205 specimens (its CNS entries are
  CNS-Medullo 154, CNS-PiloAstro 96, CNS-GBM 45, CNS-Oligo 21). TCGA never
  included meningioma either.
* Every outcome evaluation in this project is meningioma: the institutional
  cohort and the quarantined GSE189521. **Neither TCGA nor PCAWG shares a single
  patient with any evaluation cohort**, so donor overlap between them cannot
  contaminate the headline AUC by any path.
* The residual exposure is narrow: held-out *pretraining* diagnostics computed
  on PCAWG would be optimistic if a donor's WGS were in test while its SNP6 were
  in train. Those are diagnostics, not the decisive result, and the structural
  constraint removes the exposure at zero cost.
* The Xena PCAWG hub carries only ICGC ids (`SP…`/`DO…`) with no TCGA barcode or
  project code in any of its histology, survival, specimen-type, purity/ploidy
  or whitelist files. Resolving identity would require an external ICGC release
  file, i.e. real effort for a risk that is already structurally closed.

**Condition that forces revisiting:** pretraining objective N (true-SV support
prediction) *wants* these donors matched — a donor with both PCAWG WGS SV calls
and TCGA SNP6 segments is exactly the training pair for learning which array
change points are physical rearrangements. At that point the mapping becomes
required **to find pairs, not to remove them**, and must be obtained from an
official ICGC/PCAWG release file with `submitted_donor_id`.

**Alternatives:** (a) Chase the ICGC mapping now — rejected as effort spent
ahead of need. (b) Ignore the issue — rejected; the constraint is free and the
specification explicitly requires overlap not to cross splits.
**Evidence:** `tests/test_leakage.py::test_pretraining_only_corpora_are_never_evaluation`.
**Affects:** `data/registry/datasets.yaml`, split construction in Phase 2.

---

## D-010 — Acquisition runs take an exclusive lock
**Date:** 2026-08-02
**Decision:** `downloader.run()` holds an `O_EXCL` lockfile at
`data/registry/acquisition.lock` for the whole transfer. A second run refuses to
start while the owning PID is alive; a lock owned by a dead PID is reclaimed.
**Rationale:** An incident during Stage B. A background run believed dead was
still live, and a second run was started against the same targets. Both appended
to the same `.part` files at independent offsets — a silent-corruption path,
because a resumed download whose final size happens to match `Content-Length`
would pass the completeness check. The affected archives were tested with
`gzip -t` and all three passed, so nothing was lost, but that was luck rather
than design.
**Alternatives:** (a) Per-file locks — rejected as more machinery for no extra
safety at this scale. (b) Rely on operator discipline — rejected; that is what
already failed.
**Evidence:** `tests/test_acquisition.py::test_lock_blocks_a_second_live_run`
and `::test_stale_lock_from_dead_pid_is_reclaimed`.
**Affects:** `src/acquisition/downloader.py`.

---

## D-009 — TCGA payloads are fetched in parallel batches via POST /data
**Date:** 2026-08-02
**Decision:** Fetch the 21,736-file ASCAT3 inventory as tar.gz batches of 500
through `POST https://api.gdc.cancer.gov/data`, with 6 concurrent requests.
**Rationale:** Measured, not assumed. A 10-file request returns in 2.4 s but a
500-file request takes ~100 s — the GDC assembles the archive server-side, so
the cost is latency, not bandwidth (the whole ASCAT3 set is only ~190 MB).
Serial batching would take ~73 minutes. Six workers reduce this to a few minutes
without hammering the API.
**Alternatives:** (a) One request per file — 21,736 round trips, rejected.
(b) One request for all 21,736 — the GDC rejects very large id lists and there
would be no resume granularity.
**Affects:** `src/acquisition/fetch_gdc.py`.

---

## D-008 — Locked model is consumed through frozen artefacts, not re-executed
**Date:** 2026-08-02
**Decision:** The locked 62-feature model enters comparisons via its already-computed
frozen outputs (`oof_predictions_locked.csv`, `ml_breakpoint_model_dataset_2026.csv`,
`locked_model_auc_summary.csv`), not by re-running the MBPR R package.
**Rationale:** It is impossible to accidentally refit a model you never execute. This
also sidesteps the missing-R blocker. The handoff already contains honest
out-of-fold predictions for every tumour, which is exactly what a comparator needs.
**Alternatives:** Install R and re-run MBPR — rejected for now: adds a refit risk
surface for zero analytical gain, since the frozen out-of-fold scores are the
authoritative baseline.
**Evidence:** `03_ANALYSIS_READY_DATA/locked_internal/` contains per-tumour OOF
predictions and the release audit summary.
**Affects:** `src/evaluation/` (baseline adapter, Phase 3).

---

## D-007 — Institutional SEG coordinates are HALF-OPEN `[start, end)`
**Date:** 2026-08-02
**Decision:** Treat institutional segment coordinates as half-open. A pair with
`next.start == prev.end` is **contiguous**, not overlapping.
**Rationale:** Established empirically, not assumed. Grouping by
(assay sample, chromosome, track) gives **20,125** adjacent pairs with
`next.start == prev.end` and **0** pairs with `next.start < prev.end`. A shared
boundary coordinate is the file's encoding of contiguity.
**Consequence of getting it wrong:** the first implementation used a closed
convention and manufactured **18,276 phantom one-base truncations**, displacing
every affected breakpoint by 1 bp. Caught by comparing overlap counts against a
direct measurement of the raw file.
**Alternatives:** Closed `[start, end]` — refuted by the zero-true-overlap result.
**Affects:** `src/harmonization/breakpoints.py` (`_resolve_overlaps`, `tile_arms`).
**Note:** The arm reference CSV uses a *closed* convention (chr1 p ends
125,000,000; q starts 125,000,001). The 1 bp inconsistency is immaterial for arm
assignment and is deliberately left uncorrected rather than silently rewriting a
locked-model reference asset.

---

## D-006 — Copy-number and allelic events are modelled as TWO SEPARATE TRACKS
**Date:** 2026-08-02
**Decision:** Split segments into a `cn` track (CN Loss, CN Gain, Homozygous Copy
Loss, High Copy Gain) and an `allelic` track (LOH, Allelic Imbalance). Tile and
derive breakpoints within each track independently.
**Rationale:** The two are physically distinct measurements that legitimately
occupy the same coordinates. Flattening them produced 17,645 overlapping adjacent
pairs; splitting by track reduces true overlaps to zero. The file also contains
5,765 exact coordinate duplicates, which are the same region carrying both a CN
call and an allelic call — not errors.
**Alternatives:** (a) Discard LOH/AI — rejected, they are real structural
boundaries and LOH is 60% of all rows (103,029/172,091). (b) Merge into one
track — rejected, fabricates overlaps and destroys the distinction between a
copy-number change point and an allelic transition.
**Affects:** `src/harmonization/institutional.py`, `breakpoints.py`, `schema.py`.

---

## D-005 — Arms are completely tiled with explicit NEUTRAL filler segments
**Date:** 2026-08-02
**Decision:** Materialise unaltered stretches as explicit `neutral` segments so
each arm is gaplessly tiled before breakpoints are derived.
**Rationale:** The model representation is an *alternating* SEGMENT/BREAKPOINT
sequence. In an event-called file that lists only alterations, two events 40 Mb
apart would otherwise look adjacent, making a spurious "transition" out of a long
unaltered gap. Filler fraction on the institutional cohort is 0.541, i.e. most of
the genome — confirming the gaps are real and large.
**Alternatives:** Emit breakpoints only at called-segment edges — rejected, loses
the distinction between adjacency and separation, which is precisely the spatial
information the project is about.
**Affects:** `src/harmonization/breakpoints.py::tile_arms`.

---

## D-004 — Breakpoints are INTERVALS with non-zero uncertainty, never points
**Date:** 2026-08-02
**Decision:** Every breakpoint carries `breakpoint_start`, `breakpoint_end`,
`coordinate_uncertainty_bp` and `local_probe_spacing`. Uncertainty is
`max(bp_end - bp_start, local probe spacing)` and is asserted `> 0` in QC.
**Rationale:** An array cannot localise a change point more finely than its probe
spacing. Recording a point coordinate would let the model learn spurious
sub-probe precision and would make a sparse CoreExome boundary look as sharp as a
dense SNP6 one.
**Alternatives:** Store midpoints — rejected, discards the resolution covariate
required for cross-platform robustness.
**Affects:** `src/harmonization/breakpoints.py`, `src/qc/checks.py`.

---

## D-003 — Arm-edge boundaries are NOT breakpoints; centromere-spanning segments split
**Date:** 2026-08-02
**Decision:** Segments crossing the centromere are split at the arm boundary, and
boundaries coinciding with an arm edge are not emitted as breakpoints.
**Rationale:** The edge of an arm is the edge of the assayed interval, not an
observed transition. A whole-arm loss has no internal change point. Emitting one
would add a constant to every arm-level event and inflate breakpoint burden in
exactly the tumours with the most arm-level alteration — confounding the burden
control the project depends on.
**Affects:** `src/harmonization/breakpoints.py`; test
`test_centromere_spanning_segment_is_split`.

---

## D-002 — BISCUT is wrapped and pinned, never modified or vendored
**Date:** 2026-08-02
**Decision:** BISCUT lives in `external/biscut` at commit
`885110c056e3a66d41a2bdc98194f5c80e117b0c`, is gitignored, and is invoked as an
external process against a staged working directory.
**Rationale:** (a) The prompt forbids rewriting it. (b) It is **GPL-2.0**
(declared in source headers; no LICENSE file present), so copying it into `src/`
would raise copyleft obligations over our code. (c) All its parameters are
hard-coded module constants, so a staging directory is the only way to configure
it without editing source.
**Alternatives:** Fork and parameterise — rejected on both licence and
"do not modify" grounds.
**Evidence:** `docs/BISCUT_COMPATIBILITY.md`.
**Affects:** `configs/default.yaml` (`biscut:` block), `.gitignore`.

---

## D-001 — Patient-level splitting is enforced by an adversarial test, not convention
**Date:** 2026-08-02
**Decision:** All splits group by `patient_id`. The guard
`assert_no_patient_leakage` is itself tested against a deliberately bad
specimen-level split.
**Rationale:** 75 patients contribute multiple tumours, and 29 tumours were
assayed more than once. A guard that is never tested against a real failure can
silently degrade to a no-op.
**Affects:** `src/evaluation/splits.py`, `tests/test_leakage.py`.
