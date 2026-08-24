# Changelog

## external-deployment — the full ladder from a generic SEG (2026-08-23)

Makes the repository usable end-to-end by an outside group on their own data.

### Added
- **`apply/check_resolution.py`** — preflight. Any SEG will run; the question is
  whether the assay localized breakpoints finely enough for a 2.8 Mb window.
  Reports median event size against that window and returns `SUITABLE` /
  `MARGINAL` / `TOO COARSE`. Measured reference points: Agilent aCGH 0.45×
  (panel trained here), conumee on EPIC 1.80× (panel failed). Calibrated on two
  platforms — a decision aid, not a validated instrument.
- **`--clinical`** on `apply/menibreak_panel.py` — unlocks the grade
  comparators, the three `Clinical +` models and the **Restricted** scope.
  Grade is a comparator only and is barred from entering any genomic model.
  The restricted scope requires all six criteria or is refused, and a tumour
  with a missing criterion is excluded rather than assumed favourable.
- **`SCNA + enriched breakpoints`** in `apply/` — the strongest SEG-only model
  in the report (0.819 at 5 y), previously reachable only from `analyses/`.
- **`analyses/07_verify_external_path.py`** — routes the institutional cohort
  through the *external* entry point and diffs against `full_ladder.csv`. All
  9 models × 2 scopes × 3 horizons reproduce to within **0.0031** with exact
  cohort sizes and event counts.

### Fixed
- `patient_id` was required in outcomes despite the spec calling it optional.
- **Follow-up time in months returned AUC 0.085 on 4 events instead of 110,
  silently.** Implausible units are now rejected with the likely cause named.
- Tab/semicolon-separated files reported every column as missing; empty segment
  files raised `KeyError: 'end'`.
- `analyses/` fell back to a local absolute path when `MENIBREAK_DATA_ROOT` was
  unset; it now exits with instructions.
- Provenance manifests recorded an absolute path under a home directory.
- The synthetic example drew uniform 3–40 Mb blocks — coarser than real array
  data, and it failed its own resolution preflight. Event sizes now match the
  institutional cohort, and it carries `segment_mean` and an allelic track so
  the enriched model is actually exercised.

### Notes
- Test suite 53 → 88.
- The example's event *burden* is deliberately below a real cohort's; see the
  comment in `data/example/generate_example.py` for why.

## clean-python-release — fixed-panel validation, tests, figure captions (2026-08-23)

Closes the three gaps identified in the post-release audit.

### Added
- **`apply/validate_fixed_panel.py`** — validate the *published* panel on an
  external cohort with **fixed** windows and no discovery. Default scoring is an
  unweighted count of panel windows containing a breakpoint, so **zero
  parameters are fitted on the new cohort**; `--score refit` keeps the windows
  fixed and refits coefficients under patient-grouped CV. Reports per-window
  occupancy so a silent window is visible rather than hidden.
- **`reference/menibreak_panel_v1.csv`** — the frozen panel definition, derived
  reproducibly by `reference/build_panel_definition.py` and re-checked in CI.
- **`tests/` — 88 pytest tests**, including the patient-leakage guarantees
  (patients never split across folds; every tumour tested once; splits
  deterministic; the leakage detector fails a deliberately bad split), arm-tiling
  correctness, and the guards that keep breakpoint information out of the
  arm-level SCNA feature group and grade out of any model as a predictor.
- `figures/README.md` with a caption for each of the nine figures.
- Six further verified result tables (elevated-grade rescue, tail capture, panel
  size sweep, restricted transport, feature-set search).
- `requirements-dev.txt`; CI now runs the tests, both entry points, and asserts
  the panel definition is reproducible from the discovery output.

### Changed
- `RESULTS_MANIFEST` now indexes **725 result records across 11 tables** (was 437
  across 6) and documents both entry points and the panel definition.

### Documented
- **A discrepancy in the published panel.** The reproducible rule (windows
  selected in ≥40% of folds at k = 5) admits **11** windows; the report published
  **9**, additionally admitting `bp_14q_b18` (0.40) and `bp_2q_b14` (0.44) while
  keeping `bp_3p_b10` and `bp_6q_b19` at that same 0.40. The published nine are
  therefore not recoverable from any single threshold. Both sets are shipped and
  labelled (`in_published_table`); `--panel-set` selects between them.

## clean-python-release — audit pass (2026-08-23)

Documentation and metadata audit; no analysis code or result values changed.

### Changed
- Generalised the external-cohort framing: `docs/HEIDELBERG_QUICKSTART.md` is now
  `docs/QUICKSTART.md`, and site-specific wording was removed throughout.
- README caveat #1 reframed from "no external validation / methylation failure"
  to **"internal validation only"**, which is the accurate statement; the full
  measured detail of the methylation attempt is retained in `docs/LIMITATIONS.md`
  and the manifest, reframed as a resolution mismatch rather than a refutation.
- Rewrote `CITATION.cff`, `NOTICE`, `CONTRIBUTING.md` and `.gitattributes`, which
  still described the removed `mbpr` R package and told contributors to run
  `Rscript tests/smoke.R` and `./tools/audit_no_phi.sh` — neither file exists.
  CONTRIBUTING now documents the actual Python smoke test that CI runs.

### Fixed
- README stated the chr1p window is selected in "88–92%" of folds; the correct
  range is **88–96%** (0.96 / 0.88 / 0.92 at 3 / 5 / 7 years).
- `docs/THE_PANEL.md` now states that the stability table is reported at
  **k = 5** while the performance ladder is fitted at **k = 8**, and that
  `hotlist_windows.csv` contains k ∈ {3, 5, 8, 15}. Showing the two side by side
  without this note is what made the panel look like a fixed nine-window list.

### Added
- `docs/source_reports/README.md` marking those documents as historical, with an
  explicit rule that the shipped code is authoritative where prose disagrees.

## clean-python-release (2026-08-19)

Complete rebuild of the repository around the **actual recovered code** that
produced the MeniBreak report, so every reported number is backed by a runnable
script and an external cohort can apply the panel to its own SEG data.

### Added
- `src/` — the pipeline package recovered from `BreakpointsFoundationModel`
  (SEG harmonization → arm tiling → breakpoints → windows/SCNA/enriched features →
  patient-grouped evaluation), reassembled at its original import paths and
  verified to import and run.
- `apply/menibreak_panel.py` — tested entry point to run the panel on an
  external cohort's SEG + outcomes; recovers a planted chr1p window on the
  bundled synthetic data and prints a resolution-comparability QC line.
- `analyses/01–06` — the six original scripts that generated the report tables,
  with the institutional data root made configurable via `MENIBREAK_DATA_ROOT`.
- `results/reference/` — the verified result CSVs (match the report to the
  printed digit).
- `RESULTS_MANIFEST.yaml` / `.json` — machine-readable index (437 result records)
  generated from the CSVs by `build_manifest.py`.
- `docs/` — INPUT_SPEC, QUICKSTART, THE_PANEL, PROVENANCE, LIMITATIONS.
- `data/example/` — synthetic cohort + generator; `reference/` arm and hotspot
  references; `figures/` publication figures; a CI smoke test that runs the
  panel on synthetic data.

### Changed
- Provenance corrected: the report code originates in the local
  `BreakpointsFoundationModel` repository (drive export `T1:/deliverables`), not a
  ChatGPT thread. The prior "online workspace / never downloaded" claim was
  checked against every source file's metadata and rejected (see `docs/PROVENANCE.md`).

### Removed
- The previous R-package layout (`mbpr`, `renv`, R workflows, prior audits). It
  remains available in git history. This release is single-language (Python) for
  clarity.
