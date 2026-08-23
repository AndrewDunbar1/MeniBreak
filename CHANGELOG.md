# Changelog

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
