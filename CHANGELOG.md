# Changelog

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
- `docs/` — INPUT_SPEC, HEIDELBERG_QUICKSTART, THE_PANEL, PROVENANCE, LIMITATIONS.
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
