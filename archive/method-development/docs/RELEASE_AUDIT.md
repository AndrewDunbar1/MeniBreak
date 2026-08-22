# Repository release audit

Audit date: 2026-08-19

- R package build: PASS.
- `R CMD check --no-manual --no-vignettes`: PASS, status OK.
- Package installation and test suite: PASS.
- Synthetic MBPR scoring smoke test: PASS.
- Analysis and package R syntax checks: PASS.
- Frozen model checksum matches the registry and model manifest: PASS.
- Automated scan for raw SEG/IDAT/spreadsheet data, sample-level predictions,
  hard-coded home directories, and common identifier fields: PASS.
- Derivation-cohort row-level tables, clinical workbooks, raw segments, linkage
  keys, OOF predictions, and fold assignments: not included.

The repository includes only code, frozen model/reference assets, synthetic
fixtures, documentation, and aggregate analysis summaries. Automated screening
does not replace a human privacy review before future releases.
