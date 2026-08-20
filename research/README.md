# Research analysis archive

This directory preserves the analysis lineage behind MeniBreak while keeping
the supported scoring package at the repository root.

- `analysis/`: nested-validation and method-development scripts.
- `results/`: aggregate performance and feature-stability summaries only.
- `docs/`: methods, result interpretation, external-data QC, and governance.
- `config/`: example path and schema mappings.

Start with [`docs/SCRIPT_INDEX.md`](docs/SCRIPT_INDEX.md). The scripts expect
private derivation inputs that are intentionally not distributed. They are
included to make the statistical methods auditable, not to expose institutional
records or encourage refitting the frozen model on an external test cohort.
