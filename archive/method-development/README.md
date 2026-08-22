# Method-development archive

This directory preserves exploratory and internal analysis lineage. It is not
the entry point for collaborators. Start at the repository's `START_HERE.md`.

- `analysis/`: nested-validation and method-development scripts.
- `results/`: aggregate performance and feature-stability summaries only.
- `docs/`: methods, result interpretation, external-data QC, and governance.
- `config/`: example path and schema mappings.

For audit work, start with [`docs/SCRIPT_INDEX.md`](docs/SCRIPT_INDEX.md). The scripts expect
private derivation inputs that are intentionally not distributed. They are
included to make the statistical methods auditable, not to expose institutional
records or encourage refitting the frozen model on an external test cohort.
