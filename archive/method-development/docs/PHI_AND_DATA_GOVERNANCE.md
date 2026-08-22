# PHI and data-governance boundary

This transfer contains code, frozen model assets, synthetic fixtures, reference coordinates, and aggregate results only. It intentionally excludes raw institutional SEG data, clinical workbooks, outcome rows, sample-level predictions, fold assignments, patient identifiers, accession linkage files, and rendered documents that may contain identifiers.

External teams should:

- use coded research IDs in every analysis file;
- store the re-identification key in an approved clinical system, not this repository;
- keep raw inputs and generated sample-level outputs in ignored `private_input/` and `results/` directories;
- inspect HTML/JSON/CSV outputs before sharing because rare clinical combinations may still be sensitive;
- follow local IRB, DUA, institutional security, and data-retention requirements;
- never commit IDATs, SEG files, outcome manifests, or sample-level predictions to a public repository.

The included `.gitignore` is a guardrail, not a substitute for review. Before every release, run `tests/audit_no_phi.sh` and manually review the resulting manifest.
