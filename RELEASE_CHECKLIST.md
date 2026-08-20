# v0.1.0 release checklist

## Passed locally

- [x] Exact legacy reproduction across 1,273 tumors and 57 legacy features.
- [x] Prespecified legacy/V6 by 30%/50% bakeoff.
- [x] Ten repeats of patient-grouped nested five-fold validation.
- [x] Primary five-year breakpoint AUC/Brier acceptance rule.
- [x] Three-, five-, and seven-year secondary audit.
- [x] Fold-wise hotspot-window stability audit.
- [x] Source-stratified sensitivity analysis.
- [x] Frozen coefficients, feature order, imputation, calibration, percentile
      distributions, seeds, references, and SHA-256 checksums.
- [x] Assay-derived NF2 excluded; SEG-derived 22q retained.
- [x] Event-called Nexus and canonical adapters.
- [x] Explicit continuous caller profiles and coordinate conventions.
- [x] GRCh37 and GRCh38 arm-reference tests.
- [x] CSV and Parquet feature output.
- [x] Native source-tree and clean-library-installed scores identical within
      tolerance; observed maximum delta `0`.
- [x] `R CMD check --no-manual`: `Status: OK`.
- [x] No-network runtime design and local-only processing.
- [x] Apache-2.0 license, citation, model card, validation record, example data,
      and research-use-only limitations.

## Release infrastructure

- [x] Multi-architecture Dockerfile.
- [x] GitHub Actions workflow for `linux/amd64` and `linux/arm64`.
- [x] Native-versus-Docker numerical comparison script.
- [ ] Multi-architecture CI build and numerical comparison must pass after the
      source is placed in its release repository. Docker is not installed on
      the development Mac, so this is the sole unexecuted packaging gate.

Do not publish the container tag until the final CI item passes.
