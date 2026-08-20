# Contributing

Contributions that improve deterministic input adapters, QC, documentation,
tests, or external-validation support are welcome.

1. Work on a feature branch.
2. Do not commit clinical data, SEG/IDAT files, sample-level predictions,
   linkage keys, or direct/indirect identifiers.
3. Add synthetic fixtures for new formats and edge cases.
4. Preserve explicit genome-build, coordinate, and caller-profile requirements.
5. Run `R_PROFILE_USER=/dev/null Rscript tests/smoke.R` and
   `./tools/audit_no_phi.sh` before opening a pull request.
6. Document changes that affect features or scores in the model card and
   version the relevant adapter/model asset.

Changes to frozen preprocessing, feature definitions, coefficients,
calibration, or reference distributions require a new model version; they must
not silently replace `mbpr-derivation-v1`.
