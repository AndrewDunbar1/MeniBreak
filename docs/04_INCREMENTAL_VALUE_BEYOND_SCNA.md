# 4. Incremental value beyond SCNA

The primary contrast is paired: SCNA-only versus SCNA-plus-breakpoint predictions
for the same tumors, folds, repetitions, and bootstrap draws. This isolates the
information added by breakpoint architecture.

For the recovered but unverified Word-report analysis, reported broad-cohort AUCs
were 0.777 vs 0.780 at 3 years, 0.766 vs 0.797 at 5 years, and 0.767 vs 0.797 at
7 years. Restricted-cohort values were 0.796 vs 0.804, 0.766 vs 0.847, and
0.740 vs 0.812. These are document-reported internal results, not performance
guarantees for the locked software artifact or reproduced results. The report
attributes them to fold-local window selection, but its execution cannot be
reconstructed and they must not be attributed to the nine-window stability shortlist;
see `08_PROVENANCE_AND_OPEN_GAPS.md`.

External reports must identify which contrast is being shown: dynamic top-eight
selection procedure; nine-window descriptive shortlist; 18-feature breakpoint
architecture refit recipe; or frozen
`mbpr-derivation-v1` versus its frozen SCNA-only submodel.

Mixing values from these lineages in a single curve without labels is prohibited.
