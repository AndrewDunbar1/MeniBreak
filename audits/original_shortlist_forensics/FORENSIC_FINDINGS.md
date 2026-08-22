# Forensic reconstruction of the historical “eight-window” AUCs

## Question

Identify the analysis that produced the historical breakpoint-panel AUCs:

| Cohort | 3 years | 5 years | 7 years |
|---|---:|---:|---:|
| Broad | 0.702 | 0.735 | 0.752 |
| Restricted transport | 0.715 | 0.753 | 0.757 |

The restricted cohort is primary, gross-total resection, no prior or adjuvant radiotherapy, not radiation-induced, and non-angiomatous.

## Source recovered

The values come from `breakpoint_panel_word_tables_fixed.docx`, preserved in the August 10 Downloads backup. No underlying prediction file or generating script was retained with that document.

The document describes 20 equal arm-relative bins, minimum support of 20 tumors, foldwise chi-square ranking, a patient-grouped 5-fold design repeated five times, and elastic-net logistic regression with mixing parameter 0.25.

## Internal contradictions in the source document

1. The methods say **top k=8**, but the window table says **k=5**.
2. The document calls this an **eight-window panel**, but its stability table contains **nine unique windows**.
3. It says filtering leaves **745 candidates**. The described untyped-window feature class leaves 499 candidates in the retained data. Combining untyped, loss, and gain windows leaves 773 by tumor support or 750 by patient support—not 745.
4. The source does not specify the elastic-net penalty strength. A later implementation inserted scikit-learn `C=0.1`; contemporaneous R workflows instead selected `lambda` internally with `cv.glmnet`.

## Reconstruction results

All tests used the retained analysis-ready dataset (SHA-256 `e583f0f98f5b51b65e2eaf074f549f139e7b57ec5b2000cda729f21f28b760df`) and reproduced the exact horizon cohorts and comparator AUCs.

### Described procedure

The current literal scikit-learn interpretation—binary untyped windows, support ≥20, foldwise chi-square top 8, standardized elastic-net logistic regression (`l1_ratio=0.25`, `C=0.1`), five repeated patient-grouped folds—produced:

| Cohort | 3 years | 5 years | 7 years |
|---|---:|---:|---:|
| Broad | 0.724 | 0.707 | 0.741 |
| Restricted transport | 0.686 | 0.675 | 0.733 |

### Native-R interpretation

Using all untyped/loss/gain candidates with patient-level support, R-generated patient folds, foldwise chi-square top 8, and nested `cv.glmnet` with `lambda.1se`, then averaging predictions across repeats, produced:

| Cohort | 3 years | 5 years | 7 years |
|---|---:|---:|---:|
| Broad | 0.724 | 0.716 | 0.745 |
| Restricted transport | 0.703 | 0.711 | 0.757 |

This is the closest faithful reconstruction among the explicitly documented interpretations (six-value RMSE 0.0218), but it does not reproduce the 3- or 5-year restricted results.

### Broad ambiguity grid

A grid covering the retained ambiguities—untyped versus untyped/loss/gain candidates, support ≥20 versus >20, k=5 versus k=8, chi-square versus count-based F test, binary/count/log-count model inputs, four penalty strengths, stored versus regenerated patient folds, and three AUC aggregation rules—did not recover all six values. The numerically closest grid row (RMSE 0.0176) used count-valued candidates and an F test, contradicting the reported binary/chi-square method, and yielded 0.721/0.719/0.739 broad and 0.695/0.730/0.769 restricted.

### Window-stability fingerprint

The retained source claims the 5-year selector chose 19q bins 11 and 10 in 80% and 72% of folds and 2p bin 16 in 76%. Under the reported chi-square procedure, their full-data ranks are approximately 49, 102, and 59 among untyped candidates. They therefore cannot plausibly be top-five windows in most training folds.

Reconstructing candidates directly from the raw SEG using both boundaries, starts only, ends only, or removal of arm-edge boundaries did not restore the claimed rank pattern. The exact named five-window sets also fail to generate the reported AUCs; for example, the named 5-year set gives approximately 0.688 broad and 0.683 restricted under an optimistic full-data fit.

## Conclusion

The six historical AUCs are traceable to the Word report but not to a reproducible retained analysis. At least two incompatible analyses appear to have been combined in that report: one supplying the AUC table and another supplying the window-stability table/method description. The discrepancy cannot be fixed by a fold seed, k=5 versus k=8, penalty convention, AUC aggregation, or ordinary boundary-definition change.

The six values should therefore remain labeled **historical/unverified** and should not be presented as reproduced results. The defensible replacement is a single, versioned workflow with saved out-of-fold predictions, fold assignments, selected windows, and a manifest. The native-R reconstruction above is the closest faithful interpretation, but it is not an exact reproduction and should not inherit the historical values.

## Audit artifacts

- `../original_shortlist_auc_forensics.py`: ambiguity-grid harness.
- `../reconstruct_shortlist_glmnet.R`: native-R `cv.glmnet` reconstruction.
- `../raw_boundary_definition_scan.py`: raw-SEG boundary-definition and source-window rank audit.
- `full_retained_ambiguity_grid/configuration_results.csv`: ranked configuration grid.
- `glmnet_targeted_v2/glmnet_configuration_results.csv`: native-R results.
- `raw_boundary_scan/source_window_ranks.csv`: source-window ranks under raw boundary variants.
