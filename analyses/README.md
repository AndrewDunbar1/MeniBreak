# analyses/ — reproduce the report tables on the institutional cohort

These are the scripts that produced the MeniBreak report numbers. They run on the
**institutional dataset**, which is IRB-restricted and not part of this
repository. Point them at it with an environment variable:

```bash
export MENIBREAK_DATA_ROOT=/path/to/Breakpoint_Architecture_Paper_Analysis_Handoff/
python analyses/03_full_ladder.py --k 8 --repeats 5
```

Scripts `01`–`06` are the recovered originals (from `BreakpointsFoundationModel`)
that generated the shipped CSVs; the only edit is replacing the hard-coded data
root with `MENIBREAK_DATA_ROOT`, so their logic is unchanged. Note that some of
them were themselves written to complete combinations earlier analyses had left
partial — `03_full_ladder.py` says so in its own docstring — which is why its
output is the authoritative ladder. Script `07` is new in this repository.

| Script | Report content | Output (see `results/reference/`) |
|---|---|---|
| `01_who_ig_comparators.py` | WHO & integrated-grade comparator AUCs (reproduce **exactly**) | `comparators.csv` |
| `02_discover_hotlist.py` | In-fold window discovery + consensus stability windows | `hotlist_windows.csv` |
| `03_full_ladder.py` | Full model ladder, both cohorts: WHO → IG → SCNA → panel → 13/60/102 | `full_ladder.csv` |
| `04_scna_plus_panel.py` | SCNA vs SCNA+panel discrimination; matched-sensitivity operating points | `clinical_metrics_discrimination.csv`, `clinical_metrics_operating.csv` |
| `05_grade_discordance.py` | Reclassification of grade-discordant tumours | `grade_discordance.csv` |
| `06_external_validation.py` | External methylation-cohort test (negative) | `external_validation.csv` |
| `07_verify_external_path.py` | Proves `apply/` reproduces the report: rewrites the institutional cohort into the generic external format, runs it through `apply/menibreak_panel.py`, and diffs all 9 models × 2 scopes × 3 horizons against `full_ladder.csv` (largest \|Δ\| AUC 0.0031) | `_verify_external/external_path_agreement.csv` (gitignored — contains institutional data) |
| `08_fixed_panel_internal_reference.py` | **Internal reference for the FIXED panel** — what the frozen 9/11-window count achieves on the derivation cohort, scored exactly as an external cohort scores it. **Run**; result is an *optimistic anchor* — see [`../docs/WHICH_MODEL.md`](../docs/WHICH_MODEL.md). | `fixed_panel_internal_reference.csv` |
| `_hotlist_transport.py` | Helper: transport broad OOF scores to the restricted subset | — |

**To run on your own data instead**, do not use these — use
[`../apply/menibreak_panel.py`](../apply/menibreak_panel.py), which takes a
generic SEG + outcomes and needs no institutional layout. See
[`../docs/QUICKSTART.md`](../docs/QUICKSTART.md).

The shipped CSVs in `results/reference/` are the verified outputs of these
scripts and match the report to the printed precision; `RESULTS_MANIFEST.yaml`
indexes every value.
