# analyses/ — reproduce the report tables on the institutional cohort

These are the scripts that produced the MeniBreak report numbers. They run on the
**institutional dataset**, which is IRB-restricted and not part of this
repository. Point them at it with an environment variable:

```bash
export MENIBREAK_DATA_ROOT=/path/to/Breakpoint_Architecture_Paper_Analysis_Handoff/
python analyses/03_full_ladder.py --k 8 --repeats 5
```

Each script is the recovered original (from `BreakpointsFoundationModel`), with
only the hard-coded data root replaced by `MENIBREAK_DATA_ROOT`. The logic is
unchanged.

| Script | Report content | Output (see `results/reference/`) |
|---|---|---|
| `01_who_ig_comparators.py` | WHO & integrated-grade comparator AUCs (reproduce **exactly**) | `comparators.csv` |
| `02_discover_hotlist.py` | In-fold window discovery + consensus stability windows | `hotlist_windows.csv` |
| `03_full_ladder.py` | Full model ladder, both cohorts: WHO → IG → SCNA → panel → 13/60/102 | `full_ladder.csv` |
| `04_scna_plus_panel.py` | SCNA vs SCNA+panel discrimination; matched-sensitivity operating points | `clinical_metrics_discrimination.csv`, `clinical_metrics_operating.csv` |
| `05_grade_discordance.py` | Reclassification of grade-discordant tumours | `grade_discordance.csv` |
| `06_external_validation.py` | External methylation-cohort test (negative) | `external_validation.csv` |
| `_hotlist_transport.py` | Helper: transport broad OOF scores to the restricted subset | — |

**To run on your own data instead**, do not use these — use
[`../apply/menibreak_panel.py`](../apply/menibreak_panel.py), which takes a
generic SEG + outcomes and needs no institutional layout. See
[`../docs/HEIDELBERG_QUICKSTART.md`](../docs/HEIDELBERG_QUICKSTART.md).

The shipped CSVs in `results/reference/` are the verified outputs of these
scripts and match the report to the printed precision; `RESULTS_MANIFEST.yaml`
indexes every value.
