# Contributing

Contributions that improve input adapters, QC, documentation, tests, or
external-cohort support are welcome.

1. Work on a feature branch.
2. **Never commit patient data** — clinical tables, SEG/IDAT files,
   sample-level predictions, linkage keys, or direct/indirect identifiers.
   Only synthetic fixtures and aggregate results belong in this repository.
3. Add a synthetic fixture for any new input format or edge case.
4. Preserve explicit genome-build, coordinate, and event-vocabulary
   requirements. An unrecognised event label must raise, never default to
   neutral.
5. Keep window selection, imputation, and standardisation **inside training
   folds**. Any change that could leak test-fold information will be rejected.
6. Before opening a pull request, run the smoke test:

   ```bash
   pip install -r requirements.txt
   python data/example/generate_example.py --n 300 --seed 7 --outdir data/example
   python apply/menibreak_panel.py \
       --segments data/example/synthetic_segments.csv \
       --outcomes data/example/synthetic_outcomes.csv \
       --out results/example
   python build_manifest.py
   ```

   This is the same sequence CI runs (`.github/workflows/smoke-test.yml`).
7. If you change a feature definition or a modelling choice, update
   `docs/` and regenerate `RESULTS_MANIFEST.yaml` with `build_manifest.py`.

Reported results in `results/reference/` are the verified outputs of the scripts
in `analyses/`. Do not edit those CSVs by hand.
