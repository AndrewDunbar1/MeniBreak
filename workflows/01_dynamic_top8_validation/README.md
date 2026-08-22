# Dynamic top-eight internal-validation workflow

This workflow implements the literal dynamic top-eight narrative recovered from
the historical report. It is a **fold-local selection procedure**, not a fixed
eight-window score, and its fresh execution does not reproduce the report's
complete six-AUC table:

1. create binary breakpoint-presence features in 20 equal bins per chromosome arm;
2. discard candidates present in fewer than 20 tumors;
3. inside each outer training fold, rank candidates using a chi-square test;
4. retain the top eight and fit an elastic-net logistic model (`l1_ratio=0.25`, `C=0.1`);
5. score only held-out patients; and
6. repeat patient-grouped 5-fold validation five times.

Restricted-cohort results transport the honest broad-cohort out-of-fold scores
to the restricted subset. They do not refit a model with only 22–43 events.

## Run

```sh
./mbpr top8-candidates --segments cohort.seg --format auto \
  --genome-build GRCh37 --output-dir top8_candidates

python3 -m pip install -r workflows/01_dynamic_top8_validation/requirements.txt
python3 workflows/01_dynamic_top8_validation/run_validation.py \
  --features top8_candidates/breakpoint_candidates.csv \
  --outcomes outcomes.csv --output-dir top8_validation
```

This procedure requires outcome data and enough events for repeated grouped
validation. It cannot score a new outcome-free patient. The nine-window
stability shortlist remains the simpler first Heidelberg analysis; it should be
tested window by window rather than assigned the dynamic procedure's AUCs.
