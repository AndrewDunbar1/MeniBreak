# CLAUDE.md

See **[AGENTS.md](AGENTS.md)** — it contains the full instructions for working in
this repository, and applies to Claude Code equally.

Short version: you are almost certainly here to **validate the fixed nine-window
breakpoint hot list on an external cohort**. The protocol is:

1. `pip install -r requirements.txt`
2. Build two CSVs — `segments.csv` and `outcomes.csv` (schema: [docs/INPUT_SPEC.md](docs/INPUT_SPEC.md);
   `recurrence_time` **must be in years**).
3. **Resolution gate, run first:** `python apply/check_resolution.py --segments YOUR_segments.csv`
   → if `TOO COARSE`, **stop and report back**.
4. **The test:** `python apply/validate_fixed_panel.py --segments YOUR_segments.csv --outcomes YOUR_outcomes.csv --out results/validation`
5. Return `results/validation/fixed_panel_performance.csv` and `panel_window_occupancy.csv`.

A machine-parseable version is in [VALIDATION_PROTOCOL.yaml](VALIDATION_PROTOCOL.yaml).

**Do not** tune parameters to improve an AUC, edit anything under `reference/`,
or treat our derivation-cohort figures (0.690/0.772/0.785) as a target — they are
an optimistic anchor, because the windows were chosen on that cohort. A lower
number externally is expected, and a null result is a legitimate outcome worth
reporting plainly. See [docs/WHICH_MODEL.md](docs/WHICH_MODEL.md).
