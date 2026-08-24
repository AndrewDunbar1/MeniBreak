# Draft email — external validation request

**Subject:** External validation request — meningioma breakpoint panel (one command, ~30 min)

---

Hi [NAME],

Thanks again for agreeing to look at this. Below is everything you need; the
actual work is a single command once your data is in two simple CSVs.

**The claim we would like tested.** In our institutional meningioma cohort, *where*
copy-number breakpoints fall carries recurrence-risk information beyond WHO
grade, integrated grade, and arm-level copy-number alteration — and unlike those
markers, it gets *stronger* at longer horizons. We have distilled that into a
fixed nine-window "breakpoint hot list". **That hot list is the one thing we are
asking you to test.** Everything else is secondary.

Repo (I have shared access): **https://github.com/AndrewDunbar1/MeniBreak**

## What you need

Two CSVs — no clinical data required to start:

- `segments.csv` — `sample_id, chromosome, start, end, event`
- `outcomes.csv` — `sample_id, recurrence_time` (**in years**), `recurrence_event` (1/0)

Full schema, including optional columns: `docs/INPUT_SPEC.md`.

## The run

```bash
git clone https://github.com/AndrewDunbar1/MeniBreak.git && cd MeniBreak
pip install -r requirements.txt

# 1. Is your assay fine-grained enough for this to mean anything?
python apply/check_resolution.py --segments YOUR_segments.csv

# 2. The test: the fixed hot list, scored with ZERO parameters fitted on your data
python apply/validate_fixed_panel.py \
    --segments YOUR_segments.csv \
    --outcomes YOUR_outcomes.csv \
    --out results/validation
```

**Please run step 1 first and stop if it says `TOO COARSE`.** The windows are
~2.8 Mb, so a platform that localises breakpoints coarsely (methylation-derived
copy number is the usual case) cannot test this — a weak result there tells us
about the assay, not the biology. If it says `SUITABLE` or `MARGINAL`, continue.

Then send back the two files in `results/validation/`:
`fixed_panel_performance.csv` and `panel_window_occupancy.csv`. Neither contains
patient-level data.

## Easiest route: let a coding assistant set it up

If anyone on your team uses **Claude Code** or **Codex**, they can hand the whole
thing over. Both read their instructions from the repo automatically — Claude Code
picks up `CLAUDE.md`, Codex picks up `AGENTS.md`, and both are already in the repo
root. There is nothing to configure.

```bash
git clone https://github.com/AndrewDunbar1/MeniBreak.git
cd MeniBreak
claude          # or: codex
```

Then paste something like:

> Read AGENTS.md and VALIDATION_PROTOCOL.yaml. Our copy-number data is currently
> in <describe your format, e.g. "one .seg file per sample from ASCAT">, and our
> outcomes are in <file>. Convert them to the required segments.csv and
> outcomes.csv, run the resolution check, then run the primary fixed-panel
> validation and show me the results.

The data conversion is the only fiddly step and this is where an assistant helps
most. The instructions also carry the guardrails — most importantly that it must
not tune anything to chase a better AUC, and must stop and tell you if the
resolution check fails — so you get a clean result rather than an optimised one.

`VALIDATION_PROTOCOL.yaml` is the same protocol in machine-readable form if you
would rather script it.

## How to read your number

For reference, on our derivation cohort the fixed hot list scores **0.690 / 0.772
/ 0.785** at 3 / 5 / 7 years, against arm-level SCNA at 0.777 / 0.766 / 0.767.

One important caveat, so nobody over-reads it: **those are an optimistic anchor,
not an expectation.** The nine windows were chosen using our cohort's outcomes,
so scoring them there is partly circular. **A lower number in your cohort is the
expected result, not a failure.** The clean comparison is the arm-level SCNA
figure the tool prints alongside your result — that one is computed on your data
and your folds, and carries none of our bias. Details in `docs/WHICH_MODEL.md`.

## If the hot list doesn't hold up

That is a real and publishable answer, and we would rather know. If it comes out
flat, there are three progressively looser fallbacks, in this order:

1. `--score refit` — same fixed windows, but let the weights be re-fitted on your
   cohort. Tests whether the *window set* is informative even if our unweighted
   scoring is not calibrated to your population.
2. `--panel-set rule` — an 11-window variant of the list.
3. `python apply/menibreak_panel.py …` — ignore our list entirely and let the
   method rediscover windows inside your own cross-validation folds. This answers
   the broader question: does *your* cohort have recurrent prognostic breakpoint
   structure at all, and does it land anywhere near ours (we would look for
   chr1:118.75–125.00 Mb)?

Only worth doing if step 2 above is unconvincing — we would rather have a clean
negative on the primary test than a positive found after four attempts.

## Optional and secondary — the full model ladder

Entirely separate from the pass/fail question above, and only if you have the
clinical fields to hand: adding a third CSV runs **every model in our report on
your cohort at once** — the two grade comparators, arm-level SCNA, the hot list,
the combined SCNA + hot list, and the larger clinical + SCNA + breakpoint models
(9 models × broad/restricted × 3 horizons).

```bash
python apply/menibreak_panel.py \
    --segments YOUR_segments.csv --outcomes YOUR_outcomes.csv \
    --clinical YOUR_clinical.csv --out results/full_ladder
```

`clinical.csv` needs only `sample_id`; every other column is optional and the
tool scores whatever your columns support, logging what it skipped. The useful
ones are `age, mitoses, gtr, prior_rt, adj_rt` (unlock the clinical models),
`who_grade, integrated_grade` (unlock the grade comparators), and
`primary, rt_induced, angiomatous` (unlock the restricted surgical subset). Full
table in `docs/INPUT_SPEC.md` §3.

Grade is read as a **comparator only** — it never enters a genomic model as a
predictor, which is enforced in code rather than left to convention.

This is the interesting exploratory run: it shows whether breakpoints add
anything on top of *your* grading and copy-number data, not just on top of ours.
But please treat the hot-list result above as the headline either way — a clean
answer on one prespecified test is worth more than a sweep across nine.

## Two things to keep in mind

- These are **internally validated, single-institution** results. Yours is the
  first genuine external test — that is exactly why we are asking.
- Research use only; no clinical claims. Everything runs locally and nothing
  leaves your machine.

Happy to jump on a call for the data formatting step if that is easier — it is
usually the only fiddly part. And if anything errors, send me the traceback and
I will turn it around quickly.

Best,
[YOUR NAME]

---

### Notes for the sender (delete before sending)

- Confirm repo access is shared before sending.
- The full-ladder (`--clinical`) option is included but framed as secondary, so
  the primary ask stays small. Drop that section entirely if you want an even
  tighter first email.
- Their run needs no institutional data and takes minutes on a laptop.
- `docs/QUICKSTART.md` is the long-form version of this email if they want it.
- The coding-assistant section is a convenience, not a requirement — the manual
  commands above it stand alone. Cut it if the lab is unlikely to use those tools.
