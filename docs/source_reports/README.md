# source_reports/ — historical documents, kept for the record

These are the **original write-ups** from the analysis that produced the results
in this repository. They are preserved verbatim so the reasoning and the raw
methodological decisions remain inspectable.

| File | What it is |
|---|---|
| `HOTLIST_REPORT.md` | The full internal report on the breakpoint panel, including the analyses not carried into the main docs. |
| `DECISIONS.md` | Every methodological choice made during the build, with rationale. |
| `LIMITATIONS.md` | The original, exhaustive limitations inventory for the wider project (broader in scope than this repository). |
| `breakpoint_panel.tex` | The LaTeX manuscript write-up. |

## How to treat these

They are **historical**, not the maintained specification:

- Where a source report and the shipped code disagree, **the code in `analyses/`
  and `src/` is authoritative**.
- They describe a wider project (pan-cancer pretraining, TCGA generalisation,
  data acquisition) than the panel shipped here, so some sections reference
  components that are not in this repository.
- They contain known internal inconsistencies in their *method prose* — most
  notably the panel size, quoted as both top-k = 8 and top-k = 5. That tension is
  resolved in [`../THE_PANEL.md`](../THE_PANEL.md): the stability table is
  reported at k = 5 and the performance ladder is fitted at k = 8.

For the maintained, current documentation see the parent `docs/` directory:
[`../THE_PANEL.md`](../THE_PANEL.md), [`../LIMITATIONS.md`](../LIMITATIONS.md),
[`../INPUT_SPEC.md`](../INPUT_SPEC.md), [`../QUICKSTART.md`](../QUICKSTART.md),
[`../PROVENANCE.md`](../PROVENANCE.md).
