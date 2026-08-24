# Input specification

To run the panel on your cohort you provide **two CSV files**. Nothing else.

## 1. Segments — `--segments`

One row per copy-number segment call. Required columns:

| Column | Type | Meaning |
|---|---|---|
| `sample_id` | string | Tumour identifier (assay/sample key). |
| `chromosome` | string/int | `1`–`22` (or `chr1`…); X/Y are ignored. |
| `start` | integer | Segment start (bp), 1-based inclusive. |
| `end` | integer | Segment end (bp). |
| `event` | string | Copy-number call — see vocabulary below. |

Optional columns, used automatically if present:

| Column | Meaning |
|---|---|
| `patient_id` | Biological patient. **Provide this if one patient can have >1 tumour** — cross-validation folds are grouped by patient to prevent leakage. Defaults to `sample_id`. |
| `tumor_id` | Tumour/resection event, if distinct from the assay `sample_id`. |
| `segment_mean` | log2 ratio / logR (unused by the panel, retained for provenance). |
| `probe_count` | Probes supporting the segment. |
| `platform`, `caller` | Assay platform / caller name (recorded, not modelled). |
| `genome_build` | `GRCh37` or `GRCh38`; otherwise pass `--build`. |

**Only altered segments are needed.** Unaltered stretches are filled in
automatically; you do not list neutral regions.

### `event` vocabulary

Canonical labels: `loss`, `gain`, `homozygous_loss`, `high_copy_gain`, `loh`,
`allelic_imbalance`, `neutral`. Common synonyms are mapped automatically
(`del`/`deletion`/`-1` → loss; `amp`/`amplification`/`1` → gain; `homdel`/`-2` →
homozygous_loss; `dup` → gain; `cnn-loh`/`cnloh` → loh; …). An **unrecognised
label raises an error** rather than being silently treated as neutral — edit
`EVENT_SYNONYMS` in `apply/menibreak_panel.py` or recode your file.

Example (`data/example/synthetic_segments.csv`):

```csv
sample_id,patient_id,chromosome,start,end,event
TUM0001,PAT0001,1,1000000,115000000,loss
TUM0001,PAT0001,7,50000000,90000000,gain
```

## 2. Outcomes — `--outcomes`

One row per tumour. Required columns:

| Column | Type | Meaning |
|---|---|---|
| `sample_id` | string | Must match the segments `sample_id`. |
| `patient_id` | string | For patient-grouped folds. |
| `recurrence_time` | float | Years from treatment to recurrence **or** last follow-up. |
| `recurrence_event` | 0/1 | 1 = recurrence/progression observed; 0 = censored. |

The strict fixed-horizon endpoint is applied per horizon: recurrence on/before
the horizon → event; event-free follow-up through the horizon → non-event;
censored before the horizon → **excluded** (the shrinking denominator at longer
horizons is intentional and must not be imputed).

```csv
sample_id,patient_id,recurrence_time,recurrence_event
TUM0001,PAT0001,2.4,1
TUM0002,PAT0002,8.1,0
```

## 3. Clinical — `--clinical` (optional)

Supplying this unlocks the grade comparators, the three `Clinical +` models and
the **Restricted** cohort scope — the full nine-model ladder from the report.
Omit it and you still get the four SEG-only models.

Only `sample_id` is required. Every other column is optional; whichever models
your columns can support are scored, and the rest are skipped with a logged
reason. Nothing is silently imputed.

| Column | Type | Unlocks |
|---|---|---|
| `sample_id` | string, matches segments | — (required) |
| `age` | years | `Clinical +` models |
| `mitoses` | count per 10 HPF | `Clinical +` models (entered as `log1p`) |
| `gtr` | 1 = gross-total resection, 0 = subtotal | `Clinical +` models |
| `prior_rt` | 1/0 | `Clinical +` models |
| `adj_rt` | 1/0 | `Clinical +` models |
| `who_grade` | 1/2/3 | WHO grade comparator |
| `integrated_grade` | 1/2/3 | Integrated grade comparator |
| `primary` | 1 = primary presentation | Restricted scope |
| `rt_induced` | 1/0 | Restricted scope |
| `angiomatous` | 1/0 | Restricted scope |

All five of `age`, `mitoses`, `gtr`, `prior_rt`, `adj_rt` must be present for the
`Clinical +` models. **All six** restricted criteria (`primary`, `gtr`,
`prior_rt`, `adj_rt`, `rt_induced`, `angiomatous`) must be present for the
Restricted scope — a partial set would define a different subset than the
published one, so it is refused rather than approximated.

A tumour with a *missing* criterion is excluded from the restricted subset, not
assumed favourable. Treating an unknown radiotherapy history as "no
radiotherapy" would quietly enlarge the subset with exactly the patients most
likely to violate it.

`who_grade` and `integrated_grade` are read as **comparator scores only**. They
never enter any genomic model as predictors — that separation is the entire
point of the comparison, and it is enforced in code
(`src.features.clinical.DELIBERATELY_EXCLUDED`), not left to the caller.

```csv
sample_id,age,mitoses,gtr,prior_rt,adj_rt,who_grade,integrated_grade,primary,rt_induced,angiomatous
TUM0001,58,3,1,0,0,2,2,1,0,0
TUM0002,71,0,0,0,1,1,1,1,0,0
```

Restricted-scope figures use the report's transport design: models are fitted on
the broad cohort with patient-grouped out-of-fold prediction, and those scores
are then subset to eligible tumours. Nothing is refitted inside the restricted
subset — in the institutional cohort it holds only 22–43 events.

## Genome build

Coordinates are interpreted with the arm reference for `--build` (`GRCh37`
default, `GRCh38` supported; references in `reference/`). Windows are reported in
that build. The institutional panel was defined on **GRCh37**; if your data is
GRCh38 the arm-relative binning still applies, but compare windows in the same
build.

## Resolution caveat (read this)

The panel is defined at ~5 Mb window resolution and assumes a change-point
detector of roughly SNP-array/WGS sensitivity (~50–100 breakpoints/sample). On a
much coarser platform — methylation-array–derived CNA typically yields ~11
breakpoints/sample — the windows fire far less often, and a weak result reflects
the resolution mismatch rather than the panel.
`apply/menibreak_panel.py` prints your median breakpoints/sample and window
occupancy so you can check comparability before trusting the AUCs.
