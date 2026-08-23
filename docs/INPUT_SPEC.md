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
