# Input specification

## Canonical event-called schema

| Field | Type | Rule |
|---|---|---|
| `sample_id` | string | Required, non-empty, stable within a run |
| `chromosome` | string/integer | 1–22, X, or Y; `chr` prefix is normalized |
| `start` | integer | Non-negative genomic coordinate |
| `end` | integer | Must be greater than `start` |
| `event` | string | Loss, gain, homozygous/deep loss, high-copy gain/amplification, LOH, allelic imbalance, or neutral |

Coordinates are treated as zero-based half-open internally. Canonical
event-called inputs use that convention. Continuous caller profiles explicitly
declare `zero_based_half_open` or `one_based_inclusive`; one-based starts are
converted deterministically. The institutional Nexus adapter retains its native
derivation coordinates to preserve exact reproduction.

The declared genome build is mandatory. Cytoband, probe counts, clinical
variables, survival data, and assay accessions are not required for scoring.

## Built-in adapters

- `canonical_event`
- `nexus_event`: `cg`, `Chr`, `Start`, `End`, `Event`
- `cbs_log2`: `ID`, `chrom`, `loc.start`, `loc.end`, `seg.mean`
- `gistic_seg`: `sample`, `chromosome`, `start`, `end`, `segment_mean`
- `gatk_modelsegments`: `sample`, `contig`, `start`, `end`,
  `mean_log2_copy_ratio`
- `custom_mapping`: YAML column map

Continuous layouts are score-ready only when a caller profile supplies
`profile_id`, `loss_threshold`, `gain_threshold`, `validation_status`, and
`coordinate_system`.

## Outcomes for external validation

Required columns are `sample_id`, `recurrence_time` in years, and binary
`recurrence_event`. `patient_id` is optional and defaults to `sample_id`; when
present it is used as the clustered bootstrap unit.

At each horizon, recurrence on or before the horizon is an event, event-free
follow-up through the horizon is a non-event, and earlier censoring is excluded.
