# Minimal data dictionary

## Canonical SEG input

| Field | Type | Definition |
|---|---|---|
| `sample_id` | string | Deidentified tumor identifier |
| `chromosome` | string/integer | Autosomal chromosome label, normalized by adapter |
| `start` | integer | Segment start in declared coordinate convention |
| `end` | integer | Segment end in declared coordinate convention |
| `event` | categorical | Loss, gain, homozygous/deep loss, high-copy gain/amplification, LOH, or allelic imbalance |

Genome build is a required command-line argument and not a table column guessed by software.

## External outcomes

| Field | Type | Definition |
|---|---|---|
| `sample_id` | string | Must match scored SEG sample ID |
| `patient_id` | string | Deidentified patient cluster; recommended and required for repeated tumors |
| `recurrence_time` | numeric | Years from the prespecified index date to recurrence or last follow-up |
| `recurrence_event` | 0/1 | 1 recurrence; 0 censored without recurrence |

## Derived internal analysis table

Research scripts additionally expect `who_grade`, the prespecified `integrated_grade`, patient ID, recurrence fields, arm-level loss/gain indicators, breakpoint-bin features (`bp_[loss_|gain_]ARM_binNN`), architecture counts, and—for the clinical-integrated analysis—age, mitotic count, GTR, prior RT, and adjuvant RT. This derived table is not needed for locked external scoring.
