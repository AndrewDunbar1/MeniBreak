# MeniBreak

MeniBreak is deterministic, local research software that converts segment-level
copy-number data into locked meningioma SCNA-plus-breakpoint features and
recurrence predictions. The R package and command-line executable are named
`mbpr` for compatibility with the frozen model and existing manifests.

The package is research-use-only. It is not a diagnostic device and must not be
used for grading, treatment selection, radiotherapy recommendations, or patient
management. No data leave the machine, and no LLM is involved in normalization,
feature extraction, QC, or scoring.

## Release status

Version 0.1.0 contains the checksum-locked `mbpr-derivation-v1` model. Exact
legacy reproduction passed, and the prespecified 2×2 bakeoff selected portable
legacy preprocessing with a 30% arm threshold. Five-year recurrence is the
primary output; three- and seven-year estimates are secondary.

This is a research release for external validation. Absolute probabilities are
internally calibrated and provisional.

## Installation

### R

```sh
Rscript -e 'install.packages("renv"); renv::restore()'
R CMD INSTALL .
```

### Docker

```sh
docker build -t menibreak:0.1.0 .
docker run --rm menibreak:0.1.0 model-info
```

## Minimal event-called input

```text
sample_id,chromosome,start,end,event
SAMPLE_001,1,1000000,50000000,loss
SAMPLE_001,22,16000000,50000000,gain
```

The genome build (`GRCh37` or `GRCh38`) is mandatory and never guessed.
Supported events include loss, gain, homozygous/deep loss, high-copy
gain/amplification, LOH, and allelic imbalance. The latter two are retained for
QC and excluded from the loss/gain scoring features.

## Local CLI

```sh
./mbpr inspect --segments cohort.seg --format auto --genome-build GRCh37
./mbpr score --segments cohort.seg --format nexus_event \
  --genome-build GRCh37 --output-dir results
./mbpr validate --predictions results/predictions.csv \
  --outcomes outcomes.csv --output-dir validation
./mbpr model-info
./mbpr schema-template --output-dir examples
```

Continuous segment-mean formats require an explicit, versioned caller profile.
Custom thresholds are supported but are marked custom/unvalidated
preprocessing. Scoring never guesses a caller threshold or genome build.

Event-called external data are the fully supported v0.1.0 path. For continuous
segments, the caller profile must also declare whether coordinates are
zero-based half-open or one-based inclusive.

## R interface

```r
library(mbpr)
x <- normalize_segments("cohort.seg", format = "nexus_event",
                        genome_build = "GRCh37")
f <- extract_features(x, preprocess = "legacy", arm_threshold = 30)
inspect_segments("cohort.seg", genome_build = "GRCh37")
```

See the [input specification](INPUT_SPEC.md), [model card](MODEL_CARD.md), and
[external-validation guide](VALIDATION.md).

## Research analyses

The [`research/`](research/) directory documents the report-readable hotspot
model, compact six-rule hotlist, scale-free geometry features, BISCUT-inspired
peak selection, TRISCUT-inspired engineering, harmonized grade comparators, and
clinical-integrated sensitivity analyses. It contains scripts and aggregate
results only—no institutional raw data, linkage keys, or sample-level outputs.

These analyses are reproducibility resources, not alternative production
interfaces. External laboratories should use the locked `mbpr score` workflow
for primary validation.

## Reproducibility

Dependencies are pinned in `renv.lock`. The container is designed to run
without network access after it is built. `run_manifest.json` records the input,
model, reference, and software checksums for every scoring run.

Before publishing changes, run:

```sh
R_PROFILE_USER=/dev/null Rscript tests/smoke.R
./tools/audit_no_phi.sh
```

## License and citation

Code is released under Apache-2.0. Cite the software using [`CITATION.cff`](CITATION.cff).
The software is for research use only and is not intended for clinical decisions.
