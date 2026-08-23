"""Loader for the institutional (IRB-governed) master SEG file.

This module is READ-ONLY with respect to the handoff archive. It never writes
back, and it refuses to open any path listed under
`institutional.forbidden_subpaths` in the config (the restricted linkage files
that may carry MRN / surgical-pathology identifiers).

Empirical facts established by direct inspection on 2026-08-02, which drive the
design below:

* 172,091 segment rows, 1,991 unique `cg` (deidentified tumour) ids across
  2,074 unique assay `sample` ids -- so cg is NOT 1:1 with assay, and 29 cg have
  more than one assay. cg is the tumour-level key; `sample` is the assay key.
* The file carries CONTINUOUS `Probe_Median` amplitudes and `Probes` counts in
  addition to the discrete `Event` call. The locked 62-feature model consumes
  only the discrete calls, so amplitude and probe density are genuinely
  unexploited information available to the foundation model.
* Events span two physically distinct tracks that OVERLAP each other by design:
  a copy-number track (CN Loss/Gain, Homozygous Copy Loss, High Copy Gain) and
  an allelic track (LOH, Allelic Imbalance). Flattening them into one sequence
  produced 17,645 overlapping adjacent pairs; splitting by track reduces this to
  354 (CN) and 3,008 (allelic). The tracks are therefore modelled separately.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from .schema import (
    COPY_NUMBER_EVENTS,
    SEGMENT_SPEC,
    EventClass,
    map_institutional_event,
    normalise_chromosome,
)

log = logging.getLogger(__name__)

DATASET_ID = "institutional_meningioma"
GENOME_BUILD = "GRCh37"

#: Columns expected in meningioma_masterseg_ALL_2026.txt.
RAW_COLUMNS = [
    "cg", "sample", "Chr", "Start", "End",
    "Probes", "Probe_Median", "Event", "Cytoband", "source_file",
]


class Track:
    """Physically distinct segment tracks that may overlap one another."""

    COPY_NUMBER = "cn"
    ALLELIC = "allelic"


def assert_path_allowed(path: Path, forbidden_subpaths: list[str]) -> None:
    """Refuse to open restricted-linkage files carrying direct identifiers."""
    p = str(Path(path).as_posix()).lower()
    for bad in forbidden_subpaths:
        if Path(bad).as_posix().lower() in p:
            raise PermissionError(
                f"refusing to read {path}: matches forbidden subpath {bad!r}. "
                f"These files may contain MRN / surgical-pathology identifiers."
            )


def load_master_seg(
    path: str | Path,
    *,
    forbidden_subpaths: list[str] | None = None,
) -> pd.DataFrame:
    """Load the institutional master SEG into the canonical segment schema.

    Returns a frame satisfying SEGMENT_SPEC plus two extra provenance columns
    (`track`, `cytoband`) that are declared in the returned frame's attrs rather
    than smuggled through the strict schema check.
    """
    path = Path(path)
    assert_path_allowed(path, forbidden_subpaths or [])

    raw = pd.read_csv(path, sep="\t", dtype={"cg": "string", "sample": "string"})
    missing = [c for c in RAW_COLUMNS if c not in raw.columns]
    if missing:
        raise ValueError(f"{path}: master SEG missing columns {missing}")

    n_raw = len(raw)
    log.info("read %d raw segment rows from %s", n_raw, path.name)

    # Chromosome normalisation. Raises on an unrecognised contig rather than
    # dropping it, so a parsing regression cannot quietly shrink the cohort.
    chrom = raw["Chr"].map(normalise_chromosome).astype("string")

    events = raw["Event"].map(map_institutional_event)
    track = events.map(
        lambda e: Track.COPY_NUMBER if e in COPY_NUMBER_EVENTS else Track.ALLELIC
    ).astype("string")

    out = pd.DataFrame(
        {
            "dataset_id": DATASET_ID,
            # `sample` is the assay; `cg` is the deidentified tumour.
            "sample_id": raw["sample"].astype("string"),
            "tumor_id": raw["cg"].astype("string"),
            # No separate patient key exists in the SEG file itself; it is
            # supplied by outcomes.csv at link time. Left NA rather than being
            # faked from cg, because 75 patients have multiple tumours.
            "patient_id": pd.Series(pd.NA, index=raw.index, dtype="string"),
            "aliquot_id": pd.Series(pd.NA, index=raw.index, dtype="string"),
            "chromosome": chrom,
            "start": raw["Start"].astype("int64"),
            "end": raw["End"].astype("int64"),
            "genome_build": GENOME_BUILD,
            "segment_mean": raw["Probe_Median"].astype("float64"),
            # Allele-specific copy number is not produced by this caller.
            "total_copy_number": pd.Series(np.nan, index=raw.index, dtype="float64"),
            "minor_copy_number": pd.Series(np.nan, index=raw.index, dtype="float64"),
            "major_copy_number": pd.Series(np.nan, index=raw.index, dtype="float64"),
            "probe_count": raw["Probes"].astype("float64"),
            "event": pd.Series([e.value for e in events], index=raw.index, dtype="string"),
            "platform": "institutional_acgh_cna",
            "caller": "institutional_event_caller",
            "source_file": raw["source_file"].astype("string"),
            "source_coordinate_build": GENOME_BUILD,
            "qc_status": "pass",
        }
    )

    # Stable, deterministic segment key. Includes the track because CN and
    # allelic calls legitimately share identical coordinates (5,765 exact
    # coordinate duplicates exist in the source file for this reason).
    out["segment_id"] = (
        out["sample_id"].astype(str) + ":" + out["chromosome"].astype(str) + ":"
        + out["start"].astype(str) + "-" + out["end"].astype(str) + ":" + track.astype(str)
    ).astype("string")

    SEGMENT_SPEC.validate(out, strict=True)

    # Attach non-schema provenance as separate attributes-backed columns.
    out = out.assign(track=track, cytoband=raw["Cytoband"].astype("string"))
    out.attrs["extra_columns"] = ["track", "cytoband"]
    out.attrs["n_raw_rows"] = n_raw
    out.attrs["source_path"] = str(path)

    if len(out) != n_raw:
        raise AssertionError(
            f"row count changed during load: {n_raw} -> {len(out)}. "
            f"Loading must never drop segments silently."
        )
    return out


def load_outcomes(path: str | Path) -> pd.DataFrame:
    """Load the analysis-ready outcome manifest (sample/patient/time/event).

    `sample_id` in outcomes.csv is the deidentified tumour key (`cg` in the SEG
    file), NOT the assay id. The rename here is the single place that mapping is
    asserted, so it cannot drift between call sites.
    """
    df = pd.read_csv(path)
    required = {"sample_id", "patient_id", "recurrence_time", "recurrence_event"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path}: outcomes missing {sorted(missing)}")
    out = pd.DataFrame(
        {
            "dataset_id": DATASET_ID,
            "tumor_id": df["sample_id"].astype("string"),
            "patient_id": df["patient_id"].astype("string"),
            "recurrence_time": df["recurrence_time"].astype("float64"),
            "recurrence_event": df["recurrence_event"].astype("float64"),
            "provenance": str(Path(path).name),
        }
    )
    return out


def fixed_horizon_labels(
    outcomes: pd.DataFrame, horizon_years: float
) -> pd.DataFrame:
    """Apply the locked model's STRICT fixed-horizon recurrence endpoint.

    Reproduced exactly as specified in the handoff methods so that foundation
    model results are computed on the same eligible tumours as the baseline:

      * recurrence on or before the horizon      -> event (1)
      * event-free follow-up through the horizon -> non-event (0)
      * censored before the horizon, no recurrence -> EXCLUDED (dropped)

    The shrinking denominator at longer horizons is intentional and must not be
    "fixed" by imputing follow-up.
    """
    t = outcomes["recurrence_time"]
    e = outcomes["recurrence_event"]

    is_event = (e == 1) & (t <= horizon_years)
    is_nonevent = t >= horizon_years
    eligible = is_event | is_nonevent

    out = outcomes.loc[eligible].copy()
    out["label"] = is_event.loc[eligible].astype("int64")
    out["horizon_years"] = horizon_years
    return out
