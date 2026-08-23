"""Canonical data model for the breakpoint foundation model.

Three record types are defined here and nothing downstream is permitted to
invent its own column set:

    SEGMENT     - one contiguous copy-number call on one chromosome
    BREAKPOINT  - one change point between two adjacent segments
    CLINICAL    - one patient/sample outcome record

Design rules encoded here (see docs/DECISIONS.md):

1. `genome_build` is REQUIRED and never inferred. A file whose build cannot be
   established from its own header or its registry entry is rejected, not
   guessed.
2. Patient / tumour / specimen / aliquot / assay identifiers are kept as
   SEPARATE columns. They are never collapsed into one `sample_id`.
3. Array boundaries are intervals, not points. Every breakpoint carries
   `breakpoint_start`, `breakpoint_end` and `coordinate_uncertainty_bp`.
4. `boundary_class` records the physical evidence behind a change point. A
   methylation-inferred change point is never represented as equivalent to a
   WGS SV junction.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Sequence

import pandas as pd


# --------------------------------------------------------------------------
# Controlled vocabularies
# --------------------------------------------------------------------------

class BoundaryClass(str, Enum):
    """Physical evidence supporting a change point, strongest to weakest."""

    WGS_SV_SUPPORTED = "WGS_SV_SUPPORTED"
    OGM_SUPPORTED = "OGM_SUPPORTED"
    WGS_CNA_CHANGEPOINT = "WGS_CNA_CHANGEPOINT"
    SNP_ARRAY_CHANGEPOINT = "SNP_ARRAY_CHANGEPOINT"
    METHYLATION_ARRAY_CHANGEPOINT = "METHYLATION_ARRAY_CHANGEPOINT"
    UNKNOWN = "UNKNOWN"

    @property
    def credibility(self) -> int:
        """Ordinal credibility used by go/no-go acceptance test #8.

        Higher is more credible. Exposed as a feature so the model can learn to
        weight orthogonally-supported boundaries above array artefacts, and so
        that the acceptance test can assert the ordering was respected.
        """
        return {
            "WGS_SV_SUPPORTED": 5,
            "OGM_SUPPORTED": 5,
            "WGS_CNA_CHANGEPOINT": 4,
            "SNP_ARRAY_CHANGEPOINT": 3,
            "METHYLATION_ARRAY_CHANGEPOINT": 2,
            "UNKNOWN": 0,
        }[self.value]


class EventClass(str, Enum):
    """Normalised copy-number event labels.

    Deliberately mirrors the locked model's canonical event schema so the
    baseline adapter and the foundation model consume the same vocabulary.
    """

    LOSS = "loss"
    GAIN = "gain"
    HOMOZYGOUS_LOSS = "homozygous_loss"
    HIGH_COPY_GAIN = "high_copy_gain"
    LOH = "loh"
    ALLELIC_IMBALANCE = "allelic_imbalance"
    NEUTRAL = "neutral"


#: Institutional event strings -> canonical labels. Deterministic, no fallback
#: to NEUTRAL: an unrecognised label raises, so silent misclassification of a
#: new caller's vocabulary cannot happen.
INSTITUTIONAL_EVENT_MAP = {
    "CN Loss": EventClass.LOSS,
    "CN Gain": EventClass.GAIN,
    "Homozygous Copy Loss": EventClass.HOMOZYGOUS_LOSS,
    "High Copy Gain": EventClass.HIGH_COPY_GAIN,
    "LOH": EventClass.LOH,
    "Allelic Imbalance": EventClass.ALLELIC_IMBALANCE,
}

#: Events that change DNA copy number. LOH and allelic imbalance are copy-number
#: *neutral* in the general case, so they are tracked separately: they are
#: legitimate structural boundaries but do not contribute to CN burden.
COPY_NUMBER_EVENTS = frozenset(
    {
        EventClass.LOSS,
        EventClass.GAIN,
        EventClass.HOMOZYGOUS_LOSS,
        EventClass.HIGH_COPY_GAIN,
    }
)


class TransitionDirection(str, Enum):
    UP = "up"
    DOWN = "down"
    NEUTRAL = "neutral"


# --------------------------------------------------------------------------
# Column specifications
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class ColumnSpec:
    name: str
    dtype: str
    required: bool
    description: str


@dataclass(frozen=True)
class TableSpec:
    name: str
    columns: Sequence[ColumnSpec]

    @property
    def required_names(self) -> list[str]:
        return [c.name for c in self.columns if c.required]

    @property
    def all_names(self) -> list[str]:
        return [c.name for c in self.columns]

    def validate(self, df: pd.DataFrame, *, strict: bool = True) -> None:
        """Raise SchemaError if `df` does not satisfy this spec.

        `strict` additionally forbids unknown columns, which is how we stop a
        stray outcome column from leaking into an unsupervised input table
        (leakage test #4 in the acceptance criteria).
        """
        missing = [c for c in self.required_names if c not in df.columns]
        if missing:
            raise SchemaError(f"{self.name}: missing required columns {missing}")
        if strict:
            unknown = [c for c in df.columns if c not in self.all_names]
            if unknown:
                raise SchemaError(
                    f"{self.name}: unexpected columns {unknown}. Add them to the "
                    f"spec deliberately rather than passing them through."
                )


class SchemaError(ValueError):
    """Raised when a table violates the canonical data model."""


# --- A. Segment records ----------------------------------------------------

SEGMENT_SPEC = TableSpec(
    name="segments",
    columns=[
        ColumnSpec("dataset_id", "string", True, "Registry dataset_id"),
        ColumnSpec("sample_id", "string", True, "Assay-level identifier"),
        ColumnSpec("patient_id", "string", False, "Biological patient"),
        ColumnSpec("tumor_id", "string", False, "Tumour/resection event"),
        ColumnSpec("aliquot_id", "string", False, "Aliquot, where distinct"),
        ColumnSpec("chromosome", "string", True, "Normalised: '1'..'22','X','Y'"),
        ColumnSpec("start", "int64", True, "1-based inclusive start"),
        ColumnSpec("end", "int64", True, "1-based inclusive end"),
        ColumnSpec("genome_build", "string", True, "GRCh37 | GRCh38 - never guessed"),
        ColumnSpec("segment_mean", "float64", False, "log2 ratio / logR"),
        ColumnSpec("total_copy_number", "float64", False, ""),
        ColumnSpec("minor_copy_number", "float64", False, ""),
        ColumnSpec("major_copy_number", "float64", False, ""),
        ColumnSpec("probe_count", "float64", False, "Probes supporting the segment"),
        ColumnSpec("event", "string", False, "EventClass value"),
        ColumnSpec("platform", "string", True, "Assay platform"),
        ColumnSpec("caller", "string", True, "Segmentation/calling method"),
        ColumnSpec("source_file", "string", True, "Provenance: originating file"),
        ColumnSpec("source_coordinate_build", "string", True, "Build of source_file"),
        ColumnSpec("qc_status", "string", True, "pass | flag | fail"),
        ColumnSpec("segment_id", "string", True, "Stable within-sample segment key"),
    ],
)

# --- B. Breakpoint / change-point records ----------------------------------

BREAKPOINT_SPEC = TableSpec(
    name="breakpoints",
    columns=[
        ColumnSpec("dataset_id", "string", True, ""),
        ColumnSpec("sample_id", "string", True, ""),
        ColumnSpec("patient_id", "string", False, ""),
        ColumnSpec("chromosome", "string", True, ""),
        ColumnSpec("arm", "string", False, "e.g. '1p'"),
        ColumnSpec("left_segment_id", "string", False, "None at arm-start boundary"),
        ColumnSpec("right_segment_id", "string", False, "None at arm-end boundary"),
        ColumnSpec("breakpoint_start", "int64", True, "Interval start, not a point"),
        ColumnSpec("breakpoint_end", "int64", True, "Interval end, not a point"),
        ColumnSpec("coordinate_uncertainty_bp", "float64", True, "end-start, or probe spacing"),
        ColumnSpec("left_state", "string", True, "EventClass to the left"),
        ColumnSpec("right_state", "string", True, "EventClass to the right"),
        ColumnSpec("amplitude_change", "float64", False, "right segment_mean - left"),
        ColumnSpec("transition_direction", "string", True, "TransitionDirection value"),
        ColumnSpec("local_probe_spacing", "float64", False, "bp per probe, locally"),
        ColumnSpec("boundary_class", "string", True, "BoundaryClass value"),
        ColumnSpec("boundary_credibility", "int64", True, "BoundaryClass.credibility"),
        ColumnSpec("arm_relative_position", "float64", False, "0=telomere, 1=centromere"),
        ColumnSpec("distance_to_centromere", "float64", False, ""),
        ColumnSpec("distance_to_telomere", "float64", False, ""),
        ColumnSpec("is_arm_terminal", "bool", True, "Boundary coincides with arm edge"),
        ColumnSpec("caller", "string", True, ""),
        ColumnSpec("platform", "string", True, ""),
        ColumnSpec("genome_build", "string", True, ""),
        ColumnSpec("source_file", "string", True, ""),
        ColumnSpec("qc_status", "string", True, ""),
    ],
)

# --- C. Clinical / outcome records -----------------------------------------

CLINICAL_SPEC = TableSpec(
    name="clinical",
    columns=[
        ColumnSpec("dataset_id", "string", True, ""),
        ColumnSpec("patient_id", "string", True, ""),
        ColumnSpec("sample_id", "string", True, ""),
        ColumnSpec("primary_or_recurrent", "string", False, ""),
        ColumnSpec("surgery_time", "float64", False, "Relative time only - never a date"),
        ColumnSpec("recurrence_event", "float64", False, "1 = event, 0 = censored"),
        ColumnSpec("recurrence_time", "float64", False, "Years"),
        ColumnSpec("last_follow_up_time", "float64", False, "Years"),
        ColumnSpec("os_event", "float64", False, ""),
        ColumnSpec("os_time", "float64", False, ""),
        ColumnSpec("who_grade", "float64", False, ""),
        ColumnSpec("grading_system", "string", False, "e.g. 'WHO 2016', 'WHO 2021'"),
        ColumnSpec("extent_of_resection", "string", False, ""),
        ColumnSpec("prior_radiotherapy", "float64", False, ""),
        ColumnSpec("adjuvant_radiotherapy", "float64", False, ""),
        ColumnSpec("radiation_induced", "float64", False, ""),
        ColumnSpec("sex", "string", False, ""),
        ColumnSpec("age", "float64", False, ""),
        ColumnSpec("molecular_group", "string", False, ""),
        ColumnSpec("integrated_grade", "float64", False, ""),
        ColumnSpec("provenance", "string", True, "Source file for these fields"),
    ],
)

# --- D. Sample linkage -----------------------------------------------------

LINKAGE_SPEC = TableSpec(
    name="sample_linkage",
    columns=[
        ColumnSpec("dataset_id", "string", True, ""),
        ColumnSpec("patient_id", "string", True, "Biological patient"),
        ColumnSpec("tumor_id", "string", False, "Tumour/resection"),
        ColumnSpec("specimen_id", "string", False, ""),
        ColumnSpec("aliquot_id", "string", False, ""),
        ColumnSpec("assay_id", "string", True, "One physical assay run"),
        ColumnSpec("assay_type", "string", True, "snp_array | methylation | rnaseq | wgs | ogm"),
        ColumnSpec("publication", "string", False, ""),
        ColumnSpec("suspected_duplicate_of", "string", False, "Cross-dataset duplicate key"),
        ColumnSpec("linkage_confidence", "string", True, "high | medium | low"),
    ],
)


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

#: Chromosomes retained for arm-level feature construction.
AUTOSOMES = tuple(str(i) for i in range(1, 23))
ALL_CHROMOSOMES = AUTOSOMES + ("X", "Y")


def normalise_chromosome(value: object) -> str:
    """Normalise chromosome naming to bare '1'..'22','X','Y'.

    Accepts 'chr1', 'Chr1', '1', 1, '23'/'24' (which some callers use for X/Y).
    Raises on anything else rather than silently producing a junk contig, since
    a mis-parsed chromosome would scramble every downstream arm assignment.
    """
    s = str(value).strip()
    if s.lower().startswith("chr"):
        s = s[3:]
    s = s.upper()
    if s == "23":
        s = "X"
    elif s == "24":
        s = "Y"
    if s == "MT" or s == "M":
        raise ValueError("mitochondrial contig is not supported for SCNA analysis")
    if s not in ALL_CHROMOSOMES:
        raise ValueError(f"unrecognised chromosome {value!r} -> {s!r}")
    return s


def map_institutional_event(value: object) -> EventClass:
    """Map an institutional Event string to the canonical vocabulary."""
    key = str(value).strip()
    if key not in INSTITUTIONAL_EVENT_MAP:
        raise ValueError(
            f"unmapped institutional event {key!r}; add it to "
            f"INSTITUTIONAL_EVENT_MAP deliberately rather than defaulting"
        )
    return INSTITUTIONAL_EVENT_MAP[key]


ALL_SPECS = {
    "segments": SEGMENT_SPEC,
    "breakpoints": BREAKPOINT_SPEC,
    "clinical": CLINICAL_SPEC,
    "sample_linkage": LINKAGE_SPEC,
}
