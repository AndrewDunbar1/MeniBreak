"""Unit tests for the canonical data model, arm tiling and breakpoint derivation.

These run on small synthetic fixtures only. No institutional data is read, so
the suite is safe to run anywhere and contains no PHI.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.harmonization.breakpoints import derive_breakpoints, tile_arms
from src.harmonization.references import Arm, ArmReference
from src.harmonization.schema import (
    BoundaryClass,
    EventClass,
    SchemaError,
    SEGMENT_SPEC,
    map_institutional_event,
    normalise_chromosome,
)


@pytest.fixture
def arm_ref() -> ArmReference:
    """A toy genome: one chromosome, 1000 bp p-arm and 1000 bp q-arm."""
    return ArmReference(
        [
            Arm("1", "p", 0, 1000, 1000),
            Arm("1", "q", 1000, 2000, 1000),
        ],
        "GRCh37",
    )


def _segment(start: int, end: int, event: EventClass, *, track="cn",
             sample="S1", chrom="1", mean=-0.5, probes=100.0) -> dict:
    return {
        "dataset_id": "test", "sample_id": sample, "tumor_id": "T1",
        "patient_id": "P1", "aliquot_id": None, "chromosome": chrom,
        "start": start, "end": end, "genome_build": "GRCh37",
        "segment_mean": mean, "total_copy_number": None,
        "minor_copy_number": None, "major_copy_number": None,
        "probe_count": probes, "event": event.value,
        "platform": "test_array", "caller": "test_caller",
        "source_file": "test.seg", "source_coordinate_build": "GRCh37",
        "qc_status": "pass",
        "segment_id": f"{sample}:{chrom}:{start}-{end}:{track}",
        "track": track,
    }


# --------------------------------------------------------------------------
# Chromosome / event normalisation
# --------------------------------------------------------------------------

@pytest.mark.parametrize(
    "raw,expected",
    [("chr1", "1"), ("Chr1", "1"), (1, "1"), ("22", "22"),
     ("chrX", "X"), ("23", "X"), ("24", "Y")],
)
def test_normalise_chromosome(raw, expected):
    assert normalise_chromosome(raw) == expected


@pytest.mark.parametrize("bad", ["chrM", "MT", "chr99", "scaffold_7"])
def test_normalise_chromosome_rejects_junk(bad):
    """An unparseable contig must raise, never fall through to a default."""
    with pytest.raises(ValueError):
        normalise_chromosome(bad)


def test_unmapped_event_raises():
    """A new caller's vocabulary must not be silently coerced to neutral."""
    with pytest.raises(ValueError, match="unmapped"):
        map_institutional_event("Some New Event Type")


def test_schema_rejects_unknown_columns():
    """Strict validation is what stops an outcome column leaking into inputs."""
    df = pd.DataFrame([_segment(0, 100, EventClass.LOSS)])
    df["recurrence_event"] = 1          # an outcome sneaking in
    with pytest.raises(SchemaError, match="unexpected columns"):
        SEGMENT_SPEC.validate(df.drop(columns=["track"]), strict=True)


# --------------------------------------------------------------------------
# Arm geometry
# --------------------------------------------------------------------------

def test_arm_relative_position_orientation(arm_ref):
    """0.0 at the telomere and 1.0 at the centromere, on BOTH arms."""
    p, q = arm_ref.get("1p"), arm_ref.get("1q")
    # p-arm runs telomere(0) -> centromere(1000)
    assert p.relative_position(0) == pytest.approx(0.0)
    assert p.relative_position(1000) == pytest.approx(1.0)
    # q-arm is reversed: centromere(1000) -> telomere(2000). The measure is
    # distance from the TELOMERE in both cases, so the centromere is 1.0.
    assert q.relative_position(2000) == pytest.approx(0.0)
    assert q.relative_position(1000) == pytest.approx(1.0)
    assert q.relative_position(1500) == pytest.approx(0.5)


# --------------------------------------------------------------------------
# Tiling
# --------------------------------------------------------------------------

def test_tiling_is_gapless_and_complete(arm_ref):
    """A single called segment must yield a complete arm tiling with fillers."""
    seg = pd.DataFrame([_segment(200, 400, EventClass.LOSS)])
    tiles, overlaps = tile_arms(seg, arm_ref, autosomes_only=True)

    assert overlaps.empty
    p_tiles = tiles[tiles.arm == "1p"].sort_values("start")
    # neutral [0,200), loss [200,400), neutral [400,1000)
    assert list(zip(p_tiles.start, p_tiles.end)) == [(0, 200), (200, 400), (400, 1000)]
    assert list(p_tiles.is_filler) == [True, False, True]
    # gapless: each end meets the next start exactly
    assert (p_tiles.end.values[:-1] == p_tiles.start.values[1:]).all()


def test_contiguous_segments_produce_no_filler(arm_ref):
    """Shared-boundary segments are contiguous; no 1-bp filler may appear."""
    seg = pd.DataFrame([
        _segment(0, 500, EventClass.LOSS),
        _segment(500, 1000, EventClass.GAIN, mean=0.5),
    ])
    tiles, _ = tile_arms(seg, arm_ref)
    p = tiles[tiles.arm == "1p"]
    assert not p.is_filler.any(), "half-open touch was misread as a gap"
    assert len(p) == 2


def test_centromere_spanning_segment_is_split(arm_ref):
    """A segment crossing the centromere yields one piece per arm, and the
    centromere itself is NOT emitted as a breakpoint."""
    seg = pd.DataFrame([_segment(800, 1200, EventClass.LOSS)])
    tiles, _ = tile_arms(seg, arm_ref)

    called = tiles[~tiles.is_filler].sort_values("start")
    assert list(zip(called.arm, called.start, called.end)) == [
        ("1p", 800, 1000), ("1q", 1000, 1200)
    ]
    assert called.was_clipped.all()

    bps = derive_breakpoints(tiles, arm_ref, boundary_class=BoundaryClass.SNP_ARRAY_CHANGEPOINT)
    # Only the two real edges of the event (800 and 1200), not the centromere.
    assert sorted(bps.breakpoint_start.tolist()) == [800, 1200]


def test_true_overlap_is_truncated_and_reported(arm_ref):
    seg = pd.DataFrame([
        _segment(100, 500, EventClass.LOSS),
        _segment(300, 700, EventClass.GAIN, mean=0.4),
    ])
    tiles, overlaps = tile_arms(seg, arm_ref)
    assert len(overlaps) == 1
    assert overlaps.iloc[0].overlap_bp == 200
    called = tiles[~tiles.is_filler].sort_values("start")
    assert list(zip(called.start, called.end)) == [(100, 300), (300, 700)]


def test_tracks_do_not_interfere(arm_ref):
    """CN and allelic calls legitimately share coordinates and must be tiled
    independently rather than treated as mutual overlaps."""
    seg = pd.DataFrame([
        _segment(200, 400, EventClass.LOSS, track="cn"),
        _segment(200, 400, EventClass.LOH, track="allelic", mean=0.0),
    ])
    tiles, overlaps = tile_arms(seg, arm_ref)
    assert overlaps.empty
    assert set(tiles.track) == {"cn", "allelic"}


# --------------------------------------------------------------------------
# Breakpoints
# --------------------------------------------------------------------------

def test_breakpoint_uncertainty_is_never_zero(arm_ref):
    """Abutting array segments still have inter-probe uncertainty."""
    seg = pd.DataFrame([_segment(200, 400, EventClass.LOSS, probes=10.0)])
    tiles, _ = tile_arms(seg, arm_ref)
    bps = derive_breakpoints(tiles, arm_ref, boundary_class=BoundaryClass.SNP_ARRAY_CHANGEPOINT)
    assert (bps.coordinate_uncertainty_bp > 0).all()


def test_arm_edges_are_not_breakpoints(arm_ref):
    """An event covering a whole arm has no internal change point."""
    seg = pd.DataFrame([_segment(0, 1000, EventClass.LOSS)])
    tiles, _ = tile_arms(seg, arm_ref)
    bps = derive_breakpoints(tiles, arm_ref, boundary_class=BoundaryClass.SNP_ARRAY_CHANGEPOINT)
    assert bps[bps.arm == "1p"].empty


def test_boundary_credibility_ordering():
    """Acceptance test #8: orthogonally-supported boundaries outrank array ones."""
    assert (BoundaryClass.WGS_SV_SUPPORTED.credibility
            > BoundaryClass.SNP_ARRAY_CHANGEPOINT.credibility
            > BoundaryClass.METHYLATION_ARRAY_CHANGEPOINT.credibility
            > BoundaryClass.UNKNOWN.credibility)
    assert (BoundaryClass.OGM_SUPPORTED.credibility
            == BoundaryClass.WGS_SV_SUPPORTED.credibility)


def test_transition_direction_and_amplitude(arm_ref):
    seg = pd.DataFrame([_segment(200, 400, EventClass.GAIN, mean=0.8)])
    tiles, _ = tile_arms(seg, arm_ref)
    bps = derive_breakpoints(tiles, arm_ref,
                             boundary_class=BoundaryClass.SNP_ARRAY_CHANGEPOINT).sort_values("breakpoint_start")
    # neutral(0.0) -> gain(0.8) is up; gain(0.8) -> neutral(0.0) is down
    assert list(bps.transition_direction) == ["up", "down"]
    assert bps.amplitude_change.tolist() == pytest.approx([0.8, -0.8])


# --------------------------------------------------------------------------
# Public-cohort loaders
# --------------------------------------------------------------------------

# NOTE: tests covering components not shipped in this repository (TCGA/PCAWG
# harmonization, the dataset registry, pretraining-corpus separation) were removed
# when the panel was split out. See docs/source_reports/ for the wider project.
