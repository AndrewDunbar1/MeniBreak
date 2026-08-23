"""Tests for hand-engineered architecture features and the baseline protocol."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.features.architecture import (
    ARM_ALTERED_THRESHOLD,
    arm_level_features,
    breakpoint_features,
    build_feature_table,
    drop_unidentified,
    feature_groups,
)
from src.harmonization.references import Arm, ArmReference
from src.harmonization.schema import EventClass


@pytest.fixture
def arm_ref() -> ArmReference:
    return ArmReference(
        [Arm("1", "p", 0, 1000, 1000), Arm("1", "q", 1000, 2000, 1000)], "GRCh37"
    )


def _tile(start, end, event, *, sample="S1", arm="1p", filler=False) -> dict:
    return {
        "sample_id": sample, "arm": arm, "chromosome": "1",
        "start": start, "end": end, "event": event.value,
        "is_filler": filler, "track": "cn",
        "segment_id": f"{sample}:{arm}:{start}-{end}",
    }


# --------------------------------------------------------------------------
# Arm-level calls
# --------------------------------------------------------------------------

def test_arm_called_lost_only_above_threshold(arm_ref):
    """The 30% arm rule is the locked model's; it must not drift."""
    below = pd.DataFrame([_tile(0, 250, EventClass.LOSS)])       # 25% of 1p
    above = pd.DataFrame([_tile(0, 350, EventClass.LOSS)])       # 35% of 1p

    assert arm_level_features(below, arm_ref)["cnv_1p"].iloc[0] == 0
    assert arm_level_features(above, arm_ref)["cnv_1p"].iloc[0] == -1
    assert ARM_ALTERED_THRESHOLD == 0.30


def test_gain_and_loss_are_distinguished(arm_ref):
    t = pd.DataFrame([_tile(0, 400, EventClass.GAIN, arm="1q")])
    out = arm_level_features(t, arm_ref)
    assert out["cnv_1q"].iloc[0] == 1
    assert out["arm_gain_count"].iloc[0] == 1
    assert out["arm_loss_count"].iloc[0] == 0
    assert out["arm_altered_count"].iloc[0] == 1


def test_filler_segments_never_count_as_alteration(arm_ref):
    """Materialised neutral filler is an absence of signal, not a call."""
    t = pd.DataFrame([_tile(0, 900, EventClass.NEUTRAL, filler=True)])
    out = arm_level_features(t, arm_ref)
    assert out.empty or out["arm_altered_count"].iloc[0] == 0


def test_homozygous_loss_counts_toward_loss(arm_ref):
    t = pd.DataFrame([_tile(0, 400, EventClass.HOMOZYGOUS_LOSS)])
    assert arm_level_features(t, arm_ref)["cnv_1p"].iloc[0] == -1


# --------------------------------------------------------------------------
# Breakpoint features
# --------------------------------------------------------------------------

def _bp(sample, arm, rel, left=EventClass.LOSS, right=EventClass.NEUTRAL) -> dict:
    return {
        "sample_id": sample, "arm": arm, "arm_relative_position": rel,
        "left_state": left.value, "right_state": right.value,
        "coordinate_uncertainty_bp": 1000.0, "boundary_credibility": 3,
    }


def test_route_entropy_is_zero_when_concentrated_on_one_arm():
    bps = pd.DataFrame([_bp("S1", "1p", 0.5), _bp("S1", "1p", 0.6)])
    tiles = pd.DataFrame([_tile(0, 100, EventClass.LOSS)])
    out = breakpoint_features(bps, tiles)
    assert out["breakpoint_route_entropy"].iloc[0] == pytest.approx(0.0)
    assert out["max_arm_breakpoint_fraction"].iloc[0] == pytest.approx(1.0)
    assert out["breakpoint_arms_touched"].iloc[0] == 1


def test_route_entropy_rises_when_spread_across_arms():
    spread = pd.DataFrame([_bp("S1", a, 0.5) for a in ["1p", "1q", "2p", "2q"]])
    tiles = pd.DataFrame([_tile(0, 100, EventClass.LOSS)])
    out = breakpoint_features(spread, tiles)
    assert out["breakpoint_route_entropy"].iloc[0] == pytest.approx(np.log(4))
    assert out["max_arm_breakpoint_fraction"].iloc[0] == pytest.approx(0.25)


def test_positional_zones_use_telomere_to_centromere_orientation():
    """rel 0 = telomere, 1 = centromere; terciles must follow that."""
    bps = pd.DataFrame([
        _bp("S1", "1p", 0.1),      # telomeric
        _bp("S1", "1p", 0.5),      # central
        _bp("S1", "1p", 0.9),      # centromeric
    ])
    tiles = pd.DataFrame([_tile(0, 100, EventClass.LOSS)])
    out = breakpoint_features(bps, tiles)
    assert out["telomeric_breakpoint_fraction"].iloc[0] == pytest.approx(1 / 3)
    assert out["central_breakpoint_fraction"].iloc[0] == pytest.approx(1 / 3)
    assert out["centromeric_breakpoint_fraction"].iloc[0] == pytest.approx(1 / 3)


def test_neutral_only_junctions_are_not_modelled_breakpoints():
    """A junction with no CNA/AI flank carries no architecture information."""
    bps = pd.DataFrame([
        _bp("S1", "1p", 0.5, left=EventClass.NEUTRAL, right=EventClass.NEUTRAL)
    ])
    tiles = pd.DataFrame([_tile(0, 100, EventClass.LOSS)])
    assert breakpoint_features(bps, tiles).empty


# --------------------------------------------------------------------------
# Unlinkable specimens
# --------------------------------------------------------------------------

def test_null_units_are_dropped_with_a_count_not_pooled():
    """50 institutional samples have a null cg. Left in place they would
    collapse into ONE pseudo-sample pooling unrelated specimens."""
    df = pd.DataFrame({"unit": ["A", None, "B", None], "x": [1, 2, 3, 4]})
    out, n = drop_unidentified(df, "unit", what="test")
    assert n == 2
    assert list(out["unit"]) == ["A", "B"]


# --------------------------------------------------------------------------
# Feature grouping
# --------------------------------------------------------------------------

def test_feature_groups_keep_breakpoint_info_out_of_the_scna_model():
    """If breakpoint columns leaked into the SCNA-only baseline, the measured
    breakpoint increment would be erased."""
    cols = ["cnv_1p", "cnv_22q", "arm_loss_count",
            "total_breakpoint_boundaries", "breakpoint_route_entropy"]
    g = feature_groups(cols)
    assert set(g["scna"]) == {"cnv_1p", "cnv_22q", "arm_loss_count"}
    assert "total_breakpoint_boundaries" not in g["scna"]
    assert "breakpoint_route_entropy" not in g["scna"]
    assert set(g["full"]) == set(cols)
    assert not set(g["scna"]) & set(g["breakpoint"])


def test_build_feature_table_is_deterministic(arm_ref):
    tiles = pd.DataFrame([_tile(0, 400, EventClass.LOSS)])
    bps = pd.DataFrame([_bp("S1", "1p", 0.5)])
    a = build_feature_table(tiles, bps, arm_ref)
    b = build_feature_table(tiles, bps, arm_ref)
    pd.testing.assert_frame_equal(a, b)
    assert list(a.columns) == list(b.columns), "column order must be stable"


# --------------------------------------------------------------------------
# Hotspot routing (location-specific, not just burden)
# --------------------------------------------------------------------------

_WINDOWS = pd.DataFrame({
    "arm_key": ["1p", "6q"],
    "hotspot_start_norm": [0.85, 0.85],
    "hotspot_end_norm": [0.95, 0.95],
})


def _arm_calls(**kw) -> pd.DataFrame:
    """Per-sample arm call frame, e.g. _arm_calls(S1={'cnv_1p': -1})."""
    return pd.DataFrame.from_dict(kw, orient="index").fillna(0)


def test_hotspot_counted_only_inside_the_window():
    from src.features.architecture import hotspot_routing_features

    bps = pd.DataFrame([
        _bp("S1", "1p", 0.90),      # inside 0.85-0.95
        _bp("S1", "1p", 0.50),      # outside
    ])
    out = hotspot_routing_features(bps, _arm_calls(S1={"cnv_1p": -1}), _WINDOWS)
    assert out.loc["S1", "dynamic_hotspot_breakpoints_all"] == 1


def test_hotspot_requires_the_arm_to_be_lost():
    """A boundary in the window on a non-lost arm is not a routing event."""
    from src.features.architecture import hotspot_routing_features

    bps = pd.DataFrame([_bp("S1", "1p", 0.90)])
    lost = hotspot_routing_features(bps, _arm_calls(S1={"cnv_1p": -1}), _WINDOWS)
    not_lost = hotspot_routing_features(bps, _arm_calls(S1={"cnv_1p": 0}), _WINDOWS)

    assert lost.loc["S1", "dynamic_hotspot_breakpoints_all"] == 1
    assert not_lost.loc["S1", "dynamic_hotspot_breakpoints_all"] == 0


def test_non1p_variant_excludes_chr1p():
    """The non-1p features test whether routing carries signal beyond the
    single best-known locus."""
    from src.features.architecture import hotspot_routing_features

    bps = pd.DataFrame([_bp("S1", "1p", 0.90), _bp("S1", "6q", 0.90)])
    calls = _arm_calls(S1={"cnv_1p": -1, "cnv_6q": -1})
    out = hotspot_routing_features(bps, calls, _WINDOWS)

    assert out.loc["S1", "dynamic_hotspot_breakpoints_all"] == 2
    assert out.loc["S1", "dynamic_hotspot_breakpoints_non1p"] == 1
    assert out.loc["S1", "dynamic_hotspot_arms_used_all"] == 2
    assert out.loc["S1", "dynamic_hotspot_arms_used_non1p"] == 1


def test_samples_without_hotspot_routing_get_zero_not_missing():
    from src.features.architecture import hotspot_routing_features

    bps = pd.DataFrame([_bp("S1", "1p", 0.10), _bp("S2", "1p", 0.90)])
    calls = _arm_calls(S1={"cnv_1p": -1}, S2={"cnv_1p": -1})
    out = hotspot_routing_features(bps, calls, _WINDOWS)
    assert out.loc["S1", "dynamic_hotspot_breakpoints_all"] == 0
    assert out.loc["S2", "dynamic_hotspot_breakpoints_all"] == 1
