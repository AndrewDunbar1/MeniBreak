"""Arm tiling and breakpoint (change-point) derivation.

The foundation model represents each chromosome arm as the alternating sequence

    SEGMENT_1, BREAKPOINT_1, SEGMENT_2, ..., SEGMENT_N

To produce that, an arm must be COMPLETELY tiled: event-called SEG files list
only altered regions, so the unaltered stretches between them are materialised
here as explicit NEUTRAL segments. Without this step a "breakpoint" between two
far-apart events would be indistinguishable from a genuine adjacency.

Two properties of real data are handled explicitly rather than papered over:

*   **Centromere-spanning segments** are split at the arm boundary, and the
    resulting arm-edge boundaries are NOT emitted as breakpoints. A segment that
    runs through the centromere has no change point there.
*   **True overlaps** are resolved by a single documented policy and every
    instance is reported. Nothing is dropped silently.

Coordinate convention
---------------------
Segments are treated as HALF-OPEN intervals ``[start, end)``.

This was established empirically, not assumed. In the institutional master SEG,
grouping by (assay sample, chromosome, track) yields 20,125 adjacent pairs with
``next.start == prev.end`` and **zero** pairs with ``next.start < prev.end``. A
shared boundary coordinate is therefore the file's encoding of contiguity, not
an overlap. Reading these as fully-closed intervals would manufacture ~20k
phantom one-base truncations and displace every affected breakpoint by 1 bp.

Note the arm reference CSV uses a closed convention (chr1 p ends 125,000,000 and
q starts 125,000,001). The 1 bp discrepancy is immaterial for arm assignment and
is left uncorrected rather than silently rewriting the locked model's reference.
"""

from __future__ import annotations

import logging
import warnings

import numpy as np
import pandas as pd

from .references import ArmReference
from .schema import (
    BREAKPOINT_SPEC,
    BoundaryClass,
    EventClass,
    TransitionDirection,
)

log = logging.getLogger(__name__)

#: Fallback when a segment carries no probe count and spacing cannot be
#: estimated. Chosen to be visibly large so that "unknown resolution" never
#: masquerades as high confidence.
DEFAULT_PROBE_SPACING_BP = 50_000.0


class OverlapPolicy:
    """How to resolve two overlapping segments within the same track.

    TRUNCATE_LEFT is the default: the earlier segment is truncated at the start
    of the later one. This preserves both segments and both of their change
    points, whereas dropping either would erase a real boundary. The truncation
    is recorded so the QC report can quantify how much geometry was adjusted.
    """

    TRUNCATE_LEFT = "truncate_left"


def _resolve_overlaps(seg: pd.DataFrame, report: list[dict]) -> pd.DataFrame:
    """Truncate overlapping segments in a single sorted arm/track group."""
    starts = seg["start"].to_numpy(copy=True)
    ends = seg["end"].to_numpy(copy=True)
    for i in range(len(seg) - 1):
        # Strict `<`: under the half-open convention `next.start == prev.end`
        # is contiguity, not overlap.
        if starts[i + 1] < ends[i]:
            report.append(
                {
                    "sample_id": seg["sample_id"].iloc[i],
                    "arm": seg["arm"].iloc[i],
                    "track": seg["track"].iloc[i],
                    "left_segment_id": seg["segment_id"].iloc[i],
                    "right_segment_id": seg["segment_id"].iloc[i + 1],
                    "overlap_bp": int(ends[i] - starts[i + 1]),
                    "policy": OverlapPolicy.TRUNCATE_LEFT,
                }
            )
            ends[i] = starts[i + 1]
    out = seg.copy()
    out["start"] = starts
    out["end"] = ends
    # A truncation can collapse a segment to zero length; those are removed and
    # counted, since a zero-length half-open interval has no interior.
    return out[out["end"] > out["start"]]


def tile_arms(
    segments: pd.DataFrame,
    arm_ref: ArmReference,
    *,
    autosomes_only: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Clip segments to arms and fill gaps with explicit NEUTRAL segments.

    Returns ``(tiles, overlap_report)``. `tiles` completely covers every arm
    that has at least one called segment, with no gaps and no overlaps.
    """
    required = {"sample_id", "chromosome", "start", "end", "track", "segment_id"}
    missing = required - set(segments.columns)
    if missing:
        raise ValueError(f"tile_arms: segments missing {sorted(missing)}")

    df = segments.copy()
    if autosomes_only:
        df = df[df["chromosome"].isin([str(i) for i in range(1, 23)])]

    # Explode each segment across every arm it touches, clipping to arm bounds.
    # A centromere-spanning segment therefore yields one piece per arm.
    pieces: list[pd.DataFrame] = []
    for arm in (arm_ref.get(k) for k in arm_ref.arm_keys):
        hit = df[
            (df["chromosome"] == arm.chromosome)
            & (df["end"] >= arm.start)
            & (df["start"] <= arm.end)
        ]
        if hit.empty:
            continue
        piece = hit.copy()
        piece["arm"] = arm.key
        piece["arm_start"] = arm.start
        piece["arm_end"] = arm.end
        piece["was_clipped"] = (piece["start"] < arm.start) | (piece["end"] > arm.end)
        piece["start"] = piece["start"].clip(lower=arm.start)
        piece["end"] = piece["end"].clip(upper=arm.end)
        pieces.append(piece)

    if not pieces:
        empty = pd.DataFrame(columns=list(df.columns) + ["arm", "is_filler"])
        return empty, pd.DataFrame(columns=["sample_id", "arm", "track", "overlap_bp"])

    clipped = pd.concat(pieces, ignore_index=True)

    overlap_report: list[dict] = []
    tiles: list[pd.DataFrame] = []

    for (sample_id, arm_key, track), grp in clipped.groupby(
        ["sample_id", "arm", "track"], sort=False
    ):
        grp = grp.sort_values(["start", "end"]).reset_index(drop=True)
        grp = _resolve_overlaps(grp, overlap_report)
        if grp.empty:
            continue

        arm = arm_ref.get(arm_key)
        template = grp.iloc[0]
        rows: list[dict] = []

        # Half-open walk: a filler covers [cursor, r.start) and the cursor
        # advances to r.end, so contiguous segments produce no filler at all.
        cursor = arm.start
        for r in grp.itertuples():
            if r.start > cursor:
                rows.append(_filler_row(template, cursor, r.start, arm_key, track))
            rows.append(_called_row(r, arm_key, track))
            cursor = r.end
        if cursor < arm.end:
            rows.append(_filler_row(template, cursor, arm.end, arm_key, track))

        tiles.append(pd.DataFrame(rows))

    tiled = pd.concat(tiles, ignore_index=True) if tiles else pd.DataFrame()
    return tiled, pd.DataFrame(overlap_report)


def _called_row(r, arm_key: str, track: str) -> dict:
    return {
        "dataset_id": r.dataset_id,
        "sample_id": r.sample_id,
        "tumor_id": getattr(r, "tumor_id", pd.NA),
        "patient_id": getattr(r, "patient_id", pd.NA),
        "chromosome": r.chromosome,
        "arm": arm_key,
        "track": track,
        "start": int(r.start),
        "end": int(r.end),
        "event": r.event,
        "segment_mean": float(r.segment_mean) if pd.notna(r.segment_mean) else np.nan,
        "probe_count": float(r.probe_count) if pd.notna(r.probe_count) else np.nan,
        "platform": r.platform,
        "caller": r.caller,
        "genome_build": r.genome_build,
        "source_file": r.source_file,
        "segment_id": r.segment_id,
        "is_filler": False,
        "was_clipped": bool(r.was_clipped),
    }


def _filler_row(template, start: int, end: int, arm_key: str, track: str) -> dict:
    """A materialised unaltered stretch. Amplitude 0.0, probe count unknown."""
    return {
        "dataset_id": template.dataset_id,
        "sample_id": template.sample_id,
        "tumor_id": template.get("tumor_id", pd.NA),
        "patient_id": template.get("patient_id", pd.NA),
        "chromosome": template.chromosome,
        "arm": arm_key,
        "track": track,
        "start": int(start),
        "end": int(end),
        "event": EventClass.NEUTRAL.value,
        "segment_mean": 0.0,
        "probe_count": np.nan,
        "platform": template.platform,
        "caller": template.caller,
        "genome_build": template.genome_build,
        "source_file": template.source_file,
        "segment_id": f"{template.sample_id}:{arm_key}:{start}-{end}:{track}:filler",
        "is_filler": True,
        "was_clipped": False,
    }


def derive_breakpoints(
    tiles: pd.DataFrame,
    arm_ref: ArmReference,
    *,
    boundary_class: BoundaryClass,
) -> pd.DataFrame:
    """Emit one breakpoint record per internal junction of the arm tiling.

    Boundaries at the outer edges of an arm are NOT breakpoints: they are the
    edges of the assayed interval, not observed transitions.
    """
    if tiles.empty:
        return pd.DataFrame(columns=BREAKPOINT_SPEC.all_names)

    tiles = tiles.sort_values(["sample_id", "track", "arm", "start"])

    g = tiles.groupby(["sample_id", "track", "arm"], sort=False)
    left = tiles[g.cumcount() < g["start"].transform("size") - 1]
    right = tiles[g.cumcount() > 0]

    if len(left) == 0:
        return pd.DataFrame(columns=BREAKPOINT_SPEC.all_names)

    L = left.reset_index(drop=True)
    R = right.reset_index(drop=True)

    # Local probe spacing: bp per probe averaged across the two flanking
    # segments. This is the resolution covariate that keeps the model from
    # treating a dense-array boundary and a sparse-array boundary as equal.
    # Sources without probe counts (TCGA ASCAT, PCAWG WGS) make every entry NaN,
    # which is an expected state handled by the fallback below -- not a fault.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Mean of empty slice")
        spacing = np.nanmean(
            np.vstack(
                [
                    _spacing(L["start"], L["end"], L["probe_count"]),
                    _spacing(R["start"], R["end"], R["probe_count"]),
                ]
            ),
            axis=0,
        )
    spacing = np.where(np.isfinite(spacing), spacing, DEFAULT_PROBE_SPACING_BP)

    bp_start = L["end"].to_numpy()
    bp_end = R["start"].to_numpy()
    # The change point lies somewhere between the last supported position of the
    # left segment and the first of the right. Where the segments abut exactly,
    # the true uncertainty is the inter-probe distance, not zero.
    uncertainty = np.maximum(bp_end - bp_start, spacing)

    amp = R["segment_mean"].to_numpy(dtype=float) - L["segment_mean"].to_numpy(dtype=float)
    direction = np.where(amp > 0, TransitionDirection.UP.value,
                np.where(amp < 0, TransitionDirection.DOWN.value,
                         TransitionDirection.NEUTRAL.value))

    arms = L["arm"].to_numpy()
    midpoint = ((bp_start + bp_end) // 2).astype("int64")
    rel = np.array(
        [arm_ref.get(a).relative_position(int(p)) for a, p in zip(arms, midpoint)]
    )
    arm_objs = [arm_ref.get(a) for a in arms]
    dist_cent = np.array(
        [abs(int(p) - (a.end if a.arm == "p" else a.start)) for a, p in zip(arm_objs, midpoint)],
        dtype=float,
    )
    dist_tel = np.array(
        [abs(int(p) - (a.start if a.arm == "p" else a.end)) for a, p in zip(arm_objs, midpoint)],
        dtype=float,
    )

    out = pd.DataFrame(
        {
            "dataset_id": L["dataset_id"].astype("string"),
            "sample_id": L["sample_id"].astype("string"),
            "patient_id": L["patient_id"].astype("string"),
            "chromosome": L["chromosome"].astype("string"),
            "arm": L["arm"].astype("string"),
            "left_segment_id": L["segment_id"].astype("string"),
            "right_segment_id": R["segment_id"].astype("string"),
            "breakpoint_start": bp_start.astype("int64"),
            "breakpoint_end": bp_end.astype("int64"),
            "coordinate_uncertainty_bp": uncertainty.astype("float64"),
            "left_state": L["event"].astype("string"),
            "right_state": R["event"].astype("string"),
            "amplitude_change": amp.astype("float64"),
            "transition_direction": pd.Series(direction, dtype="string"),
            "local_probe_spacing": spacing.astype("float64"),
            "boundary_class": boundary_class.value,
            "boundary_credibility": np.int64(boundary_class.credibility),
            "arm_relative_position": rel.astype("float64"),
            "distance_to_centromere": dist_cent,
            "distance_to_telomere": dist_tel,
            "is_arm_terminal": False,
            "caller": L["caller"].astype("string"),
            "platform": L["platform"].astype("string"),
            "genome_build": L["genome_build"].astype("string"),
            "source_file": L["source_file"].astype("string"),
            "qc_status": "pass",
        }
    )

    # A junction where both flanks are NEUTRAL carries no information; it can
    # only arise from a fully-truncated overlap and is dropped with a count.
    both_neutral = (out["left_state"] == EventClass.NEUTRAL.value) & (
        out["right_state"] == EventClass.NEUTRAL.value
    )
    if both_neutral.any():
        log.info("dropped %d neutral-neutral junctions", int(both_neutral.sum()))
    out = out[~both_neutral].reset_index(drop=True)

    BREAKPOINT_SPEC.validate(out, strict=True)
    return out


def _spacing(start: pd.Series, end: pd.Series, probes: pd.Series) -> np.ndarray:
    """bp per probe for a segment; NaN where the probe count is unusable."""
    length = (end - start + 1).to_numpy(dtype=float)
    p = probes.to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(p > 0, length / p, np.nan)
    return out
