"""Enriched breakpoint features that exploit information the locked model discards.

The locked 62-feature model reduces each tumour to discrete event calls above a
30% arm threshold, and its breakpoint family is counts, dispersion and four
frozen hotspot windows. Direct inspection of the institutional master SEG shows
three substantial sources of information it never touches:

1. **Continuous amplitude.** `Probe_Median` is a real log-ratio per segment. The
   locked model sees only "loss"/"gain". Amplitude carries clonality and
   focality information: a -0.9 deletion in 30% of an arm is a different event
   from a -0.15 shallow loss over the same span.

2. **The allelic track.** LOH and allelic imbalance are 60% of all segment rows
   (103,029 + 22,776 of 172,091) and are excluded from the locked breakpoint
   features entirely. Copy-neutral LOH is a real structural transition that a
   total-copy-number view cannot see.

3. **Fine positional structure.** Routing is summarised by four hotspot counts.
   Phase 3 measured that *where* a breakpoint lands carries more prognostic
   signal than every global summary combined (D-012), so a finer positional
   representation is the obvious place to look for more.

4. **Probe density.** Available per segment and never used.

Every feature here is computed from this project's canonical tables, so it can
be recomputed inside a training fold. Nothing depends on the frozen windows.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..harmonization.schema import COPY_NUMBER_EVENTS, EventClass

_CN = frozenset(e.value for e in COPY_NUMBER_EVENTS)
_LOSS = frozenset({EventClass.LOSS.value, EventClass.HOMOZYGOUS_LOSS.value})
_GAIN = frozenset({EventClass.GAIN.value, EventClass.HIGH_COPY_GAIN.value})
_ALLELIC = frozenset({EventClass.LOH.value, EventClass.ALLELIC_IMBALANCE.value})

#: Positional bins per arm. 10 keeps the feature count tractable (390 columns)
#: while still resolving a hotspot an order of magnitude finer than the
#: telomeric/central/centromeric terciles the locked model uses.
N_POSITION_BINS = 10


def amplitude_features(
    tiles: pd.DataFrame, bps: pd.DataFrame, *, sample_col: str = "sample_id"
) -> pd.DataFrame:
    """Continuous-amplitude descriptors, ignored entirely by the locked model."""
    called = tiles[~tiles["is_filler"]].copy()
    called = called[called["segment_mean"].notna()]
    if called.empty:
        return pd.DataFrame()

    called["abs_amp"] = called["segment_mean"].abs()
    called["length"] = called["end"] - called["start"]

    g = called.groupby(sample_col)
    out = pd.DataFrame(index=sorted(called[sample_col].unique()))
    out.index.name = sample_col

    out["amp_mean_abs"] = g["abs_amp"].mean()
    out["amp_max_abs"] = g["abs_amp"].max()
    out["amp_sd"] = g["segment_mean"].std().fillna(0.0)
    out["amp_p90_abs"] = g["abs_amp"].quantile(0.9)
    # Length-weighted mean amplitude approximates genome-average ploidy shift.
    out["amp_weighted_mean"] = (
        called.assign(w=called["segment_mean"] * called["length"])
        .groupby(sample_col)["w"].sum()
        / g["length"].sum().replace(0, np.nan)
    )
    # Deep events are clonal/homozygous; shallow ones are subclonal or noise.
    out["n_deep_loss"] = called[called["segment_mean"] <= -0.5].groupby(sample_col).size()
    out["n_deep_gain"] = called[called["segment_mean"] >= 0.5].groupby(sample_col).size()
    out["n_shallow"] = called[called["abs_amp"] < 0.2].groupby(sample_col).size()
    out["deep_fraction"] = out["n_deep_loss"].fillna(0) + out["n_deep_gain"].fillna(0)
    out["deep_fraction"] = out["deep_fraction"] / g.size().replace(0, np.nan)

    # Amplitude CHANGE across breakpoints -- the step size of each transition.
    if not bps.empty and "amplitude_change" in bps.columns:
        b = bps[bps[sample_col].notna()].copy()
        b["abs_step"] = b["amplitude_change"].abs()
        bg = b.groupby(sample_col)
        out["bp_step_mean"] = bg["abs_step"].mean()
        out["bp_step_max"] = bg["abs_step"].max()
        out["bp_step_sd"] = bg["abs_step"].std()
        out["n_large_steps"] = b[b["abs_step"] >= 0.5].groupby(sample_col).size()

    return out.fillna(0.0)


def allelic_features(
    tiles: pd.DataFrame, bps: pd.DataFrame, *, sample_col: str = "sample_id"
) -> pd.DataFrame:
    """LOH / allelic-imbalance architecture -- 60% of the data, unused by the locked model."""
    called = tiles[~tiles["is_filler"]]
    allelic = called[called["event"].isin(_ALLELIC)]
    if allelic.empty:
        return pd.DataFrame()

    a = allelic.copy()
    a["length"] = a["end"] - a["start"]
    g = a.groupby(sample_col)

    out = pd.DataFrame(index=sorted(called[sample_col].unique()))
    out.index.name = sample_col
    out["loh_segment_count"] = a[a["event"] == EventClass.LOH.value].groupby(sample_col).size()
    out["ai_segment_count"] = (
        a[a["event"] == EventClass.ALLELIC_IMBALANCE.value].groupby(sample_col).size()
    )
    out["allelic_total_length"] = g["length"].sum()
    out["allelic_arms_touched"] = g["arm"].nunique()
    out["allelic_median_length"] = g["length"].median()

    # Copy-neutral LOH: allelic transition with no copy-number change. Invisible
    # to the locked model by construction.
    cn = called[called["event"].isin(_CN)]
    cn_arms = cn.groupby(sample_col)["arm"].apply(set) if not cn.empty else pd.Series(dtype=object)
    loh_arms = a.groupby(sample_col)["arm"].apply(set)
    cnn = {}
    for s, arms in loh_arms.items():
        cnn[s] = len(arms - cn_arms.get(s, set()))
    out["copy_neutral_loh_arms"] = pd.Series(cnn)

    # Allelic-track breakpoints, which the locked feature set excludes.
    if not bps.empty:
        ab = bps[bps["left_state"].isin(_ALLELIC) | bps["right_state"].isin(_ALLELIC)]
        if not ab.empty:
            out["allelic_breakpoints"] = ab.groupby(sample_col).size()
            out["allelic_bp_mean_pos"] = ab.groupby(sample_col)["arm_relative_position"].mean()
            out["allelic_bp_sd_pos"] = ab.groupby(sample_col)["arm_relative_position"].std()

    return out.fillna(0.0)


def positional_bin_features(
    bps: pd.DataFrame,
    *,
    sample_col: str = "sample_id",
    n_bins: int = N_POSITION_BINS,
    arms: list[str] | None = None,
) -> pd.DataFrame:
    """Breakpoint counts per (arm, positional bin).

    A data-driven alternative to the four frozen hotspot windows: the model is
    given the whole positional distribution and left to find its own structure.
    Derived purely from this project's breakpoint table, so it is recomputable
    inside a training fold and carries none of the frozen windows' optimism.
    """
    if bps.empty:
        return pd.DataFrame()

    b = bps[bps[sample_col].notna()].copy()
    b = b[b["arm_relative_position"].between(0, 1)]
    if b.empty:
        return pd.DataFrame()

    b["bin"] = np.minimum((b["arm_relative_position"] * n_bins).astype(int), n_bins - 1)
    b["key"] = "pbin_" + b["arm"].astype(str) + "_" + b["bin"].astype(str)

    wide = (
        b.groupby([sample_col, "key"]).size().unstack(fill_value=0)
    )
    if arms is not None:
        keep = [f"pbin_{a}_{i}" for a in arms for i in range(n_bins)]
        wide = wide.reindex(columns=keep, fill_value=0)
    return wide.astype("float64")


def density_features(
    tiles: pd.DataFrame, bps: pd.DataFrame, *, sample_col: str = "sample_id"
) -> pd.DataFrame:
    """Local clustering of breakpoints -- fragmentation beyond a raw count."""
    if bps.empty:
        return pd.DataFrame()
    b = bps[bps[sample_col].notna()].copy()
    out = pd.DataFrame(index=sorted(b[sample_col].unique()))
    out.index.name = sample_col

    # Inter-breakpoint spacing within an arm: small gaps mean clustered events
    # (chromothripsis-like), large gaps mean isolated ones.
    gaps = []
    for (s, arm), grp in b.groupby([sample_col, "arm"], sort=False):
        if len(grp) < 2:
            continue
        pos = np.sort(grp["breakpoint_start"].to_numpy(dtype=float))
        d = np.diff(pos)
        gaps.append(pd.DataFrame({sample_col: s, "gap": d}))
    if gaps:
        gd = pd.concat(gaps, ignore_index=True)
        gg = gd.groupby(sample_col)["gap"]
        out["bp_gap_median"] = gg.median()
        out["bp_gap_min"] = gg.min()
        out["n_close_pairs_1mb"] = gd[gd["gap"] < 1e6].groupby(sample_col).size()
        out["n_close_pairs_100kb"] = gd[gd["gap"] < 1e5].groupby(sample_col).size()

    # Maximum breakpoints within any single arm -- a focal-catastrophe signal.
    per_arm = b.groupby([sample_col, "arm"]).size()
    out["max_bp_single_arm"] = per_arm.groupby(level=0).max()
    out["mean_bp_per_touched_arm"] = per_arm.groupby(level=0).mean()

    return out.fillna(0.0)


def build_enriched(
    tiles: pd.DataFrame,
    bps: pd.DataFrame,
    *,
    sample_col: str = "sample_id",
    include_bins: bool = True,
) -> pd.DataFrame:
    """All enriched families joined into one per-sample table."""
    # 50 institutional assay samples carry a null cg and cannot be linked to a
    # tumour (17 blood-derived `BL-*`, 33 unassigned surgical accessions). Left
    # in, they collapse into a single NA pseudo-sample pooling 50 unrelated
    # specimens -- the same hazard guarded in features/architecture.py.
    from .architecture import drop_unidentified

    tiles, _ = drop_unidentified(tiles, sample_col, what="enriched tiles")
    bps, _ = drop_unidentified(bps, sample_col, what="enriched breakpoints")

    parts = [
        amplitude_features(tiles, bps, sample_col=sample_col),
        allelic_features(tiles, bps, sample_col=sample_col),
        density_features(tiles, bps, sample_col=sample_col),
    ]
    if include_bins:
        parts.append(positional_bin_features(bps, sample_col=sample_col))

    parts = [p for p in parts if not p.empty]
    if not parts:
        return pd.DataFrame()

    out = parts[0]
    for p in parts[1:]:
        out = out.join(p, how="outer")
    return out.fillna(0.0)
