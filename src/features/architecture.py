"""Hand-engineered arm-level and breakpoint-architecture features.

These reproduce the *families* of features in the locked 62-feature model, but
are computed independently from this project's canonical breakpoint tables. They
serve three purposes:

1. a **baseline** the foundation model must beat;
2. a **positive control** — if this pipeline is faithful, an elastic net on
   these features should land near the locked model's published AUCs
   (0.753 SCNA-only / 0.788 SCNA+breakpoint at 5 years);
3. an **interpretability reference** for whatever the route embedding learns.

This is NOT the locked model and must never be presented as it. The locked model
stays frozen and is compared against via its own out-of-fold predictions (D-008).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..harmonization.references import ArmReference
from ..harmonization.schema import COPY_NUMBER_EVENTS, EventClass

#: Arm coverage at which an arm is called altered. 0.30 is the locked model's
#: selected threshold (legacy preprocessing, 30% arm rule). Fixed here so the
#: baseline is comparable; it is NOT re-tuned.
ARM_ALTERED_THRESHOLD = 0.30

#: Events counted as CNA segments for burden features. LOH/AI are tracked
#: separately because they are copy-number neutral.
_CN_EVENT_VALUES = frozenset(e.value for e in COPY_NUMBER_EVENTS)
_LOSS_VALUES = frozenset({EventClass.LOSS.value, EventClass.HOMOZYGOUS_LOSS.value})
_GAIN_VALUES = frozenset({EventClass.GAIN.value, EventClass.HIGH_COPY_GAIN.value})


def arm_level_features(
    tiles: pd.DataFrame,
    arm_ref: ArmReference,
    *,
    threshold: float = ARM_ALTERED_THRESHOLD,
    sample_col: str = "sample_id",
) -> pd.DataFrame:
    """Per-sample arm-level loss/gain calls plus altered-arm burden.

    An arm is called lost or gained when the unioned coverage of the
    corresponding events reaches `threshold` of the arm's length.
    """
    called = tiles[~tiles["is_filler"]].copy()
    called["length"] = called["end"] - called["start"]
    called["is_loss"] = called["event"].isin(_LOSS_VALUES)
    called["is_gain"] = called["event"].isin(_GAIN_VALUES)

    arm_len = {k: arm_ref.get(k).length for k in arm_ref.arm_keys}

    rows: list[dict] = []
    for (sample, arm), g in called.groupby([sample_col, "arm"], sort=False):
        L = arm_len.get(arm)
        if not L:
            continue
        rows.append({
            sample_col: sample,
            "arm": arm,
            "loss_frac": g.loc[g["is_loss"], "length"].sum() / L,
            "gain_frac": g.loc[g["is_gain"], "length"].sum() / L,
        })
    if not rows:
        return pd.DataFrame()

    long = pd.DataFrame(rows)
    long["cnv"] = np.where(
        long["loss_frac"] >= threshold, -1,
        np.where(long["gain_frac"] >= threshold, 1, 0),
    )

    wide = long.pivot(index=sample_col, columns="arm", values="cnv")
    # Autosomal arms only, in a stable order so feature columns are reproducible.
    arms = sorted(
        [a for a in arm_ref.arm_keys if a[:-1].isdigit()],
        key=lambda a: (int(a[:-1]), a[-1]),
    )
    wide = wide.reindex(columns=arms).fillna(0).astype("int64")
    wide.columns = [f"cnv_{c}" for c in wide.columns]

    wide["arm_loss_count"] = (wide < 0).sum(axis=1)
    wide["arm_gain_count"] = (wide > 0).sum(axis=1)
    wide["arm_altered_count"] = wide["arm_loss_count"] + wide["arm_gain_count"]
    return wide


def breakpoint_features(
    bps: pd.DataFrame,
    tiles: pd.DataFrame,
    *,
    sample_col: str = "sample_id",
    small_segment_bp: int = 5_000_000,
) -> pd.DataFrame:
    """Breakpoint burden, spatial distribution and fragmentation features.

    Mirrors the locked model's breakpoint family: total boundaries, arms
    touched, distributional entropy across arms, maximum single-arm
    concentration, telomeric/central/centromeric fractions, arm-relative
    position moments, segment counts and small-segment burden.
    """
    if bps.empty:
        return pd.DataFrame()

    # The locked model counts boundaries of loss/gain/allelic-imbalance
    # segments; a junction qualifies if either flank is such an event.
    modelled = frozenset(_CN_EVENT_VALUES | {EventClass.ALLELIC_IMBALANCE.value})
    b = bps[bps["left_state"].isin(modelled) | bps["right_state"].isin(modelled)].copy()
    if b.empty:
        return pd.DataFrame()

    rel = b["arm_relative_position"]
    # 0 = telomere, 1 = centromere, so terciles read telomeric -> centromeric.
    b["zone"] = np.where(rel < 1 / 3, "tel", np.where(rel < 2 / 3, "central", "cent"))

    out = pd.DataFrame(index=sorted(b[sample_col].unique()))
    out.index.name = sample_col

    g = b.groupby(sample_col)
    out["total_breakpoint_boundaries"] = g.size()
    out["breakpoint_arms_touched"] = g["arm"].nunique()
    out["mean_rel_tel_to_cent"] = g["arm_relative_position"].mean()
    out["sd_rel_tel_to_cent"] = g["arm_relative_position"].std().fillna(0.0)
    # NOTE: `boundary_credibility` and `coordinate_uncertainty_bp` are NOT
    # architecture features. `boundary_credibility` is a deterministic function
    # of `boundary_class`, which is a deterministic function of platform, so
    # including it hands the platform label straight to the model. The Phase 3
    # negative control caught exactly this: it ranked as the single strongest
    # platform discriminator (|coef| 2.11).
    #
    # Coordinate uncertainty is likewise an assay-resolution property. Both are
    # exposed as NUISANCE covariates (see models/tokenizer.py) so they can be
    # conditioned on or adversarially removed -- never as free predictors.

    # Distribution of boundaries across arms: entropy and concentration.
    per_arm = b.groupby([sample_col, "arm"]).size().rename("n").reset_index()
    tot = per_arm.groupby(sample_col)["n"].transform("sum")
    per_arm["p"] = per_arm["n"] / tot
    per_arm["plogp"] = -per_arm["p"] * np.log(per_arm["p"].clip(lower=1e-12))
    out["breakpoint_route_entropy"] = per_arm.groupby(sample_col)["plogp"].sum()
    out["max_arm_breakpoint_fraction"] = per_arm.groupby(sample_col)["p"].max()

    zone = (
        b.groupby([sample_col, "zone"]).size().unstack(fill_value=0)
        .reindex(columns=["tel", "central", "cent"], fill_value=0)
    )
    zsum = zone.sum(axis=1).replace(0, np.nan)
    out["telomeric_breakpoint_fraction"] = zone["tel"] / zsum
    out["central_breakpoint_fraction"] = zone["central"] / zsum
    out["centromeric_breakpoint_fraction"] = zone["cent"] / zsum

    # Segment-level fragmentation, from the called (non-filler) tiles.
    called = tiles[~tiles["is_filler"]].copy()
    called["length"] = called["end"] - called["start"]
    cn = called[called["event"].isin(_CN_EVENT_VALUES)]
    cg = cn.groupby(sample_col)
    out["cna_segment_count"] = cg.size()
    out["loss_segment_count"] = (
        cn[cn["event"].isin(_LOSS_VALUES)].groupby(sample_col).size()
    )
    out["gain_segment_count"] = (
        cn[cn["event"].isin(_GAIN_VALUES)].groupby(sample_col).size()
    )
    out["median_segment_length"] = cg["length"].median()
    out["small_segment_count_lt5mb"] = (
        cn[cn["length"] < small_segment_bp].groupby(sample_col).size()
    )

    return out.fillna(0.0)


def drop_unidentified(
    df: pd.DataFrame, sample_col: str, *, what: str
) -> tuple[pd.DataFrame, int]:
    """Remove rows whose analysis unit is null, returning the count dropped.

    The institutional master SEG contains 50 assay samples with a null `cg`
    (17 named ``BL-*``, apparently blood/germline, and 33 named ``S##-####``,
    surgical accessions never assigned a deidentified id). They cannot be linked
    to an outcome, so they cannot enter a labelled analysis — but they remain
    perfectly usable for self-supervised pretraining, which needs no label.

    They are removed here with an explicit count rather than being allowed to
    become a single NA pseudo-sample that silently pools 50 unrelated specimens.
    """
    if sample_col not in df.columns:
        return df, 0
    bad = df[sample_col].isna()
    n = int(bad.sum())
    if n:
        import logging
        logging.getLogger(__name__).info(
            "%s: dropped %d rows with null %s (unlinkable specimens; "
            "still usable for unsupervised pretraining)", what, n, sample_col
        )
    return df[~bad], n


def hotspot_routing_features(
    bps: pd.DataFrame,
    arm_calls: pd.DataFrame,
    hotspot_windows: pd.DataFrame,
    *,
    sample_col: str = "sample_id",
) -> pd.DataFrame:
    """Location-specific routing: breakpoints inside frozen hotspot windows.

    This is the feature family that distinguishes *where* a breakpoint lands
    from *how many* there are. Following the locked model's rule, a hotspot
    breakpoint is counted only where the arm carries the required LOSS call AND
    the boundary falls inside that arm's frozen arm-relative window.

    The `non1p` variants exclude chr1p, whose 112-123 Mb hotspot dominates this
    disease. They test whether routing carries signal beyond the single
    best-known locus.

    IMPORTANT: the frozen windows are used here to reproduce the locked BASELINE
    and as a positive control. They must never be used as a training label for
    the foundation model, and any newly discovered window must be derived inside
    training folds only (specification section 10).
    """
    if bps.empty or hotspot_windows.empty:
        return pd.DataFrame()

    win = hotspot_windows.set_index("arm_key")
    b = bps[bps["arm"].isin(win.index)].copy()
    if b.empty:
        return pd.DataFrame()

    b["_lo"] = b["arm"].map(win["hotspot_start_norm"])
    b["_hi"] = b["arm"].map(win["hotspot_end_norm"])
    in_window = (b["arm_relative_position"] >= b["_lo"]) & (
        b["arm_relative_position"] <= b["_hi"]
    )

    # The arm must actually be lost for a routing event to be meaningful.
    lost = {
        (s, a.replace("cnv_", ""))
        for s, row in arm_calls.iterrows()
        for a, v in row.items()
        if a.startswith("cnv_") and v < 0
    }
    has_loss = pd.Series(
        [(s, a) in lost for s, a in zip(b[sample_col], b["arm"])], index=b.index
    )

    hit = b[in_window & has_loss]
    idx = sorted(bps[sample_col].dropna().unique())
    out = pd.DataFrame(index=idx)
    out.index.name = sample_col

    g = hit.groupby(sample_col)
    out["dynamic_hotspot_breakpoints_all"] = g.size()
    out["dynamic_hotspot_arms_used_all"] = g["arm"].nunique()

    non1p = hit[hit["arm"] != "1p"]
    g2 = non1p.groupby(sample_col)
    out["dynamic_hotspot_breakpoints_non1p"] = g2.size()
    out["dynamic_hotspot_arms_used_non1p"] = g2["arm"].nunique()

    return out.fillna(0.0)


def build_feature_table(
    tiles: pd.DataFrame,
    bps: pd.DataFrame,
    arm_ref: ArmReference,
    *,
    sample_col: str = "sample_id",
    hotspot_windows: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Join arm-level, breakpoint and (optionally) hotspot-routing features."""
    tiles, _ = drop_unidentified(tiles, sample_col, what="tiles")
    bps, _ = drop_unidentified(bps, sample_col, what="breakpoints")

    arm = arm_level_features(tiles, arm_ref, sample_col=sample_col)
    bpf = breakpoint_features(bps, tiles, sample_col=sample_col)
    out = arm.join(bpf, how="outer")

    if hotspot_windows is not None:
        hot = hotspot_routing_features(
            bps, arm, hotspot_windows, sample_col=sample_col
        )
        if not hot.empty:
            out = out.join(hot, how="outer")

    out = out.fillna(0.0)
    out.index.name = sample_col
    return out


#: Column groups used to define the nested baselines. Kept explicit so an
#: evaluation cannot accidentally give the SCNA-only model breakpoint
#: information, which would erase the very increment being measured.
def feature_groups(columns) -> dict[str, list[str]]:
    cols = list(columns)
    scna = [c for c in cols if c.startswith("cnv_") or c.startswith("arm_")]
    breakpoint = [c for c in cols if c not in scna]
    return {"scna": scna, "breakpoint": breakpoint, "full": scna + breakpoint}
