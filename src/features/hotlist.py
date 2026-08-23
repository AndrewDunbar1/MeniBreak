"""Discovery of a parsimonious breakpoint 'hot list'.

Goal: identify a small number (<10) of specific breakpoint windows whose
presence adds prognostic information over WHO grade, integrated grade and
arm-level SCNA.

The methodological trap here is circularity. The locked model's four frozen
hotspot windows were derived on this same cohort, so any AUC that uses them is
optimistic — its own stability audit found the exact windows reproducible in only
38.3% of training folds (median Jaccard 0.60). Reusing them, or selecting new
windows on the full data and then scoring on the full data, measures how well the
cohort fits itself.

Everything here therefore obeys one rule:

    **Window selection happens strictly inside the training folds.**

Two outputs follow, and they answer different questions:

* an **honest AUC**, where each fold selects its own windows from its own
  training rows and is scored on rows it never saw;
* a **consensus hot list**, the windows that recur across folds, reported with
  the selection frequency that makes their reliability explicit.

A window appearing in 95% of folds is a finding. One appearing in 40% is a
coincidence, and the frozen windows' own audit suggests that is common.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

log = logging.getLogger(__name__)


@dataclass
class HotWindow:
    """One candidate breakpoint window."""

    arm: str
    bin_index: int
    n_bins: int

    @property
    def key(self) -> str:
        return f"bp_{self.arm}_b{self.bin_index}"

    @property
    def start_norm(self) -> float:
        return self.bin_index / self.n_bins

    @property
    def end_norm(self) -> float:
        return (self.bin_index + 1) / self.n_bins

    def __str__(self) -> str:
        return (f"{self.arm}:{self.start_norm:.2f}-{self.end_norm:.2f} "
                f"(telomere=0, centromere=1)")


def build_window_matrix(
    bps: pd.DataFrame,
    *,
    sample_col: str = "sample_id",
    n_bins: int = 20,
    binary: bool = True,
    min_prevalence: int = 20,
) -> pd.DataFrame:
    """Per-sample indicator (or count) for every (arm x positional bin) window.

    `binary` gives presence/absence, which is what a clinically usable panel
    would record. `min_prevalence` drops windows too rare to estimate, which
    also keeps the multiple-testing burden honest.
    """
    b = bps[bps[sample_col].notna()].copy()
    b = b[b["arm_relative_position"].between(0, 1)]
    if b.empty:
        return pd.DataFrame()

    b["bin"] = np.minimum((b["arm_relative_position"] * n_bins).astype(int), n_bins - 1)
    b["key"] = "bp_" + b["arm"].astype(str) + "_b" + b["bin"].astype(str)

    wide = b.groupby([sample_col, "key"]).size().unstack(fill_value=0)
    if binary:
        wide = (wide > 0).astype("float64")

    keep = wide.columns[(wide > 0).sum(axis=0) >= min_prevalence]
    dropped = wide.shape[1] - len(keep)
    if dropped:
        log.info("dropped %d windows below prevalence %d; %d remain",
                 dropped, min_prevalence, len(keep))
    return wide[keep]


def rank_windows(
    X: pd.DataFrame, y: np.ndarray, *, method: str = "logrank_chi2"
) -> pd.Series:
    """Rank windows by univariate association with the outcome.

    Uses a chi-square test of independence for binary windows -- cheap enough to
    run inside every fold, which is the whole point. Returns a Series of scores
    indexed by window, descending.
    """
    scores = {}
    for col in X.columns:
        v = X[col].to_numpy()
        present = v > 0
        if present.sum() < 5 or (~present).sum() < 5:
            continue
        table = np.array([
            [int(((y == 1) & present).sum()), int(((y == 0) & present).sum())],
            [int(((y == 1) & ~present).sum()), int(((y == 0) & ~present).sum())],
        ])
        if table.min() < 1:
            continue
        chi2, p, _, _ = stats.chi2_contingency(table)
        scores[col] = chi2
    return pd.Series(scores).sort_values(ascending=False)


def select_in_fold(
    X_train: pd.DataFrame, y_train: np.ndarray, *, k: int
) -> list[str]:
    """Top-k windows chosen from TRAINING rows only."""
    ranked = rank_windows(X_train, y_train)
    return list(ranked.head(k).index)


def consensus_hotlist(
    selections: list[list[str]], *, n_bins: int = 20
) -> pd.DataFrame:
    """Aggregate per-fold selections into a stability-ranked hot list.

    `selection_frequency` is the fraction of folds that chose a window. This is
    the number that says whether a window is a finding or a coincidence, and it
    is reported alongside every candidate rather than hidden.
    """
    if not selections:
        return pd.DataFrame()
    n = len(selections)
    counts: dict[str, int] = {}
    for sel in selections:
        for w in sel:
            counts[w] = counts.get(w, 0) + 1

    rows = []
    for key, c in sorted(counts.items(), key=lambda kv: -kv[1]):
        arm, _, bin_part = key[3:].rpartition("_b")
        try:
            bi = int(bin_part)
        except ValueError:
            continue
        w = HotWindow(arm=arm, bin_index=bi, n_bins=n_bins)
        rows.append({
            "window": key,
            "arm": arm,
            "bin": bi,
            "start_norm": w.start_norm,
            "end_norm": w.end_norm,
            "n_folds_selected": c,
            "selection_frequency": c / n,
            "description": str(w),
        })
    return pd.DataFrame(rows)


def windows_to_genomic_coordinates(
    hotlist: pd.DataFrame, arm_ref, *, n_bins: int = 20
) -> pd.DataFrame:
    """Translate arm-relative windows into genomic Mb for interpretation.

    Arm-relative position runs 0 at the telomere to 1 at the centromere, so a
    p-arm window maps forward from the arm start and a q-arm window maps
    backward from the arm end.
    """
    out = []
    for r in hotlist.itertuples():
        try:
            arm = arm_ref.get(r.arm)
        except KeyError:
            continue
        if arm.arm == "p":
            g0 = arm.start + r.start_norm * arm.length
            g1 = arm.start + r.end_norm * arm.length
        else:
            g0 = arm.end - r.end_norm * arm.length
            g1 = arm.end - r.start_norm * arm.length
        out.append({
            **r._asdict(),
            "chrom": arm.chromosome,
            "genomic_start_mb": round(min(g0, g1) / 1e6, 2),
            "genomic_end_mb": round(max(g0, g1) / 1e6, 2),
        })
    return pd.DataFrame(out).drop(columns=["Index"], errors="ignore")
