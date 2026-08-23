"""Deterministic, patient-grouped data splitting.

Every split in this project is at the PATIENT level. The institutional cohort
contains 75 patients with more than one tumour record, so a specimen-level split
would train and test on different resections from the same person -- the exact
leak the locked model's nested design was built to avoid.

All functions here are pure and seeded; the same seed always yields the same
partition, and the resulting assignment is written to disk so a run can be
audited after the fact.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def _patient_order(patients: np.ndarray, seed: int) -> np.ndarray:
    """Deterministically shuffle unique patients.

    Uses a seeded Generator rather than hashing the id strings, so the split is
    reproducible across Python versions (string hashing is salted per process).
    """
    rng = np.random.default_rng(seed)
    uniq = np.array(sorted(pd.unique(patients)))
    rng.shuffle(uniq)
    return uniq


def patient_grouped_folds(
    df: pd.DataFrame,
    *,
    n_splits: int = 5,
    seed: int = 0,
    patient_col: str = "patient_id",
) -> list[tuple[np.ndarray, np.ndarray]]:
    """K-fold split in which all rows of a patient stay in one fold.

    Returns a list of ``(train_positional_idx, test_positional_idx)`` pairs.
    Indices are POSITIONAL (usable with ``.iloc``) so callers cannot accidentally
    mix them with a non-default DataFrame index.
    """
    if patient_col not in df.columns:
        raise ValueError(f"missing patient column {patient_col!r}")
    if df[patient_col].isna().any():
        raise ValueError(
            "patient_id contains nulls; a null patient cannot be grouped and "
            "would silently be treated as one shared pseudo-patient"
        )

    patients = df[patient_col].to_numpy()
    order = _patient_order(patients, seed)
    # Round-robin patients into folds: keeps fold sizes balanced even when
    # patients contribute different numbers of tumours.
    fold_of_patient = {p: i % n_splits for i, p in enumerate(order)}
    fold_id = np.array([fold_of_patient[p] for p in patients])

    positions = np.arange(len(df))
    folds = []
    for k in range(n_splits):
        test = positions[fold_id == k]
        train = positions[fold_id != k]
        folds.append((train, test))
    return folds


def assert_no_patient_leakage(
    df: pd.DataFrame,
    train_idx,
    test_idx,
    *,
    patient_col: str = "patient_id",
    duplicate_col: str | None = "suspected_duplicate_of",
) -> None:
    """Raise AssertionError if a patient or suspected duplicate spans the split."""
    train_p = set(df.iloc[list(train_idx)][patient_col])
    test_p = set(df.iloc[list(test_idx)][patient_col])
    shared = train_p & test_p
    if shared:
        raise AssertionError(
            f"patient leakage: {len(shared)} patient(s) in both splits, "
            f"e.g. {sorted(shared)[:5]}"
        )

    if duplicate_col and duplicate_col in df.columns:
        tr = set(df.iloc[list(train_idx)][duplicate_col].dropna())
        te = set(df.iloc[list(test_idx)][duplicate_col].dropna())
        dup_shared = tr & te
        if dup_shared:
            raise AssertionError(
                f"suspected-duplicate leakage across splits: {sorted(dup_shared)[:5]}"
            )


def record_split(
    df: pd.DataFrame,
    folds: list[tuple[np.ndarray, np.ndarray]],
    out_path: str | Path,
    *,
    seed: int,
    patient_col: str = "patient_id",
) -> dict:
    """Persist a split assignment plus a checksum, for after-the-fact audit."""
    assignment = np.empty(len(df), dtype=int)
    for k, (_, test_idx) in enumerate(folds):
        assignment[test_idx] = k

    table = pd.DataFrame(
        {patient_col: df[patient_col].to_numpy(), "fold": assignment}
    )
    digest = hashlib.sha256(
        table.to_csv(index=False).encode("utf-8")
    ).hexdigest()

    meta = {
        "seed": seed,
        "n_splits": len(folds),
        "n_rows": len(df),
        "n_patients": int(df[patient_col].nunique()),
        "sha256": digest,
    }
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_path.with_suffix(".csv"), index=False)
    out_path.with_suffix(".json").write_text(json.dumps(meta, indent=2))
    return meta
