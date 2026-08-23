"""Leakage and quarantine guardrails.

These tests are adversarial by design: each one is written to FAIL if a
plausible mistake has been made. They implement section 12 of the project
specification and go/no-go acceptance tests 1 and 2.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.evaluation.splits import assert_no_patient_leakage, patient_grouped_folds

CONFIGS = ROOT / "configs"


# --------------------------------------------------------------------------
# Quarantine
# --------------------------------------------------------------------------

def _registry() -> dict:
    return yaml.safe_load(REGISTRY.read_text())


@pytest.fixture
def nested_cohort() -> pd.DataFrame:
    """A cohort where several patients contribute more than one tumour.

    This mirrors the institutional structure (75 patients with repeat
    resections), which is precisely where specimen-level splitting leaks.
    """
    rows = []
    for p in range(40):
        for t in range(1 + (p % 3)):        # 1-3 tumours per patient
            rows.append({"patient_id": f"P{p:03d}", "tumor_id": f"P{p:03d}_T{t}"})
    return pd.DataFrame(rows)


def test_patient_grouped_folds_keep_patients_together(nested_cohort):
    folds = patient_grouped_folds(nested_cohort, n_splits=5, seed=1)
    for train_idx, test_idx in folds:
        train_p = set(nested_cohort.iloc[train_idx].patient_id)
        test_p = set(nested_cohort.iloc[test_idx].patient_id)
        assert not (train_p & test_p), "a patient appeared in both train and test"


def test_every_tumor_is_tested_exactly_once(nested_cohort):
    folds = patient_grouped_folds(nested_cohort, n_splits=5, seed=1)
    seen: list[int] = []
    for _, test_idx in folds:
        seen.extend(test_idx.tolist())
    assert sorted(seen) == list(range(len(nested_cohort)))


def test_splits_are_deterministic_under_seed(nested_cohort):
    a = patient_grouped_folds(nested_cohort, n_splits=5, seed=7)
    b = patient_grouped_folds(nested_cohort, n_splits=5, seed=7)
    for (_, ta), (_, tb) in zip(a, b):
        assert ta.tolist() == tb.tolist()


def test_leakage_detector_catches_a_deliberately_bad_split(nested_cohort):
    """The detector must fail on a naive specimen-level split.

    Without this test the guard could silently be a no-op.
    """
    n = len(nested_cohort)
    train_idx = list(range(0, n, 2))
    test_idx = list(range(1, n, 2))
    with pytest.raises(AssertionError, match="patient"):
        assert_no_patient_leakage(nested_cohort, train_idx, test_idx)


def test_leakage_detector_passes_a_good_split(nested_cohort):
    folds = patient_grouped_folds(nested_cohort, n_splits=5, seed=3)
    for train_idx, test_idx in folds:
        assert_no_patient_leakage(nested_cohort, train_idx, test_idx)


# --------------------------------------------------------------------------
# Outcome columns must not reach unsupervised inputs
# --------------------------------------------------------------------------

OUTCOME_TOKENS = (
    "recur", "rfs", "os_event", "os_time", "survival", "event_5y",
    "time_to_death", "deceased", "label",
)


def test_segment_schema_declares_no_outcome_columns():
    """Acceptance test: outcome information cannot enter unsupervised features."""
    from src.harmonization.schema import BREAKPOINT_SPEC, SEGMENT_SPEC

    for spec in (SEGMENT_SPEC, BREAKPOINT_SPEC):
        for col in spec.all_names:
            assert not any(tok in col.lower() for tok in OUTCOME_TOKENS), (
                f"{spec.name}.{col} looks like an outcome field"
            )


# --------------------------------------------------------------------------
# TCGA / PCAWG structural containment (D-011)
# --------------------------------------------------------------------------

#: Corpora that may contain the same donor as each other. Overlap between them
#: is neutralised by keeping them all on the same side of every split, so they
#: must never carry an evaluation role.
PRETRAINING_ONLY = {"tcga_snp6_ascat", "pcawg_consensus_cn"}

# NOTE: tests covering components not shipped in this repository (TCGA/PCAWG
# harmonization, the dataset registry, pretraining-corpus separation) were removed
# when the panel was split out. See docs/source_reports/ for the wider project.
