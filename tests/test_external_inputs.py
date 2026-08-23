"""Tests for the external-cohort entry point.

These cover the ways a real collaborator's SEG export differs from the tidy
example file. Each case here was found by stress-testing the deployment path;
three of them were live defects.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from apply.menibreak_panel import adapt_segments, load_outcomes_table


@pytest.fixture
def outcomes() -> pd.DataFrame:
    return pd.DataFrame({
        "sample_id": [f"T{i}" for i in range(20)],
        "patient_id": [f"P{i}" for i in range(20)],
        "recurrence_time": [1.0, 2.5, 3.0, 4.5, 5.0, 6.0, 7.5, 8.0, 2.0, 3.5] * 2,
        "recurrence_event": [1, 0] * 10,
    })


@pytest.fixture
def segments() -> pd.DataFrame:
    return pd.DataFrame({
        "sample_id": ["T0", "T0", "T1"],
        "chromosome": ["1", "22", "6"],
        "start": [1_000_000, 20_000_000, 60_000_000],
        "end": [115_000_000, 45_000_000, 90_000_000],
        "event": ["loss", "loss", "gain"],
    })


# --------------------------------------------------------------------------
# patient_id is optional (INPUT_SPEC.md says so; the loader used to require it)
# --------------------------------------------------------------------------

def test_outcomes_without_patient_id_defaults_to_sample_id(tmp_path, outcomes):
    p = tmp_path / "o.csv"
    outcomes.drop(columns=["patient_id"]).to_csv(p, index=False)
    oc = load_outcomes_table(p)
    assert (oc["patient_id"] == oc["tumor_id"]).all()


def test_outcomes_with_patient_id_is_respected(tmp_path, outcomes):
    p = tmp_path / "o.csv"
    outcomes.to_csv(p, index=False)
    oc = load_outcomes_table(p)
    assert oc["patient_id"].tolist() == outcomes["patient_id"].tolist()


# --------------------------------------------------------------------------
# Follow-up time units. Passing months does not crash -- it silently returns a
# meaningless AUC on a handful of events, which is worse than failing.
# --------------------------------------------------------------------------

def test_months_instead_of_years_is_rejected(tmp_path, outcomes):
    p = tmp_path / "o.csv"
    o = outcomes.copy()
    o["recurrence_time"] = o["recurrence_time"] * 12
    o.to_csv(p, index=False)
    with pytest.raises(ValueError, match="implausible in YEARS"):
        load_outcomes_table(p)


def test_days_instead_of_years_is_rejected(tmp_path, outcomes):
    p = tmp_path / "o.csv"
    o = outcomes.copy()
    o["recurrence_time"] = o["recurrence_time"] * 365.25
    o.to_csv(p, index=False)
    with pytest.raises(ValueError, match="implausible in YEARS"):
        load_outcomes_table(p)


def test_plausible_years_are_accepted(tmp_path, outcomes):
    p = tmp_path / "o.csv"
    outcomes.to_csv(p, index=False)
    assert len(load_outcomes_table(p)) == len(outcomes)


def test_non_binary_event_is_rejected(tmp_path, outcomes):
    p = tmp_path / "o.csv"
    o = outcomes.copy()
    o.loc[0, "recurrence_event"] = 2
    o.to_csv(p, index=False)
    with pytest.raises(ValueError, match="must be 0/1"):
        load_outcomes_table(p)


def test_negative_time_is_rejected(tmp_path, outcomes):
    p = tmp_path / "o.csv"
    o = outcomes.copy()
    o.loc[0, "recurrence_time"] = -1.0
    o.to_csv(p, index=False)
    with pytest.raises(ValueError, match="negative"):
        load_outcomes_table(p)


def test_missing_required_outcome_column_names_the_column(tmp_path, outcomes):
    p = tmp_path / "o.csv"
    outcomes.drop(columns=["recurrence_event"]).to_csv(p, index=False)
    with pytest.raises(ValueError, match="recurrence_event"):
        load_outcomes_table(p)


# --------------------------------------------------------------------------
# Segment file shape problems
# --------------------------------------------------------------------------

def test_tab_separated_file_names_the_delimiter(tmp_path, segments):
    p = tmp_path / "s.csv"
    segments.to_csv(p, sep="\t", index=False)
    with pytest.raises(ValueError, match="tab"):
        adapt_segments(p)


def test_semicolon_separated_file_names_the_delimiter(tmp_path, segments):
    p = tmp_path / "s.csv"
    segments.to_csv(p, sep=";", index=False)
    with pytest.raises(ValueError, match="semicolon"):
        adapt_segments(p)


def test_empty_segments_file_says_so(tmp_path, segments):
    p = tmp_path / "s.csv"
    segments.iloc[:0].to_csv(p, index=False)
    with pytest.raises(ValueError, match="no rows"):
        adapt_segments(p)


def test_missing_segment_column_lists_what_was_found(tmp_path, segments):
    p = tmp_path / "s.csv"
    segments.drop(columns=["event"]).to_csv(p, index=False)
    with pytest.raises(ValueError, match="Found:"):
        adapt_segments(p)


# --------------------------------------------------------------------------
# Formats a real cohort will actually arrive in
# --------------------------------------------------------------------------

def test_chr_prefixed_chromosomes(tmp_path, segments):
    p = tmp_path / "s.csv"
    s = segments.copy()
    s["chromosome"] = "chr" + s["chromosome"].astype(str)
    s.to_csv(p, index=False)
    assert len(adapt_segments(p)) == len(segments)


@pytest.mark.parametrize("vocab", [
    {"loss": "Deletion", "gain": "Amplification"},
    {"loss": "del", "gain": "amp"},
    {"loss": -1, "gain": 1},
])
def test_vendor_event_vocabularies(tmp_path, segments, vocab):
    p = tmp_path / "s.csv"
    s = segments.copy()
    s["event"] = s["event"].map(vocab)
    s.to_csv(p, index=False)
    assert len(adapt_segments(p)) == len(segments)


def test_unknown_event_label_is_rejected(tmp_path, segments):
    p = tmp_path / "s.csv"
    s = segments.copy()
    s.loc[0, "event"] = "totally_unknown_call"
    s.to_csv(p, index=False)
    with pytest.raises(ValueError, match="unrecognised event label"):
        adapt_segments(p)
