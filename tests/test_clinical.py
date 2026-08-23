"""Tests for the clinical-integrated comparator's covariate set.

The locked project has two genomic benchmarks and conflating them would
misstate every downstream claim:

    SEG-only            0.7792 / 0.7880 / 0.7798  (genomic features alone)
    clinical-integrated 0.8464 / 0.8442 / 0.8418  (+ five clinical covariates)

These tests pin the covariate set so it cannot drift.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.features.clinical import (
    CLINICAL_FEATURES,
    CLINICAL_SOURCE_COLUMNS,
    DELIBERATELY_EXCLUDED,
    assert_no_excluded_covariates,
    load_clinical_covariates,
)


def test_exactly_five_clinical_covariates():
    """Age, log(1+mitoses), GTR, prior RT, adjuvant RT -- and nothing else."""
    assert CLINICAL_FEATURES == ["age", "log1p_mitoses", "gtr", "prior_rt", "adj_rt"]
    assert len(CLINICAL_FEATURES) == 5


def test_sex_is_not_a_covariate():
    """Sex is NOT in the locked clinical-integrated model. It is a plausible
    thing to add by mistake, so the omission is asserted explicitly."""
    assert "sex" not in CLINICAL_FEATURES
    assert "sex" in DELIBERATELY_EXCLUDED
    assert "sex" not in CLINICAL_SOURCE_COLUMNS


@pytest.mark.parametrize("banned", ["who_grade", "integrated_grade"])
def test_grade_is_a_comparator_not_a_predictor(banned):
    """WHO and integrated grade are comparator SCORES. Using either as a
    predictor would inflate the genomic models and destroy the comparison."""
    assert banned in DELIBERATELY_EXCLUDED
    with pytest.raises(AssertionError, match=banned):
        assert_no_excluded_covariates(["cnv_1p", banned])


def test_guard_allows_a_legitimate_feature_matrix():
    assert_no_excluded_covariates(
        ["cnv_1p", "arm_loss_count", "total_breakpoint_boundaries"] + CLINICAL_FEATURES
    )


def test_guard_is_case_insensitive():
    with pytest.raises(AssertionError):
        assert_no_excluded_covariates(["CNV_1P", "WHO_GRADE"])


def test_mitoses_are_log1p_transformed(tmp_path):
    """The specification says log(1 + mitotic count); a raw count would give a
    different model, and log(0) would be -inf."""
    csv = tmp_path / "clin.csv"
    pd.DataFrame({
        "cg": ["T1", "T2", "T3"],
        "age": [50.0, 60.0, 70.0],
        "mitoses": [0.0, 3.0, 20.0],
        "gtr_num": [1, 0, 1],
        "prior_rt_num": [0, 1, 0],
        "adj_rt_num": [0, 0, 1],
    }).to_csv(csv, index=False)

    out = load_clinical_covariates(csv)
    assert np.isfinite(out["log1p_mitoses"]).all()
    assert out.loc["T1", "log1p_mitoses"] == pytest.approx(0.0)       # log1p(0)
    assert out.loc["T2", "log1p_mitoses"] == pytest.approx(np.log(4))
    assert out.loc["T3", "log1p_mitoses"] == pytest.approx(np.log(21))


def test_only_declared_columns_are_read(tmp_path):
    """An outcome column adjacent in the source file must not come along."""
    csv = tmp_path / "clin.csv"
    pd.DataFrame({
        "cg": ["T1"], "age": [50.0], "mitoses": [1.0],
        "gtr_num": [1], "prior_rt_num": [0], "adj_rt_num": [0],
        "rfs_event": [1], "who_grade": [3], "sex": ["F"],   # must not appear
    }).to_csv(csv, index=False)

    out = load_clinical_covariates(csv)
    assert list(out.columns) == CLINICAL_FEATURES
    for leaked in ("rfs_event", "who_grade", "sex"):
        assert leaked not in out.columns
    assert_no_excluded_covariates(out.columns)


def test_duplicate_tumour_ids_are_collapsed(tmp_path):
    csv = tmp_path / "clin.csv"
    pd.DataFrame({
        "cg": ["T1", "T1"], "age": [50.0, 50.0], "mitoses": [1.0, 1.0],
        "gtr_num": [1, 1], "prior_rt_num": [0, 0], "adj_rt_num": [0, 0],
    }).to_csv(csv, index=False)
    assert len(load_clinical_covariates(csv)) == 1
