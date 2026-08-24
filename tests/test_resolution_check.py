"""Tests for the resolution preflight and the realism of the shipped example.

The panel scores breakpoints in ~2.8 Mb windows. Whether an incoming SEG file
can support that depends on the assay behind it, not on the file format. These
tests pin the preflight's behaviour at both ends and guard the example data
against drifting back to arm-scale blocks.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def _seg(tmp_path: Path, lengths_mb: list[float], n_samples: int = 30) -> Path:
    """A SEG file whose events have the requested sizes."""
    rng = np.random.default_rng(0)
    rows = []
    for s in range(n_samples):
        for Lmb in lengths_mb:
            c = int(rng.integers(1, 23))
            a = int(rng.integers(1_000_000, 40_000_000))
            rows.append([f"T{s}", str(c), a, a + int(Lmb * 1e6), "loss"])
    p = tmp_path / "seg.csv"
    pd.DataFrame(rows, columns=["sample_id", "chromosome", "start", "end", "event"]
                 ).to_csv(p, index=False)
    return p


def _run(seg: Path) -> tuple[int, str]:
    r = subprocess.run(
        [sys.executable, str(ROOT / "apply" / "check_resolution.py"),
         "--segments", str(seg)],
        capture_output=True, text=True, cwd=ROOT)
    return r.returncode, r.stdout


def test_fine_grained_data_is_suitable(tmp_path):
    code, out = _run(_seg(tmp_path, [0.3, 0.6, 1.0, 1.5, 2.0]))
    assert "SUITABLE" in out
    assert code == 0


def test_coarse_data_is_flagged_too_coarse(tmp_path):
    code, out = _run(_seg(tmp_path, [20.0, 30.0, 45.0, 60.0]))
    assert "TOO COARSE" in out
    assert code == 2


def test_marginal_data_is_flagged_marginal(tmp_path):
    # ~1.0-1.5x the 2.77 Mb window
    code, out = _run(_seg(tmp_path, [3.2, 3.5, 3.8]))
    assert "MARGINAL" in out
    assert code == 0


def test_verdict_reports_the_ratio(tmp_path):
    _, out = _run(_seg(tmp_path, [1.0, 1.2, 1.4]))
    assert "median event / window" in out
    assert "x" in out


# --------------------------------------------------------------------------
# The shipped example must resemble the platform the panel was trained on.
# Targets come from the institutional aCGH cohort (reports/platform_resolution
# .csv in the development tree): median event 1.26 Mb, 47% under 1 Mb,
# 15 events per sample. A toy set of uniform 3-40 Mb blocks reads as TOO COARSE
# and misleads anyone evaluating the pipeline.
# --------------------------------------------------------------------------

@pytest.fixture(scope="module")
def example_segments(tmp_path_factory) -> pd.DataFrame:
    d = tmp_path_factory.mktemp("example")
    subprocess.run(
        [sys.executable, str(ROOT / "data" / "example" / "generate_example.py"),
         "--n", "250", "--seed", "7", "--outdir", str(d)],
        capture_output=True, text=True, cwd=ROOT, check=True)
    return pd.read_csv(d / "synthetic_segments.csv")


def test_example_event_size_matches_acgh(example_segments):
    med_mb = ((example_segments["end"] - example_segments["start"]).median()) / 1e6
    assert 0.7 <= med_mb <= 2.5, f"median event {med_mb:.2f} Mb; aCGH is 1.26 Mb"


def test_example_has_focal_events(example_segments):
    frac = ((example_segments["end"] - example_segments["start"]) < 1e6).mean()
    assert frac > 0.30, f"only {frac:.0%} sub-Mb events; aCGH is 47%"


def test_example_event_burden_matches_acgh(example_segments):
    per = example_segments.groupby("sample_id").size().median()
    assert 8 <= per <= 25, f"{per} events/sample; aCGH median is 15"


def test_example_reads_as_suitable(tmp_path, example_segments):
    p = tmp_path / "ex.csv"
    example_segments.to_csv(p, index=False)
    code, out = _run(p)
    assert "SUITABLE" in out, out
    assert code == 0


def test_example_segments_are_well_formed(example_segments):
    assert (example_segments["end"] > example_segments["start"]).all()
    assert (example_segments["start"] > 0).all()


# --------------------------------------------------------------------------
# The enriched model is the strongest SEG-only model in the report. It needs
# `segment_mean` (amplitude block) and LOH/allelic rows (allelic block), both
# optional in the input spec -- so the example must carry them or it cannot
# exercise the model at all.
# --------------------------------------------------------------------------

def test_example_carries_segment_mean(example_segments):
    assert "segment_mean" in example_segments.columns
    assert example_segments["segment_mean"].notna().all()


def test_example_carries_an_allelic_track(example_segments):
    frac = example_segments["event"].isin(["loh", "allelic_imbalance"]).mean()
    assert frac > 0.4, f"only {frac:.0%} allelic rows; institutional is 63.5%"


def test_enriched_model_is_scored_on_the_example(tmp_path, example_segments):
    """The external ladder must include the enriched model, not just the panel."""
    import subprocess
    d = tmp_path / "run"
    seg = tmp_path / "s.csv"
    example_segments.to_csv(seg, index=False)
    oc = pd.DataFrame({
        "sample_id": sorted(example_segments["sample_id"].unique()),
    })
    rng = np.random.default_rng(1)
    oc["patient_id"] = oc["sample_id"]
    oc["recurrence_time"] = rng.uniform(0.5, 9.0, len(oc)).round(2)
    oc["recurrence_event"] = rng.integers(0, 2, len(oc))
    op = tmp_path / "o.csv"
    oc.to_csv(op, index=False)

    r = subprocess.run(
        [sys.executable, str(ROOT / "apply" / "menibreak_panel.py"),
         "--segments", str(seg), "--outcomes", str(op),
         "--repeats", "1", "--out", str(d)],
        capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stderr[-2000:]
    perf = pd.read_csv(d / "performance.csv")
    assert "SCNA + enriched breakpoints" in set(perf["model"])


def test_missing_segment_mean_warns_but_still_runs(tmp_path, example_segments):
    import subprocess
    seg = tmp_path / "s.csv"
    example_segments.drop(columns=["segment_mean"]).to_csv(seg, index=False)
    r = subprocess.run(
        [sys.executable, str(ROOT / "apply" / "check_resolution.py"),
         "--segments", str(seg)], capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0


# --------------------------------------------------------------------------
# Clinical input: unlocks the grade comparators, the Clinical+ models and the
# Restricted scope. Partial input must degrade to a documented subset of the
# ladder rather than silently redefining the restricted cohort.
# --------------------------------------------------------------------------

def _clinical_for(seg: pd.DataFrame, *, full: bool = True) -> pd.DataFrame:
    ids = sorted(seg["sample_id"].unique())
    rng = np.random.default_rng(5)
    c = pd.DataFrame({"sample_id": ids})
    c["age"] = rng.integers(30, 85, len(ids))
    c["mitoses"] = rng.integers(0, 8, len(ids))
    c["gtr"] = rng.integers(0, 2, len(ids))
    c["prior_rt"] = 0
    c["adj_rt"] = rng.integers(0, 2, len(ids))
    c["who_grade"] = rng.integers(1, 4, len(ids))
    c["integrated_grade"] = rng.integers(1, 4, len(ids))
    if full:
        c["primary"] = 1
        c["rt_induced"] = 0
        c["angiomatous"] = 0
    return c


def _run_panel(tmp_path, seg_df, clin_df=None, repeats="1"):
    import subprocess
    seg = tmp_path / "s.csv"; seg_df.to_csv(seg, index=False)
    ids = sorted(seg_df["sample_id"].unique())
    rng = np.random.default_rng(2)
    oc = pd.DataFrame({"sample_id": ids, "patient_id": ids,
                       "recurrence_time": rng.uniform(0.5, 9.0, len(ids)).round(2),
                       "recurrence_event": rng.integers(0, 2, len(ids))})
    op = tmp_path / "o.csv"; oc.to_csv(op, index=False)
    cmd = [sys.executable, str(ROOT / "apply" / "menibreak_panel.py"),
           "--segments", str(seg), "--outcomes", str(op),
           "--repeats", repeats, "--out", str(tmp_path / "run")]
    if clin_df is not None:
        cp = tmp_path / "c.csv"; clin_df.to_csv(cp, index=False)
        cmd += ["--clinical", str(cp)]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stderr[-2000:]
    return pd.read_csv(tmp_path / "run" / "performance.csv"), r.stderr


def test_clinical_unlocks_the_full_ladder(tmp_path, example_segments):
    perf, _ = _run_panel(tmp_path, example_segments, _clinical_for(example_segments))
    models = set(perf["model"])
    for m in ("WHO grade", "Integrated grade", "Clinical + 8-window panel",
              "Clinical + SCNA + 8-window panel",
              "Clinical + SCNA + enriched breakpoints"):
        assert m in models, f"{m} missing from {sorted(models)}"


def test_clinical_unlocks_the_restricted_scope(tmp_path, example_segments):
    perf, _ = _run_panel(tmp_path, example_segments, _clinical_for(example_segments))
    assert "Restricted" in set(perf["scope"])


def test_without_clinical_only_seg_models_and_broad(tmp_path, example_segments):
    perf, _ = _run_panel(tmp_path, example_segments, None)
    assert set(perf["scope"]) == {"Broad"}
    assert not any(m.startswith("Clinical") for m in perf["model"])


def test_partial_restricted_criteria_refuses_the_scope(tmp_path, example_segments):
    """Five of six criteria must not silently produce a different subset."""
    clin = _clinical_for(example_segments, full=True).drop(columns=["angiomatous"])
    perf, err = _run_panel(tmp_path, example_segments, clin)
    assert set(perf["scope"]) == {"Broad"}
    assert "angiomatous" in err


def test_grade_is_a_comparator_not_a_predictor(tmp_path, example_segments):
    """who_grade / integrated_grade must never enter a genomic model."""
    from src.features.clinical import DELIBERATELY_EXCLUDED
    assert "who_grade" in DELIBERATELY_EXCLUDED
    assert "integrated_grade" in DELIBERATELY_EXCLUDED
