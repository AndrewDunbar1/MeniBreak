#!/usr/bin/env python3
"""Outcome-linked validation of the documented dynamic top-eight procedure.

This is not a frozen eight-window score. In every outer training fold, binary
breakpoint windows are ranked by a chi-square association test and the top
eight are fitted with the documented elastic-net logistic model. Held-out
patients never contribute to selection, scaling, or model fitting.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency
from sklearn.feature_selection import chi2
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler


HORIZONS = (3, 5, 7)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fixed_endpoint(time: pd.Series, event: pd.Series, horizon: int) -> np.ndarray:
    t = pd.to_numeric(time, errors="coerce").to_numpy()
    e = pd.to_numeric(event, errors="coerce").to_numpy()
    return np.where((e == 1) & (t <= horizon), 1, np.where(t >= horizon, 0, np.nan))


def direct_auc(y: np.ndarray, score: pd.Series) -> float:
    x = pd.to_numeric(score, errors="coerce").to_numpy()
    ok = np.isfinite(x) & np.isfinite(y)
    return float(roc_auc_score(y[ok].astype(int), x[ok])) if len(np.unique(y[ok])) == 2 else np.nan


def candidate_columns(frame: pd.DataFrame, mode: str) -> list[str]:
    any_pattern = re.compile(r"^bp_[0-9]+[pq]_bin[0-9]{2}$")
    typed_pattern = re.compile(r"^bp_(loss|gain)_[0-9]+[pq]_bin[0-9]{2}$")
    any_columns = [name for name in frame.columns if any_pattern.match(name)]
    if mode == "any":
        return any_columns
    return any_columns + [name for name in frame.columns if typed_pattern.match(name)]


def restricted_mask(frame: pd.DataFrame) -> np.ndarray | None:
    required = ("primary", "gtr", "prior_rt", "adj_rt", "rt_induced", "angiomatous")
    if not all(name in frame.columns for name in required):
        return None
    numeric = {name: pd.to_numeric(frame[name], errors="coerce").to_numpy() for name in required}
    return (
        (numeric["primary"] == 1)
        & (numeric["gtr"] == 1)
        & (numeric["prior_rt"] == 0)
        & (numeric["adj_rt"] == 0)
        & (numeric["rt_induced"] == 0)
        & (numeric["angiomatous"] == 0)
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", required=True, type=Path)
    parser.add_argument("--outcomes", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--sample-id-column", default="sample_id")
    parser.add_argument("--patient-id-column", default="patient_id")
    parser.add_argument("--time-column", default="recurrence_time")
    parser.add_argument("--event-column", default="recurrence_event")
    parser.add_argument("--candidate-mode", choices=("any", "any-loss-gain"), default="any")
    parser.add_argument("--minimum-support", type=int, default=20)
    parser.add_argument("--top-k", type=int, default=8)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260803)
    parser.add_argument("--write-sample-predictions", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    features = pd.read_csv(args.features, low_memory=False)
    outcomes = pd.read_csv(args.outcomes, low_memory=False)
    if args.sample_id_column not in features or args.sample_id_column not in outcomes:
        raise ValueError(f"Both files must contain {args.sample_id_column!r}")
    if features[args.sample_id_column].duplicated().any() or outcomes[args.sample_id_column].duplicated().any():
        raise ValueError("Each input must contain one row per sample ID")
    data = features.merge(outcomes, on=args.sample_id_column, how="inner", suffixes=("", "_outcome"))
    if not len(data):
        raise ValueError("Feature and outcome files have no matching sample IDs")
    if args.patient_id_column not in data:
        data[args.patient_id_column] = data[args.sample_id_column]
    if "qc_status" in data:
        data = data[~data.qc_status.astype(str).str.startswith("FAIL_")].copy()

    columns = candidate_columns(data, args.candidate_mode)
    if not columns:
        raise ValueError("No bp_<arm>_binNN breakpoint candidate columns were found")
    matrix = (data[columns].apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy() > 0).astype(np.uint8)
    support = matrix.sum(axis=0)
    eligible = support >= args.minimum_support
    columns = list(np.asarray(columns)[eligible])
    matrix = matrix[:, eligible]
    if len(columns) < args.top_k:
        raise ValueError(f"Only {len(columns)} candidates meet minimum support; top-k is {args.top_k}")

    metric_rows: list[dict] = []
    selection_rows: list[dict] = []
    prediction_rows: list[dict] = []
    restriction = restricted_mask(data)

    for horizon in HORIZONS:
        endpoint = fixed_endpoint(data[args.time_column], data[args.event_column], horizon)
        known = np.isfinite(endpoint)
        y = endpoint[known].astype(int)
        x = matrix[known]
        groups = data.loc[known, args.patient_id_column].astype(str).to_numpy()
        samples = data.loc[known, args.sample_id_column].astype(str).to_numpy()
        restricted = restriction[known] if restriction is not None else None
        repeat_predictions: list[np.ndarray] = []
        selected = Counter()

        for repetition in range(1, args.repeats + 1):
            splitter = StratifiedGroupKFold(
                n_splits=args.folds, shuffle=True, random_state=args.seed + horizon * 1000 + repetition
            )
            prediction = np.full(len(y), np.nan)
            for fold, (train, test) in enumerate(splitter.split(x, y, groups), start=1):
                scores, _ = chi2(x[train], y[train])
                ranking = np.argsort(-np.nan_to_num(scores, nan=-1.0), kind="stable")
                chosen = ranking[: args.top_k]
                chosen_names = [columns[index] for index in chosen]
                selected.update(chosen_names)
                for name in chosen_names:
                    selection_rows.append(
                        {"horizon_years": horizon, "repetition": repetition, "fold": fold, "feature": name}
                    )
                scaler = StandardScaler().fit(x[train][:, chosen])
                train_x = scaler.transform(x[train][:, chosen])
                test_x = scaler.transform(x[test][:, chosen])
                model = LogisticRegression(
                    penalty="elasticnet", solver="saga", l1_ratio=0.25, C=0.1,
                    max_iter=20000, random_state=args.seed + repetition
                )
                model.fit(train_x, y[train])
                prediction[test] = model.predict_proba(test_x)[:, 1]
            repeat_predictions.append(prediction)
            for cohort, mask in (("Broad", np.ones(len(y), dtype=bool)), ("Restricted transport", restricted)):
                if mask is None or len(np.unique(y[mask])) < 2:
                    continue
                metric_rows.append(
                    {"cohort": cohort, "horizon_years": horizon, "repetition": repetition,
                     "model": "Dynamic top-eight breakpoint procedure", "n": int(mask.sum()),
                     "events": int(y[mask].sum()), "auc": float(roc_auc_score(y[mask], prediction[mask]))}
                )
            if args.write_sample_predictions:
                prediction_rows.extend(
                    {"sample_id": sample, "horizon_years": horizon, "repetition": repetition,
                     "endpoint": int(outcome), "dynamic_top8_score": float(score)}
                    for sample, outcome, score in zip(samples, y, prediction)
                )

        for label, name in (("Recorded WHO grade", "who_grade"),
                            ("Prespecified integrated grade", "integrated_grade")):
            if name not in data:
                continue
            score = data.loc[known, name]
            for cohort, mask in (("Broad", np.ones(len(y), dtype=bool)), ("Restricted", restricted)):
                if mask is None:
                    continue
                metric_rows.append(
                    {"cohort": cohort, "horizon_years": horizon, "repetition": 0, "model": label,
                     "n": int(mask.sum()), "events": int(y[mask].sum()), "auc": direct_auc(y[mask], score.iloc[np.where(mask)[0]])}
                )

    metrics = pd.DataFrame(metric_rows)
    summary = (
        metrics.groupby(["cohort", "horizon_years", "model"], dropna=False)
        .agg(n=("n", "first"), events=("events", "first"), repetitions=("auc", "count"),
             mean_auc=("auc", "mean"), sd_auc=("auc", "std"), min_auc=("auc", "min"), max_auc=("auc", "max"))
        .reset_index()
    )
    selections = pd.DataFrame(selection_rows)
    frequency = (
        selections.groupby(["horizon_years", "feature"]).size().rename("selected_folds").reset_index()
    )
    frequency["total_folds"] = args.repeats * args.folds
    frequency["selection_frequency"] = frequency.selected_folds / frequency.total_folds
    frequency = frequency.sort_values(["horizon_years", "selection_frequency", "feature"], ascending=[True, False, True])

    metrics.to_csv(args.output_dir / "auc_by_repetition.csv", index=False)
    summary.to_csv(args.output_dir / "auc_summary.csv", index=False)
    selections.to_csv(args.output_dir / "selected_features_by_fold.csv", index=False)
    frequency.to_csv(args.output_dir / "selection_frequency.csv", index=False)
    if args.write_sample_predictions:
        pd.DataFrame(prediction_rows).to_csv(args.output_dir / "sample_predictions.csv", index=False)
    manifest = {
        "workflow": "dynamic-top8-v1", "features_sha256": sha256(args.features),
        "outcomes_sha256": sha256(args.outcomes), "matched_samples": int(len(data)),
        "candidate_mode": args.candidate_mode, "candidate_features_before_support_filter": int(len(eligible)),
        "candidate_features_after_support_filter": int(len(columns)), "minimum_support": args.minimum_support,
        "top_k": args.top_k, "repeats": args.repeats, "folds": args.folds, "seed": args.seed,
        "selection": "chi-square ranking inside every outer training fold",
        "model": {"family": "binomial", "penalty": "elasticnet", "l1_ratio": 0.25,
                  "C": 0.1, "standardize": True},
        "restricted_estimate": "transport of broad-cohort out-of-fold predictions; no restricted refit",
        "sample_predictions_written": bool(args.write_sample_predictions),
    }
    (args.output_dir / "validation_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
