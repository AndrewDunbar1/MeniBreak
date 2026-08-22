#!/usr/bin/env python3
"""Forensic reconstruction of the historical six-AUC shortlist fingerprint.

This script deliberately varies only ambiguities found in retained artifacts:
candidate family, support rule, k=5 versus k=8, binary versus count-valued
features, the stored cv_fold versus reconstructed patient folds, and standard
elastic-net strength. It ranks configurations against all six reported AUCs.
"""

from __future__ import annotations

import argparse
import itertools
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.feature_selection import chi2, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler


TARGET = {(3, "broad"): .702, (5, "broad"): .735, (7, "broad"): .752,
          (3, "restricted"): .715, (5, "restricted"): .753, (7, "restricted"): .757}


def endpoint(df: pd.DataFrame, horizon: int) -> np.ndarray:
    t = pd.to_numeric(df.rfs_time_years, errors="coerce").to_numpy()
    e = pd.to_numeric(df.rfs_event, errors="coerce").to_numpy()
    return np.where((e == 1) & (t <= horizon), 1., np.where(t >= horizon, 0., np.nan))


def restricted(df: pd.DataFrame) -> np.ndarray:
    n = lambda c: pd.to_numeric(df[c], errors="coerce").to_numpy()
    return ((n("primary") == 1) & (n("gtr") == 1) & (n("prior_rt") == 0) &
            (n("adj_rt") == 0) & (n("rt_induced") == 0) & (n("angiomatous") == 0))


def candidate_columns(df: pd.DataFrame, family: str) -> list[str]:
    any_re = re.compile(r"^bp_[0-9]+[pq]_bin[0-9]{2}$")
    typed_re = re.compile(r"^bp_(loss|gain)_[0-9]+[pq]_bin[0-9]{2}$")
    any_cols = [c for c in df if any_re.match(c)]
    typed_cols = [c for c in df if typed_re.match(c)]
    if family == "any":
        return any_cols
    if family == "typed":
        return typed_cols
    return any_cols + typed_cols


def cyclic_patient_folds(groups: np.ndarray, seed: int) -> np.ndarray:
    ids = np.unique(groups.astype(str))
    ids = np.random.default_rng(seed).permutation(ids)
    mapping = {patient: i % 5 for i, patient in enumerate(ids)}
    return np.asarray([mapping[str(g)] for g in groups])


def auc(y: np.ndarray, score: np.ndarray, mask: np.ndarray) -> float:
    return float(roc_auc_score(y[mask], score[mask]))


def evaluate(df: pd.DataFrame, raw: np.ndarray, names: list[str], config: dict) -> tuple[dict, list[dict]]:
    binary = (raw > 0).astype(np.uint8)
    support = binary.sum(axis=0)
    keep = support >= config["support"] if config["support_op"] == "ge" else support > config["support"]
    raw, binary = raw[:, keep], binary[:, keep]
    kept_names = np.asarray(names)[keep]
    restriction = restricted(df)
    result, selections = {}, []
    for horizon in (3, 5, 7):
        ep = endpoint(df, horizon)
        known = np.isfinite(ep)
        y = ep[known].astype(int)
        xb, xr = binary[known], raw[known]
        groups = df.loc[known, "patient_id"].astype(str).to_numpy()
        stored = pd.to_numeric(df.loc[known, "cv_fold"], errors="coerce").astype(int).to_numpy() - 1
        restricted_h = restriction[known]
        repetitions = 1 if config["fold_scheme"] == "stored" else config["repeats"]
        predictions = []
        fold_aucs = {"broad": [], "restricted": []}
        for rep in range(repetitions):
            folds = stored if config["fold_scheme"] == "stored" else cyclic_patient_folds(groups, config["seed"] + rep)
            pred = np.full(len(y), np.nan)
            for fold in range(5):
                train, test = folds != fold, folds == fold
                select_x = xb[train] if config["selector_input"] == "binary" else xr[train]
                if config["selector"] == "chi2":
                    scores = chi2(select_x, y[train])[0]
                else:
                    scores = f_classif(select_x, y[train])[0]
                order = np.argsort(-np.nan_to_num(scores, nan=-np.inf), kind="stable")
                chosen = order[:config["k"]]
                selections.extend({"horizon": horizon, "repetition": rep + 1, "fold": fold + 1,
                                   "feature": kept_names[j], "configuration": config["id"]} for j in chosen)
                if config["model_input"] == "binary":
                    x = xb[:, chosen].astype(float)
                elif config["model_input"] == "log1p":
                    x = np.log1p(xr[:, chosen])
                else:
                    x = xr[:, chosen].astype(float)
                scaler = StandardScaler().fit(x[train])
                model = LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=.25,
                                           C=config["C"], max_iter=10000,
                                           random_state=config["seed"] + rep, tol=1e-4)
                model.fit(scaler.transform(x[train]), y[train])
                pred[test] = model.predict_proba(scaler.transform(x[test]))[:, 1]
                fold_aucs["broad"].append(auc(y, pred, test))
                rm = test & restricted_h
                if len(np.unique(y[rm])) == 2:
                    fold_aucs["restricted"].append(auc(y, pred, rm))
            predictions.append(pred)
        p = np.vstack(predictions)
        for cohort, mask in (("broad", np.ones(len(y), bool)), ("restricted", restricted_h)):
            repeat_aucs = [auc(y, row, mask) for row in p]
            result[(horizon, cohort, "mean_repeat")] = float(np.mean(repeat_aucs))
            result[(horizon, cohort, "averaged_prediction")] = auc(y, p.mean(axis=0), mask)
            result[(horizon, cohort, "mean_fold")] = float(np.mean(fold_aucs[cohort]))
    return result, selections


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.data, low_memory=False)
    grid = itertools.product(
        ("any", "all"), ((20, "ge"), (20, "gt")), (5, 8),
        (("chi2", "binary"), ("f_classif", "count")),
        ("binary", "count", "log1p"), (.03, .1, .3, 1.0),
        ("stored", "cyclic"),
    )
    rows, all_selections = [], []
    config_number = 0
    cache = {}
    for family, (support, op), k, (selector, selector_input), model_input, C, fold_scheme in grid:
        if args.quick and fold_scheme != "stored":
            continue
        if family not in cache:
            names = candidate_columns(df, family)
            cache[family] = (df[names].apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy(), names)
        config_number += 1
        config = dict(id=f"c{config_number:04d}", family=family, support=support, support_op=op, k=k,
                      selector=selector, selector_input=selector_input, model_input=model_input, C=C,
                      fold_scheme=fold_scheme, repeats=5, seed=20260803)
        values, selections = evaluate(df, *cache[family], config)
        all_selections.extend(selections)
        for aggregation in ("mean_repeat", "averaged_prediction", "mean_fold"):
            observed = {(h, c): values[(h, c, aggregation)] for h, c in TARGET}
            error = np.asarray([observed[key] - TARGET[key] for key in TARGET])
            rows.append(config | {"aggregation": aggregation, "rmse": float(np.sqrt(np.mean(error**2))),
                                  "max_abs_error": float(np.max(np.abs(error))),
                                  **{f"auc_{h}y_{c}": v for (h, c), v in observed.items()}})
    results = pd.DataFrame(rows).sort_values(["rmse", "max_abs_error"])
    results.to_csv(args.output_dir / "configuration_results.csv", index=False)
    best_ids = set(results.head(20).id)
    pd.DataFrame([x for x in all_selections if x["configuration"] in best_ids]).to_csv(
        args.output_dir / "top_configuration_selections.csv", index=False)
    (args.output_dir / "target.json").write_text(json.dumps({f"{h}y_{c}": v for (h, c), v in TARGET.items()}, indent=2) + "\n")
    print(results.head(25).to_string(index=False))


if __name__ == "__main__":
    main()
