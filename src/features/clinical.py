"""Clinical covariates for the clinical-integrated comparator.

The locked project has TWO genomic benchmarks and they must not be confused:

* **SEG-only locked model** — genomic features alone.
  AUC 0.7792 / 0.7880 / 0.7798 at 3 / 5 / 7 years.
* **Clinical-integrated extension** — the same genomic features PLUS five
  clinical covariates. AUC 0.8464 / 0.8442 / 0.8418.

The clinical-integrated figures are the higher bar and the realistic starting
point for any attempt at 0.90, so the foundation model must ultimately be
compared against them, not only against the SEG-only model.

The five covariates, exactly as specified in the handoff methods:

    age                       continuous
    log(1 + mitotic count)
    gross-total resection     indicator
    prior radiotherapy        indicator
    adjuvant/post-op RT       indicator

**Sex is NOT among them**, nor are WHO grade or integrated grade — those were
comparator scores, never predictors inside the genomic models. Adding any of
them would silently change what is being measured.

Source: `03_ANALYSIS_READY_DATA/ml_breakpoint_model_dataset_2026.csv`, the
deidentified analysis-ready table. The restricted linkage workbook is never
touched.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

#: Column in the analysis-ready table -> canonical covariate name.
CLINICAL_SOURCE_COLUMNS = {
    "age": "age",
    "mitoses": "mitoses",
    "gtr_num": "gtr",
    "prior_rt_num": "prior_rt",
    "adj_rt_num": "adj_rt",
}

#: The five covariates the clinical-integrated model actually uses.
CLINICAL_FEATURES = ["age", "log1p_mitoses", "gtr", "prior_rt", "adj_rt"]

#: Deliberately excluded. Kept as an explicit list so the omission is visible
#: and testable rather than an accident of which columns happened to be read.
DELIBERATELY_EXCLUDED = {
    "sex": "not a predictor in the locked clinical-integrated model",
    "who_grade": "comparator score, never a predictor inside the genomic models",
    "integrated_grade": "comparator score, never a predictor",
    "ki67": "not part of the specified covariate set",
    "brain_invasion": "not part of the specified covariate set",
    "nf2_loss": "assay-derived; not in the portable scoring core",
    "cdkn2ab_loss": "not part of the specified covariate set",
}


def load_clinical_covariates(
    path: str | Path, *, key: str = "cg"
) -> pd.DataFrame:
    """Load the five clinical covariates, indexed by deidentified tumour id.

    Only the required columns are read. Nothing else from the analysis-ready
    table enters the model, so an outcome or grade column cannot leak in by
    being adjacent in the file.
    """
    usecols = [key] + list(CLINICAL_SOURCE_COLUMNS)
    df = pd.read_csv(path, usecols=lambda c: c in usecols, low_memory=False)

    missing = [c for c in usecols if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: clinical source missing {missing}")

    out = pd.DataFrame({"tumor_id": df[key].astype("string")})
    out["age"] = pd.to_numeric(df["age"], errors="coerce")
    # log1p compresses the long right tail of mitotic counts, as specified.
    out["log1p_mitoses"] = np.log1p(pd.to_numeric(df["mitoses"], errors="coerce"))
    out["gtr"] = pd.to_numeric(df["gtr_num"], errors="coerce")
    out["prior_rt"] = pd.to_numeric(df["prior_rt_num"], errors="coerce")
    out["adj_rt"] = pd.to_numeric(df["adj_rt_num"], errors="coerce")

    out = out.dropna(subset=["tumor_id"]).drop_duplicates("tumor_id")
    log.info(
        "clinical covariates for %d tumours; missing rates: %s",
        len(out),
        {c: round(float(out[c].isna().mean()), 3) for c in CLINICAL_FEATURES},
    )
    return out.set_index("tumor_id")


#: The locked model's restricted clinical subset, reproduced verbatim from the
#: handoff methods (section 4.3). All six conditions must hold.
RESTRICTED_CRITERIA = {
    "primary": 1.0,        # primary presentation
    "gtr": 1.0,            # gross-total resection
    "prior_rt": 0.0,       # no prior radiotherapy
    "adj_rt": 0.0,         # no adjuvant / post-operative radiotherapy
    "rt_induced": 0.0,     # not radiation-induced
    "angiomatous": 0.0,    # non-angiomatous histology
}


def load_restricted_flag(path: str | Path, *, key: str = "cg") -> pd.DataFrame:
    """Per-tumour boolean for the locked model's restricted subset.

    Rows with a MISSING indicator are excluded rather than assumed favourable.
    Treating an unknown radiotherapy history as "no radiotherapy" would quietly
    enlarge the subset with exactly the patients most likely to violate it.

    The locked model applies broad-cohort models to this subset **without
    subset-specific refitting**, so it tests transport into a cleaner clinical
    setting rather than a separately trained restricted model.
    """
    usecols = [key] + list(RESTRICTED_CRITERIA)
    df = pd.read_csv(path, usecols=lambda c: c in usecols, low_memory=False)

    missing = [c for c in usecols if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: restricted criteria missing {missing}")

    ok = pd.Series(True, index=df.index)
    for col, want in RESTRICTED_CRITERIA.items():
        v = pd.to_numeric(df[col], errors="coerce")
        ok &= v.eq(want)          # NaN != want, so missing values drop out

    out = pd.DataFrame({
        "tumor_id": df[key].astype("string"),
        "is_restricted": ok.to_numpy(),
    }).dropna(subset=["tumor_id"]).drop_duplicates("tumor_id")

    log.info("restricted subset: %d of %d tumours meet all %d criteria",
             int(out["is_restricted"].sum()), len(out), len(RESTRICTED_CRITERIA))
    return out.set_index("tumor_id")


def assert_no_excluded_covariates(columns) -> None:
    """Fail loudly if a forbidden covariate reached the feature matrix.

    Grade in particular is a comparator; letting it in as a predictor would
    inflate the genomic models and destroy the comparison it is meant to serve.
    """
    lowered = {str(c).lower() for c in columns}
    for banned in DELIBERATELY_EXCLUDED:
        if banned in lowered:
            raise AssertionError(
                f"{banned!r} must not be a predictor: "
                f"{DELIBERATELY_EXCLUDED[banned]}"
            )
