\name{mbpr-api}
\alias{inspect_segments}
\alias{normalize_segments}
\alias{extract_features}
\alias{freeze_selected_model}
\alias{score_segments}
\alias{validate_external}
\alias{model_info}
\alias{run_bakeoff}
\alias{audit_selected_model}
\title{Core mbpr research workflow}
\usage{
inspect_segments(segments, format = "auto", genome_build,
  mapping_file = NULL, caller_profile = NULL)
normalize_segments(segments, format = "auto", genome_build,
  mapping_file = NULL, caller_profile = NULL, allow_incompatible = FALSE)
extract_features(segments, genome_build = NULL,
  preprocess = c("legacy", "v6"), arm_threshold = 50, dgv_path = NULL)
score_segments(segments, format = "auto", genome_build, output_dir = NULL,
  mapping_file = NULL, caller_profile = NULL, features_format = "csv")
validate_external(predictions, outcomes, output_dir = NULL,
  bootstrap_reps = 1000, seed = 20260703)
model_info(as_json = FALSE)
run_bakeoff(segments, outcomes, genome_build = "GRCh37",
  v6_segments = NULL, dgv_path = NULL,
  output_dir, repeats = 10, folds = 5, seed = 20260703,
  legacy_reproduction_verified = FALSE, lock_model = FALSE)
audit_selected_model(segments, outcomes, output_dir,
  genome_build = "GRCh37", repeats = 10, folds = 5, seed = 20260703)
freeze_selected_model(segments, outcomes, oof_predictions,
  legacy_audit_json, selection_json, release_audit_json, output_path,
  genome_build = "GRCh37", seed = 20260703)
}
\arguments{
  \item{segments}{A segment data frame or path to a delimited segment file.}
  \item{format}{A built-in adapter name, \code{custom_mapping}, or \code{auto}.}
  \item{genome_build}{Explicitly \code{GRCh37} or \code{GRCh38}.}
  \item{mapping_file}{Optional custom column-mapping YAML path.}
  \item{caller_profile}{Explicit caller-profile YAML for continuous segments.}
  \item{allow_incompatible}{Internal legacy-reproduction override. External scoring must leave this false.}
  \item{preprocess}{Legacy or V6 segment-cleaning lineage.}
  \item{arm_threshold}{Percent arm coverage required for an arm call: 30 or 50.}
  \item{dgv_path}{Optional Zarrei DGV reference path used by V6 preprocessing.}
  \item{v6_segments}{Optional protocol-produced V6 segment file for the derivation bakeoff.}
  \item{output_dir}{Directory for deterministic output artifacts.}
  \item{features_format}{\code{csv} or optional \code{parquet}.}
  \item{predictions}{Prediction data frame or file path.}
  \item{outcomes}{Outcome data frame or file path.}
  \item{oof_predictions}{Repeated out-of-fold predictions from the release audit.}
  \item{legacy_audit_json}{Passed legacy-reproduction audit JSON path.}
  \item{selection_json}{Passed prespecified candidate-selection JSON path.}
  \item{release_audit_json}{Completed selected-model release-audit JSON path.}
  \item{output_path}{Destination path for the frozen RDS artifact.}
  \item{bootstrap_reps}{Number of patient-clustered bootstrap replicates.}
  \item{seed}{Fixed random seed.}
  \item{as_json}{Return model information as JSON text.}
  \item{repeats}{Number of repeated outer validations.}
  \item{folds}{Number of patient-grouped outer folds.}
  \item{legacy_reproduction_verified}{Logical legacy-control gate.}
  \item{lock_model}{Whether to write an artifact after all gates pass.}
}
\value{
An inspection list, canonical segment table, feature table, prediction table,
validation result, model metadata, or bakeoff result depending on the function.
}
\description{
Functions implementing deterministic inspection, normalization, feature
extraction, fail-closed scoring, outcome-linked validation, and the prespecified
derivation bakeoff.
}
