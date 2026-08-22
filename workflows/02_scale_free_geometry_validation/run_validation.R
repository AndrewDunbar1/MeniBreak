#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(glmnet)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript run_validation.R geometry_features.csv outcomes.csv output_directory")
}
feature_path <- args[[1]]
outcome_path <- args[[2]]
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

features <- read_csv(feature_path, show_col_types = FALSE)
outcomes <- read_csv(outcome_path, show_col_types = FALSE)
required_outcomes <- c("sample_id", "recurrence_time", "recurrence_event")
if (!all(required_outcomes %in% names(outcomes))) stop("Missing required outcome columns.")
if (!"patient_id" %in% names(outcomes)) outcomes$patient_id <- outcomes$sample_id
if ("qc_status" %in% names(features)) features <- filter(features, qc_status == "PASS")
geometry_names <- setdiff(names(features), c("sample_id", "qc_status"))
if (length(geometry_names) != 56L) stop("Expected exactly 56 scale-free geometry features.")
data <- inner_join(outcomes, features, by = "sample_id")

endpoint_at <- function(time, event, horizon) {
  ifelse(event == 1 & time <= horizon, 1,
         ifelse(time >= horizon, 0, NA_real_))
}
auc_safe <- function(y, score) {
  ok <- is.finite(y) & is.finite(score); y <- y[ok]; score <- score[ok]
  positives <- sum(y == 1); negatives <- sum(y == 0)
  if (!positives || !negatives) return(NA_real_)
  (sum(rank(score, ties.method = "average")[y == 1]) - positives * (positives + 1) / 2) /
    (positives * negatives)
}
group_folds <- function(groups, k, seed) {
  set.seed(seed)
  ids <- sample(unique(as.character(groups)))
  map <- setNames(rep(seq_len(k), length.out = length(ids)), ids)
  unname(map[as.character(groups)])
}
fit_predict <- function(training, testing, seed) {
  x_training <- as.matrix(training[, geometry_names, drop = FALSE])
  x_testing <- as.matrix(testing[, geometry_names, drop = FALSE])
  storage.mode(x_training) <- storage.mode(x_testing) <- "double"
  usable <- apply(x_training, 2, function(x) length(unique(x[is.finite(x)])) > 1L)
  x_training <- x_training[, usable, drop = FALSE]
  x_testing <- x_testing[, usable, drop = FALSE]
  if (!ncol(x_training)) return(rep(mean(training$endpoint), nrow(testing)))
  medians <- apply(x_training, 2, function(x) median(x[is.finite(x)], na.rm = TRUE))
  medians[!is.finite(medians)] <- 0
  for (j in seq_len(ncol(x_training))) {
    x_training[!is.finite(x_training[, j]), j] <- medians[[j]]
    x_testing[!is.finite(x_testing[, j]), j] <- medians[[j]]
  }
  inner <- group_folds(training$patient_id,
                       min(5L, length(unique(training$patient_id))), seed)
  fit <- cv.glmnet(x_training, training$endpoint, family = "binomial", alpha = 0.25,
                   foldid = inner, type.measure = "deviance", standardize = TRUE)
  drop(predict(fit, x_testing, s = "lambda.1se", type = "link"))
}

predictions <- list()
for (horizon in c(3L, 5L, 7L)) {
  eligible <- data %>%
    mutate(endpoint = endpoint_at(recurrence_time, recurrence_event, horizon)) %>%
    filter(!is.na(endpoint))
  if (n_distinct(eligible$endpoint) != 2L) stop("Endpoint lacks both classes at horizon ", horizon)
  for (repetition in seq_len(10L)) {
    outer <- group_folds(eligible$patient_id, 5L, 20260803L + horizon * 10000L + repetition)
    score <- rep(NA_real_, nrow(eligible))
    for (fold in seq_len(5L)) {
      score[outer == fold] <- fit_predict(
        eligible[outer != fold, , drop = FALSE], eligible[outer == fold, , drop = FALSE],
        20260806L + horizon * 100000L + repetition * 1000L + fold
      )
    }
    predictions[[length(predictions) + 1L]] <- tibble(
      sample_id = eligible$sample_id, patient_id = eligible$patient_id,
      horizon_years = horizon, repetition = repetition,
      endpoint = eligible$endpoint, geometry_score = score
    )
  }
}
prediction_table <- bind_rows(predictions)
performance_by_repetition <- prediction_table %>%
  group_by(horizon_years, repetition) %>%
  summarise(n = n(), events = sum(endpoint), auc = auc_safe(endpoint, geometry_score),
            .groups = "drop")
summary <- performance_by_repetition %>%
  group_by(horizon_years) %>%
  summarise(n = first(n), events = first(events), mean_auc = mean(auc),
            sd_across_repetitions = sd(auc), .groups = "drop")
write_csv(prediction_table, file.path(output_dir, "geometry_oof_predictions.csv"))
write_csv(performance_by_repetition, file.path(output_dir, "geometry_auc_by_repetition.csv"))
write_csv(summary, file.path(output_dir, "geometry_auc_summary.csv"))
writeLines(c(
  "model_id: scale-free-geometry-v1",
  "features: exact 56-feature scale-free geometry definition",
  "validation: 10 repeats of patient-grouped nested 5-fold validation",
  "estimator: binomial cv.glmnet; alpha=0.25; lambda.1se; inner patient grouping",
  "aggregation: mean of repetition-level AUCs"
), file.path(output_dir, "run_specification.txt"))
print(summary)
