#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(dplyr)
  library(glmnet)
  library(jsonlite)
  library(readr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript run_validation.R HOTLIST_FEATURES.csv OUTCOMES.csv OUTPUT_DIR")
}
feature_path <- args[[1]]
outcome_path <- args[[2]]
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

seed0 <- 20260803L
repetitions <- 10L
outer_folds <- 5L
alpha <- 0.25

auc_safe <- function(y, score) {
  ok <- is.finite(y) & is.finite(score)
  y <- y[ok]; score <- score[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(score, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

fixed_endpoint <- function(time, event, horizon) {
  ifelse(event == 1 & time <= horizon, 1,
         ifelse(time >= horizon, 0, NA_real_))
}

group_folds <- function(groups, k, seed) {
  set.seed(seed)
  ids <- sample(unique(as.character(groups)))
  map <- setNames(rep(seq_len(k), length.out = length(ids)), ids)
  unname(map[as.character(groups)])
}

fit_predict <- function(train, test, columns, seed) {
  xtr <- as.matrix(train[, columns, drop = FALSE])
  xte <- as.matrix(test[, columns, drop = FALSE])
  storage.mode(xtr) <- storage.mode(xte) <- "double"
  med <- apply(xtr, 2, function(z) median(z[is.finite(z)], na.rm = TRUE))
  med[!is.finite(med)] <- 0
  for (j in seq_along(columns)) {
    xtr[!is.finite(xtr[, j]), j] <- med[[j]]
    xte[!is.finite(xte[, j]), j] <- med[[j]]
  }
  variable <- apply(xtr, 2, function(z) length(unique(z)) > 1L)
  xtr <- xtr[, variable, drop = FALSE]
  xte <- xte[, variable, drop = FALSE]
  if (!ncol(xtr)) return(rep(NA_real_, nrow(test)))
  k <- min(5L, length(unique(train$patient_id)))
  if (k < 3L) return(rep(NA_real_, nrow(test)))
  inner <- group_folds(train$patient_id, k, seed)
  fit <- try(cv.glmnet(
    xtr, train$endpoint, family = "binomial", alpha = alpha,
    foldid = inner, type.measure = "deviance", standardize = TRUE
  ), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, nrow(test)))
  drop(predict(fit, xte, s = "lambda.1se", type = "link"))
}

features <- fread(feature_path, data.table = FALSE, check.names = FALSE) %>% as_tibble()
outcomes <- fread(outcome_path, data.table = FALSE, check.names = FALSE) %>% as_tibble()
required <- c("sample_id", "recurrence_time", "recurrence_event")
if (!all(required %in% names(outcomes))) {
  stop("Outcomes must contain: ", paste(required, collapse = ", "))
}
if (!"patient_id" %in% names(outcomes)) outcomes$patient_id <- outcomes$sample_id
if (anyDuplicated(features$sample_id)) stop("hotlist_features.csv must have one row per sample_id")
if (anyDuplicated(outcomes$sample_id)) stop("outcomes.csv must have one row per sample_id")

hotlist_columns <- grep("^hotlist_.+_count$", names(features), value = TRUE)
if (length(hotlist_columns) != 9L) {
  stop("Expected exactly nine hotlist count columns; found ", length(hotlist_columns))
}
d_join <- inner_join(features, outcomes, by = "sample_id")
if (!nrow(d_join)) stop("Feature and outcome files have no matching sample_id values")
qc_excluded <- if ("qc_status" %in% names(d_join)) sum(startsWith(d_join$qc_status, "FAIL_")) else 0L
d0 <- if ("qc_status" %in% names(d_join)) {
  d_join %>% filter(!startsWith(qc_status, "FAIL_"))
} else d_join
if (!nrow(d0)) stop("All matched samples failed hotlist QC")

clinical_raw <- c("age", "mitoses", "gtr", "prior_rt", "adj_rt")
has_clinical <- all(clinical_raw %in% names(d0))
if (has_clinical) d0$log1p_mitoses <- log1p(as.numeric(d0$mitoses))
clinical_columns <- c("age", "log1p_mitoses", "gtr", "prior_rt", "adj_rt")

# Secondary sensitivity only: fit the nine-window stability union as a fixed
# feature set. This is not the dynamic top-eight procedure that generated the
# historical panel AUCs.
horizon_features <- list(
  `3` = hotlist_columns,
  `5` = hotlist_columns,
  `7` = hotlist_columns
)

prediction_rows <- list()
for (h in c(3, 5, 7)) {
  active_hotlist <- horizon_features[[as.character(h)]]
  d <- d0 %>% mutate(endpoint = fixed_endpoint(recurrence_time, recurrence_event, h)) %>%
    filter(!is.na(endpoint))
  if (length(unique(d$endpoint)) < 2L) next
  k <- min(outer_folds, length(unique(d$patient_id)))
  if (k < 3L) next
  for (r in seq_len(repetitions)) {
    outer <- group_folds(d$patient_id, k, seed0 + h * 1000L + r)
    panel_lp <- rep(NA_real_, nrow(d))
    clinical_panel_lp <- rep(NA_real_, nrow(d))
    for (fold in seq_len(k)) {
      tr <- d[outer != fold, , drop = FALSE]
      te <- d[outer == fold, , drop = FALSE]
      panel_lp[outer == fold] <- fit_predict(
        tr, te, active_hotlist, seed0 + h * 100000L + r * 100L + fold
      )
      if (has_clinical) {
        clinical_panel_lp[outer == fold] <- fit_predict(
          tr, te, c(clinical_columns, active_hotlist),
          seed0 + h * 200000L + r * 100L + fold
        )
      }
    }
    prediction_rows[[length(prediction_rows) + 1L]] <- tibble(
      horizon_years = h, repetition = r, sample_id = d$sample_id,
      patient_id = d$patient_id, endpoint = d$endpoint,
      fixed_hotlist = panel_lp,
      clinical_plus_hotlist = clinical_panel_lp,
      who_grade = if ("who_grade" %in% names(d)) as.numeric(d$who_grade) else NA_real_,
      integrated_grade = if ("integrated_grade" %in% names(d)) as.numeric(d$integrated_grade) else NA_real_
    )
  }
}
pred <- bind_rows(prediction_rows)
if (!nrow(pred)) stop("No horizon had enough eligible events and non-events")

model_columns <- c(
  "Nine-window stability-union sensitivity" = "fixed_hotlist",
  "Clinical plus nine-window sensitivity" = "clinical_plus_hotlist",
  "Recorded WHO grade" = "who_grade",
  "Prespecified integrated grade" = "integrated_grade"
)
by_rep <- bind_rows(lapply(names(model_columns), function(label) {
  column <- model_columns[[label]]
  pred %>% group_by(horizon_years, repetition) %>% summarise(
    model = label, n = n(), events = sum(endpoint),
    auc = auc_safe(endpoint, .data[[column]]), .groups = "drop"
  )
})) %>% filter(is.finite(auc))
summary <- by_rep %>% group_by(horizon_years, model) %>% summarise(
  n = first(n), events = first(events), repetitions = n(),
  mean_auc = mean(auc), sd_auc = sd(auc), min_auc = min(auc), max_auc = max(auc),
  .groups = "drop"
)

association_rows <- list()
for (h in c(3, 5, 7)) {
  d <- d0 %>% mutate(endpoint = fixed_endpoint(recurrence_time, recurrence_event, h)) %>%
    filter(!is.na(endpoint))
  for (column in hotlist_columns) {
    z <- as.numeric(d[[column]])
    zsd <- sd(z, na.rm = TRUE)
    fit <- if (is.finite(zsd) && zsd > 0) try(glm(endpoint ~ I((z - mean(z, na.rm = TRUE)) / zsd),
      data = d, family = binomial()), silent = TRUE) else NULL
    if (is.null(fit) || inherits(fit, "try-error")) {
      estimate <- low <- high <- p <- NA_real_
    } else {
      sm <- summary(fit)$coefficients[2, ]
      estimate <- exp(sm[["Estimate"]])
      low <- exp(sm[["Estimate"]] - 1.96 * sm[["Std. Error"]])
      high <- exp(sm[["Estimate"]] + 1.96 * sm[["Std. Error"]])
      p <- sm[["Pr(>|z|)"]]
    }
    association_rows[[length(association_rows) + 1L]] <- tibble(
      horizon_years = h, feature = column, n = nrow(d), events = sum(d$endpoint),
      prevalence = mean(z > 0, na.rm = TRUE), mean_boundary_count = mean(z, na.rm = TRUE),
      odds_ratio_per_sd = estimate, ci_low = low, ci_high = high, p_value = p
    )
  }
}

write_csv(summary, file.path(output_dir, "panel_auc_summary.csv"))
write_csv(by_rep, file.path(output_dir, "panel_auc_by_repetition.csv"))
write_csv(bind_rows(association_rows), file.path(output_dir, "individual_window_associations.csv"))
write_json(list(
  workflow = "hotlist-stability-v1 secondary fixed-union sensitivity",
  timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  feature_file_sha256 = digest(file = feature_path, algo = "sha256", serialize = FALSE),
  outcome_file_sha256 = digest(file = outcome_path, algo = "sha256", serialize = FALSE),
  feature_union = hotlist_columns, horizon_features = horizon_features,
  horizons_years = c(3, 5, 7),
  repetitions = repetitions, outer_folds = outer_folds, inner_folds = 5,
  alpha = alpha, lambda = "lambda.1se", seed = seed0,
  integrated_grade_rule = "recorded prespecified label only; no proxy reconstruction",
  matched_samples_before_qc = nrow(d_join), qc_failed_samples_excluded = qc_excluded,
  sample_level_predictions_written = FALSE
), file.path(output_dir, "validation_manifest.json"), auto_unbox = TRUE, pretty = TRUE)

print(summary)
