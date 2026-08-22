#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(openxlsx)
  library(glmnet)
  library(ggplot2)
})

root <- Sys.getenv("MENIBREAK_PROJECT_ROOT")
if (!nzchar(root)) stop("Set MENIBREAK_PROJECT_ROOT to the analysis project root; see config/path_config.example.sh")
adir <- file.path(root, "analysis_folders", "Confirmatory_Comparator_Harmonization_2026")
outcomes_path <- file.path(root, "software", "mbpr", "analysis", "bakeoff_2x2", "outcomes.csv")
features_path <- file.path(adir, "legacy30_features.csv")
clinical_path <- file.path(root, "analysis_folders", "v6_processed_rerun_July6",
                           "clinical_arm_calls_merged.xlsx")
dir.create(adir, recursive = TRUE, showWarnings = FALSE)

master_seed <- 20260729L
repeats <- 10L
outer_folds <- 5L
inner_folds <- 5L
bootstrap_reps <- 2000L
alpha <- 0.25

auc_safe <- function(y, score) {
  ok <- is.finite(y) & is.finite(score)
  y <- y[ok]
  score <- score[ok]
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(score, ties.method = "average")[y == 1]) -
     n1 * (n1 + 1) / 2) / (n1 * n0)
}

average_precision <- function(y, score) {
  ok <- is.finite(y) & is.finite(score)
  y <- y[ok]
  score <- score[ok]
  if (!sum(y == 1) || !sum(y == 0)) return(NA_real_)
  ord <- order(score, decreasing = TRUE)
  yy <- y[ord]
  mean(cumsum(yy)[yy == 1] / which(yy == 1))
}

fixed_horizon_endpoint <- function(time, event, horizon) {
  ifelse(event == 1 & time <= horizon, 1,
         ifelse(time >= horizon, 0, NA_real_))
}

group_folds <- function(groups, k, seed) {
  set.seed(seed)
  ids <- sample(unique(as.character(groups)))
  map <- setNames(rep(seq_len(k), length.out = length(ids)), ids)
  unname(map[as.character(groups)])
}

scna_feature_names <- function(features) {
  grep("^cnv_|^scna_(loss_count|gain_count|altered_arm_count)$",
       names(features), value = TRUE)
}

full_feature_names <- function(features) {
  c(
    scna_feature_names(features),
    intersect(
      c("total_breakpoint_boundaries", "breakpoint_arms_touched",
        "breakpoint_route_entropy", "max_arm_breakpoint_fraction",
        "telomeric_breakpoint_fraction", "central_breakpoint_fraction",
        "centromeric_breakpoint_fraction", "mean_rel_tel_to_cent",
        "sd_rel_tel_to_cent", "cna_segment_count", "loss_segment_count",
        "gain_segment_count", "median_segment_length",
        "small_segment_count_lt5mb", "dynamic_hotspot_breakpoints_all",
        "dynamic_hotspot_breakpoints_non1p", "dynamic_hotspot_arms_used_all",
        "dynamic_hotspot_arms_used_non1p"),
      names(features)
    )
  )
}

grade_matrix <- function(x, prefix) {
  x <- as.integer(x)
  ans <- cbind(as.numeric(x == 2L), as.numeric(x == 3L))
  colnames(ans) <- paste0(prefix, c("_2", "_3"))
  ans
}

numeric_matrix <- function(data, vars) {
  x <- as.matrix(data[, vars, drop = FALSE])
  storage.mode(x) <- "double"
  x
}

impute_train_test <- function(xtr, xte) {
  med <- apply(xtr, 2, function(z) median(z[is.finite(z)], na.rm = TRUE))
  med[!is.finite(med)] <- 0
  for (j in seq_len(ncol(xtr))) {
    xtr[!is.finite(xtr[, j]), j] <- med[[j]]
    xte[!is.finite(xte[, j]), j] <- med[[j]]
  }
  list(train = xtr, test = xte)
}

fit_predict <- function(xtr, xte, y, groups, foldid) {
  imp <- impute_train_test(xtr, xte)
  fit <- cv.glmnet(
    imp$train, y, family = "binomial", alpha = alpha,
    foldid = foldid, type.measure = "deviance", standardize = TRUE
  )
  raw <- drop(predict(fit, imp$test, s = "lambda.1se", type = "link"))
  prob <- drop(predict(fit, imp$test, s = "lambda.1se", type = "response"))
  list(raw = raw, probability = prob, lambda = fit$lambda.1se)
}

outcomes <- read_csv(outcomes_path, show_col_types = FALSE)
features <- read_csv(features_path, show_col_types = FALSE)
clinical <- read.xlsx(clinical_path) %>%
  transmute(
    sample_id = trimws(as.character(cg)),
    who_grade = as.integer(who_grade),
    integrated_grade_30 = as.integer(IG30)
  ) %>%
  distinct(sample_id, .keep_all = TRUE)

stopifnot(!anyDuplicated(outcomes$sample_id), !anyDuplicated(features$sample_id),
          !anyDuplicated(clinical$sample_id))

base <- outcomes %>%
  inner_join(features, by = "sample_id") %>%
  inner_join(clinical, by = "sample_id") %>%
  filter(who_grade %in% 1:3, integrated_grade_30 %in% 1:3)

scna_vars <- scna_feature_names(features)
full_vars <- full_feature_names(features)
stopifnot(length(scna_vars) == 44L, length(full_vars) == 62L)

cohort_rows <- list()
prediction_rows <- list()
metric_rows <- list()
fold_rows <- list()
lambda_rows <- list()

for (h in c(3L, 5L, 7L)) {
  d <- base %>%
    mutate(endpoint = fixed_horizon_endpoint(recurrence_time, recurrence_event, h)) %>%
    filter(!is.na(endpoint)) %>%
    arrange(patient_id, sample_id)

  cohort_rows[[as.character(h)]] <- tibble(
    horizon_years = h,
    tumors = nrow(d),
    patients = n_distinct(d$patient_id),
    events = sum(d$endpoint),
    nonevents = sum(d$endpoint == 0),
    who_grade_1 = sum(d$who_grade == 1),
    who_grade_2 = sum(d$who_grade == 2),
    who_grade_3 = sum(d$who_grade == 3),
    integrated_grade_1 = sum(d$integrated_grade_30 == 1),
    integrated_grade_2 = sum(d$integrated_grade_30 == 2),
    integrated_grade_3 = sum(d$integrated_grade_30 == 3)
  )

  for (r in seq_len(repeats)) {
    outer_seed <- master_seed + h * 10000L + r
    outer <- group_folds(d$patient_id, outer_folds, outer_seed)
    pred <- list(
      WHO = list(raw = rep(NA_real_, nrow(d)), probability = rep(NA_real_, nrow(d))),
      `Integrated grade` = list(raw = rep(NA_real_, nrow(d)), probability = rep(NA_real_, nrow(d))),
      `SCNA alone` = list(raw = rep(NA_real_, nrow(d)), probability = rep(NA_real_, nrow(d))),
      `SCNA + breakpoint` = list(raw = rep(NA_real_, nrow(d)), probability = rep(NA_real_, nrow(d)))
    )

    fold_rows[[length(fold_rows) + 1L]] <- tibble(
      sample_id = d$sample_id, patient_id = d$patient_id,
      horizon_years = h, repetition = r, outer_fold = outer,
      outer_seed = outer_seed
    )

    for (fold in seq_len(outer_folds)) {
      train_idx <- which(outer != fold)
      test_idx <- which(outer == fold)
      inner_seed <- master_seed + h * 100000L + r * 1000L + fold
      inner <- group_folds(d$patient_id[train_idx], inner_folds, inner_seed)

      matrices <- list(
        WHO = list(train = grade_matrix(d$who_grade[train_idx], "who"),
                   test = grade_matrix(d$who_grade[test_idx], "who")),
        `Integrated grade` = list(
          train = grade_matrix(d$integrated_grade_30[train_idx], "ig"),
          test = grade_matrix(d$integrated_grade_30[test_idx], "ig")),
        `SCNA alone` = list(train = numeric_matrix(d[train_idx, ], scna_vars),
                            test = numeric_matrix(d[test_idx, ], scna_vars)),
        `SCNA + breakpoint` = list(train = numeric_matrix(d[train_idx, ], full_vars),
                                   test = numeric_matrix(d[test_idx, ], full_vars))
      )

      for (model_name in names(matrices)) {
        fit <- fit_predict(
          matrices[[model_name]]$train, matrices[[model_name]]$test,
          d$endpoint[train_idx], d$patient_id[train_idx], inner
        )
        pred[[model_name]]$raw[test_idx] <- fit$raw
        pred[[model_name]]$probability[test_idx] <- fit$probability
        lambda_rows[[length(lambda_rows) + 1L]] <- tibble(
          horizon_years = h, repetition = r, outer_fold = fold,
          model = model_name, inner_seed = inner_seed, lambda_1se = fit$lambda
        )
      }
    }

    for (model_name in names(pred)) {
      metric_rows[[length(metric_rows) + 1L]] <- tibble(
        horizon_years = h, repetition = r, model = model_name,
        n = nrow(d), events = sum(d$endpoint),
        auc = auc_safe(d$endpoint, pred[[model_name]]$raw),
        average_precision = average_precision(d$endpoint, pred[[model_name]]$raw),
        brier = mean((d$endpoint - pred[[model_name]]$probability)^2)
      )
      prediction_rows[[length(prediction_rows) + 1L]] <- tibble(
        sample_id = d$sample_id, patient_id = d$patient_id,
        horizon_years = h, repetition = r, endpoint = d$endpoint,
        who_grade = d$who_grade,
        integrated_grade_30 = d$integrated_grade_30,
        model = model_name, raw_score = pred[[model_name]]$raw,
        probability = pred[[model_name]]$probability
      )
    }
  }
}

cohorts <- bind_rows(cohort_rows)
predictions <- bind_rows(prediction_rows)
repeat_metrics <- bind_rows(metric_rows)
fold_assignments <- bind_rows(fold_rows)
lambdas <- bind_rows(lambda_rows)

averaged <- predictions %>%
  group_by(sample_id, patient_id, horizon_years, endpoint, who_grade,
           integrated_grade_30, model) %>%
  summarise(raw_score = mean(raw_score), probability = mean(probability),
            .groups = "drop")

bootstrap_rows <- list()
for (h in c(3L, 5L, 7L)) {
  dh <- predictions %>% filter(horizon_years == h)
  patients <- sort(unique(dh$patient_id))
  set.seed(master_seed + h * 1000000L)
  for (b in seq_len(bootstrap_reps)) {
    sampled <- sample(patients, length(patients), replace = TRUE)
    weights <- table(sampled)
    ids <- rep(names(weights), as.integer(weights))
    pieces <- lapply(seq_along(ids), function(i) {
      x <- dh[as.character(dh$patient_id) == ids[[i]], , drop = FALSE]
      x$cluster_copy <- i
      x
    })
    db <- bind_rows(pieces)
    aucs <- db %>%
      group_by(model, repetition) %>%
      summarise(repetition_auc = auc_safe(endpoint, raw_score),
                .groups = "drop_last") %>%
      summarise(auc = mean(repetition_auc, na.rm = TRUE), .groups = "drop")
    for (model in aucs$model) {
      bootstrap_rows[[length(bootstrap_rows) + 1L]] <- tibble(
        horizon_years = h, bootstrap_draw = b,
        model = model, auc = aucs$auc[aucs$model == model]
      )
    }
  }
}
bootstrap_auc <- bind_rows(bootstrap_rows)

model_summary <- repeat_metrics %>%
  group_by(horizon_years, model) %>%
  summarise(
    tumors = first(n), events = first(events),
    mean_auc = mean(auc), sd_auc_across_repeats = sd(auc),
    mean_average_precision = mean(average_precision),
    mean_brier = mean(brier), .groups = "drop"
  ) %>%
  left_join(
    bootstrap_auc %>%
      group_by(horizon_years, model) %>%
      summarise(bootstrap_auc = median(auc),
                auc_ci_low = quantile(auc, 0.025, na.rm = TRUE),
                auc_ci_high = quantile(auc, 0.975, na.rm = TRUE),
                .groups = "drop"),
    by = c("horizon_years", "model")
  )

paired <- bootstrap_auc %>%
  select(horizon_years, bootstrap_draw, model, auc) %>%
  pivot_wider(names_from = model, values_from = auc)

comparators <- c("WHO", "Integrated grade", "SCNA alone")
delta_summary <- bind_rows(lapply(comparators, function(comp) {
  paired %>%
    transmute(horizon_years, bootstrap_draw, comparator = comp,
              delta_auc = .data[["SCNA + breakpoint"]] - .data[[comp]]) %>%
    group_by(horizon_years, comparator) %>%
    summarise(
      ci_low = quantile(delta_auc, 0.025, na.rm = TRUE),
      ci_high = quantile(delta_auc, 0.975, na.rm = TRUE),
      bootstrap_two_sided_p = 2 * min(mean(delta_auc <= 0, na.rm = TRUE),
                                      mean(delta_auc >= 0, na.rm = TRUE)),
      delta_auc = median(delta_auc, na.rm = TRUE),
      .groups = "drop"
    )
}))

write_csv(cohorts, file.path(adir, "common_cohort_by_horizon.csv"))
write_csv(fold_assignments, file.path(adir, "shared_outer_fold_assignments.csv"))
write_csv(lambdas, file.path(adir, "nested_cv_lambda_audit.csv"))
write_csv(predictions, file.path(adir, "all_oof_predictions.csv"))
write_csv(averaged, file.path(adir, "averaged_oof_predictions.csv"))
write_csv(repeat_metrics, file.path(adir, "repeat_metrics.csv"))
write_csv(bootstrap_auc, file.path(adir, "shared_patient_bootstrap_auc.csv"))
write_csv(model_summary, file.path(adir, "confirmatory_model_summary.csv"))
write_csv(delta_summary, file.path(adir, "confirmatory_paired_auc_deltas.csv"))

palette <- c(
  "WHO" = "#7F7F7F",
  "Integrated grade" = "#4C78A8",
  "SCNA alone" = "#E69F00",
  "SCNA + breakpoint" = "#B2182B"
)
model_summary$model <- factor(
  model_summary$model,
  levels = c("WHO", "Integrated grade", "SCNA alone", "SCNA + breakpoint")
)

p <- ggplot(model_summary,
            aes(x = factor(horizon_years), y = mean_auc,
                color = model, group = model)) +
  geom_hline(yintercept = 0.5, color = "#D9D9D9", linewidth = 0.5) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.6) +
  geom_errorbar(aes(ymin = auc_ci_low, ymax = auc_ci_high),
                width = 0.08, linewidth = 0.65) +
  scale_color_manual(values = palette) +
  scale_y_continuous(limits = c(0.45, 0.90), breaks = seq(0.5, 0.9, 0.1)) +
  labs(x = "Recurrence horizon (years)", y = "Cross-validated AUC",
       color = NULL,
       title = "Harmonized internal recurrence discrimination",
       subtitle = "Identical tumors, patient-grouped nested resampling, and paired patient bootstrap") +
  theme_classic(base_family = "sans", base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "#444444"),
    legend.position = "bottom",
    legend.text = element_text(size = 9)
  )

ggsave(file.path(adir, "confirmatory_harmonized_auc.png"), p,
       width = 7.6, height = 5.2, dpi = 400, bg = "white")
ggsave(file.path(adir, "confirmatory_harmonized_auc.pdf"), p,
       width = 7.6, height = 5.2, device = "pdf")

manifest <- tibble(
  item = c("analysis_timestamp", "master_seed", "outer_repeats",
           "outer_folds", "inner_folds", "bootstrap_reps", "alpha",
           "outcomes_path", "features_path", "clinical_path",
           "outcomes_sha256", "features_sha256", "clinical_sha256"),
  value = c(
    format(Sys.time(), tz = "UTC", usetz = TRUE),
    master_seed, repeats, outer_folds, inner_folds, bootstrap_reps, alpha,
    outcomes_path, features_path, clinical_path,
    digest::digest(file = outcomes_path, algo = "sha256"),
    digest::digest(file = features_path, algo = "sha256"),
    digest::digest(file = clinical_path, algo = "sha256")
  )
)
write_csv(manifest, file.path(adir, "run_manifest.csv"))

cat("Confirmatory analysis complete.\n")
print(cohorts)
print(model_summary)
print(delta_summary)
