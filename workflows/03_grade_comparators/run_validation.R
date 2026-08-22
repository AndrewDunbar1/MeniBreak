#!/usr/bin/env Rscript

suppressPackageStartupMessages({ library(dplyr); library(readr) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript run_validation.R outcomes.csv output_directory")
data <- read_csv(args[[1]], show_col_types = FALSE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
required <- c("sample_id", "recurrence_time", "recurrence_event", "who_grade", "integrated_grade")
if (!all(required %in% names(data))) stop("Missing: ", paste(setdiff(required, names(data)), collapse = ", "))

endpoint_at <- function(time, event, horizon) ifelse(
  event == 1 & time <= horizon, 1, ifelse(time >= horizon, 0, NA_real_)
)
auc_safe <- function(y, score) {
  ok <- is.finite(y) & is.finite(score); y <- y[ok]; score <- score[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(score, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
cohorts <- list(Broad = rep(TRUE, nrow(data)))
if ("restricted_eligible" %in% names(data)) {
  cohorts[["Restricted"]] <- as.logical(data$restricted_eligible)
}
results <- bind_rows(lapply(names(cohorts), function(cohort) {
  bind_rows(lapply(c(3L, 5L, 7L), function(horizon) {
    d <- data %>% mutate(endpoint = endpoint_at(recurrence_time, recurrence_event, horizon)) %>%
      filter(.env$cohorts[[cohort]], !is.na(endpoint), who_grade %in% 1:3, integrated_grade %in% 1:3)
    bind_rows(
      tibble(cohort = cohort, horizon_years = horizon, model = "Recorded WHO grade",
             n = nrow(d), events = sum(d$endpoint), auc = auc_safe(d$endpoint, d$who_grade)),
      tibble(cohort = cohort, horizon_years = horizon, model = "Prespecified integrated grade",
             n = nrow(d), events = sum(d$endpoint), auc = auc_safe(d$endpoint, d$integrated_grade))
    )
  }))
}))
write_csv(results, file.path(output_dir, "recorded_grade_auc_summary.csv"))
print(results)
