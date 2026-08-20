#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(readr); library(dplyr); library(tidyr); library(glmnet)})
root <- Sys.getenv("MENIBREAK_PROJECT_ROOT")
if (!nzchar(root)) stop("Set MENIBREAK_PROJECT_ROOT to the analysis project root; see config/path_config.example.sh")
out_dir <- file.path(root, "analysis_folders", "Breakpoint_Geometry_Route_2026")
raw_path <- file.path(root, "analysis_folders", "ml_breakpoint_survival_stratification_2026", "ml_breakpoint_model_dataset_2026.csv")
outcome_path <- file.path(root, "software", "mbpr", "analysis", "bakeoff_2x2", "outcomes.csv")
geometry_path <- file.path(out_dir, "geometry_route_features.csv")
stopifnot(all(file.exists(c(raw_path, outcome_path, geometry_path))))

repeats <- 10L; seed0 <- 20260806L
auc_safe <- function(y, s) {
  ok <- is.finite(y) & is.finite(s); y <- y[ok]; s <- s[ok]; n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(s, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
fixed_endpoint <- function(time, event, h) ifelse(event == 1 & time <= h, 1, ifelse(time >= h, 0, NA_real_))
group_folds <- function(groups, k, seed) {
  set.seed(seed); ids <- sample(unique(as.character(groups))); map <- setNames(rep(seq_len(k), length.out = length(ids)), ids)
  unname(map[as.character(groups)])
}
fit_predict <- function(xtr, ytr, groups, xte, seed) {
  usable <- apply(xtr, 2, function(z) length(unique(z[is.finite(z)])) > 1)
  xtr <- xtr[, usable, drop = FALSE]; xte <- xte[, usable, drop = FALSE]
  foldid <- group_folds(groups, min(5L, length(unique(groups))), seed)
  f <- cv.glmnet(xtr, ytr, family = "binomial", alpha = .25, foldid = foldid,
                 type.measure = "deviance", standardize = TRUE)
  drop(predict(f, xte, s = "lambda.1se", type = "link"))
}

raw <- read_csv(raw_path, show_col_types = FALSE)
outcomes <- read_csv(outcome_path, show_col_types = FALSE)
geo <- read_csv(geometry_path, show_col_types = FALSE)
# Recover the exact fixed scale-free set by definition. The saved table stores
# topology columns interleaved by event class, so positional slicing is unsafe.
topology_suffix <- "__(arms_touched|chromosomes_touched|both_arms_chromosomes|occupied_bins|multi_break_arms|route_clusters)$"
complexity_names <- c("log_total_breakpoints", "log_cna_segments", "log_loss_segments", "log_gain_segments",
                      "log_focal_segments", "log_median_segment_length", "focal_segment_fraction", "loss_segment_fraction")
shape_names <- setdiff(names(geo), c("sample_id", grep(topology_suffix, names(geo), value = TRUE), complexity_names))
stopifnot(length(shape_names) == 56L)
Xall <- as.matrix(geo[, shape_names]); storage.mode(Xall) <- "double"; Xall[!is.finite(Xall)] <- 0
base <- raw %>% transmute(sample_id = cg, patient_id, who_grade = as.numeric(who_grade), integrated_grade = as.numeric(integrated_grade)) %>%
  inner_join(outcomes, by = c("sample_id", "patient_id"))

all_pred <- list()
for (h in c(3, 5, 7)) {
  dat <- base %>% mutate(endpoint = fixed_endpoint(recurrence_time, recurrence_event, h)) %>%
    filter(!is.na(endpoint), !is.na(who_grade), !is.na(integrated_grade))
  idx <- match(dat$sample_id, geo$sample_id); X <- Xall[idx, , drop = FALSE]
  Xwho <- cbind(WHO2 = as.integer(dat$who_grade == 2), WHO3 = as.integer(dat$who_grade == 3))
  Xig <- cbind(IG2 = as.integer(dat$integrated_grade == 2), IG3 = as.integer(dat$integrated_grade == 3))
  for (rep in seq_len(repeats)) {
    message("Scale-free geometry horizon ", h, " years, repetition ", rep, "/", repeats)
    outer <- group_folds(dat$patient_id, 5, 20260803L + rep)
    for (fold in 1:5) {
      tr <- which(outer != fold); te <- which(outer == fold)
      all_pred[[length(all_pred) + 1L]] <- tibble(
        sample_id = dat$sample_id[te], patient_id = dat$patient_id[te], horizon = h, repetition = rep, fold = fold,
        endpoint = dat$endpoint[te], who_grade = dat$who_grade[te], integrated_grade = dat$integrated_grade[te],
        geometry_alone = fit_predict(X[tr, ], dat$endpoint[tr], dat$patient_id[tr], X[te, ], seed0 + h * 10000 + rep * 100 + fold),
        who_geometry = fit_predict(cbind(Xwho[tr, ], X[tr, ]), dat$endpoint[tr], dat$patient_id[tr],
                                   cbind(Xwho[te, ], X[te, ]), seed0 + 100000 + h * 10000 + rep * 100 + fold),
        ig_geometry = fit_predict(cbind(Xig[tr, ], X[tr, ]), dat$endpoint[tr], dat$patient_id[tr],
                                  cbind(Xig[te, ], X[te, ]), seed0 + 200000 + h * 10000 + rep * 100 + fold)
      )
    }
  }
}
pred <- bind_rows(all_pred)
write_csv(pred, file.path(out_dir, "scale_free_geometry_multihorizon_oof_predictions.csv"))

model_map <- c("WHO grade" = "who_grade", "WHO + scale-free geometry" = "who_geometry",
               "Integrated grade" = "integrated_grade", "IG + scale-free geometry" = "ig_geometry",
               "Scale-free geometry alone" = "geometry_alone")
per_rep <- pred %>% group_by(horizon, repetition) %>% group_modify(~tibble(
  model = names(model_map), auc = vapply(model_map, function(v) auc_safe(.x$endpoint, .x[[v]]), numeric(1)))) %>% ungroup()

set.seed(seed0 + 9000L)
boot <- list()
for (h in c(3, 5, 7)) {
  dh <- pred %>% filter(horizon == h); pats <- unique(dh$patient_id)
  for (b in 1:2000) {
    wt <- table(sample(pats, length(pats), replace = TRUE))
    vals <- lapply(seq_len(repeats), function(rep) {
      d <- dh[dh$repetition == rep, ]; times <- as.integer(wt[match(d$patient_id, names(wt))]); times[is.na(times)] <- 0L
      ii <- rep(seq_len(nrow(d)), times)
      vapply(model_map, function(v) auc_safe(d$endpoint[ii], d[[v]][ii]), numeric(1))
    })
    a <- rowMeans(do.call(cbind, vals), na.rm = TRUE)
    boot[[length(boot) + 1L]] <- tibble(horizon = h, draw = b, model = names(a), auc = unname(a))
  }
}
boot <- bind_rows(boot)
summary <- per_rep %>% group_by(horizon, model) %>% summarise(auc = mean(auc), .groups = "drop") %>%
  left_join(boot %>% group_by(horizon, model) %>% summarise(low = quantile(auc, .025), high = quantile(auc, .975), .groups = "drop"),
            by = c("horizon", "model")) %>%
  left_join(pred %>% distinct(horizon, sample_id, .keep_all = TRUE) %>% group_by(horizon) %>%
              summarise(n = n(), events = sum(endpoint), .groups = "drop"), by = "horizon")
wide <- boot %>% pivot_wider(names_from = model, values_from = auc)
gains <- bind_rows(
  wide %>% transmute(horizon, draw, comparator = "WHO grade", gain = `WHO + scale-free geometry` - `WHO grade`),
  wide %>% transmute(horizon, draw, comparator = "Integrated grade", gain = `IG + scale-free geometry` - `Integrated grade`)
) %>% group_by(horizon, comparator) %>% summarise(delta_auc = mean(gain), low = quantile(gain, .025), high = quantile(gain, .975),
                                                   bootstrap_two_sided_p = 2 * min(mean(gain <= 0), mean(gain >= 0)), .groups = "drop")
write_csv(per_rep, file.path(out_dir, "scale_free_geometry_multihorizon_performance_by_repetition.csv"))
write_csv(summary, file.path(out_dir, "scale_free_geometry_multihorizon_summary.csv"))
write_csv(gains, file.path(out_dir, "scale_free_geometry_multihorizon_paired_gains.csv"))
