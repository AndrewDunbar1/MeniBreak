#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(glmnet)
})

root <- Sys.getenv("MENIBREAK_PROJECT_ROOT")
if (!nzchar(root)) stop("Set MENIBREAK_PROJECT_ROOT to the analysis project root; see config/path_config.example.sh")
out_dir <- file.path(root, "analysis_folders", "BISCUT_Breakpoint_Selection_Bakeoff_2026")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seg_path <- Sys.getenv("MENIBREAK_SEGMENTS")
if (!nzchar(seg_path)) stop("Set MENIBREAK_SEGMENTS to the event-called SEG file")
feature_path <- file.path(root, "analysis_folders", "ml_breakpoint_survival_stratification_2026",
                          "ml_breakpoint_model_dataset_2026.csv")
outcome_path <- file.path(root, "software", "mbpr", "analysis", "bakeoff_2x2", "outcomes.csv")
arm_path <- file.path(root, "software", "mbpr", "inst", "extdata", "arm_reference_GRCh37.csv")
biscut_root <- file.path(root, "software", "vendor", "BISCUT-py3")
tel_background_path <- file.path(biscut_root, "docs", "PANCAN_tels_200605.txt")
cent_background_path <- file.path(biscut_root, "docs", "PANCAN_cents_200605.txt")

required_paths <- c(seg_path, feature_path, outcome_path, arm_path,
                    tel_background_path, cent_background_path)
stopifnot(all(file.exists(required_paths)))

args <- commandArgs(trailingOnly = TRUE)
repeats <- if (length(args) >= 1) as.integer(args[[1]]) else 10L
bootstrap_draws <- if (length(args) >= 2) as.integer(args[[2]]) else 2000L
peak_bootstrap_draws <- if (length(args) >= 3) as.integer(args[[3]]) else 100L
background_mode <- if (length(args) >= 4) args[[4]] else "tcga"
if (!background_mode %in% c("tcga", "lineage")) stop("background_mode must be tcga or lineage")
output_tag <- if (background_mode == "tcga") "" else paste0("_", background_mode)
resume <- identical(Sys.getenv("BISCUT_RESUME", "1"), "1")
master_seed <- 20260803L

message("Configuration: repeats=", repeats, ", metric bootstrap=", bootstrap_draws,
        ", peak bootstrap=", peak_bootstrap_draws, ", background=", background_mode)

auc_safe <- function(y, score) {
  ok <- is.finite(y) & is.finite(score)
  y <- y[ok]; score <- score[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(score, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
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

fit_enet_predict <- function(train, test, feature_names, seed) {
  feature_names <- unique(feature_names)
  usable <- vapply(feature_names, function(v) {
    z <- suppressWarnings(as.numeric(train[[v]]))
    length(unique(z[is.finite(z)])) > 1L
  }, logical(1))
  feature_names <- feature_names[usable]
  if (!length(feature_names)) stop("No variable predictors remain")
  xtr <- as.matrix(train[, feature_names, drop = FALSE])
  xte <- as.matrix(test[, feature_names, drop = FALSE])
  storage.mode(xtr) <- storage.mode(xte) <- "double"
  med <- apply(xtr, 2, function(z) median(z[is.finite(z)], na.rm = TRUE))
  med[!is.finite(med)] <- 0
  for (j in seq_along(feature_names)) {
    xtr[!is.finite(xtr[, j]), j] <- med[[j]]
    xte[!is.finite(xte[, j]), j] <- med[[j]]
  }
  inner_k <- min(5L, n_distinct(train$patient_id))
  foldid <- group_folds(train$patient_id, inner_k, seed)
  fit <- cv.glmnet(xtr, train$endpoint, family = "binomial", alpha = 0.25,
                   foldid = foldid, type.measure = "deviance", standardize = TRUE)
  drop(predict(fit, xte, s = "lambda.1se", type = "link"))
}

union_length <- function(lo, hi) {
  z <- tibble(lo = pmax(0, lo), hi = pmin(1, hi)) %>%
    filter(is.finite(lo), is.finite(hi), hi > lo) %>% arrange(lo, hi)
  if (!nrow(z)) return(0)
  total <- 0; cur_lo <- z$lo[[1]]; cur_hi <- z$hi[[1]]
  if (nrow(z) > 1) for (i in 2:nrow(z)) {
    if (z$lo[[i]] <= cur_hi) cur_hi <- max(cur_hi, z$hi[[i]]) else {
      total <- total + cur_hi - cur_lo
      cur_lo <- z$lo[[i]]; cur_hi <- z$hi[[i]]
    }
  }
  total + cur_hi - cur_lo
}

build_partial_scna_events <- function() {
  cache <- file.path(out_dir, "biscut_partial_scna_events.csv")
  if (resume && file.exists(cache)) return(read_csv(cache, show_col_types = FALSE))
  message("Building terminal partial-SCNA event table")
  seg <- as_tibble(fread(seg_path, showProgress = FALSE)) %>%
    transmute(
      sample_id = trimws(as.character(cg)),
      chromosome = if_else(grepl("^chr", Chr), as.character(Chr), paste0("chr", Chr)),
      seg_start = as.numeric(Start), seg_end = as.numeric(End),
      event = case_when(
        Event %in% c("CN Loss", "Homozygous Copy Loss") ~ "del",
        Event %in% c("CN Gain", "High Copy Gain") ~ "amp",
        TRUE ~ NA_character_
      )
    ) %>% filter(!is.na(sample_id), nzchar(sample_id), !is.na(event), seg_end > seg_start)
  arm <- read_csv(arm_path, show_col_types = FALSE) %>%
    transmute(chromosome, arm, arm_start = as.numeric(start), arm_end = as.numeric(end),
              arm_length = as.numeric(arm_length), arm_key = paste0(sub("chr", "", chromosome), arm))
  ov <- seg %>% inner_join(arm, by = "chromosome", relationship = "many-to-many") %>%
    mutate(lo_bp = pmax(seg_start, arm_start), hi_bp = pmin(seg_end, arm_end)) %>%
    filter(hi_bp > lo_bp) %>%
    mutate(
      rel_lo = if_else(arm == "p", (lo_bp - arm_start) / arm_length,
                       (arm_end - hi_bp) / arm_length),
      rel_hi = if_else(arm == "p", (hi_bp - arm_start) / arm_length,
                       (arm_end - lo_bp) / arm_length),
      rel_lo = pmax(0, pmin(rel_lo, 1)), rel_hi = pmax(0, pmin(rel_hi, 1))
    ) %>%
    filter(!arm_key %in% c("13p", "14p", "15p", "21p", "22p"))

  terminal_tolerance <- 0.02
  minimum_altered_density <- 0.50
  events <- ov %>% group_by(sample_id, arm_key, event) %>% group_modify(~{
    min_rel <- min(.x$rel_lo); max_rel <- max(.x$rel_hi)
    span <- max_rel - min_rel
    density <- if (span > 0) union_length(.x$rel_lo, .x$rel_hi) / span else 0
    out <- list()
    if (density >= minimum_altered_density && min_rel <= terminal_tolerance) {
      out[[length(out) + 1L]] <- tibble(boundary = "tel", fraction = max_rel,
                                         altered_density = density)
    }
    if (density >= minimum_altered_density && max_rel >= 1 - terminal_tolerance) {
      out[[length(out) + 1L]] <- tibble(boundary = "cent", fraction = 1 - min_rel,
                                         altered_density = density)
    }
    if (!length(out)) tibble(boundary = character(), fraction = numeric(),
                             altered_density = numeric()) else bind_rows(out)
  }) %>% ungroup() %>%
    filter(fraction >= 0.001, fraction <= 0.999) %>%
    distinct(sample_id, arm_key, event, boundary, .keep_all = TRUE)
  write_csv(events, cache)
  events
}

fit_beta_mle <- function(x) {
  x <- x[is.finite(x) & x > 0.001 & x < 0.999]
  m <- mean(x); v <- var(x)
  common <- max(m * (1 - m) / v - 1, 0.1)
  start <- log(c(max(m * common, 0.1), max((1 - m) * common, 0.1)))
  sx <- sum(log(x)); s1 <- sum(log1p(-x)); n <- length(x)
  objective <- function(lp) {
    a <- exp(lp[[1]]); b <- exp(lp[[2]])
    -((a - 1) * sx + (b - 1) * s1 - n * lbeta(a, b))
  }
  fit <- optim(start, objective, method = "BFGS")
  setNames(exp(fit$par), c("alpha", "beta"))
}

load_background <- function() {
  cache <- file.path(out_dir, "biscut_tcga_background_beta.csv")
  if (resume && file.exists(cache)) return(read_csv(cache, show_col_types = FALSE))
  message("Fitting official TCGA BISCUT background beta distributions")
  tel <- fread(tel_background_path, select = c("percent", "amp_del"), showProgress = FALSE) %>%
    as_tibble() %>% mutate(boundary = "tel")
  cent <- fread(cent_background_path, select = c("percent", "amp_del"), showProgress = FALSE) %>%
    as_tibble() %>% mutate(boundary = "cent")
  bg <- bind_rows(tel, cent) %>% filter(percent > 0.001, percent < 0.999,
                                        amp_del %in% c("amp", "del")) %>%
    group_by(boundary, event = amp_del) %>% group_modify(~{
      fit <- fit_beta_mle(.x$percent)
      tibble(alpha = fit[["alpha"]], beta = fit[["beta"]], n_background = nrow(.x))
    }) %>% ungroup()
  write_csv(bg, cache)
  bg
}

fit_lineage_background <- function(events, train_samples) {
  events %>% filter(sample_id %in% train_samples) %>%
    group_by(boundary, event) %>% group_modify(~{
      x <- .x$fraction[is.finite(.x$fraction) & .x$fraction > 0.001 & .x$fraction < 0.999]
      if (length(x) < 20L) return(tibble(alpha = NA_real_, beta = NA_real_, n_background = length(x)))
      fit <- fit_beta_mle(x)
      tibble(alpha = fit[["alpha"]], beta = fit[["beta"]], n_background = length(x))
    }) %>% ungroup()
}

linearize_ties <- function(x) {
  x <- sort(pmax(0.001, pmin(0.999, x)))
  runs <- rle(x)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  duplicated_runs <- which(runs$lengths > 1L)
  if (!length(duplicated_runs)) return(x)
  for (j in duplicated_runs) {
    idx <- starts[[j]]:ends[[j]]; last <- ends[[j]]
    next_value <- if (last < length(x)) x[[last + 1L]] else 0.999
    if (next_value <= x[[idx[[1]]]]) next_value <- min(0.999, x[[idx[[1]]]] + 1e-4)
    spread <- seq(x[[idx[[1]]]], next_value, length.out = length(idx) + 1L)
    x[idx] <- spread[seq_along(idx)]
  }
  sort(x)
}

peak_from_distribution <- function(x, alpha, beta, low = 0.001, high = 0.999) {
  x <- linearize_ties(x[x >= low & x <= high])
  n <- length(x)
  if (n < 4) return(NULL)
  denom <- pbeta(high, alpha, beta) - pbeta(low, alpha, beta)
  bg_cdf <- (pbeta(x, alpha, beta) - pbeta(low, alpha, beta)) / denom
  divergence <- bg_cdf - seq_len(n) / n
  positive <- abs(max(divergence)) > abs(min(divergence))
  idx <- if (positive) which.max(divergence) else which.min(divergence)
  cdf_fun <- function(q) pmax(0, pmin(1,
    (pbeta(q, alpha, beta) - pbeta(low, alpha, beta)) / denom))
  ks <- suppressWarnings(ks.test(x, cdf_fun, exact = FALSE))
  list(x = x, peak = x[[idx]], selection = if (positive) "positive" else "negative",
       ks_stat = unname(ks$statistic), p_value = ks$p.value)
}

bootstrap_peak_interval <- function(x, alpha, beta, low, high, B, seed) {
  if (B < 1) return(c(max(low, min(x)), min(high, max(x))))
  set.seed(seed)
  peaks <- replicate(B, {
    z <- sample(x, length(x), replace = TRUE)
    p <- peak_from_distribution(z, alpha, beta, low, high)
    if (is.null(p)) NA_real_ else p$peak
  })
  qs <- quantile(peaks[is.finite(peaks)], c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  if (!all(is.finite(qs))) qs <- range(x)
  if (diff(qs) < 0.025) qs <- mean(qs) + c(-0.0125, 0.0125)
  pmax(low, pmin(high, qs))
}

discover_biscut_peaks <- function(events, train_samples, background, seed) {
  d <- events %>% filter(sample_id %in% train_samples)
  if (!nrow(d)) return(tibble())
  strata <- d %>% group_by(arm_key, event, boundary) %>%
    filter(n() >= 4) %>% group_split()
  tests <- vector("list", length(strata))
  for (i in seq_along(strata)) {
    s <- strata[[i]]
    bg <- background %>% filter(event == first(s$event), boundary == first(s$boundary))
    if (nrow(bg) != 1) next
    p <- peak_from_distribution(s$fraction, bg$alpha, bg$beta)
    if (is.null(p)) next
    tests[[i]] <- tibble(
      arm_key = first(s$arm_key), event = first(s$event), boundary = first(s$boundary),
      n_events = length(p$x), peak = p$peak, selection = p$selection,
      ks_stat = p$ks_stat, p_value = p$p_value,
      alpha = bg$alpha, beta = bg$beta,
      fractions = list(p$x)
    )
  }
  tests <- bind_rows(tests)
  if (!nrow(tests)) return(tibble())
  tests <- tests %>% mutate(q_value = p.adjust(p_value, method = "BY")) %>%
    filter(q_value <= 0.05)
  if (!nrow(tests)) return(tibble())
  bounds <- lapply(seq_len(nrow(tests)), function(i) {
    bootstrap_peak_interval(tests$fractions[[i]], tests$alpha[[i]], tests$beta[[i]],
                            0.001, 0.999, peak_bootstrap_draws, seed + i)
  })
  tests %>% mutate(
    lower = vapply(bounds, `[[`, numeric(1), 1),
    upper = vapply(bounds, `[[`, numeric(1), 2),
    lower = pmin(lower, peak), upper = pmax(upper, peak),
    significance_score = ks_stat * -log10(pmax(q_value, 1e-300)),
    peak_id = sprintf("%s_%s_%s_%s", arm_key, event, boundary, row_number())
  ) %>% select(-fractions, -alpha, -beta)
}

make_biscut_features <- function(events, sample_ids, peaks) {
  samples <- tibble(sample_id = unique(sample_ids))
  summary <- samples %>% mutate(
    biscut_endpoint_hits = 0, biscut_cover_hits = 0,
    biscut_endpoint_arms = 0, biscut_weighted_endpoint_score = 0,
    biscut_positive_cover_hits = 0, biscut_negative_cover_hits = 0
  )
  if (!nrow(peaks)) return(summary)
  endpoint_arms <- setNames(vector("list", nrow(summary)), summary$sample_id)
  for (i in seq_len(nrow(peaks))) {
    pk <- peaks[i, ]
    ev <- events %>% filter(arm_key == pk$arm_key, event == pk$event,
                            boundary == pk$boundary, sample_id %in% sample_ids)
    endpoint_ids <- unique(ev$sample_id[ev$fraction >= pk$lower & ev$fraction <= pk$upper])
    cover_ids <- unique(ev$sample_id[ev$fraction >= pk$peak])
    ep_col <- paste0("biscut_endpoint_", sprintf("%03d", i))
    cv_col <- paste0("biscut_covers_", sprintf("%03d", i))
    summary[[ep_col]] <- as.integer(summary$sample_id %in% endpoint_ids)
    summary[[cv_col]] <- as.integer(summary$sample_id %in% cover_ids)
    summary$biscut_endpoint_hits <- summary$biscut_endpoint_hits + summary[[ep_col]]
    summary$biscut_cover_hits <- summary$biscut_cover_hits + summary[[cv_col]]
    summary$biscut_weighted_endpoint_score <- summary$biscut_weighted_endpoint_score +
      summary[[ep_col]] * pk$significance_score
    if (pk$selection == "positive") {
      summary$biscut_positive_cover_hits <- summary$biscut_positive_cover_hits + summary[[cv_col]]
    } else {
      summary$biscut_negative_cover_hits <- summary$biscut_negative_cover_hits + summary[[cv_col]]
    }
    for (sid in endpoint_ids) endpoint_arms[[sid]] <- unique(c(endpoint_arms[[sid]], pk$arm_key))
  }
  summary$biscut_endpoint_arms <- vapply(summary$sample_id,
                                         function(s) length(endpoint_arms[[s]]), integer(1))
  summary
}

raw <- read_csv(feature_path, show_col_types = FALSE)
outcomes <- read_csv(outcome_path, show_col_types = FALSE)
partial_events <- build_partial_scna_events()
background <- load_background()

scna_features <- grep("^cnv_|^scna_(loss_count|gain_count|altered_arm_count)$",
                      names(raw), value = TRUE)
breakpoint_candidates <- c(
  "total_breakpoint_boundaries", "breakpoint_arms_touched",
  "breakpoint_route_entropy", "max_arm_breakpoint_fraction",
  "telomeric_breakpoint_fraction", "central_breakpoint_fraction",
  "centromeric_breakpoint_fraction", "mean_rel_tel_to_cent",
  "sd_rel_tel_to_cent", "cna_segment_count", "loss_segment_count",
  "gain_segment_count", "median_segment_length", "small_segment_count_lt5mb",
  "dynamic_hotspot_breakpoints_all", "dynamic_hotspot_breakpoints_non1p",
  "dynamic_hotspot_arms_used_all", "dynamic_hotspot_arms_used_non1p"
)
breakpoint_features <- intersect(breakpoint_candidates, names(raw))
dynamic_hotspot_features <- grep("^dynamic_hotspot_", breakpoint_features, value = TRUE)
generic_breakpoint_features <- setdiff(breakpoint_features, dynamic_hotspot_features)
clinical_features <- c("age", "log1p_mitoses", "gtr", "prior_rt", "adj_rt")

dat <- raw %>% transmute(
  sample_id = cg,
  age = as.numeric(age), log1p_mitoses = log1p(as.numeric(mitoses)),
  gtr = as.numeric(gtr), prior_rt = as.numeric(prior_rt), adj_rt = as.numeric(adj_rt),
  across(all_of(c(scna_features, breakpoint_features)))
) %>% inner_join(outcomes, by = "sample_id")

baseline_features <- c(scna_features, clinical_features)
current_features <- c(scna_features, clinical_features, breakpoint_features)

prediction_cache <- file.path(out_dir, sprintf("biscut_oof_predictions%s_%drepeats.rds", output_tag, repeats))
peak_cache <- file.path(out_dir, sprintf("biscut_fold_peaks%s_%drepeats.rds", output_tag, repeats))
if (resume && file.exists(prediction_cache) && file.exists(peak_cache)) {
  pred <- readRDS(prediction_cache)
  peak_audit <- readRDS(peak_cache)
} else {
  pred_list <- list(); peak_list <- list()
  for (h in c(3, 5, 7)) {
    d <- dat %>% mutate(endpoint = fixed_horizon_endpoint(recurrence_time, recurrence_event, h)) %>%
      filter(!is.na(endpoint))
    for (r in seq_len(repeats)) {
      message("Horizon ", h, ", repetition ", r, "/", repeats)
      outer <- group_folds(d$patient_id, 5L, 20260703L + h * 10000L + r)
      p_base <- p_current <- p_replace <- p_hybrid <- rep(NA_real_, nrow(d))
      for (fold in seq_len(5L)) {
        tr <- d[outer != fold, , drop = FALSE]
        te <- d[outer == fold, , drop = FALSE]
        fold_background <- if (background_mode == "lineage") {
          fit_lineage_background(partial_events, tr$sample_id)
        } else background
        peaks <- discover_biscut_peaks(partial_events, tr$sample_id, fold_background,
                                       master_seed + h * 100000L + r * 1000L + fold)
        if (nrow(peaks)) {
          peak_list[[length(peak_list) + 1L]] <- peaks %>% mutate(
            horizon_years = h, repetition = r, fold = fold,
            background_mode = background_mode,
            n_train = nrow(tr), events_train = sum(tr$endpoint)
          )
        } else {
          peak_list[[length(peak_list) + 1L]] <- tibble(
            horizon_years = h, repetition = r, fold = fold,
            background_mode = background_mode,
            n_train = nrow(tr), events_train = sum(tr$endpoint)
          )
        }
        bf <- make_biscut_features(partial_events, c(tr$sample_id, te$sample_id), peaks)
        trb <- tr %>% left_join(bf, by = "sample_id")
        teb <- te %>% left_join(bf, by = "sample_id")
        bcols <- grep("^biscut_", names(trb), value = TRUE)
        replacement_features <- c(baseline_features, generic_breakpoint_features, bcols)
        hybrid_features <- c(current_features, bcols)
        p_base[outer == fold] <- fit_enet_predict(
          tr, te, baseline_features,
          20260703L + h * 100000L + r * 1000L + fold
        )
        p_current[outer == fold] <- fit_enet_predict(
          tr, te, current_features,
          20260703L + h * 200000L + r * 1000L + fold
        )
        p_replace[outer == fold] <- fit_enet_predict(
          trb, teb, replacement_features,
          20260703L + h * 300000L + r * 1000L + fold
        )
        p_hybrid[outer == fold] <- fit_enet_predict(
          trb, teb, hybrid_features,
          20260703L + h * 400000L + r * 1000L + fold
        )
      }
      pred_list[[length(pred_list) + 1L]] <- tibble(
        sample_id = d$sample_id, patient_id = d$patient_id,
        horizon_years = h, repetition = r, endpoint = d$endpoint,
        clinical_scna = p_base, current_breakpoints = p_current,
        biscut_replacement = p_replace, biscut_hybrid = p_hybrid
      )
      pred <- bind_rows(pred_list); peak_audit <- bind_rows(peak_list)
      saveRDS(pred, prediction_cache, version = 3)
      saveRDS(peak_audit, peak_cache, version = 3)
    }
  }
}

write_csv(pred, file.path(out_dir, paste0("biscut_oof_predictions", output_tag, ".csv")))
write_csv(peak_audit, file.path(out_dir, paste0("biscut_fold_peak_audit", output_tag, ".csv")))

model_map <- c(
  "Clinical + SCNA" = "clinical_scna",
  "Current clinical + SCNA + breakpoint architecture" = "current_breakpoints",
  "BISCUT replacement" = "biscut_replacement",
  "Current + BISCUT hybrid" = "biscut_hybrid"
)

bootstrap_horizon <- function(d, B, seed) {
  set.seed(seed)
  patients <- unique(d$patient_id); reps <- sort(unique(d$repetition))
  draws <- vector("list", B)
  for (b in seq_len(B)) {
    wt <- table(sample(patients, length(patients), replace = TRUE))
    per_rep <- lapply(reps, function(r) {
      dr <- d[d$repetition == r, , drop = FALSE]
      times <- as.integer(wt[match(dr$patient_id, names(wt))]); times[is.na(times)] <- 0L
      idx <- rep(seq_len(nrow(dr)), times = times)
      vapply(model_map, function(v) auc_safe(dr$endpoint[idx], dr[[v]][idx]), numeric(1))
    })
    aucs <- rowMeans(do.call(cbind, per_rep), na.rm = TRUE)
    draws[[b]] <- tibble(draw = b, model = names(model_map), auc = aucs,
                         delta_vs_baseline = aucs - aucs[["Clinical + SCNA"]],
                         delta_vs_current = aucs - aucs[["Current clinical + SCNA + breakpoint architecture"]])
  }
  bind_rows(draws)
}

point <- bind_rows(lapply(c(3, 5, 7), function(h) {
  d <- filter(pred, horizon_years == h)
  bind_rows(lapply(names(model_map), function(m) {
    per_rep <- d %>% group_by(repetition) %>% summarise(
      n = n(), events = sum(endpoint), auc = auc_safe(endpoint, .data[[model_map[[m]]]]), .groups = "drop")
    tibble(horizon_years = h, model = m, n = first(per_rep$n), events = first(per_rep$events),
           auc = mean(per_rep$auc))
  }))
}))

boots <- bind_rows(lapply(c(3, 5, 7), function(h) {
  bootstrap_horizon(filter(pred, horizon_years == h), bootstrap_draws, master_seed + h) %>%
    mutate(horizon_years = h)
}))

summary <- point %>% left_join(
  boots %>% group_by(horizon_years, model) %>% summarise(
    low = quantile(auc, 0.025), high = quantile(auc, 0.975),
    delta_vs_baseline_low = quantile(delta_vs_baseline, 0.025),
    delta_vs_baseline_high = quantile(delta_vs_baseline, 0.975),
    delta_vs_baseline = median(delta_vs_baseline),
    delta_vs_current_low = quantile(delta_vs_current, 0.025),
    delta_vs_current_high = quantile(delta_vs_current, 0.975),
    p_vs_current = pmin(1, 2 * min(mean(delta_vs_current <= 0), mean(delta_vs_current >= 0))),
    delta_vs_current = median(delta_vs_current),
    .groups = "drop"), by = c("horizon_years", "model"))
write_csv(summary, file.path(out_dir, paste0("biscut_model_performance_summary", output_tag, ".csv")))

peak_stability <- peak_audit %>% filter(!is.na(arm_key)) %>%
  group_by(arm_key, event, boundary, selection) %>% summarise(
    folds_selected = n(), total_folds = repeats * 5L * 3L,
    selection_fraction = folds_selected / total_folds,
    median_peak = median(peak), median_lower = median(lower), median_upper = median(upper),
    median_q = median(q_value), median_significance_score = median(significance_score),
    .groups = "drop") %>% arrange(desc(selection_fraction), median_q)
write_csv(peak_stability, file.path(out_dir, paste0("biscut_peak_stability_summary", output_tag, ".csv")))

fold_peak_counts <- peak_audit %>% group_by(horizon_years, repetition, fold) %>%
  summarise(n_peaks = sum(!is.na(arm_key)), .groups = "drop")
write_csv(fold_peak_counts, file.path(out_dir, paste0("biscut_fold_peak_counts", output_tag, ".csv")))

cols <- c(
  "Clinical + SCNA" = "#0072B5",
  "Current clinical + SCNA + breakpoint architecture" = "#BC3C29",
  "BISCUT replacement" = "#7876B1",
  "Current + BISCUT hybrid" = "#E18727"
)
shapes <- c(
  "Clinical + SCNA" = 15,
  "Current clinical + SCNA + breakpoint architecture" = 18,
  "BISCUT replacement" = 17,
  "Current + BISCUT hybrid" = 16
)
display_names <- c(
  "Clinical + SCNA" = "Clinical + SCNA",
  "Current clinical + SCNA + breakpoint architecture" = "Current breakpoint model",
  "BISCUT replacement" = "BISCUT replacement",
  "Current + BISCUT hybrid" = "Current + BISCUT"
)
theme_pub <- theme_classic(base_family = "sans", base_size = 12.2) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10.3, color = "grey35"),
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 8.8), panel.grid.major.y = element_line(color = "grey90", linewidth = .35),
        plot.margin = margin(10, 12, 10, 10))

plot_summary <- summary %>% mutate(model = factor(model, levels = names(cols)))
p1 <- ggplot(plot_summary, aes(horizon_years, auc, color = model, shape = model, group = model)) +
  geom_line(linewidth = .85) + geom_point(size = 2.4) +
  geom_errorbar(aes(ymin = low, ymax = high), width = .09, linewidth = .5) +
  scale_color_manual(values = cols, labels = display_names) +
  scale_shape_manual(values = shapes, labels = display_names) +
  scale_x_continuous(breaks = c(3, 5, 7)) +
  scale_y_continuous(limits = c(.75, .90), breaks = seq(.75, .90, .025)) +
  labs(title = "A  Recurrence discrimination", x = "Recurrence horizon (years)", y = "Cross-validated AUC") +
  theme_pub + theme(legend.position = "none")

delta_base <- plot_summary %>% filter(model != "Clinical + SCNA")
p2 <- ggplot(delta_base, aes(horizon_years, delta_vs_baseline, color = model, shape = model, group = model)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60") + geom_line(linewidth = .8) +
  geom_point(size = 2.3) +
  geom_errorbar(aes(ymin = delta_vs_baseline_low, ymax = delta_vs_baseline_high), width = .09, linewidth = .5) +
  scale_color_manual(values = cols, labels = display_names) +
  scale_shape_manual(values = shapes, labels = display_names) +
  scale_x_continuous(breaks = c(3, 5, 7)) +
  labs(title = "B  Breakpoint benefit", x = "Recurrence horizon (years)",
       y = "AUC difference vs clinical + SCNA") + theme_pub + theme(legend.position = "none")

delta_current <- plot_summary %>% filter(model %in% c("BISCUT replacement", "Current + BISCUT hybrid"))
p3 <- ggplot(delta_current, aes(horizon_years, delta_vs_current, color = model, shape = model, group = model)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60") + geom_line(linewidth = .8) +
  geom_point(size = 2.3) +
  geom_errorbar(aes(ymin = delta_vs_current_low, ymax = delta_vs_current_high), width = .09, linewidth = .5) +
  scale_color_manual(values = cols, labels = display_names) +
  scale_shape_manual(values = shapes, labels = display_names) +
  scale_x_continuous(breaks = c(3, 5, 7)) +
  labs(title = "C  Direct comparison with current breakpoint model",
       x = "Recurrence horizon (years)", y = "AUC difference vs current model") +
  theme_pub + theme(legend.position = "none")

fig <- (p1 | p2) / p3 + plot_layout(heights = c(1, .72)) +
  plot_annotation(
    title = paste0("BISCUT-guided breakpoint selection bakeoff (", background_mode, " background)"),
    subtitle = paste0(
      "Candidate peaks were learned from outer-training tumors only; all models used identical patient-grouped nested folds\n",
      "Blue square: Clinical + SCNA   |   Red diamond: Current breakpoint model   |   Purple triangle: BISCUT replacement   |   Orange circle: Current + BISCUT"
    ),
    caption = "BISCUT replacement swaps the four frozen hotspot summaries for fold-derived BISCUT peak features; hybrid adds BISCUT features to the current architecture.",
    theme = theme(plot.title = element_text(family = "sans", face = "bold", size = 16),
                  plot.subtitle = element_text(family = "sans", size = 10.2, color = "grey35"),
                  plot.caption = element_text(family = "sans", size = 8.8, color = "grey40", hjust = 0),
                  plot.margin = margin(12, 14, 10, 14)))

ggsave(file.path(out_dir, paste0("figure4_biscut_breakpoint_selection_bakeoff", output_tag, ".png")), fig,
       width = 10.2, height = 8.0, dpi = 400, bg = "white")
ggsave(file.path(out_dir, paste0("figure4_biscut_breakpoint_selection_bakeoff", output_tag, ".pdf")), fig,
       width = 10.2, height = 8.0)

write_lines(c(
  "BISCUT-guided breakpoint-selection bakeoff",
  "",
  "Implementation boundary:",
  paste0("- Background mode: ", background_mode, ". TCGA uses official pan-cancer distributions; lineage fits pooled institutional partial-SCNA lengths inside each outer-training fold."),
  "- Event-called loss/gain segments were converted into telomere- or centromere-bounded partial-SCNA lengths on GRCh37 arms.",
  "- Within each outer training fold, arm/event/boundary length distributions were compared with the matched TCGA beta background.",
  "- The maximal signed CDF divergence defined the BISCUT peak; Benjamini-Yekutieli q <= 0.05 defined selection.",
  "- A fold-specific bootstrap interval around each peak created endpoint and peak-coverage features.",
  "- No recurrence outcomes and no outer-test tumors participated in BISCUT peak discovery.",
  "- This is an event-called SEG adaptation of BISCUT, not a claim that the original TCGA continuous-log2 pipeline ran unchanged.",
  "",
  paste0("Repeated patient-grouped nested validation: ", repeats, " repeats x 5 outer folds; alpha=0.25; lambda.1se."),
  paste0("Patient-bootstrap draws for performance intervals: ", bootstrap_draws, "."),
  paste0("Bootstrap draws for fold-specific peak intervals: ", peak_bootstrap_draws, ".")
), file.path(out_dir, paste0("METHODS_biscut_breakpoint_selection", output_tag, ".md")))

print(summary)
print(head(peak_stability, 20))
