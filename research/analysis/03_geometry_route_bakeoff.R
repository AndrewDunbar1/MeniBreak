#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(glmnet)
  library(ggplot2)
})

root <- Sys.getenv("MENIBREAK_PROJECT_ROOT")
if (!nzchar(root)) stop("Set MENIBREAK_PROJECT_ROOT to the analysis project root; see config/path_config.example.sh")
out_dir <- file.path(root, "analysis_folders", "Breakpoint_Geometry_Route_2026")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
feature_path <- file.path(root, "analysis_folders", "ml_breakpoint_survival_stratification_2026",
                          "ml_breakpoint_model_dataset_2026.csv")
outcome_path <- file.path(root, "software", "mbpr", "analysis", "bakeoff_2x2", "outcomes.csv")
reference_path <- file.path(root, "analysis_folders", "Breakpoint_Report_Hotlist_2026", "nested_oof_predictions.csv")
stopifnot(all(file.exists(c(feature_path, outcome_path, reference_path))))

repeats <- 10L; outer_k <- 5L; bootstrap_draws <- 2000L; seed0 <- 20260805L
arms <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
          "7p","7q","8p","8q","9p","9q","10p","10q","11p","11q","12p","12q",
          "13q","14q","15q","16p","16q","17p","17q","18p","18q","19p","19q",
          "20p","20q","21q","22q")
bins <- 1:20

auc_safe <- function(y, s) {
  ok <- is.finite(y) & is.finite(s); y <- y[ok]; s <- s[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(s, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
fixed_endpoint <- function(time, event, h = 5) ifelse(event == 1 & time <= h, 1, ifelse(time >= h, 0, NA_real_))
group_folds <- function(groups, k, seed) {
  set.seed(seed); ids <- sample(unique(as.character(groups)))
  map <- setNames(rep(seq_len(k), length.out = length(ids)), ids)
  unname(map[as.character(groups)])
}
entropy_rows <- function(p) {
  z <- p; z[z <= 0] <- 1
  -rowSums(ifelse(p > 0, p * log(z), 0))
}
safe_div <- function(a, b) { z <- a / b; z[!is.finite(z)] <- 0; z }

raw <- read_csv(feature_path, show_col_types = FALSE)
outcomes <- read_csv(outcome_path, show_col_types = FALSE)
n <- nrow(raw)

make_array <- function(event = c("any", "loss", "gain")) {
  event <- match.arg(event)
  prefix <- if (event == "any") "bp_" else paste0("bp_", event, "_")
  a <- array(0, dim = c(n, length(arms), length(bins)),
             dimnames = list(NULL, arms, sprintf("%02d", bins)))
  for (j in seq_along(arms)) for (b in bins) {
    nm <- paste0(prefix, arms[[j]], "_bin", sprintf("%02d", b))
    if (nm %in% names(raw)) a[, j, b] <- as.numeric(raw[[nm]])
  }
  a[!is.finite(a)] <- 0
  a
}

arrays <- list(any = make_array("any"), loss = make_array("loss"), gain = make_array("gain"))

route_features_one <- function(a, tag) {
  n <- dim(a)[1]; arm_counts <- apply(a, c(1, 2), sum); bin_counts <- apply(a, c(1, 3), sum)
  total <- rowSums(arm_counts); arm_p <- safe_div(arm_counts, total); bin_p <- safe_div(bin_counts, total)
  touched <- arm_counts > 0; occupied <- a > 0
  arm_entropy <- entropy_rows(arm_p)
  radial_entropy <- entropy_rows(bin_p)
  pos <- (bins - .5) / 20
  radial_mean <- rowSums(sweep(bin_p, 2, pos, `*`))
  radial_var <- rowSums(sweep(bin_p, 2, pos^2, `*`)) - radial_mean^2
  radial_sd <- sqrt(pmax(radial_var, 0))
  chromosomes <- sub("[pq]$", "", arms)
  chr_touched <- sapply(unique(chromosomes), function(chr) rowSums(touched[, chromosomes == chr, drop = FALSE]) > 0)
  if (is.vector(chr_touched)) chr_touched <- matrix(chr_touched, ncol = 1)
  both_arms <- sapply(unique(chromosomes), function(chr) {
    ii <- which(chromosomes == chr)
    if (length(ii) == 2) rowSums(touched[, ii, drop = FALSE]) == 2 else rep(FALSE, n)
  })
  if (is.vector(both_arms)) both_arms <- matrix(both_arms, ncol = 1)
  clusters <- matrix(0, n, length(arms)); spans <- matrix(0, n, length(arms)); isolated <- matrix(0, n, length(arms))
  for (i in seq_len(n)) for (j in seq_along(arms)) {
    z <- occupied[i, j, ]
    if (!any(z)) next
    rr <- rle(z); clusters[i, j] <- sum(rr$values)
    idx <- which(z); spans[i, j] <- (max(idx) - min(idx) + 1) / 20
    run_lengths <- rr$lengths[rr$values]
    isolated[i, j] <- sum(run_lengths == 1)
  }
  arm_n <- rowSums(touched)
  occupied_n <- apply(occupied, 1, sum)
  shape <- tibble(
    !!paste0(tag, "__arm_entropy_norm") := safe_div(arm_entropy, log(pmax(arm_n, 2))),
    !!paste0(tag, "__effective_arm_fraction") := safe_div(exp(arm_entropy), pmax(arm_n, 1)),
    !!paste0(tag, "__arm_hhi") := rowSums(arm_p^2),
    !!paste0(tag, "__max_arm_fraction") := apply(arm_p, 1, max),
    !!paste0(tag, "__p_arm_fraction") := rowSums(arm_p[, grepl("p$", arms), drop = FALSE]),
    !!paste0(tag, "__radial_mean") := radial_mean,
    !!paste0(tag, "__radial_sd") := radial_sd,
    !!paste0(tag, "__radial_entropy_norm") := radial_entropy / log(20),
    !!paste0(tag, "__radial_hhi") := rowSums(bin_p^2),
    !!paste0(tag, "__telomeric_fraction") := rowSums(bin_p[, 1:4, drop = FALSE]),
    !!paste0(tag, "__outer_mid_fraction") := rowSums(bin_p[, 5:8, drop = FALSE]),
    !!paste0(tag, "__middle_fraction") := rowSums(bin_p[, 9:12, drop = FALSE]),
    !!paste0(tag, "__inner_mid_fraction") := rowSums(bin_p[, 13:16, drop = FALSE]),
    !!paste0(tag, "__centromeric_fraction") := rowSums(bin_p[, 17:20, drop = FALSE]),
    !!paste0(tag, "__cluster_per_occupied_bin") := safe_div(rowSums(clusters), occupied_n),
    !!paste0(tag, "__isolated_cluster_fraction") := safe_div(rowSums(isolated), rowSums(clusters)),
    !!paste0(tag, "__mean_within_arm_span") := safe_div(rowSums(spans), arm_n)
  )
  topology <- bind_cols(shape, tibble(
    !!paste0(tag, "__arms_touched") := arm_n,
    !!paste0(tag, "__chromosomes_touched") := rowSums(chr_touched),
    !!paste0(tag, "__both_arms_chromosomes") := rowSums(both_arms),
    !!paste0(tag, "__occupied_bins") := occupied_n,
    !!paste0(tag, "__multi_break_arms") := rowSums(arm_counts >= 2),
    !!paste0(tag, "__route_clusters") := rowSums(clusters)
  ))
  arm_fraction <- as_tibble(arm_p, .name_repair = "minimal")
  names(arm_fraction) <- paste0(tag, "__arm_fraction__", arms)
  list(shape = shape, topology = topology, arm_fraction = arm_fraction,
       arm_counts = arm_counts, bin_p = bin_p, total = total, touched = touched)
}

r <- lapply(names(arrays), function(tag) route_features_one(arrays[[tag]], tag)); names(r) <- names(arrays)

js_rows <- function(p, q) {
  m <- (p + q) / 2
  kl <- function(a, b) rowSums(ifelse(a > 0, a * log(a / pmax(b, 1e-12)), 0))
  (kl(p, m) + kl(q, m)) / 2
}
loss_arm_p <- safe_div(r$loss$arm_counts, r$loss$total)
gain_arm_p <- safe_div(r$gain$arm_counts, r$gain$total)
cross_route <- tibble(
  loss_boundary_fraction = safe_div(r$loss$total, r$loss$total + r$gain$total),
  loss_gain_arm_js = js_rows(loss_arm_p, gain_arm_p),
  loss_gain_radial_js = js_rows(r$loss$bin_p, r$gain$bin_p),
  loss_gain_arm_jaccard = safe_div(rowSums(r$loss$touched & r$gain$touched), rowSums(r$loss$touched | r$gain$touched)),
  loss_gain_arm_cosine = safe_div(rowSums(loss_arm_p * gain_arm_p),
                                  sqrt(rowSums(loss_arm_p^2) * rowSums(gain_arm_p^2)))
)

shape <- bind_cols(r$any$shape, r$loss$shape, r$gain$shape, cross_route)
topology <- bind_cols(r$any$topology, r$loss$topology, r$gain$topology, cross_route)
complexity <- bind_cols(topology, tibble(
  log_total_breakpoints = log1p(as.numeric(raw$total_breakpoint_boundaries)),
  log_cna_segments = log1p(as.numeric(raw$cna_segment_count)),
  log_loss_segments = log1p(as.numeric(raw$loss_segment_count)),
  log_gain_segments = log1p(as.numeric(raw$gain_segment_count)),
  log_focal_segments = log1p(as.numeric(raw$small_segment_count_lt5mb)),
  log_median_segment_length = log1p(as.numeric(raw$median_segment_length)),
  focal_segment_fraction = safe_div(as.numeric(raw$small_segment_count_lt5mb), as.numeric(raw$cna_segment_count)),
  loss_segment_fraction = safe_div(as.numeric(raw$loss_segment_count),
                                   as.numeric(raw$loss_segment_count) + as.numeric(raw$gain_segment_count))
))
coarse_map <- bind_cols(complexity, r$any$arm_fraction, r$loss$arm_fraction, r$gain$arm_fraction)

clean_matrix <- function(x) {
  x <- as.matrix(x); storage.mode(x) <- "double"; x[!is.finite(x)] <- 0
  x
}
feature_sets <- list(shape = clean_matrix(shape), topology = clean_matrix(topology),
                     geometry_complexity = clean_matrix(complexity), coarse_route_map = clean_matrix(coarse_map))
write_csv(bind_cols(tibble(sample_id = raw$cg), complexity), file.path(out_dir, "geometry_route_features.csv"))
write_csv(tibble(feature_set = names(feature_sets), n_features = vapply(feature_sets, ncol, integer(1))),
          file.path(out_dir, "feature_set_manifest.csv"))

dat <- raw %>% transmute(sample_id = cg, patient_id, who_grade = as.numeric(who_grade),
                         integrated_grade = as.numeric(integrated_grade)) %>%
  inner_join(outcomes, by = c("sample_id", "patient_id")) %>%
  mutate(endpoint = fixed_endpoint(recurrence_time, recurrence_event, 5)) %>%
  filter(!is.na(endpoint), !is.na(who_grade), !is.na(integrated_grade))
idx <- match(dat$sample_id, raw$cg)
feature_sets <- lapply(feature_sets, function(x) x[idx, , drop = FALSE])
Xwho <- cbind(WHO2 = as.integer(dat$who_grade == 2), WHO3 = as.integer(dat$who_grade == 3))
Xig <- cbind(IG2 = as.integer(dat$integrated_grade == 2), IG3 = as.integer(dat$integrated_grade == 3))

fit_predict <- function(xtr, ytr, groups, xte, seed) {
  usable <- apply(xtr, 2, function(z) length(unique(z[is.finite(z)])) > 1)
  xtr <- xtr[, usable, drop = FALSE]; xte <- xte[, usable, drop = FALSE]
  foldid <- group_folds(groups, min(5L, length(unique(groups))), seed)
  fit <- cv.glmnet(xtr, ytr, family = "binomial", alpha = .25, foldid = foldid,
                   type.measure = "deviance", standardize = TRUE)
  drop(predict(fit, xte, s = "lambda.1se", type = "link"))
}

pred_rows <- list()
for (rep in seq_len(repeats)) {
  message("Geometry route repetition ", rep, "/", repeats)
  outer <- group_folds(dat$patient_id, outer_k, 20260803L + rep)
  for (fold in seq_len(outer_k)) {
    tr <- which(outer != fold); te <- which(outer == fold)
    row <- tibble(sample_id = dat$sample_id[te], patient_id = dat$patient_id[te], repetition = rep,
                  fold = fold, endpoint = dat$endpoint[te], who_grade = dat$who_grade[te],
                  integrated_grade = dat$integrated_grade[te])
    for (set_name in names(feature_sets)) {
      x <- feature_sets[[set_name]]
      row[[paste0("who_", set_name)]] <- fit_predict(cbind(Xwho[tr, ], x[tr, ]), dat$endpoint[tr], dat$patient_id[tr],
                                                       cbind(Xwho[te, ], x[te, ]), seed0 + rep * 1000 + fold * 10)
      row[[paste0("ig_", set_name)]] <- fit_predict(cbind(Xig[tr, ], x[tr, ]), dat$endpoint[tr], dat$patient_id[tr],
                                                      cbind(Xig[te, ], x[te, ]), seed0 + 100000 + rep * 1000 + fold * 10)
      row[[paste0("only_", set_name)]] <- fit_predict(x[tr, ], dat$endpoint[tr], dat$patient_id[tr], x[te, ],
                                                       seed0 + 200000 + rep * 1000 + fold * 10)
    }
    pred_rows[[length(pred_rows) + 1L]] <- row
  }
}
pred <- bind_rows(pred_rows)
reference <- read_csv(reference_path, show_col_types = FALSE) %>% select(sample_id, repetition, full_seg_lp)
pred <- pred %>% left_join(reference, by = c("sample_id", "repetition"))
write_csv(pred, file.path(out_dir, "geometry_route_nested_oof_predictions.csv"))

model_map <- c("WHO grade" = "who_grade", "Integrated grade" = "integrated_grade",
               "WHO + scale-free geometry" = "who_shape", "IG + scale-free geometry" = "ig_shape",
               "WHO + route topology" = "who_topology", "IG + route topology" = "ig_topology",
               "WHO + geometry/complexity" = "who_geometry_complexity", "IG + geometry/complexity" = "ig_geometry_complexity",
               "WHO + coarse route map" = "who_coarse_route_map", "IG + coarse route map" = "ig_coarse_route_map",
               "Geometry/complexity alone" = "only_geometry_complexity",
               "Full SEG model (reference)" = "full_seg_lp")
per_rep <- lapply(seq_len(repeats), function(rep) {
  d <- pred[pred$repetition == rep, ]
  tibble(repetition = rep, model = names(model_map),
         auc = vapply(model_map, function(v) auc_safe(d$endpoint, d[[v]]), numeric(1)))
}) %>% bind_rows()

set.seed(seed0 + 9000L); pats <- unique(pred$patient_id)
boots <- vector("list", bootstrap_draws)
for (b in seq_len(bootstrap_draws)) {
  wt <- table(sample(pats, length(pats), replace = TRUE))
  vals <- lapply(seq_len(repeats), function(rep) {
    d <- pred[pred$repetition == rep, ]; times <- as.integer(wt[match(d$patient_id, names(wt))]); times[is.na(times)] <- 0L
    ii <- rep(seq_len(nrow(d)), times)
    vapply(model_map, function(v) auc_safe(d$endpoint[ii], d[[v]][ii]), numeric(1))
  })
  a <- rowMeans(do.call(cbind, vals), na.rm = TRUE)
  boots[[b]] <- tibble(draw = b, model = names(a), auc = unname(a))
}
boots <- bind_rows(boots)
summary <- per_rep %>% group_by(model) %>% summarise(auc = mean(auc), .groups = "drop") %>%
  left_join(boots %>% group_by(model) %>% summarise(low = quantile(auc, .025), high = quantile(auc, .975), .groups = "drop"), by = "model")
wide <- boots %>% pivot_wider(names_from = model, values_from = auc)
comparisons <- bind_rows(lapply(c("scale-free geometry", "route topology", "geometry/complexity", "coarse route map"), function(suffix) {
  who_aug <- paste("WHO +", suffix); ig_aug <- paste("IG +", suffix)
  bind_rows(
    tibble(comparator = "WHO grade", augmentation = suffix,
           gain = wide[[who_aug]] - wide[["WHO grade"]]),
    tibble(comparator = "Integrated grade", augmentation = suffix,
           gain = wide[[ig_aug]] - wide[["Integrated grade"]])
  )
})) %>% group_by(comparator, augmentation) %>% summarise(
  delta_auc = mean(gain), low = quantile(gain, .025), high = quantile(gain, .975),
  bootstrap_two_sided_p = 2 * min(mean(gain <= 0), mean(gain >= 0)), .groups = "drop")
write_csv(per_rep, file.path(out_dir, "geometry_route_performance_by_repetition.csv"))
write_csv(boots, file.path(out_dir, "geometry_route_patient_bootstrap.csv"))
write_csv(summary, file.path(out_dir, "geometry_route_performance_summary.csv"))
write_csv(comparisons, file.path(out_dir, "geometry_route_paired_gains.csv"))

plot_models <- c("WHO grade", "WHO + scale-free geometry", "WHO + route topology", "WHO + geometry/complexity",
                 "Integrated grade", "IG + scale-free geometry", "IG + route topology", "IG + geometry/complexity")
pd <- summary %>% filter(model %in% plot_models) %>%
  mutate(model = factor(model, levels = rev(plot_models)), family = ifelse(grepl("WHO", model), "WHO", "Integrated grade"))
p <- ggplot(pd, aes(auc, model, color = family)) +
  geom_errorbar(aes(xmin = low, xmax = high), orientation = "y", width = .16, linewidth = .6) +
  geom_point(size = 3.2) + geom_text(aes(label = sprintf("%.3f", auc)), color = "black", nudge_x = .007, hjust = 0, size = 3.3) +
  scale_color_manual(values = c("WHO" = "#0072B5", "Integrated grade" = "#7876B1")) +
  scale_x_continuous(limits = c(.60, .86), breaks = seq(.60, .85, .05)) +
  labs(title = "Breakpoint geometry adds route information beyond grade",
       subtitle = "No thresholded arm-level SCNA calls; strict 5-year recurrence, repeated patient-grouped nested validation",
       x = "Cross-validated AUC", y = NULL, color = NULL) +
  theme_classic(base_family = "sans", base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(color = "grey35", size = 9.5),
        axis.text = element_text(color = "black"), panel.grid.major.x = element_line(color = "grey90", linewidth = .35),
        legend.position = "bottom", plot.margin = margin(8, 18, 8, 8))
ggsave(file.path(out_dir, "figure_geometry_route_vs_grade.png"), p, width = 9.4, height = 6.2, dpi = 400, bg = "white")
ggsave(file.path(out_dir, "figure_geometry_route_vs_grade.pdf"), p, width = 9.4, height = 6.2)

writeLines(c(
  "Breakpoint geometry/route bakeoff",
  "Primary endpoint: strict 5-year recurrence.",
  "All geometry features are fixed transformations of arm-relative breakpoint-bin occupancy and counts.",
  "No thresholded arm-level loss/gain calls, SCNA burden counts, clinical treatment variables, or mitotic count enter these models.",
  "Scale-free geometry captures normalized dispersion, concentration, radial position, and loss/gain route asymmetry.",
  "Route topology adds numbers of arms/chromosomes traversed, bilateral chromosome involvement, occupied bins, clusters, and multi-break arms.",
  "Geometry/complexity additionally adds continuous breakpoint and segment fragmentation measures.",
  "All coefficients and penalties were estimated inside patient-grouped training folds; held-out tumors contributed only to evaluation."
), file.path(out_dir, "METHODS_AND_FEATURE_DEFINITIONS.txt"))
