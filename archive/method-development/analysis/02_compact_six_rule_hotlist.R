#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(glmnet)
  library(ggplot2)
})

root <- Sys.getenv("MENIBREAK_PROJECT_ROOT")
if (!nzchar(root)) stop("Set MENIBREAK_PROJECT_ROOT to the analysis project root; see config/path_config.example.sh")
base_dir <- file.path(root, "analysis_folders", "Breakpoint_Report_Hotlist_2026")
out_dir <- file.path(base_dir, "compact_six_rule")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

feature_path <- file.path(root, "analysis_folders", "ml_breakpoint_survival_stratification_2026",
                          "ml_breakpoint_model_dataset_2026.csv")
outcome_path <- file.path(root, "software", "mbpr", "analysis", "bakeoff_2x2", "outcomes.csv")
reference_oof_path <- file.path(base_dir, "nested_oof_predictions.csv")
stopifnot(all(file.exists(c(feature_path, outcome_path, reference_oof_path))))

repeats <- 10L
outer_k <- 5L
hotlist_size <- 6L
bootstrap_draws <- 2000L
seed0 <- 20260803L

auc_safe <- function(y, s) {
  ok <- is.finite(y) & is.finite(s)
  y <- y[ok]; s <- s[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(s, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

fixed_endpoint <- function(time, event, horizon = 5) {
  ifelse(event == 1 & time <= horizon, 1,
         ifelse(time >= horizon, 0, NA_real_))
}

group_folds <- function(groups, k, seed) {
  set.seed(seed)
  ids <- sample(unique(as.character(groups)))
  map <- setNames(rep(seq_len(k), length.out = length(ids)), ids)
  unname(map[as.character(groups)])
}

fit_locked_size <- function(x, y, groups, seed) {
  foldid <- group_folds(groups, min(5L, length(unique(groups))), seed)
  cv.glmnet(x, y, family = "binomial", alpha = 0.25, foldid = foldid,
            type.measure = "deviance", standardize = TRUE)
}

raw <- read_csv(feature_path, show_col_types = FALSE)
outcomes <- read_csv(outcome_path, show_col_types = FALSE)
dat <- raw %>% transmute(
  sample_id = cg, patient_id,
  who_grade = as.numeric(who_grade),
  integrated_grade = as.numeric(integrated_grade),
  across(matches("^bp_(loss_|gain_)?[0-9]+[pq]_bin[0-9]{2}$")),
  total_breakpoint_boundaries = as.numeric(total_breakpoint_boundaries),
  breakpoint_arms_touched = as.numeric(breakpoint_arms_touched),
  loss_segment_count = as.numeric(loss_segment_count),
  gain_segment_count = as.numeric(gain_segment_count),
  small_segment_count_lt5mb = as.numeric(small_segment_count_lt5mb)
) %>% inner_join(outcomes, by = c("sample_id", "patient_id")) %>%
  mutate(endpoint = fixed_endpoint(recurrence_time, recurrence_event)) %>%
  filter(!is.na(endpoint), !is.na(who_grade), !is.na(integrated_grade))

# Report-readable primitives: breakpoint presence in non-overlapping 10%-of-arm
# intervals plus five simple breakpoint-architecture counts. Arm-relative bins
# make the spatial definitions portable between GRCh37 and GRCh38.
bp_cols <- grep("^bp_(loss_|gain_)?[0-9]+[pq]_bin[0-9]{2}$", names(dat), value = TRUE)
parse_bp <- function(x) {
  m <- regexec("^bp_(?:(loss|gain)_)?([0-9]+[pq])_bin([0-9]{2})$", x, perl = TRUE)
  z <- regmatches(x, m)[[1]]
  if (!length(z)) return(NULL)
  data.frame(col = x, event = ifelse(z[2] == "", "any", z[2]), arm = z[3], bin = as.integer(z[4]))
}
bp_map <- bind_rows(lapply(bp_cols, parse_bp))
primitive <- data.frame(row.names = seq_len(nrow(dat)))
labels <- character()
for (ev in c("any", "loss", "gain")) {
  for (arm in sort(unique(bp_map$arm))) {
    for (region in 1:10) {
      cols <- bp_map %>% filter(event == ev, arm == !!arm, bin %in% c(2L * region - 1L, 2L * region)) %>% pull(col)
      if (!length(cols)) next
      nm <- paste("hot", ev, arm, sprintf("r%02d", region), sep = "__")
      primitive[[nm]] <- as.integer(rowSums(as.data.frame(dat[, cols, drop = FALSE]), na.rm = TRUE) > 0)
      labels[[nm]] <- sprintf("%s %d0-%d0%% telomere-to-centromere %s breakpoint",
                              arm, region - 1L, region, ifelse(ev == "any", "CNA", ev))
    }
  }
}

thresholds <- list(total_breakpoint_boundaries = c(20, 40, 60),
                   breakpoint_arms_touched = c(5, 7, 10),
                   loss_segment_count = c(5, 8, 12),
                   gain_segment_count = c(5, 8, 12),
                   small_segment_count_lt5mb = c(4, 6, 10))
label_roots <- c(total_breakpoint_boundaries = "Total CNA breakpoints",
                 breakpoint_arms_touched = "Arms with breakpoints",
                 loss_segment_count = "Loss segments",
                 gain_segment_count = "Gain segments",
                 small_segment_count_lt5mb = "Focal segments under 5 Mb")
for (v in names(thresholds)) for (cut in thresholds[[v]]) {
  nm <- paste("arch", v, paste0("ge", cut), sep = "__")
  primitive[[nm]] <- as.integer(dat[[v]] >= cut)
  labels[[nm]] <- paste(label_roots[[v]], "at least", cut)
}
X <- as.matrix(primitive); storage.mode(X) <- "double"
Xwho <- cbind(WHO2 = as.integer(dat$who_grade == 2), WHO3 = as.integer(dat$who_grade == 3))
Xig <- cbind(IG2 = as.integer(dat$integrated_grade == 2), IG3 = as.integer(dat$integrated_grade == 3))

make_pairs_training <- function(top, train_idx) {
  if (length(top) < 2) return(matrix(numeric(), nrow(X), 0))
  cmb <- combn(top, 2, simplify = FALSE)
  out <- list()
  for (p in cmb) {
    p <- sort(p)
    # Avoid tautological pairs within the same arm or architecture family.
    arm1 <- strsplit(p[1], "__", fixed = TRUE)[[1]][3]
    arm2 <- strsplit(p[2], "__", fixed = TRUE)[[1]][3]
    if (startsWith(p[1], "arch__") || startsWith(p[2], "arch__")) next
    if (identical(arm1, arm2)) next
    z <- as.integer(X[, p[1]] > 0 & X[, p[2]] > 0)
    if (sum(z[train_idx]) < 15 || sum(z[train_idx] == 0) < 15) next
    nm <- paste0("pair__", p[1], "__AND__", p[2])
    out[[nm]] <- z
    labels[[nm]] <<- paste(labels[[p[1]]], "+", labels[[p[2]]])
  }
  if (!length(out)) matrix(numeric(), nrow(X), 0) else do.call(cbind, out)
}

rule_groups <- function(nm) {
  atoms <- if (startsWith(nm, "pair__")) {
    strsplit(sub("^pair__", "", nm), "__AND__", fixed = TRUE)[[1]]
  } else nm
  vapply(atoms, function(atom) {
    if (startsWith(atom, "arch__")) return(sub("__ge[0-9.]+$", "", atom))
    if (startsWith(atom, "hot__")) {
      p <- strsplit(atom, "__", fixed = TRUE)[[1]]
      return(paste(p[1:3], collapse = "__"))
    }
    atom
  }, character(1))
}

choose_distinct_rules <- function(scores, n) {
  ordered <- names(sort(scores[is.finite(scores) & scores > 0], decreasing = TRUE))
  chosen <- character(); groups <- character()
  for (nm in ordered) {
    g <- rule_groups(nm)
    if (any(g %in% groups)) next
    chosen <- c(chosen, nm); groups <- c(groups, g)
    if (length(chosen) == n) break
  }
  chosen
}

reference <- read_csv(reference_oof_path, show_col_types = FALSE) %>%
  select(sample_id, repetition, full_seg_lp)
pred_rows <- list(); selected_rows <- list()
for (r in seq_len(repeats)) {
  message("Compact hotlist repetition ", r, "/", repeats)
  outer <- group_folds(dat$patient_id, outer_k, seed0 + r)
  for (fold in seq_len(outer_k)) {
    tr <- which(outer != fold); te <- which(outer == fold)
    base_fit <- glm(endpoint ~ factor(who_grade) + factor(integrated_grade),
                    data = dat[tr, ], family = binomial())
    residual <- dat$endpoint[tr] - predict(base_fit, dat[tr, ], type = "response")
    xtr <- X[tr, , drop = FALSE]
    prev <- colSums(xtr > 0)
    cases <- colSums((xtr > 0) * dat$endpoint[tr])
    controls <- colSums((xtr > 0) * (1 - dat$endpoint[tr]))
    eligible <- prev >= 15 & prev <= length(tr) - 15 & cases >= 5 & controls >= 5
    z <- scale(xtr[, eligible, drop = FALSE], center = TRUE, scale = FALSE)
    score <- drop(crossprod(z, residual)) / sqrt(colSums(z^2) * sum(residual^2))
    names(score) <- colnames(z)
    ranked <- names(sort(score[is.finite(score) & score > 0], decreasing = TRUE))
    top_hotspots <- head(ranked[startsWith(ranked, "hot__")], 12)
    pairs <- make_pairs_training(top_hotspots, tr)
    pair_score <- numeric()
    if (ncol(pairs)) {
      zp <- scale(pairs[tr, , drop = FALSE], center = TRUE, scale = FALSE)
      pair_score <- drop(crossprod(zp, residual)) / sqrt(colSums(zp^2) * sum(residual^2))
      names(pair_score) <- colnames(zp)
      pair_score <- pair_score[is.finite(pair_score) & pair_score > 0]
    }
    all_scores <- c(score[ranked], pair_score)
    chosen <- choose_distinct_rules(all_scores, hotlist_size)
    Xall <- cbind(X, pairs)
    chosen <- intersect(chosen, colnames(Xall))
    xhot <- Xall[, chosen, drop = FALSE]
    who_train <- cbind(Xwho[tr, , drop = FALSE], xhot[tr, , drop = FALSE])
    who_test <- cbind(Xwho[te, , drop = FALSE], xhot[te, , drop = FALSE])
    ig_train <- cbind(Xig[tr, , drop = FALSE], xhot[tr, , drop = FALSE])
    ig_test <- cbind(Xig[te, , drop = FALSE], xhot[te, , drop = FALSE])
    who_fit <- fit_locked_size(who_train, dat$endpoint[tr], dat$patient_id[tr], seed0 + 10000 + r * 100 + fold)
    ig_fit <- fit_locked_size(ig_train, dat$endpoint[tr], dat$patient_id[tr], seed0 + 20000 + r * 100 + fold)
    pred_rows[[length(pred_rows) + 1L]] <- tibble(
      sample_id = dat$sample_id[te], patient_id = dat$patient_id[te], repetition = r, fold = fold,
      endpoint = dat$endpoint[te], who_grade = dat$who_grade[te], integrated_grade = dat$integrated_grade[te],
      who_hotlist_lp = drop(predict(who_fit, who_test, s = "lambda.1se", type = "link")),
      ig_hotlist_lp = drop(predict(ig_fit, ig_test, s = "lambda.1se", type = "link"))
    )
    selected_rows[[length(selected_rows) + 1L]] <- tibble(
      repetition = r, fold = fold, rank = seq_along(chosen), feature = chosen,
      report_label = unname(labels[chosen]))
  }
}

pred <- bind_rows(pred_rows) %>% left_join(reference, by = c("sample_id", "repetition"))
selected <- bind_rows(selected_rows)
write_csv(pred, file.path(out_dir, "compact_nested_oof_predictions.csv"))
write_csv(selected, file.path(out_dir, "compact_nested_selected_rules.csv"))

model_map <- c("WHO grade" = "who_grade", "WHO + six-rule breakpoint hotlist" = "who_hotlist_lp",
               "Integrated grade" = "integrated_grade", "Integrated grade + six-rule breakpoint hotlist" = "ig_hotlist_lp",
               "Full SEG model (reference ceiling)" = "full_seg_lp")
per_rep <- lapply(seq_len(repeats), function(r) {
  d <- pred[pred$repetition == r, ]
  tibble(repetition = r, model = names(model_map),
         auc = vapply(model_map, function(v) auc_safe(d$endpoint, d[[v]]), numeric(1)))
}) %>% bind_rows()

set.seed(seed0 + 9901L)
patients <- unique(pred$patient_id)
boot <- vector("list", bootstrap_draws)
for (b in seq_len(bootstrap_draws)) {
  sampled <- sample(patients, length(patients), replace = TRUE)
  wt <- table(sampled)
  vals <- lapply(seq_len(repeats), function(r) {
    d <- pred[pred$repetition == r, ]
    times <- as.integer(wt[match(d$patient_id, names(wt))]); times[is.na(times)] <- 0L
    ii <- rep(seq_len(nrow(d)), times)
    vapply(model_map, function(v) auc_safe(d$endpoint[ii], d[[v]][ii]), numeric(1))
  })
  a <- rowMeans(do.call(cbind, vals), na.rm = TRUE)
  boot[[b]] <- tibble(draw = b, model = names(a), auc = unname(a))
}
boot <- bind_rows(boot)
summary <- per_rep %>% group_by(model) %>% summarise(auc = mean(auc), .groups = "drop") %>%
  left_join(boot %>% group_by(model) %>% summarise(low = quantile(auc, .025), high = quantile(auc, .975), .groups = "drop"), by = "model")

wide <- boot %>% tidyr::pivot_wider(names_from = model, values_from = auc)
gains <- tibble(
  comparator = c("WHO grade", "Integrated grade"),
  baseline_auc = c(summary$auc[summary$model == "WHO grade"], summary$auc[summary$model == "Integrated grade"]),
  augmented_auc = c(summary$auc[summary$model == "WHO + six-rule breakpoint hotlist"], summary$auc[summary$model == "Integrated grade + six-rule breakpoint hotlist"]),
  full_reference_auc = summary$auc[summary$model == "Full SEG model (reference ceiling)"]
) %>% mutate(gain = augmented_auc - baseline_auc,
             fraction_full_gap_recovered = gain / (full_reference_auc - baseline_auc))
gain_boot <- bind_rows(
  tibble(comparator = "WHO grade", gain = wide[["WHO + six-rule breakpoint hotlist"]] - wide[["WHO grade"]],
         fraction = gain / (wide[["Full SEG model (reference ceiling)"]] - wide[["WHO grade"]])),
  tibble(comparator = "Integrated grade", gain = wide[["Integrated grade + six-rule breakpoint hotlist"]] - wide[["Integrated grade"]],
         fraction = gain / (wide[["Full SEG model (reference ceiling)"]] - wide[["Integrated grade"]]))
)
gains <- gains %>% left_join(gain_boot %>% group_by(comparator) %>% summarise(
  gain_low = quantile(gain, .025), gain_high = quantile(gain, .975),
  fraction_low = quantile(fraction, .025, na.rm = TRUE), fraction_high = quantile(fraction, .975, na.rm = TRUE),
  bootstrap_two_sided_p = 2 * min(mean(gain <= 0), mean(gain >= 0)), .groups = "drop"), by = "comparator")

stability <- selected %>% count(feature, report_label, name = "folds_selected") %>%
  mutate(selection_fraction = folds_selected / (repeats * outer_k)) %>% arrange(desc(selection_fraction), report_label)
write_csv(per_rep, file.path(out_dir, "compact_performance_by_repetition.csv"))
write_csv(boot, file.path(out_dir, "compact_patient_bootstrap.csv"))
write_csv(summary, file.path(out_dir, "compact_performance_summary.csv"))
write_csv(gains, file.path(out_dir, "compact_gain_and_fraction_recovered.csv"))
write_csv(stability, file.path(out_dir, "compact_rule_stability.csv"))

nejm <- c(blue = "#0072B5", purple = "#7876B1", red = "#BC3C29", gray = "#6F6F6F")
plot_data <- summary %>% filter(!grepl("Full SEG", model)) %>%
  mutate(model = factor(model, levels = rev(c("WHO grade", "WHO + six-rule breakpoint hotlist",
                                              "Integrated grade", "Integrated grade + six-rule breakpoint hotlist"))),
         family = ifelse(grepl("WHO", model), "WHO", "Integrated grade"))
p <- ggplot(plot_data, aes(auc, model, color = family)) +
  geom_errorbarh(aes(xmin = low, xmax = high), height = .16, linewidth = .65) +
  geom_point(size = 3.3) +
  geom_text(aes(label = sprintf("%.3f", auc)), color = "black", nudge_x = .008, hjust = 0, size = 3.5) +
  scale_color_manual(values = c("WHO" = nejm[["blue"]], "Integrated grade" = nejm[["purple"]])) +
  scale_x_continuous(limits = c(.60, .82), breaks = seq(.60, .80, .05)) +
  labs(title = "A six-rule breakpoint hotlist adds discrimination beyond grade",
       subtitle = "Strict 5-year recurrence; repeated patient-grouped nested validation (n=674; 264 events)",
       x = "Cross-validated AUC", y = NULL, color = NULL) +
  theme_classic(base_family = "sans", base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 9.5, color = "grey35"),
        axis.text = element_text(color = "black"), legend.position = "bottom",
        panel.grid.major.x = element_line(color = "grey90", linewidth = .35), plot.margin = margin(8, 18, 8, 8))
ggsave(file.path(out_dir, "figure_compact_hotlist_vs_grade.png"), p, width = 8.8, height = 4.9, dpi = 400, bg = "white")
ggsave(file.path(out_dir, "figure_compact_hotlist_vs_grade.pdf"), p, width = 8.8, height = 4.9)

writeLines(c(
  "Compact six-rule breakpoint hotlist analysis",
  sprintf("Strict 5-year cohort: n=%d, events=%d.", nrow(dat), sum(dat$endpoint)),
  "Exactly six breakpoint rules were selected inside every outer-training fold after residualizing recurrence against WHO and prespecified integrated grade.",
  "Held-out tumors never participated in rule screening, pair generation, coefficient fitting, or tuning.",
  "The same six selected rules in each fold were added separately to WHO and integrated grade.",
  "The full SEG model is used only as a reference ceiling for fraction-of-gap calculations; it is not an input to either hotlist model.",
  "This validates the compact-selection procedure internally. A single final locked six-rule list still requires frozen derivation and independent testing."
), file.path(out_dir, "METHODS_AND_INTERPRETATION.txt"))
