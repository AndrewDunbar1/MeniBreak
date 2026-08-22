#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(glmnet)
  library(survival)
  library(ggplot2)
  library(patchwork)
})

root <- Sys.getenv("MENIBREAK_PROJECT_ROOT")
if (!nzchar(root)) stop("Set MENIBREAK_PROJECT_ROOT to the analysis project root; see config/path_config.example.sh")
out_dir <- file.path(root, "analysis_folders", "Breakpoint_Report_Hotlist_2026")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

feature_path <- file.path(root, "analysis_folders", "ml_breakpoint_survival_stratification_2026",
                          "ml_breakpoint_model_dataset_2026.csv")
outcome_path <- file.path(root, "software", "mbpr", "analysis", "bakeoff_2x2", "outcomes.csv")
arm_path <- file.path(root, "software", "mbpr", "inst", "extdata", "arm_reference_GRCh37.csv")
cytoband_path <- file.path(out_dir, "cytoBand_hg19.txt.gz")
stopifnot(all(file.exists(c(feature_path, outcome_path, arm_path, cytoband_path))))

args <- commandArgs(trailingOnly = TRUE)
repeats <- if (length(args) >= 1) as.integer(args[[1]]) else 10L
bootstrap_draws <- if (length(args) >= 2) as.integer(args[[2]]) else 2000L
seed0 <- 20260803L

auc_safe <- function(y, s) {
  ok <- is.finite(y) & is.finite(s); y <- y[ok]; s <- s[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(s, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
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

median_impute <- function(xtr, xte) {
  med <- apply(xtr, 2, function(z) median(z[is.finite(z)], na.rm = TRUE))
  med[!is.finite(med)] <- 0
  for (j in seq_len(ncol(xtr))) {
    xtr[!is.finite(xtr[, j]), j] <- med[[j]]
    xte[!is.finite(xte[, j]), j] <- med[[j]]
  }
  list(train = xtr, test = xte, median = med)
}

fit_cv <- function(xtr, ytr, groups, seed, alpha = 0.25, family = "binomial") {
  usable <- apply(xtr, 2, function(z) length(unique(z[is.finite(z)])) > 1L)
  xtr <- xtr[, usable, drop = FALSE]
  foldid <- group_folds(groups, min(5L, length(unique(groups))), seed)
  fit <- cv.glmnet(xtr, ytr, family = family, alpha = alpha, foldid = foldid,
                   type.measure = if (family == "binomial") "deviance" else "deviance",
                   standardize = TRUE)
  list(fit = fit, columns = colnames(xtr))
}

predict_fit <- function(obj, x, type = "link") {
  x <- x[, obj$columns, drop = FALSE]
  drop(predict(obj$fit, x, s = "lambda.1se", type = type))
}

coef_fit <- function(obj) {
  z <- as.matrix(coef(obj$fit, s = "lambda.1se"))[, 1]
  z[z != 0 & names(z) != "(Intercept)"]
}

arm_ref <- read_csv(arm_path, show_col_types = FALSE) %>%
  mutate(arm_key = paste0(sub("chr", "", chromosome), arm))
cyto <- fread(cytoband_path, header = FALSE, col.names = c("chromosome", "start", "end", "band", "stain")) %>%
  as_tibble()

region_annotation <- function(arm_key, region) {
  a <- arm_ref %>% filter(.data$arm_key == !!arm_key)
  if (!nrow(a)) return(tibble(chr = NA_character_, start_bp = NA_real_, end_bp = NA_real_, cytoband = NA_character_))
  rel_lo <- (region - 1) / 10; rel_hi <- region / 10
  if (a$arm[[1]] == "p") {
    lo <- a$start[[1]] + rel_lo * a$arm_length[[1]]
    hi <- a$start[[1]] + rel_hi * a$arm_length[[1]]
  } else {
    lo <- a$end[[1]] - rel_hi * a$arm_length[[1]]
    hi <- a$end[[1]] - rel_lo * a$arm_length[[1]]
  }
  bands <- cyto %>% filter(chromosome == a$chromosome[[1]], end > lo, start < hi) %>% pull(band)
  band_label <- if (!length(bands)) NA_character_ else if (length(bands) == 1) bands else paste0(first(bands), "-", last(bands))
  tibble(chr = a$chromosome[[1]], start_bp = round(lo), end_bp = round(hi), cytoband = band_label)
}

raw <- read_csv(feature_path, show_col_types = FALSE)
outcomes <- read_csv(outcome_path, show_col_types = FALSE)

scna_features <- grep("^cnv_|^scna_(loss_count|gain_count|altered_arm_count)$", names(raw), value = TRUE)
clinical_features <- c("age", "log1p_mitoses", "gtr", "prior_rt", "adj_rt")
current_breakpoint_candidates <- c(
  "total_breakpoint_boundaries", "breakpoint_arms_touched", "breakpoint_route_entropy",
  "max_arm_breakpoint_fraction", "telomeric_breakpoint_fraction", "central_breakpoint_fraction",
  "centromeric_breakpoint_fraction", "mean_rel_tel_to_cent", "sd_rel_tel_to_cent",
  "cna_segment_count", "loss_segment_count", "gain_segment_count", "median_segment_length",
  "small_segment_count_lt5mb", "dynamic_hotspot_breakpoints_all",
  "dynamic_hotspot_breakpoints_non1p", "dynamic_hotspot_arms_used_all",
  "dynamic_hotspot_arms_used_non1p"
)
current_breakpoint_features <- intersect(current_breakpoint_candidates, names(raw))

dat <- raw %>% transmute(
  sample_id = cg, patient_id,
  age = as.numeric(age), log1p_mitoses = log1p(as.numeric(mitoses)),
  gtr = as.numeric(gtr), prior_rt = as.numeric(prior_rt), adj_rt = as.numeric(adj_rt),
  primary = as.numeric(primary), who_grade = as.numeric(who_grade), integrated_grade = as.numeric(integrated_grade),
  across(all_of(scna_features)),
  across(all_of(current_breakpoint_features)),
  across(matches("^bp_(loss_|gain_)?[0-9]+[pq]_bin[0-9]{2}$")),
  total_breakpoint_boundaries = as.numeric(total_breakpoint_boundaries),
  breakpoint_arms_touched = as.numeric(breakpoint_arms_touched),
  loss_segment_count = as.numeric(loss_segment_count),
  gain_segment_count = as.numeric(gain_segment_count),
  small_segment_count_lt5mb = as.numeric(small_segment_count_lt5mb)
) %>% inner_join(outcomes, by = c("sample_id", "patient_id"))

# Build fixed, report-readable primitive features: presence of a breakpoint in
# a non-overlapping 10%-of-arm interval, stratified as any/loss/gain.
primitive <- tibble(sample_id = dat$sample_id)
meta <- list()
arms <- arm_ref$arm_key
for (ev in c("any", "loss", "gain")) {
  prefix <- if (ev == "any") "bp_" else paste0("bp_", ev, "_")
  for (arm in arms) for (region in 1:10) {
    bins <- c(2 * region - 1, 2 * region)
    cols <- paste0(prefix, arm, "_bin", sprintf("%02d", bins))
    cols <- intersect(cols, names(dat))
    if (!length(cols)) next
    nm <- paste("hot", ev, arm, sprintf("r%02d", region), sep = "__")
    primitive[[nm]] <- as.integer(rowSums(as.data.frame(dat[, cols, drop = FALSE]), na.rm = TRUE) > 0)
    ann <- region_annotation(arm, region)
    direction <- if (ev == "any") "Any CNA" else tools::toTitleCase(ev)
    meta[[length(meta) + 1L]] <- tibble(
      feature = nm, feature_type = "hotspot", event_type = ev, arm_key = arm,
      region = region, start_fraction_tel_to_cent = (region - 1) / 10,
      end_fraction_tel_to_cent = region / 10,
      chr = ann$chr, start_bp = ann$start_bp, end_bp = ann$end_bp, cytoband = ann$cytoband,
      report_label = paste0(arm, " ", ann$cytoband, " ", ifelse(ev == "any", "CNA", tolower(direction)), " breakpoint")
    )
  }
}

# Simple global patterns that can be abstracted from a report without a SEG file.
threshold_specs <- list(
  total_breakpoint_boundaries = c(10, 20, 30, 40, 60),
  breakpoint_arms_touched = c(3, 5, 7, 10),
  loss_segment_count = c(3, 5, 8, 12),
  gain_segment_count = c(3, 5, 8, 12),
  small_segment_count_lt5mb = c(2, 4, 6, 10)
)
for (v in names(threshold_specs)) for (cut in threshold_specs[[v]]) {
  nm <- paste("arch", v, paste0("ge", cut), sep = "__")
  primitive[[nm]] <- as.integer(dat[[v]] >= cut)
  label_root <- c(
    total_breakpoint_boundaries = "Total CNA breakpoints",
    breakpoint_arms_touched = "Chromosome arms with breakpoints",
    loss_segment_count = "Loss segments", gain_segment_count = "Gain segments",
    small_segment_count_lt5mb = "Focal (<5 Mb) segments"
  )[[v]]
  meta[[length(meta) + 1L]] <- tibble(
    feature = nm, feature_type = "architecture", event_type = NA_character_, arm_key = NA_character_,
    region = NA_integer_, start_fraction_tel_to_cent = NA_real_, end_fraction_tel_to_cent = NA_real_,
    chr = NA_character_, start_bp = NA_real_, end_bp = NA_real_, cytoband = NA_character_,
    report_label = paste0(label_root, " >= ", cut)
  )
}
meta <- bind_rows(meta)
candidate_names <- setdiff(names(primitive), "sample_id")
Xcand <- as.matrix(primitive[, candidate_names, drop = FALSE]); storage.mode(Xcand) <- "double"

Xwho <- cbind(who_grade_2 = as.numeric(dat$who_grade == 2), who_grade_3 = as.numeric(dat$who_grade == 3))
Xig <- cbind(ig_grade_2 = as.numeric(dat$integrated_grade == 2), ig_grade_3 = as.numeric(dat$integrated_grade == 3))
Xboth <- cbind(Xwho, Xig)
current_scna_names <- c(scna_features, clinical_features)
current_full_names <- c(scna_features, clinical_features, current_breakpoint_features)
Xcurrent_scna <- as.matrix(dat[, current_scna_names, drop = FALSE]); storage.mode(Xcurrent_scna) <- "double"
Xcurrent_full <- as.matrix(dat[, current_full_names, drop = FALSE]); storage.mode(Xcurrent_full) <- "double"

make_pairs <- function(top_names, x, metadata, min_n = 15L) {
  if (length(top_names) < 2) return(list(x = matrix(numeric(), nrow(x), 0), meta = tibble()))
  cmb <- combn(top_names, 2, simplify = FALSE)
  out <- list(); mm <- list()
  for (pair in cmb) {
    pair <- sort(pair)
    m1 <- metadata %>% filter(feature == pair[[1]])
    m2 <- metadata %>% filter(feature == pair[[2]])
    if (nrow(m1) && nrow(m2) && !is.na(m1$arm_key) && !is.na(m2$arm_key) && m1$arm_key == m2$arm_key) next
    z <- as.integer(x[, pair[[1]]] > 0 & x[, pair[[2]]] > 0)
    if (sum(z) < min_n || sum(z == 0) < min_n) next
    nm <- paste0("pair__", pair[[1]], "__AND__", pair[[2]])
    out[[nm]] <- z
    mm[[length(mm) + 1L]] <- tibble(
      feature = nm, feature_type = "pair", event_type = NA_character_, arm_key = NA_character_,
      region = NA_integer_, start_fraction_tel_to_cent = NA_real_, end_fraction_tel_to_cent = NA_real_,
      chr = NA_character_, start_bp = NA_real_, end_bp = NA_real_, cytoband = NA_character_,
      report_label = paste0(m1$report_label, " + ", m2$report_label)
    )
  }
  if (!length(out)) return(list(x = matrix(numeric(), nrow(x), 0), meta = tibble()))
  list(x = do.call(cbind, out), meta = bind_rows(mm))
}

primary <- dat %>% mutate(endpoint = fixed_endpoint(recurrence_time, recurrence_event, 5)) %>%
  filter(!is.na(endpoint))
idx_primary <- match(primary$sample_id, dat$sample_id)
Xc_primary <- Xcand[idx_primary, , drop = FALSE]
Xwho_primary <- Xwho[idx_primary, , drop = FALSE]
Xig_primary <- Xig[idx_primary, , drop = FALSE]
Xboth_primary <- Xboth[idx_primary, , drop = FALSE]
Xcurrent_scna_primary <- Xcurrent_scna[idx_primary, , drop = FALSE]
Xcurrent_full_primary <- Xcurrent_full[idx_primary, , drop = FALSE]

pred_list <- list(); selection_list <- list(); screen_list <- list(); contribution_list <- list()
for (r in seq_len(repeats)) {
  message("Primary 5-year nested stability repetition ", r, "/", repeats)
  outer <- group_folds(primary$patient_id, 5L, seed0 + r)
  for (fold in 1:5) {
    tri <- which(outer != fold); tei <- which(outer == fold)
    grade_train <- primary[tri, , drop = FALSE]
    grade_fit <- glm(endpoint ~ factor(who_grade) + factor(integrated_grade),
                     data = grade_train, family = binomial())
    grade_both_lp_train <- drop(predict(grade_fit, grade_train, type = "link"))
    base_prob_train <- plogis(grade_both_lp_train)
    residual <- primary$endpoint[tri] - base_prob_train

    xc_tr <- Xc_primary[tri, , drop = FALSE]
    prev <- colSums(xc_tr > 0)
    case_n <- colSums((xc_tr > 0) * primary$endpoint[tri])
    ctrl_n <- colSums((xc_tr > 0) * (1 - primary$endpoint[tri]))
    eligible <- prev >= 15 & prev <= length(tri) - 15 & case_n >= 5 & ctrl_n >= 5
    score <- rep(-Inf, ncol(xc_tr)); names(score) <- colnames(xc_tr)
    if (any(eligible)) {
      z <- scale(xc_tr[, eligible, drop = FALSE], center = TRUE, scale = FALSE)
      denom <- sqrt(colSums(z^2) * sum(residual^2))
      score[eligible] <- drop(crossprod(z, residual)) / denom
    }
    ranked <- names(sort(score[is.finite(score) & score > 0], decreasing = TRUE))
    top_primitive <- head(ranked, 30)
    # Recursive combinations are restricted to spatial hotspots. Combining
    # nested burden thresholds (for example >=20 AND >=30 breakpoints) creates
    # tautological rules rather than a biologically interpretable pattern.
    hotspot_ranked <- ranked[ranked %in% meta$feature[meta$feature_type == "hotspot"]]
    top_for_pairs <- head(hotspot_ranked, 12)
    pairs <- make_pairs(top_for_pairs, Xc_primary, meta, min_n = 15L)
    pair_scores <- numeric()
    if (ncol(pairs$x)) {
      xp <- pairs$x[tri, , drop = FALSE]
      zp <- scale(xp, center = TRUE, scale = FALSE)
      pair_scores <- drop(crossprod(zp, residual)) / sqrt(colSums(zp^2) * sum(residual^2))
      names(pair_scores) <- colnames(xp)
      pair_scores <- sort(pair_scores[is.finite(pair_scores) & pair_scores > 0], decreasing = TRUE)
    }
    top_pairs <- head(names(pair_scores), 30)
    fold_candidates <- c(top_primitive, top_pairs)
    xhot <- cbind(Xc_primary[, top_primitive, drop = FALSE],
                  if (length(top_pairs)) pairs$x[, top_pairs, drop = FALSE] else NULL)
    x_who_hot <- cbind(Xwho_primary, xhot)
    x_ig_hot <- cbind(Xig_primary, xhot)
    x_both_hot <- cbind(Xboth_primary, xhot)
    imp_hot <- median_impute(xhot[tri, , drop = FALSE], xhot[tei, , drop = FALSE])
    imp_who <- median_impute(x_who_hot[tri, , drop = FALSE], x_who_hot[tei, , drop = FALSE])
    imp_ig <- median_impute(x_ig_hot[tri, , drop = FALSE], x_ig_hot[tei, , drop = FALSE])
    imp_both <- median_impute(x_both_hot[tri, , drop = FALSE], x_both_hot[tei, , drop = FALSE])
    imp_scna <- median_impute(Xcurrent_scna_primary[tri, , drop = FALSE], Xcurrent_scna_primary[tei, , drop = FALSE])
    imp_full <- median_impute(Xcurrent_full_primary[tri, , drop = FALSE], Xcurrent_full_primary[tei, , drop = FALSE])
    hot_obj <- fit_cv(imp_hot$train, primary$endpoint[tri], primary$patient_id[tri],
                      seed0 + 10000L + r * 100L + fold)
    who_hot_obj <- fit_cv(imp_who$train, primary$endpoint[tri], primary$patient_id[tri],
                          seed0 + 20000L + r * 100L + fold)
    ig_hot_obj <- fit_cv(imp_ig$train, primary$endpoint[tri], primary$patient_id[tri],
                         seed0 + 30000L + r * 100L + fold)
    both_hot_obj <- fit_cv(imp_both$train, primary$endpoint[tri], primary$patient_id[tri],
                           seed0 + 40000L + r * 100L + fold)
    scna_obj <- fit_cv(imp_scna$train, primary$endpoint[tri], primary$patient_id[tri],
                       seed0 + 50000L + r * 100L + fold)
    full_obj <- fit_cv(imp_full$train, primary$endpoint[tri], primary$patient_id[tri],
                       seed0 + 60000L + r * 100L + fold)
    grade_both_test <- primary[tei, ]
    grade_both_lp_test <- drop(predict(grade_fit, grade_both_test, type = "link"))
    pred_list[[length(pred_list) + 1L]] <- tibble(
      sample_id = primary$sample_id[tei], patient_id = primary$patient_id[tei],
      repetition = r, fold = fold, endpoint = primary$endpoint[tei],
      recurrence_time = primary$recurrence_time[tei], recurrence_event = primary$recurrence_event[tei],
      who_grade = primary$who_grade[tei], integrated_grade = primary$integrated_grade[tei],
      grade_both_lp = grade_both_lp_test,
      hotlist_lp = predict_fit(hot_obj, imp_hot$test),
      who_hotlist_lp = predict_fit(who_hot_obj, imp_who$test),
      ig_hotlist_lp = predict_fit(ig_hot_obj, imp_ig$test),
      both_hotlist_lp = predict_fit(both_hot_obj, imp_both$test),
      clinical_scna_lp = predict_fit(scna_obj, imp_scna$test),
      full_seg_lp = predict_fit(full_obj, imp_full$test)
    )
    co <- coef_fit(both_hot_obj)
    selected <- intersect(names(co), fold_candidates)
    if (length(selected)) {
      selection_list[[length(selection_list) + 1L]] <- tibble(
        repetition = r, fold = fold, feature = selected, coefficient = unname(co[selected])
      )
      x_test_named <- imp_both$test[, selected, drop = FALSE]
      contrib <- sweep(x_test_named, 2, co[selected], `*`)
      contribution_list[[length(contribution_list) + 1L]] <- as_tibble(contrib) %>%
        mutate(sample_id = primary$sample_id[tei], repetition = r, fold = fold, .before = 1) %>%
        pivot_longer(-c(sample_id, repetition, fold), names_to = "feature", values_to = "contribution")
    }
    screen_list[[length(screen_list) + 1L]] <- tibble(
      repetition = r, fold = fold, feature = names(score)[is.finite(score)],
      screen_score = score[is.finite(score)], screened_top30 = names(score)[is.finite(score)] %in% top_primitive
    )
  }
}

pred <- bind_rows(pred_list)
selected <- bind_rows(selection_list)
screens <- bind_rows(screen_list)
contributions <- bind_rows(contribution_list)
write_csv(pred, file.path(out_dir, "nested_oof_predictions.csv"))
write_csv(selected, file.path(out_dir, "nested_selected_features.csv"))
write_csv(screens, file.path(out_dir, "nested_screening_audit.csv"))
write_csv(contributions, file.path(out_dir, "nested_case_contributions.csv"))

# Attach pair metadata discovered in at least one fold.
pair_features <- unique(selected$feature[startsWith(selected$feature, "pair__")])
pair_meta <- lapply(pair_features, function(nm) {
  pieces <- strsplit(sub("^pair__", "", nm), "__AND__", fixed = FALSE)[[1]]
  m1 <- meta %>% filter(feature == pieces[[1]]); m2 <- meta %>% filter(feature == pieces[[2]])
  tibble(feature = nm, feature_type = "pair", event_type = NA_character_, arm_key = NA_character_,
         region = NA_integer_, start_fraction_tel_to_cent = NA_real_, end_fraction_tel_to_cent = NA_real_,
         chr = NA_character_, start_bp = NA_real_, end_bp = NA_real_, cytoband = NA_character_,
         report_label = paste0(m1$report_label, " + ", m2$report_label))
}) %>% bind_rows()
all_meta <- bind_rows(meta, pair_meta) %>% distinct(feature, .keep_all = TRUE)

stability <- selected %>% group_by(feature) %>% summarise(
  folds_selected = n(), selection_fraction = n() / (repeats * 5),
  positive_fraction = mean(coefficient > 0), median_coefficient = median(coefficient),
  .groups = "drop"
) %>% full_join(
  screens %>% group_by(feature) %>% summarise(
    mean_screen_score = mean(screen_score), top30_fraction = mean(screened_top30), .groups = "drop"
  ), by = "feature"
) %>% left_join(all_meta, by = "feature")

# Build matrices for every shortlisted feature so adjusted survival associations
# and report-level prevalences use a single fixed definition.
get_feature_vector <- function(nm) {
  if (nm %in% colnames(Xcand)) return(Xcand[, nm])
  if (startsWith(nm, "pair__")) {
    pieces <- strsplit(sub("^pair__", "", nm), "__AND__", fixed = FALSE)[[1]]
    return(as.integer(get_feature_vector(pieces[[1]]) > 0 & get_feature_vector(pieces[[2]]) > 0))
  }
  rep(NA_real_, nrow(dat))
}

shortlist <- stability %>% filter(selection_fraction > 0 | top30_fraction >= 0.10) %>% pull(feature) %>% unique()

# Descriptive survival associations are adjusted for both grade systems. No
# arm-level SCNA features enter hotlist discovery or this association model.
surv_ok <- is.finite(dat$recurrence_time) & is.finite(dat$recurrence_event) & dat$recurrence_time > 0

association_rows <- lapply(shortlist, function(nm) {
  z_all <- get_feature_vector(nm)
  z <- z_all[surv_ok]
  dd <- tibble(time = dat$recurrence_time[surv_ok], event = dat$recurrence_event[surv_ok],
               patient_id = dat$patient_id[surv_ok], marker = z,
               who_grade = factor(dat$who_grade[surv_ok]),
               integrated_grade = factor(dat$integrated_grade[surv_ok]))
  if (sum(z > 0, na.rm = TRUE) < 10 || length(unique(z[is.finite(z)])) < 2) return(NULL)
  fit <- try(coxph(Surv(time, event) ~ marker + who_grade + integrated_grade + cluster(patient_id),
                   data = dd, ties = "efron", singular.ok = TRUE), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  sm <- summary(fit); rr <- sm$coefficients["marker", ]; ci <- sm$conf.int["marker", ]
  p3 <- fixed_endpoint(dat$recurrence_time, dat$recurrence_event, 3)
  p5 <- fixed_endpoint(dat$recurrence_time, dat$recurrence_event, 5)
  p7 <- fixed_endpoint(dat$recurrence_time, dat$recurrence_event, 7)
  eligible3 <- !is.na(p3); eligible7 <- !is.na(p7)
  eligible5 <- !is.na(p5)
  pos3 <- eligible3 & z_all > 0; neg3 <- eligible3 & z_all == 0
  pos5 <- eligible5 & z_all > 0; neg5 <- eligible5 & z_all == 0
  pos7 <- eligible7 & z_all > 0; neg7 <- eligible7 & z_all == 0
  rec_times <- dat$recurrence_time[dat$recurrence_event == 1 & z_all > 0]
  tibble(feature = nm, prevalence = mean(z_all > 0, na.rm = TRUE), n_positive = sum(z_all > 0, na.rm = TRUE),
         events_positive = sum(dat$recurrence_event == 1 & z_all > 0, na.rm = TRUE),
         adjustment_model = "WHO grade + prespecified integrated grade",
         adjusted_hr = unname(ci[["exp(coef)"]]), hr_low = unname(ci[["lower .95"]]),
         hr_high = unname(ci[["upper .95"]]), p_value = unname(rr[["Pr(>|z|)"]]),
         recurrence_3y_positive = mean(p3[pos3]), recurrence_3y_negative = mean(p3[neg3]),
         recurrence_5y_positive = mean(p5[pos5]), recurrence_5y_negative = mean(p5[neg5]),
         risk_difference_5y = mean(p5[pos5]) - mean(p5[neg5]),
         recurrence_7y_positive = mean(p7[pos7]), recurrence_7y_negative = mean(p7[neg7]),
         median_recurrence_years_positive = if (length(rec_times)) median(rec_times) else NA_real_)
}) %>% bind_rows() %>% mutate(q_value = p.adjust(p_value, method = "BH"))

hotlist <- stability %>% left_join(association_rows, by = "feature") %>%
  mutate(
    evidence_tier = case_when(
      selection_fraction >= 0.40 & positive_fraction >= 0.90 & adjusted_hr >= 1.5 & q_value <= 0.10 ~ "Tier 1: high-signal",
      selection_fraction >= 0.20 & positive_fraction >= 0.80 & adjusted_hr >= 1.3 & q_value <= 0.20 ~ "Tier 2: supported",
      selection_fraction >= 0.10 & positive_fraction >= 0.70 & adjusted_hr > 1 ~ "Tier 3: exploratory",
      TRUE ~ "Not retained"
    ),
    gradient_score = selection_fraction * pmax(log(adjusted_hr), 0) * pmax(risk_difference_5y, 0)
  ) %>% arrange(desc(gradient_score), desc(selection_fraction))

write_csv(all_meta, file.path(out_dir, "candidate_feature_dictionary.csv"))
write_csv(stability, file.path(out_dir, "nested_feature_stability.csv"))
write_csv(association_rows, file.path(out_dir, "adjusted_survival_associations.csv"))
write_csv(hotlist, file.path(out_dir, "breakpoint_hotlist_ranked.csv"))

# Collapse correlated report rules: retain one cutoff per architecture family,
# one 10%-of-arm region per event/arm, and individually stable cross-arm pairs.
hotlist_reduced <- hotlist %>%
  mutate(collapse_group = case_when(
    feature_type == "architecture" ~ sub("__ge[0-9.]+$", "", feature),
    feature_type == "hotspot" ~ paste(feature_type, event_type, arm_key, sep = "__"),
    TRUE ~ feature
  )) %>%
  group_by(collapse_group) %>% arrange(desc(gradient_score), desc(selection_fraction), .by_group = TRUE) %>%
  slice_head(n = 1) %>% ungroup() %>%
  filter(evidence_tier != "Not retained") %>% arrange(desc(gradient_score), desc(selection_fraction))
write_csv(hotlist_reduced, file.path(out_dir, "breakpoint_hotlist_reduced_reportable.csv"))

# Performance of the report-friendly models against WHO and integrated grade;
# the full SEG model is a reference ceiling rather than the hotlist comparator.
model_map <- c(
  "WHO grade" = "who_grade",
  "Integrated grade" = "integrated_grade",
  "Breakpoint hotlist alone" = "hotlist_lp",
  "WHO + breakpoint hotlist" = "who_hotlist_lp",
  "Integrated grade + breakpoint hotlist" = "ig_hotlist_lp",
  "WHO + integrated grade + breakpoint hotlist" = "both_hotlist_lp",
  "Clinical + SCNA" = "clinical_scna_lp",
  "Full BiLab SEG model" = "full_seg_lp"
)
per_rep <- pred %>% group_by(repetition) %>% summarise(
  across(all_of(unname(model_map)), ~auc_safe(endpoint, .x), .names = "auc__{.col}"), .groups = "drop")
write_csv(per_rep, file.path(out_dir, "nested_performance_by_repetition.csv"))

set.seed(seed0 + 901L)
patients <- unique(pred$patient_id)
boot_rows <- vector("list", bootstrap_draws)
for (b in seq_len(bootstrap_draws)) {
  wt <- table(sample(patients, length(patients), replace = TRUE))
  per_r <- lapply(sort(unique(pred$repetition)), function(r) {
    d <- pred[pred$repetition == r, ]
    times <- as.integer(wt[match(d$patient_id, names(wt))]); times[is.na(times)] <- 0
    ii <- rep(seq_len(nrow(d)), times)
    vapply(model_map, function(v) auc_safe(d$endpoint[ii], d[[v]][ii]), numeric(1))
  })
  aucs <- rowMeans(do.call(cbind, per_r), na.rm = TRUE)
  who_gain <- aucs[["WHO + breakpoint hotlist"]] - aucs[["WHO grade"]]
  ig_gain <- aucs[["Integrated grade + breakpoint hotlist"]] - aucs[["Integrated grade"]]
  who_ceiling <- aucs[["Full BiLab SEG model"]] - aucs[["WHO grade"]]
  ig_ceiling <- aucs[["Full BiLab SEG model"]] - aucs[["Integrated grade"]]
  boot_rows[[b]] <- tibble(
    draw = b, model = names(aucs), auc = unname(aucs),
    who_gain = who_gain, ig_gain = ig_gain,
    who_fraction_recovered = ifelse(who_ceiling > 0, who_gain / who_ceiling, NA_real_),
    ig_fraction_recovered = ifelse(ig_ceiling > 0, ig_gain / ig_ceiling, NA_real_)
  )
}
boots <- bind_rows(boot_rows)
point_auc <- vapply(model_map, function(v) mean(per_rep[[paste0("auc__", v)]]), numeric(1))
performance <- tibble(model = names(point_auc), auc = unname(point_auc)) %>%
  left_join(boots %>% group_by(model) %>% summarise(
    auc_low = quantile(auc, .025), auc_high = quantile(auc, .975), .groups = "drop"), by = "model") %>%
  mutate(horizon_years = 5, n = nrow(primary), events = sum(primary$endpoint))
recovery <- tibble(
  comparator = c("WHO grade", "Integrated grade"),
  grade_auc = c(point_auc[["WHO grade"]], point_auc[["Integrated grade"]]),
  grade_plus_hotlist_auc = c(point_auc[["WHO + breakpoint hotlist"]], point_auc[["Integrated grade + breakpoint hotlist"]]),
  full_seg_auc = point_auc[["Full BiLab SEG model"]]
) %>% mutate(
  hotlist_gain = grade_plus_hotlist_auc - grade_auc,
  full_gap = full_seg_auc - grade_auc,
  fraction_of_full_gap_recovered = hotlist_gain / full_gap,
  gain_low = c(quantile(unique(boots$who_gain), .025), quantile(unique(boots$ig_gain), .025)),
  gain_high = c(quantile(unique(boots$who_gain), .975), quantile(unique(boots$ig_gain), .975)),
  fraction_low = c(quantile(unique(boots$who_fraction_recovered), .025, na.rm = TRUE),
                   quantile(unique(boots$ig_fraction_recovered), .025, na.rm = TRUE)),
  fraction_high = c(quantile(unique(boots$who_fraction_recovered), .975, na.rm = TRUE),
                    quantile(unique(boots$ig_fraction_recovered), .975, na.rm = TRUE))
)
write_csv(performance, file.path(out_dir, "nested_performance_summary.csv"))
write_csv(recovery, file.path(out_dir, "hotlist_fraction_of_full_model_benefit.csv"))
write_csv(boots, file.path(out_dir, "nested_performance_patient_bootstrap.csv"))

# Reverse analysis: favorable-treatment early recurrences whose risk ranking was
# increased by the report-friendly breakpoint layer.
pred_ranked <- pred %>% group_by(repetition) %>% mutate(
  baseline_pct = percent_rank(grade_both_lp), hotlist_pct = percent_rank(both_hotlist_lp),
  delta_pct = hotlist_pct - baseline_pct) %>% ungroup() %>%
  group_by(sample_id, patient_id) %>% summarise(
    endpoint = first(endpoint), recurrence_time = first(recurrence_time), recurrence_event = first(recurrence_event),
    baseline_pct = mean(baseline_pct), hotlist_pct = mean(hotlist_pct), delta_pct = mean(delta_pct), .groups = "drop") %>%
  left_join(dat %>% select(sample_id, gtr, prior_rt, adj_rt, primary, who_grade, integrated_grade), by = "sample_id")

reverse_cases <- pred_ranked %>% filter(recurrence_event == 1, recurrence_time <= 3,
                                        gtr == 1, prior_rt == 0, adj_rt == 0,
                                        baseline_pct <= 0.50, delta_pct > 0) %>%
  arrange(desc(delta_pct))

case_contrib <- contributions %>% group_by(sample_id, feature) %>% summarise(
  mean_contribution = mean(contribution), positive_repeat_fraction = mean(contribution > 0), .groups = "drop") %>%
  filter(sample_id %in% reverse_cases$sample_id, mean_contribution > 0) %>%
  left_join(all_meta %>% select(feature, report_label), by = "feature") %>%
  group_by(sample_id) %>% arrange(desc(mean_contribution), .by_group = TRUE) %>%
  mutate(contribution_rank = row_number()) %>% filter(contribution_rank <= 5) %>% ungroup()

reverse_summary <- case_contrib %>% group_by(feature, report_label) %>% summarise(
  rescued_cases_with_contribution = n_distinct(sample_id),
  mean_positive_contribution = mean(mean_contribution), .groups = "drop") %>%
  arrange(desc(rescued_cases_with_contribution), desc(mean_positive_contribution))

write_csv(reverse_cases, file.path(out_dir, "reverse_early_recurrence_cases.csv"))
write_csv(case_contrib, file.path(out_dir, "reverse_case_top_contributors.csv"))
write_csv(reverse_summary, file.path(out_dir, "reverse_early_recurrence_feature_summary.csv"))

# Publication-style figures.
nejm <- c(red = "#BC3C29", orange = "#E18727", blue = "#0072B5", purple = "#7876B1",
          gray = "#6F6F6F", light = "#E7E7E7")
theme_pub <- theme_classic(base_family = "sans", base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(color = "grey35", size = 9.5),
        axis.text = element_text(color = "black"), panel.grid.major.x = element_line(color = "grey90", linewidth = .35),
        plot.margin = margin(8, 12, 8, 8))

perf_order <- c("WHO grade", "WHO + breakpoint hotlist", "Integrated grade",
                "Integrated grade + breakpoint hotlist", "Breakpoint hotlist alone",
                "Full BiLab SEG model")
perf_plot <- performance %>% filter(model %in% perf_order) %>%
  mutate(model = factor(model, levels = rev(perf_order)),
         family = case_when(
           grepl("WHO", model) ~ "WHO",
           grepl("Integrated", model) ~ "Integrated grade",
           model == "Full BiLab SEG model" ~ "Full SEG reference",
           TRUE ~ "Hotlist alone"
         ))
p0 <- ggplot(perf_plot, aes(auc, model, color = family)) +
  geom_errorbarh(aes(xmin = auc_low, xmax = auc_high), height = .18, linewidth = .55) +
  geom_point(size = 3.2) +
  geom_text(aes(label = sprintf("%.3f", auc)), nudge_x = .006, hjust = 0, size = 3.1, color = "black") +
  scale_color_manual(values = c("WHO" = nejm[["blue"]], "Integrated grade" = nejm[["purple"]],
                                "Hotlist alone" = nejm[["orange"]], "Full SEG reference" = nejm[["red"]])) +
  scale_x_continuous(limits = c(.58, .89), breaks = seq(.60, .85, .05)) +
  labs(title = "Report-readable breakpoint hotlist versus grade",
       subtitle = "Patient-grouped nested validation; strict 5-year recurrence endpoint",
       x = "Cross-validated AUC", y = NULL, color = NULL) + theme_pub +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "figure_hotlist_vs_who_integrated_grade.png"), p0, width = 8.8, height = 5.5, dpi = 400, bg = "white")
ggsave(file.path(out_dir, "figure_hotlist_vs_who_integrated_grade.pdf"), p0, width = 8.8, height = 5.5)

top_plot <- hotlist_reduced %>% filter(is.finite(adjusted_hr)) %>% slice_head(n = 12)
if (nrow(top_plot) < 6) top_plot <- hotlist %>% filter(is.finite(adjusted_hr), adjusted_hr > 1) %>%
  distinct(report_label, .keep_all = TRUE) %>% slice_head(n = 12)
top_plot <- top_plot %>% mutate(label = stringr::str_wrap(report_label, 38),
                                label = factor(label, levels = rev(label)))
p1 <- ggplot(top_plot, aes(adjusted_hr, label, color = evidence_tier)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey60") +
  geom_errorbarh(aes(xmin = hr_low, xmax = hr_high), height = .18, linewidth = .55) +
  geom_point(size = 3.1, shape = 16) +
  scale_color_manual(values = c("Tier 1: high-signal" = nejm[["red"]], "Tier 2: supported" = nejm[["orange"]],
                                "Tier 3: exploratory" = nejm[["purple"]], "Not retained" = nejm[["gray"]])) +
  labs(title = "A  Breakpoint hotlist", subtitle = "Adjusted for WHO and prespecified integrated grade; no arm-level SCNA model",
       x = "Adjusted recurrence hazard ratio", y = NULL, color = NULL) + theme_pub

hot_features <- hotlist_reduced %>% slice_head(n = 8) %>% pull(feature)
if (length(hot_features) < 3) hot_features <- hotlist %>% filter(adjusted_hr > 1) %>% slice_head(n = 8) %>% pull(feature)
hot_mat <- sapply(hot_features, get_feature_vector)
if (is.null(dim(hot_mat))) hot_mat <- matrix(hot_mat, ncol = 1, dimnames = list(NULL, hot_features))
hot_burden <- rowSums(hot_mat > 0)
p5all <- fixed_endpoint(dat$recurrence_time, dat$recurrence_event, 5)
grad <- tibble(endpoint = p5all, burden = hot_burden) %>% filter(!is.na(endpoint)) %>%
  mutate(group = case_when(burden == 0 ~ "0 patterns", burden == 1 ~ "1 pattern", burden >= 2 ~ "2+ patterns"),
         group = factor(group, levels = c("0 patterns", "1 pattern", "2+ patterns"))) %>%
  group_by(group) %>% summarise(n = n(), events = sum(endpoint), rate = mean(endpoint),
                                se = sqrt(rate * (1 - rate) / n), .groups = "drop")
p2 <- ggplot(grad, aes(group, rate, fill = group)) + geom_col(width = .68) +
  geom_errorbar(aes(ymin = pmax(0, rate - 1.96 * se), ymax = pmin(1, rate + 1.96 * se)), width = .14) +
  geom_text(aes(label = paste0(events, "/", n)), vjust = -0.45, size = 3.2) +
  scale_fill_manual(values = c(nejm[["blue"]], nejm[["orange"]], nejm[["red"]])) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, .14))) +
  labs(title = "B  Outcome gradient", subtitle = "Strict 5-year recurrence endpoint",
       x = "Number of retained hotlist patterns", y = "Five-year recurrence", fill = NULL) +
  theme_pub + theme(legend.position = "none", panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color = "grey90"))

fig1 <- p1 | p2
ggsave(file.path(out_dir, "figure_breakpoint_report_hotlist.png"), fig1, width = 11.2, height = 6.2, dpi = 400, bg = "white")
ggsave(file.path(out_dir, "figure_breakpoint_report_hotlist.pdf"), fig1, width = 11.2, height = 6.2)

rev_plot <- reverse_summary %>% slice_head(n = 10) %>%
  mutate(label = factor(str_wrap(report_label, 42), levels = rev(str_wrap(report_label, 42))))
p3 <- ggplot(rev_plot, aes(rescued_cases_with_contribution, label)) +
  geom_segment(aes(x = 0, xend = rescued_cases_with_contribution, yend = label), color = nejm[["light"]], linewidth = 2.5) +
  geom_point(color = nejm[["red"]], size = 3) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(title = "Reverse attribution in unexpectedly early recurrences",
       subtitle = paste0("GTR, no prior/adjuvant RT, recurrence <=3 years, and baseline risk below the cohort median; n = ", nrow(reverse_cases)),
       x = "Cases with a positive out-of-fold contribution", y = NULL) + theme_pub
ggsave(file.path(out_dir, "figure_reverse_early_recurrence_attribution.png"), p3, width = 8.6, height = 5.8, dpi = 400, bg = "white")
ggsave(file.path(out_dir, "figure_reverse_early_recurrence_attribution.pdf"), p3, width = 8.6, height = 5.8)

write_lines(c(
  "Breakpoint report hotlist analysis", "",
  paste0("Primary strict 5-year cohort: n=", nrow(primary), ", events=", sum(primary$endpoint), "."),
  paste0("Nested stability selection: ", repeats, " repetitions x 5 patient-grouped outer folds."),
  "Within each outer training fold, breakpoint candidates were screened after fitting WHO and prespecified integrated grade; no arm-level SCNA features entered discovery.",
  "The 30 highest positive residual-association primitives were retained; pairwise AND patterns were generated recursively from the top 12 and filtered for prevalence.",
  "An elastic-net augmented model was tuned in patient-grouped inner folds. Held-out tumors did not participate in screening, pair generation, or tuning.",
  "Hotspots are fixed 10%-of-arm intervals formed from adjacent normalized 5% bins and mapped to UCSC hg19 cytobands.",
  "Adjusted hazard ratios include WHO and prespecified integrated grade plus one candidate with patient-clustered standard errors; they are descriptive after nested discovery.",
  "Tier 1 requires >=40% nested selection, >=90% positive coefficients, adjusted HR >=1.5, and BH q<=0.10.",
  "Tier 2 requires >=20% nested selection, >=80% positive coefficients, adjusted HR >=1.3, and BH q<=0.20.",
  "Tier 3 requires >=10% nested selection, >=70% positive coefficients, and adjusted HR >1.",
  "These are internally validated research candidates, not a clinical reporting standard."
), file.path(out_dir, "METHODS_AND_THRESHOLDS.md"))

print(performance)
print(recovery)
print(hotlist %>% filter(evidence_tier != "Not retained") %>%
        select(evidence_tier, report_label, selection_fraction, adjusted_hr, hr_low, hr_high, q_value,
               recurrence_5y_positive, recurrence_5y_negative) %>% head(20))
print(reverse_summary %>% head(15))
