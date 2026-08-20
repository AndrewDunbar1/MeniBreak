#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(glmnet)
  library(ggplot2)
  library(patchwork)
})

root <- Sys.getenv("MENIBREAK_PROJECT_ROOT")
if (!nzchar(root)) stop("Set MENIBREAK_PROJECT_ROOT to the analysis project root; see config/path_config.example.sh")
out_dir <- file.path(root, "analysis_folders", "TRISCUT_Inspired_Breakpoint_Hotlist_2026")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dataset_path <- file.path(root, "analysis_folders", "ml_breakpoint_survival_stratification_2026",
                          "ml_breakpoint_model_dataset_2026.csv")
outcome_path <- file.path(root, "software", "mbpr", "analysis", "bakeoff_2x2", "outcomes.csv")
seg_path <- Sys.getenv("MENIBREAK_SEGMENTS")
if (!nzchar(seg_path)) stop("Set MENIBREAK_SEGMENTS to the event-called SEG file")
arm_path <- file.path(root, "software", "mbpr", "inst", "extdata", "arm_reference_GRCh37.csv")
old_path <- file.path(root, "analysis_folders", "Breakpoint_Report_Hotlist_2026", "compact_six_rule",
                      "compact_nested_oof_predictions.csv")
old_selected_path <- file.path(root, "analysis_folders", "Breakpoint_Report_Hotlist_2026", "compact_six_rule",
                               "compact_nested_selected_rules.csv")
stopifnot(all(file.exists(c(dataset_path, outcome_path, seg_path, arm_path, old_path, old_selected_path))))

repeats <- 10L
outer_k <- 5L
hotlist_size <- 6L
bootstrap_draws <- 2000L
seed0 <- 20260803L

auc_safe <- function(y, s) {
  ok <- is.finite(y) & is.finite(s); y <- y[ok]; s <- s[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(s, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

ap_safe <- function(y, s) {
  ok <- is.finite(y) & is.finite(s); y <- y[ok]; s <- s[ok]
  if (!sum(y == 1)) return(NA_real_)
  ord <- order(s, decreasing = TRUE)
  yy <- y[ord]
  mean(cumsum(yy)[yy == 1] / which(yy == 1))
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

fit_cv <- function(xtr, ytr, groups, seed) {
  keep <- apply(xtr, 2, function(z) length(unique(z[is.finite(z)])) > 1L)
  xtr <- xtr[, keep, drop = FALSE]
  foldid <- group_folds(groups, min(5L, length(unique(groups))), seed)
  fit <- cv.glmnet(xtr, ytr, family = "binomial", alpha = 0.25,
                   foldid = foldid, type.measure = "deviance",
                   standardize = TRUE)
  list(fit = fit, columns = colnames(xtr))
}

predict_cv <- function(obj, x, type = "link") {
  drop(predict(obj$fit, x[, obj$columns, drop = FALSE],
               s = "lambda.1se", type = type))
}

coef_cv <- function(obj) {
  z <- as.matrix(coef(obj$fit, s = "lambda.1se"))[, 1]
  z[names(z) != "(Intercept)" & z != 0]
}

safe_glm_residual <- function(y, grade, selected = NULL) {
  x <- cbind(WHO2 = as.numeric(grade$who_grade == 2),
             WHO3 = as.numeric(grade$who_grade == 3),
             IG2 = as.numeric(grade$integrated_grade == 2),
             IG3 = as.numeric(grade$integrated_grade == 3))
  if (!is.null(selected) && ncol(selected)) x <- cbind(x, selected)
  fit <- try(glm.fit(cbind(1, x), y, family = binomial()), silent = TRUE)
  if (inherits(fit, "try-error") || any(!is.finite(fit$fitted.values))) {
    return(y - mean(y))
  }
  y - pmin(pmax(fit$fitted.values, 1e-5), 1 - 1e-5)
}

residualize_matrix <- function(x, z) {
  if (!ncol(x)) return(x)
  z <- as.matrix(z)
  z[!is.finite(z)] <- 0
  q <- qr(z)
  qr.resid(q, x)
}

message("Loading deidentified modeling table")
raw <- read_csv(dataset_path, show_col_types = FALSE)
outcomes <- read_csv(outcome_path, show_col_types = FALSE)
dat <- raw %>% transmute(
  sample_id = str_to_upper(str_trim(as.character(cg))), patient_id,
  who_grade = as.numeric(who_grade), integrated_grade = as.numeric(integrated_grade),
  total_breakpoint_boundaries = as.numeric(total_breakpoint_boundaries),
  breakpoint_arms_touched = as.numeric(breakpoint_arms_touched),
  breakpoint_route_entropy = as.numeric(breakpoint_route_entropy),
  max_arm_breakpoint_fraction = as.numeric(max_arm_breakpoint_fraction),
  small_segment_count_lt5mb = as.numeric(small_segment_count_lt5mb),
  loss_segment_count = as.numeric(loss_segment_count),
  gain_segment_count = as.numeric(gain_segment_count)
) %>% inner_join(outcomes, by = c("sample_id", "patient_id")) %>%
  mutate(endpoint = fixed_endpoint(recurrence_time, recurrence_event, 5)) %>%
  filter(!is.na(endpoint), !is.na(who_grade), !is.na(integrated_grade))

message("Building deterministic event-call geometry candidates")
allowed_arms <- c(
  "1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q",
  "8p","8q","9p","9q","10p","10q","11p","11q","12p","12q","13q","14q","15q",
  "16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21q","22q"
)
arm <- fread(arm_path) %>% as_tibble() %>% transmute(
  chromosome = sub("^chr", "", chromosome), arm = tolower(arm),
  arm_key = paste0(chromosome, arm), arm_start = as.numeric(start),
  arm_end = as.numeric(end), arm_length = as.numeric(arm_length)
) %>% filter(arm_key %in% allowed_arms)

seg <- fread(seg_path, select = c("cg", "Chr", "Start", "End", "Event"))
seg[, `:=`(
  cg = toupper(trimws(as.character(cg))),
  chromosome = sub("^chr", "", as.character(Chr)),
  Start = as.numeric(Start), End = as.numeric(End),
  event = fcase(
    Event %chin% c("CN Loss", "Homozygous Copy Loss"), "loss",
    Event %chin% c("CN Gain", "High Copy Gain"), "gain",
    Event == "LOH", "loh",
    Event == "Allelic Imbalance", "ai",
    default = NA_character_
  )
)]
seg <- seg[!is.na(event) & cg %chin% dat$sample_id & is.finite(Start) & is.finite(End) & End > Start]
seg[, segment_id := .I]

# CNA is the union of loss and gain. LOH and allelic imbalance remain distinct;
# they are not treated as absolute or phased copy-number states.
seg_specific <- seg[, .(cg, chromosome, Start, End, event, segment_id)]
seg_cna <- seg[event %chin% c("loss", "gain"),
               .(cg, chromosome, Start, End, event = "cna", segment_id)]
seg_ev <- rbindlist(list(seg_specific, seg_cna), use.names = TRUE)

arm_dt <- as.data.table(arm)
bound <- melt(seg_ev, id.vars = c("cg", "chromosome", "event", "segment_id"),
              measure.vars = c("Start", "End"), variable.name = "boundary", value.name = "pos")
bound <- merge(bound, arm_dt, by = "chromosome", allow.cartesian = TRUE)
bound <- bound[pos >= arm_start & pos <= arm_end]
bound <- bound[!((arm == "p" & pos == arm_start) | (arm == "q" & pos == arm_end))]
bound[, rel := fifelse(arm == "p", (pos - arm_start) / arm_length,
                       (arm_end - pos) / arm_length)]
bound[, rel := pmin(pmax(rel, 0), .999999)]
bound[, region10 := pmin(10L, floor(rel * 10) + 1L)]
bound[, hot_feature := sprintf("hot__%s__%s__r%02d", event, arm_key, region10)]

hot_counts <- dcast(bound, cg ~ hot_feature, fun.aggregate = length, value.var = "pos", fill = 0)
hot_names <- setdiff(names(hot_counts), "cg")
hot_binary <- copy(hot_counts)
for (nm in hot_names) set(hot_binary, j = nm, value = as.integer(hot_binary[[nm]] > 0))
cluster_binary <- hot_counts[, c(list(cg = cg), lapply(.SD, function(z) as.integer(z >= 2L))), .SDcols = hot_names]
setnames(cluster_binary, hot_names, sub("^hot__", "cluster__", hot_names))

# Complete same-arm segments yield simple route features between coarse 20%-of-arm
# locations. These represent terminal, internal, and centromere-anchored routes.
route <- merge(seg_ev, arm_dt, by = "chromosome", allow.cartesian = TRUE)
route <- route[Start >= arm_start & End <= arm_end]
route[, rel_start := fifelse(arm == "p", (Start - arm_start) / arm_length,
                             (arm_end - Start) / arm_length)]
route[, rel_end := fifelse(arm == "p", (End - arm_start) / arm_length,
                           (arm_end - End) / arm_length)]
route[, `:=`(rel_lo = pmin(rel_start, rel_end), rel_hi = pmax(rel_start, rel_end))]
route[, `:=`(rlo = pmin(5L, floor(pmin(pmax(rel_lo, 0), .999999) * 5) + 1L),
             rhi = pmin(5L, floor(pmin(pmax(rel_hi, 0), .999999) * 5) + 1L))]
route[, route_feature := sprintf("route__%s__%s__r%02d_r%02d", event, arm_key, rlo, rhi)]
route_wide <- dcast(unique(route[, .(cg, route_feature)]), cg ~ route_feature,
                    fun.aggregate = length, value.var = "route_feature", fill = 0)
route_names <- setdiff(names(route_wide), "cg")
for (nm in route_names) set(route_wide, j = nm, value = as.integer(route_wide[[nm]] > 0))

route[, focal_feature := sprintf("focal__%s__%s", event, arm_key)]
focal <- route[(End - Start) < 5e6]
focal_wide <- dcast(unique(focal[, .(cg, focal_feature)]), cg ~ focal_feature,
                    fun.aggregate = length, value.var = "focal_feature", fill = 0)
focal_names <- setdiff(names(focal_wide), "cg")
for (nm in focal_names) set(focal_wide, j = nm, value = as.integer(focal_wide[[nm]] > 0))

route[, internal_feature := sprintf("internal__%s__%s", event, arm_key)]
internal <- route[rel_lo > .10 & rel_hi < .90]
internal_wide <- dcast(unique(internal[, .(cg, internal_feature)]), cg ~ internal_feature,
                       fun.aggregate = length, value.var = "internal_feature", fill = 0)
internal_names <- setdiff(names(internal_wide), "cg")
for (nm in internal_names) set(internal_wide, j = nm, value = as.integer(internal_wide[[nm]] > 0))

candidate_dt <- Reduce(function(x, y) merge(x, y, by = "cg", all = TRUE),
                       list(hot_binary, cluster_binary, route_wide, focal_wide, internal_wide))
candidate_dt[is.na(candidate_dt)] <- 0

# Fixed report-readable architecture rules. Thresholds match or refine the
# previously audited hotlist; no outcome-derived cut point is created.
arch <- dat %>% transmute(
  cg = sample_id,
  arch__total_breakpoints__ge20 = as.integer(total_breakpoint_boundaries >= 20),
  arch__total_breakpoints__ge40 = as.integer(total_breakpoint_boundaries >= 40),
  arch__arms_touched__ge5 = as.integer(breakpoint_arms_touched >= 5),
  arch__arms_touched__ge7 = as.integer(breakpoint_arms_touched >= 7),
  arch__focal_segments__ge6 = as.integer(small_segment_count_lt5mb >= 6),
  arch__focal_segments__ge10 = as.integer(small_segment_count_lt5mb >= 10),
  arch__route_entropy__ge2 = as.integer(breakpoint_route_entropy >= 2),
  arch__arm_concentration__ge50pct = as.integer(max_arm_breakpoint_fraction >= .50)
)
candidate_dt <- merge(candidate_dt, as.data.table(arch), by = "cg", all.y = TRUE)
candidate_dt[is.na(candidate_dt)] <- 0
candidate_dt <- candidate_dt[match(dat$sample_id, cg)]
stopifnot(identical(candidate_dt$cg, dat$sample_id))
X <- as.matrix(candidate_dt[, setdiff(names(candidate_dt), "cg"), with = FALSE])
storage.mode(X) <- "double"

parse_meta <- function(nm) {
  p <- strsplit(nm, "__", fixed = TRUE)[[1]]
  type <- p[[1]]
  event_label <- function(x) c(cna = "CNA", loss = "loss", gain = "gain",
                               loh = "LOH", ai = "allelic-imbalance")[[x]]
  if (type %in% c("hot", "cluster")) {
    event <- p[[2]]; arm_key <- p[[3]]; loc <- p[[4]]
    group <- paste("region", arm_key, loc, sep = "__")
    region <- as.integer(sub("^r", "", loc))
    label <- sprintf("%s %d–%d%% %s breakpoint%s", arm_key, (region - 1L) * 10L, region * 10L,
                     event_label(event), ifelse(type == "cluster", " cluster (≥2)", ""))
  } else if (type == "route") {
    event <- p[[2]]; arm_key <- p[[3]]; loc <- p[[4]]
    group <- paste("route", arm_key, loc, sep = "__")
    rr <- as.integer(sub("^r", "", strsplit(loc, "_", fixed = TRUE)[[1]]))
    label <- sprintf("%s %s route: %d–%d%% ↔ %d–%d%%", arm_key, event_label(event),
                     (rr[[1]] - 1L) * 20L, rr[[1]] * 20L,
                     (rr[[2]] - 1L) * 20L, rr[[2]] * 20L)
  } else if (type %in% c("focal", "internal")) {
    event <- p[[2]]; arm_key <- p[[3]]
    group <- paste(type, arm_key, sep = "__")
    label <- sprintf("%s %s %s segment", arm_key, event_label(event),
                     ifelse(type == "focal", "focal (<5 Mb)", "internal"))
  } else {
    event <- "architecture"; arm_key <- NA_character_
    group <- paste(p[1:min(2, length(p))], collapse = "__")
    label <- c(
      arch__total_breakpoints__ge20 = "≥20 total CNA breakpoints",
      arch__total_breakpoints__ge40 = "≥40 total CNA breakpoints",
      arch__arms_touched__ge5 = "Breakpoints on ≥5 chromosome arms",
      arch__arms_touched__ge7 = "Breakpoints on ≥7 chromosome arms",
      arch__focal_segments__ge6 = "≥6 focal (<5 Mb) CNA segments",
      arch__focal_segments__ge10 = "≥10 focal (<5 Mb) CNA segments",
      arch__route_entropy__ge2 = "Breakpoint-route entropy ≥2",
      arch__arm_concentration__ge50pct = "≥50% of breakpoints on one arm"
    )[[nm]]
  }
  tibble(feature = nm, feature_type = type, event = event, arm_key = arm_key,
         rule_group = group, report_label = label)
}
meta <- bind_rows(lapply(colnames(X), parse_meta))

# Opportunity variables used only for within-training background correction.
total_ctx <- bound[, .(n_boundary = .N), by = .(cg, event)]
total_wide <- dcast(total_ctx, cg ~ event, value.var = "n_boundary", fill = 0)
setnames(total_wide, setdiff(names(total_wide), "cg"),
         paste0("ctx_total_", setdiff(names(total_wide), "cg")))
arm_ctx <- bound[, .(n_boundary = .N), by = .(cg, event, arm_key)]
arm_ctx[, key := paste0("ctx_arm_", event, "_", arm_key)]
arm_wide <- dcast(arm_ctx, cg ~ key, value.var = "n_boundary", fill = 0)
ctx <- merge(total_wide, arm_wide, by = "cg", all = TRUE)
ctx[is.na(ctx)] <- 0
ctx <- merge(data.table(cg = dat$sample_id), ctx, by = "cg", all.x = TRUE, sort = FALSE)
ctx[is.na(ctx)] <- 0
ctx <- ctx[match(dat$sample_id, cg)]

background_residuals <- function(xtr, metadata, context_tr, use_background = TRUE) {
  if (!use_background) return(scale(xtr, center = TRUE, scale = FALSE))
  out <- matrix(0, nrow(xtr), ncol(xtr), dimnames = dimnames(xtr))
  groups <- split(seq_len(nrow(metadata)), paste(metadata$event, metadata$arm_key, sep = "||"))
  global <- if ("ctx_total_cna" %in% names(context_tr)) context_tr$ctx_total_cna else rep(0, nrow(context_tr))
  for (idx in groups) {
    ev <- metadata$event[idx[[1]]]
    ar <- metadata$arm_key[idx[[1]]]
    if (ev == "architecture") {
      out[, idx] <- scale(xtr[, idx, drop = FALSE], center = TRUE, scale = FALSE)
      next
    }
    total_nm <- paste0("ctx_total_", ev)
    arm_nm <- paste0("ctx_arm_", ev, "_", ar)
    ev_total <- if (total_nm %in% names(context_tr)) context_tr[[total_nm]] else rep(0, nrow(context_tr))
    arm_total <- if (arm_nm %in% names(context_tr)) context_tr[[arm_nm]] else rep(0, nrow(context_tr))
    z <- cbind(1, log1p(global), log1p(ev_total), log1p(arm_total), as.numeric(arm_total > 0))
    out[, idx] <- residualize_matrix(xtr[, idx, drop = FALSE], z)
  }
  out
}

feature_atoms <- function(nm, metadata) {
  if (startsWith(nm, "pair__")) return(strsplit(sub("^pair__", "", nm), "__AND__", fixed = TRUE)[[1]])
  nm
}

feature_groups <- function(nm, metadata) {
  atoms <- feature_atoms(nm, metadata)
  metadata$rule_group[match(atoms, metadata$feature)]
}

select_rules <- function(train_idx, event_set, use_background, k, variant,
                         exclude_arm_spanning_routes = FALSE) {
  md <- meta
  permitted <- md$event %in% c(event_set, "architecture")
  if (exclude_arm_spanning_routes) {
    permitted <- permitted & !grepl("^route__.*__r01_r05$", md$feature)
  }
  xc <- X[train_idx, permitted, drop = FALSE]
  md <- md[permitted, , drop = FALSE]
  y <- dat$endpoint[train_idx]
  prev <- colSums(xc > 0)
  case_n <- colSums((xc > 0) * y)
  ctrl_n <- colSums((xc > 0) * (1 - y))
  eligible <- prev >= 15 & prev <= length(train_idx) - 15 & case_n >= 5 & ctrl_n >= 5
  xc <- xc[, eligible, drop = FALSE]
  md <- md[eligible, , drop = FALSE]
  bg <- background_residuals(xc, md, as.data.frame(ctx[train_idx]), use_background)
  yres0 <- safe_glm_residual(y, dat[train_idx, c("who_grade", "integrated_grade")])
  denom <- sqrt(colSums(bg^2) * sum(yres0^2))
  score0 <- drop(crossprod(bg, yres0)) / denom
  names(score0) <- colnames(bg)
  score0[!is.finite(score0)] <- -Inf

  # Candidate cross-arm interactions are generated only from the strongest
  # training-fold spatial primitives. Their score tests synergy beyond both
  # constituent rules and global CNA-breakpoint burden.
  spatial <- md$feature[md$feature_type %in% c("hot", "cluster", "route", "focal", "internal")]
  top_spatial <- head(names(sort(score0[names(score0) %in% spatial & score0 > 0], decreasing = TRUE)), 18)
  pair_x <- list(); pair_meta <- list(); pair_bg <- list()
  if (length(top_spatial) >= 2) {
    for (pp in combn(top_spatial, 2, simplify = FALSE)) {
      m1 <- md[match(pp[[1]], md$feature), ]; m2 <- md[match(pp[[2]], md$feature), ]
      if (is.na(m1$arm_key) || is.na(m2$arm_key) || m1$arm_key == m2$arm_key) next
      zall <- as.integer(X[, pp[[1]]] > 0 & X[, pp[[2]]] > 0)
      z <- zall[train_idx]
      if (sum(z) < 15 || sum(z == 0) < 15 || sum(z * y) < 5 || sum(z * (1 - y)) < 5) next
      nm <- paste0("pair__", pp[[1]], "__AND__", pp[[2]])
      global <- if ("ctx_total_cna" %in% names(ctx)) ctx$ctx_total_cna[train_idx] else rep(0, length(train_idx))
      rz <- if (use_background) residualize_matrix(matrix(z, ncol = 1),
              cbind(1, X[train_idx, pp[[1]]], X[train_idx, pp[[2]]], log1p(global)))[, 1] else z - mean(z)
      pair_x[[nm]] <- zall
      pair_bg[[nm]] <- rz
      pair_meta[[nm]] <- tibble(feature = nm, feature_type = "pair", event = "pair",
                                arm_key = NA_character_, rule_group = nm,
                                report_label = paste(m1$report_label, "+", m2$report_label))
    }
  }
  if (length(pair_x)) {
    px <- do.call(cbind, pair_x); pbg <- do.call(cbind, pair_bg)
    pmd <- bind_rows(pair_meta)
    xc_all <- cbind(X[, colnames(xc), drop = FALSE], px)
    bg_all <- cbind(bg, pbg)
    md_all <- bind_rows(md, pmd)
  } else {
    xc_all <- X[, colnames(xc), drop = FALSE]
    bg_all <- bg
    md_all <- md
  }

  selected <- character(); selected_groups <- character(); scores_chosen <- numeric()
  yres <- yres0
  remaining <- colnames(bg_all)
  for (step in seq_len(k)) {
    scores <- drop(crossprod(bg_all[, remaining, drop = FALSE], yres)) /
      sqrt(colSums(bg_all[, remaining, drop = FALSE]^2) * sum(yres^2))
    names(scores) <- remaining
    scores[!is.finite(scores)] <- -Inf
    ordered <- names(sort(scores, decreasing = TRUE))
    pick <- NA_character_
    for (nm in ordered) {
      if (!is.finite(scores[[nm]]) || scores[[nm]] <= 0) next
      atoms <- feature_atoms(nm, md_all)
      groups <- if (startsWith(nm, "pair__")) {
        md_all$rule_group[match(atoms, md_all$feature)]
      } else md_all$rule_group[match(nm, md_all$feature)]
      if (any(groups %in% selected_groups)) next
      if (length(selected)) {
        zz <- xc_all[train_idx, c(selected, nm), drop = FALSE]
        cors <- suppressWarnings(cor(zz[, nm], zz[, selected, drop = FALSE]))
        if (any(abs(cors) > .85, na.rm = TRUE)) next
      }
      pick <- nm; selected_groups <- c(selected_groups, groups); break
    }
    if (is.na(pick)) break
    selected <- c(selected, pick); scores_chosen <- c(scores_chosen, scores[[pick]])
    remaining <- setdiff(remaining, pick)
    yres <- safe_glm_residual(
      y, dat[train_idx, c("who_grade", "integrated_grade")],
      xc_all[train_idx, selected, drop = FALSE]
    )
  }
  list(features = selected, scores = scores_chosen, x = xc_all,
       meta = md_all, variant = variant)
}

variants <- tibble(
  variant = c("triscut_full_background", "triscut_cna_background",
              "triscut_full_no_background", "triscut_cna_no_background",
              "triscut_cna_no_background_no_armspan",
              "triscut_full_no_background_no_armspan",
              "triscut_full_no_background_k4", "triscut_full_no_background_k8"),
  event_set = list(c("cna", "loss", "gain", "loh", "ai"),
                   c("cna", "loss", "gain"),
                   c("cna", "loss", "gain", "loh", "ai"),
                   c("cna", "loss", "gain"),
                   c("cna", "loss", "gain"),
                   c("cna", "loss", "gain", "loh", "ai"),
                   c("cna", "loss", "gain", "loh", "ai"),
                   c("cna", "loss", "gain", "loh", "ai")),
  use_background = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  rule_count = c(6L, 6L, 6L, 6L, 6L, 6L, 4L, 8L),
  exclude_arm_spanning_routes = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE)
)

old <- read_csv(old_path, show_col_types = FALSE) %>% transmute(
  sample_id = str_to_upper(sample_id), patient_id, repetition,
  old_who_lp = who_hotlist_lp, old_ig_lp = ig_hotlist_lp,
  full_seg_lp
)

pred_rows <- list(); selected_rows <- list()
for (r in seq_len(repeats)) {
  message("Nested repetition ", r, "/", repeats)
  outer <- group_folds(dat$patient_id, outer_k, seed0 + r)
  for (fold in seq_len(outer_k)) {
    tr <- which(outer != fold); te <- which(outer == fold)
    who_base <- glm(endpoint ~ factor(who_grade), data = dat[tr, ], family = binomial())
    ig_base <- glm(endpoint ~ factor(integrated_grade), data = dat[tr, ], family = binomial())
    fold_pred <- tibble(
      sample_id = dat$sample_id[te], patient_id = dat$patient_id[te],
      repetition = r, fold = fold, endpoint = dat$endpoint[te],
      who_grade = dat$who_grade[te], integrated_grade = dat$integrated_grade[te],
      who_base_lp = drop(predict(who_base, dat[te, ], type = "link")),
      ig_base_lp = drop(predict(ig_base, dat[te, ], type = "link"))
    )
    for (v in seq_len(nrow(variants))) {
      variant <- variants$variant[[v]]
      sel <- select_rules(tr, variants$event_set[[v]], variants$use_background[[v]],
                          variants$rule_count[[v]], variant,
                          variants$exclude_arm_spanning_routes[[v]])
      chosen <- sel$features
      if (length(chosen) < variants$rule_count[[v]]) warning(variant, " selected only ", length(chosen), " rules")
      xhot <- sel$x[, chosen, drop = FALSE]
      Xwho <- cbind(WHO2 = as.numeric(dat$who_grade == 2),
                    WHO3 = as.numeric(dat$who_grade == 3), xhot)
      Xig <- cbind(IG2 = as.numeric(dat$integrated_grade == 2),
                   IG3 = as.numeric(dat$integrated_grade == 3), xhot)
      who_fit <- fit_cv(Xwho[tr, , drop = FALSE], dat$endpoint[tr], dat$patient_id[tr],
                        seed0 + 100000L * v + r * 100L + fold)
      ig_fit <- fit_cv(Xig[tr, , drop = FALSE], dat$endpoint[tr], dat$patient_id[tr],
                       seed0 + 200000L * v + r * 100L + fold)
      fold_pred[[paste0(variant, "__who_lp")]] <- predict_cv(who_fit, Xwho[te, , drop = FALSE])
      fold_pred[[paste0(variant, "__ig_lp")]] <- predict_cv(ig_fit, Xig[te, , drop = FALSE])
      cw <- coef_cv(who_fit); ci <- coef_cv(ig_fit)
      sm <- sel$meta[match(chosen, sel$meta$feature), ]
      selected_rows[[length(selected_rows) + 1L]] <- tibble(
        repetition = r, fold = fold, variant = variant,
        rank = seq_along(chosen), feature = chosen,
        screen_score = sel$scores,
        report_label = sm$report_label,
        feature_type = sm$feature_type,
        event = sm$event,
        who_coefficient = unname(cw[chosen]),
        ig_coefficient = unname(ci[chosen])
      )
    }
    pred_rows[[length(pred_rows) + 1L]] <- fold_pred
  }
}

pred <- bind_rows(pred_rows) %>% left_join(old, by = c("sample_id", "patient_id", "repetition"))
selected <- bind_rows(selected_rows)
write_csv(pred, file.path(out_dir, "nested_oof_predictions.csv"))
write_csv(selected, file.path(out_dir, "nested_selected_rules.csv"))

model_spec <- tribble(
  ~model, ~score_col, ~prob_lp_col,
  "WHO grade", "who_grade", "who_base_lp",
  "WHO + current six-rule hotlist", "old_who_lp", "old_who_lp",
  "WHO + TRISCUT-inspired hotlist", "triscut_full_background__who_lp", "triscut_full_background__who_lp",
  "WHO + TRISCUT-inspired CNA-only ablation", "triscut_cna_background__who_lp", "triscut_cna_background__who_lp",
  "WHO + TRISCUT-inspired no-background ablation", "triscut_full_no_background__who_lp", "triscut_full_no_background__who_lp",
  "WHO + iterative CNA-only geometry hotlist", "triscut_cna_no_background__who_lp", "triscut_cna_no_background__who_lp",
  "WHO + CNA geometry hotlist excluding arm-spanning routes", "triscut_cna_no_background_no_armspan__who_lp", "triscut_cna_no_background_no_armspan__who_lp",
  "WHO + all-event geometry hotlist excluding arm-spanning routes", "triscut_full_no_background_no_armspan__who_lp", "triscut_full_no_background_no_armspan__who_lp",
  "WHO + iterative geometry hotlist (4 rules)", "triscut_full_no_background_k4__who_lp", "triscut_full_no_background_k4__who_lp",
  "WHO + iterative geometry hotlist (8 rules)", "triscut_full_no_background_k8__who_lp", "triscut_full_no_background_k8__who_lp",
  "Integrated grade", "integrated_grade", "ig_base_lp",
  "Integrated grade + current six-rule hotlist", "old_ig_lp", "old_ig_lp",
  "Integrated grade + TRISCUT-inspired hotlist", "triscut_full_background__ig_lp", "triscut_full_background__ig_lp",
  "Integrated grade + TRISCUT-inspired CNA-only ablation", "triscut_cna_background__ig_lp", "triscut_cna_background__ig_lp",
  "Integrated grade + TRISCUT-inspired no-background ablation", "triscut_full_no_background__ig_lp", "triscut_full_no_background__ig_lp",
  "Integrated grade + iterative CNA-only geometry hotlist", "triscut_cna_no_background__ig_lp", "triscut_cna_no_background__ig_lp",
  "Integrated grade + CNA geometry hotlist excluding arm-spanning routes", "triscut_cna_no_background_no_armspan__ig_lp", "triscut_cna_no_background_no_armspan__ig_lp",
  "Integrated grade + all-event geometry hotlist excluding arm-spanning routes", "triscut_full_no_background_no_armspan__ig_lp", "triscut_full_no_background_no_armspan__ig_lp",
  "Integrated grade + iterative geometry hotlist (4 rules)", "triscut_full_no_background_k4__ig_lp", "triscut_full_no_background_k4__ig_lp",
  "Integrated grade + iterative geometry hotlist (8 rules)", "triscut_full_no_background_k8__ig_lp", "triscut_full_no_background_k8__ig_lp",
  "Full clinical + SCNA + breakpoint model", "full_seg_lp", "full_seg_lp"
)

per_rep <- lapply(seq_len(repeats), function(r) {
  d <- pred[pred$repetition == r, ]
  bind_rows(lapply(seq_len(nrow(model_spec)), function(j) {
    m <- model_spec$model[[j]]
    s <- d[[model_spec$score_col[[j]]]]
    prob <- plogis(d[[model_spec$prob_lp_col[[j]]]])
    tibble(repetition = r, model = m,
           auc = auc_safe(d$endpoint, s),
           average_precision = ap_safe(d$endpoint, s),
           brier = mean((prob - d$endpoint)^2),
           top20_event_capture = {
             top <- order(s, decreasing = TRUE)[seq_len(ceiling(.20 * nrow(d)))]
             sum(d$endpoint[top]) / sum(d$endpoint)
           },
           top20_precision = {
             top <- order(s, decreasing = TRUE)[seq_len(ceiling(.20 * nrow(d)))]
             mean(d$endpoint[top])
           })
  }))
}) %>% bind_rows()
write_csv(per_rep, file.path(out_dir, "performance_by_repetition.csv"))

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
    setNames(vapply(model_spec$score_col,
                    function(v) auc_safe(d$endpoint[ii], d[[v]][ii]), numeric(1)),
             model_spec$model)
  })
  a <- rowMeans(do.call(cbind, vals), na.rm = TRUE)
  boot[[b]] <- tibble(draw = b, model = names(a), auc = unname(a))
}
boot <- bind_rows(boot)
write_csv(boot, file.path(out_dir, "patient_bootstrap_auc.csv"))

summary <- per_rep %>% group_by(model) %>% summarise(
  auc_sd_repetition = sd(auc), auc = mean(auc),
  average_precision = mean(average_precision), brier = mean(brier),
  top20_event_capture = mean(top20_event_capture), top20_precision = mean(top20_precision),
  .groups = "drop"
) %>% left_join(boot %>% group_by(model) %>% summarise(
  auc_low = quantile(auc, .025), auc_high = quantile(auc, .975), .groups = "drop"
), by = "model")
write_csv(summary, file.path(out_dir, "performance_summary.csv"))

boot_wide <- boot %>% pivot_wider(names_from = model, values_from = auc)
comparisons <- tribble(
  ~comparison, ~new_model, ~reference_model,
  "WHO: TRISCUT-inspired versus WHO alone", "WHO + TRISCUT-inspired hotlist", "WHO grade",
  "WHO: TRISCUT-inspired versus current hotlist", "WHO + TRISCUT-inspired hotlist", "WHO + current six-rule hotlist",
  "IG: TRISCUT-inspired versus IG alone", "Integrated grade + TRISCUT-inspired hotlist", "Integrated grade",
  "IG: TRISCUT-inspired versus current hotlist", "Integrated grade + TRISCUT-inspired hotlist", "Integrated grade + current six-rule hotlist",
  "WHO: all-event versus CNA-only", "WHO + TRISCUT-inspired hotlist", "WHO + TRISCUT-inspired CNA-only ablation",
  "IG: all-event versus CNA-only", "Integrated grade + TRISCUT-inspired hotlist", "Integrated grade + TRISCUT-inspired CNA-only ablation",
  "WHO: background correction versus none", "WHO + TRISCUT-inspired hotlist", "WHO + TRISCUT-inspired no-background ablation",
  "IG: background correction versus none", "Integrated grade + TRISCUT-inspired hotlist", "Integrated grade + TRISCUT-inspired no-background ablation"
  ,"WHO: iterative geometry versus current hotlist", "WHO + TRISCUT-inspired no-background ablation", "WHO + current six-rule hotlist"
  ,"IG: iterative geometry versus current hotlist", "Integrated grade + TRISCUT-inspired no-background ablation", "Integrated grade + current six-rule hotlist"
  ,"WHO: iterative CNA-only geometry versus current hotlist", "WHO + iterative CNA-only geometry hotlist", "WHO + current six-rule hotlist"
  ,"IG: iterative CNA-only geometry versus current hotlist", "Integrated grade + iterative CNA-only geometry hotlist", "Integrated grade + current six-rule hotlist"
  ,"WHO: exclude arm-spanning routes versus CNA geometry", "WHO + CNA geometry hotlist excluding arm-spanning routes", "WHO + iterative CNA-only geometry hotlist"
  ,"IG: exclude arm-spanning routes versus CNA geometry", "Integrated grade + CNA geometry hotlist excluding arm-spanning routes", "Integrated grade + iterative CNA-only geometry hotlist"
  ,"WHO: arm-span-excluded CNA geometry versus current hotlist", "WHO + CNA geometry hotlist excluding arm-spanning routes", "WHO + current six-rule hotlist"
  ,"IG: arm-span-excluded CNA geometry versus current hotlist", "Integrated grade + CNA geometry hotlist excluding arm-spanning routes", "Integrated grade + current six-rule hotlist"
  ,"WHO: exclude arm-spanning routes versus all-event geometry", "WHO + all-event geometry hotlist excluding arm-spanning routes", "WHO + TRISCUT-inspired no-background ablation"
  ,"IG: exclude arm-spanning routes versus all-event geometry", "Integrated grade + all-event geometry hotlist excluding arm-spanning routes", "Integrated grade + TRISCUT-inspired no-background ablation"
  ,"WHO: arm-span-excluded all-event geometry versus current hotlist", "WHO + all-event geometry hotlist excluding arm-spanning routes", "WHO + current six-rule hotlist"
  ,"IG: arm-span-excluded all-event geometry versus current hotlist", "Integrated grade + all-event geometry hotlist excluding arm-spanning routes", "Integrated grade + current six-rule hotlist"
  ,"WHO: all-event versus CNA-only without background correction", "WHO + TRISCUT-inspired no-background ablation", "WHO + iterative CNA-only geometry hotlist"
  ,"IG: all-event versus CNA-only without background correction", "Integrated grade + TRISCUT-inspired no-background ablation", "Integrated grade + iterative CNA-only geometry hotlist"
  ,"WHO: 4 versus 6 iterative geometry rules", "WHO + iterative geometry hotlist (4 rules)", "WHO + TRISCUT-inspired no-background ablation"
  ,"WHO: 8 versus 6 iterative geometry rules", "WHO + iterative geometry hotlist (8 rules)", "WHO + TRISCUT-inspired no-background ablation"
  ,"IG: 4 versus 6 iterative geometry rules", "Integrated grade + iterative geometry hotlist (4 rules)", "Integrated grade + TRISCUT-inspired no-background ablation"
  ,"IG: 8 versus 6 iterative geometry rules", "Integrated grade + iterative geometry hotlist (8 rules)", "Integrated grade + TRISCUT-inspired no-background ablation"
)
comparisons <- comparisons %>% rowwise() %>% mutate(
  delta_auc = summary$auc[match(new_model, summary$model)] - summary$auc[match(reference_model, summary$model)],
  delta_low = quantile(boot_wide[[new_model]] - boot_wide[[reference_model]], .025),
  delta_high = quantile(boot_wide[[new_model]] - boot_wide[[reference_model]], .975),
  bootstrap_two_sided_p = 2 * min(mean((boot_wide[[new_model]] - boot_wide[[reference_model]]) <= 0),
                                  mean((boot_wide[[new_model]] - boot_wide[[reference_model]]) >= 0)),
  positive_repetitions = sum(
    per_rep$auc[match(paste(seq_len(repeats), new_model), paste(per_rep$repetition, per_rep$model))] >
      per_rep$auc[match(paste(seq_len(repeats), reference_model), paste(per_rep$repetition, per_rep$model))]
  )
) %>% ungroup()
write_csv(comparisons, file.path(out_dir, "paired_auc_comparisons.csv"))

stability <- selected %>% group_by(variant, feature, report_label, feature_type, event) %>% summarise(
  folds_selected = n(), selection_fraction = n() / (repeats * outer_k),
  median_rank = median(rank), median_screen_score = median(screen_score),
  who_positive_fraction = mean(who_coefficient > 0, na.rm = TRUE),
  ig_positive_fraction = mean(ig_coefficient > 0, na.rm = TRUE),
  .groups = "drop"
) %>% arrange(variant, desc(selection_fraction), median_rank)
write_csv(stability, file.path(out_dir, "rule_stability.csv"))

jaccard_rows <- list()
for (v in unique(selected$variant)) {
  sets <- split(selected$feature[selected$variant == v],
                interaction(selected$repetition[selected$variant == v], selected$fold[selected$variant == v]))
  cmb <- combn(seq_along(sets), 2)
  vals <- apply(cmb, 2, function(ii) length(intersect(sets[[ii[1]]], sets[[ii[2]]])) /
                  length(union(sets[[ii[1]]], sets[[ii[2]]])))
  jaccard_rows[[v]] <- tibble(variant = v, median_pairwise_jaccard = median(vals),
                              mean_pairwise_jaccard = mean(vals),
                              q25 = quantile(vals, .25), q75 = quantile(vals, .75),
                              unique_rules = n_distinct(unlist(sets)))
}
jaccard <- bind_rows(jaccard_rows)
old_selected <- read_csv(old_selected_path, show_col_types = FALSE)
old_sets <- split(old_selected$feature, interaction(old_selected$repetition, old_selected$fold))
old_cmb <- combn(seq_along(old_sets), 2)
old_vals <- apply(old_cmb, 2, function(ii) length(intersect(old_sets[[ii[1]]], old_sets[[ii[2]]])) /
                    length(union(old_sets[[ii[1]]], old_sets[[ii[2]]])))
jaccard <- bind_rows(jaccard, tibble(
  variant = "current_six_rule_hotlist",
  median_pairwise_jaccard = median(old_vals), mean_pairwise_jaccard = mean(old_vals),
  q25 = quantile(old_vals, .25), q75 = quantile(old_vals, .75),
  unique_rules = n_distinct(old_selected$feature)
))
write_csv(jaccard, file.path(out_dir, "selection_robustness.csv"))

# Source/platform sensitivity is descriptive because the older remapped source
# is much smaller. Only aggregate metrics are emitted; no sample identifiers.
source_map <- fread(seg_path, select = c("cg", "source_file")) %>% as_tibble() %>%
  transmute(sample_id = str_to_upper(str_trim(as.character(cg))), source_file) %>% distinct()
source_models <- c(
  "WHO + current six-rule hotlist" = "old_who_lp",
  "WHO + iterative geometry hotlist" = "triscut_full_no_background__who_lp",
  "Integrated grade + current six-rule hotlist" = "old_ig_lp",
  "Integrated grade + iterative geometry hotlist" = "triscut_full_no_background__ig_lp"
)
source_rep <- pred %>% left_join(source_map, by = "sample_id") %>%
  group_by(source_file, repetition) %>% group_modify(~ bind_rows(lapply(names(source_models), function(m) {
    tibble(model = m, auc = auc_safe(.x$endpoint, .x[[source_models[[m]]]]),
           n = nrow(.x), events = sum(.x$endpoint))
  }))) %>% ungroup()
source_summary <- source_rep %>% group_by(source_file, model) %>% summarise(
  n = first(n), events = first(events), auc_sd_repetition = sd(auc), auc = mean(auc), .groups = "drop"
)
write_csv(source_summary, file.path(out_dir, "source_sensitivity.csv"))

nejm <- c(blue = "#0072B5", purple = "#7876B1", red = "#BC3C29",
          orange = "#E18727", gray = "#6F6F6F", light = "#D9DEE3")
theme_pub <- theme_classic(base_family = "sans", base_size = 10.5) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey35"),
        axis.text = element_text(color = "black"),
        panel.grid.major.x = element_line(color = "grey90", linewidth = .35),
        legend.position = "bottom", plot.margin = margin(8, 10, 8, 8))

main_models <- c("WHO grade", "WHO + current six-rule hotlist", "WHO + TRISCUT-inspired hotlist",
                 "WHO + TRISCUT-inspired no-background ablation",
                 "Integrated grade", "Integrated grade + current six-rule hotlist",
                 "Integrated grade + TRISCUT-inspired hotlist",
                 "Integrated grade + TRISCUT-inspired no-background ablation")
pdat <- summary %>% filter(model %in% main_models) %>% mutate(
  model = factor(model, levels = rev(main_models)),
  family = ifelse(grepl("^WHO", as.character(model)), "WHO", "Integrated grade"),
  method = case_when(grepl("no-background", as.character(model)) ~ "Iterative geometry",
                     grepl("TRISCUT", as.character(model)) ~ "Background-adjusted",
                     grepl("current", as.character(model)) ~ "Current hotlist",
                     TRUE ~ "Grade alone")
)
p1 <- ggplot(pdat, aes(auc, model, color = family, shape = method)) +
  geom_errorbarh(aes(xmin = auc_low, xmax = auc_high), height = .15, linewidth = .55) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.3f", auc)), color = "black", nudge_x = .006,
            hjust = 0, size = 3) +
  scale_color_manual(values = c("WHO" = nejm[["blue"]], "Integrated grade" = nejm[["purple"]])) +
  scale_shape_manual(values = c("Grade alone" = 16, "Current hotlist" = 15,
                                "Background-adjusted" = 17, "Iterative geometry" = 18)) +
  scale_x_continuous(limits = c(.62, .82), breaks = seq(.65, .80, .05)) +
  labs(title = "A  Held-out discrimination", subtitle = "Strict 5-year endpoint; identical tumors and folds",
       x = "Area under the ROC curve", y = NULL, color = NULL, shape = NULL) + theme_pub

cdat <- comparisons %>% filter(comparison %in% c(
  "WHO: TRISCUT-inspired versus current hotlist", "IG: TRISCUT-inspired versus current hotlist",
  "WHO: iterative geometry versus current hotlist", "IG: iterative geometry versus current hotlist")) %>% mutate(
  comparator = ifelse(grepl("^WHO", comparison), "WHO", "Integrated grade"),
  method = ifelse(grepl("iterative geometry", comparison), "Iterative geometry", "Background-adjusted"),
  row_label = paste(comparator, method, sep = " — "),
  row_label = factor(row_label, levels = rev(c("WHO — Background-adjusted", "WHO — Iterative geometry",
                                               "Integrated grade — Background-adjusted",
                                               "Integrated grade — Iterative geometry")))
)
p2 <- ggplot(cdat, aes(delta_auc, row_label, color = comparator, shape = method)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
  geom_errorbarh(aes(xmin = delta_low, xmax = delta_high), height = .16, linewidth = .65) +
  geom_point(size = 3.2) +
  geom_text(aes(label = sprintf("%+.3f", delta_auc)), color = "black", nudge_x = .003,
            hjust = 0, size = 3.1) +
  scale_color_manual(values = c("WHO" = nejm[["blue"]], "Integrated grade" = nejm[["purple"]]), guide = "none") +
  scale_shape_manual(values = c("Background-adjusted" = 17, "Iterative geometry" = 18), guide = "none") +
  labs(title = "B  Change versus current hotlist", subtitle = "Paired patient-bootstrap 95% intervals",
       x = "Difference in area under the ROC curve", y = NULL) + theme_pub

sdat <- stability %>% filter(variant == "triscut_full_no_background") %>%
  slice_max(selection_fraction, n = 10, with_ties = FALSE) %>%
  mutate(report_label = factor(report_label, levels = rev(report_label)))
p3 <- ggplot(sdat, aes(selection_fraction, report_label, fill = event)) +
  geom_col(width = .68) +
  scale_fill_manual(values = c(cna = nejm[["blue"]], loss = nejm[["red"]], gain = nejm[["orange"]],
                               loh = nejm[["purple"]], ai = "#4DBBD5", architecture = nejm[["gray"]],
                               pair = "#8491B4"), na.value = nejm[["gray"]]) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(title = "C  Iterative-geometry rule stability", subtitle = "Top rules across 50 outer training fits",
       x = "Training folds selecting rule", y = NULL, fill = "Event class") + theme_pub

fig <- (p1 | p2) / p3 + plot_layout(heights = c(1.1, 1.15), guides = "keep") +
  plot_annotation(title = "TRISCUT-inspired breakpoint hotlist: performance and robustness",
                  subtitle = "Six-rule background-adjusted and iterative-geometry variants; not a faithful TRISCUT implementation",
                  theme = theme(plot.title = element_text(family = "sans", face = "bold", size = 16),
                                plot.subtitle = element_text(family = "sans", size = 10, color = "grey35"),
                                plot.margin = margin(12, 14, 10, 14)))
ggsave(file.path(out_dir, "figure_triscut_inspired_hotlist_comparison.png"), fig,
       width = 11.5, height = 10.0, dpi = 400, bg = "white")
ggsave(file.path(out_dir, "figure_triscut_inspired_hotlist_comparison.pdf"), fig,
       width = 11.5, height = 10.0)

write_lines(c(
  "TRISCUT-inspired breakpoint hotlist analysis",
  "",
  sprintf("Primary strict 5-year cohort: n=%d; events=%d.", nrow(dat), sum(dat$endpoint)),
  sprintf("Validation: %d repetitions x %d patient-grouped outer folds; fixed %d-rule hotlist.", repeats, outer_k, hotlist_size),
  "All breakpoint candidate screening, opportunity/background correction, interaction generation, iterative selection, and coefficient tuning occurred inside each outer training fold.",
  "Candidate events were CNA loss, CNA gain, LOH, and allelic imbalance. LOH/AI were treated as categorical event-call geometry only, not as phased or absolute copy number.",
  "Background correction residualized each candidate against training-fold total CNA breakpoint burden, event-specific breakpoint burden, and event-by-arm breakpoint opportunity.",
  "Iterative selection refit the WHO-plus-integrated-grade outcome residual after every selected rule and prohibited highly correlated or spatially redundant rules.",
  "The primary comparison used six rules to match the existing compact hotlist. CNA-only and no-background variants were engineering ablations defined before the primary run was inspected.",
  "Four- and eight-rule variants were added afterward as exploratory rule-count sensitivity analyses and are not primary results.",
  "This analysis borrows TRISCUT concepts but is not faithful TRISCUT because allele-specific absolute copy number, haplotype phasing, WGD-aware histories, Hi-C partial-event backgrounds, and TRISCUT's generative fitness model were unavailable."
), file.path(out_dir, "METHODS.txt"))

message("Completed: ", out_dir)
