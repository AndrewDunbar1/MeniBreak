scna_feature_names <- function(features) {
  grep(
    "^cnv_|^scna_(loss_count|gain_count|altered_arm_count)$",
    names(features), value = TRUE
  )
}

breakpoint_feature_names <- function(features) {
  c(
    scna_feature_names(features),
    intersect(
      c(
        "total_breakpoint_boundaries", "breakpoint_arms_touched",
        "breakpoint_route_entropy", "max_arm_breakpoint_fraction",
        "telomeric_breakpoint_fraction", "central_breakpoint_fraction",
        "centromeric_breakpoint_fraction", "mean_rel_tel_to_cent",
        "sd_rel_tel_to_cent", "cna_segment_count", "loss_segment_count",
        "gain_segment_count", "median_segment_length",
        "small_segment_count_lt5mb", "dynamic_hotspot_breakpoints_all",
        "dynamic_hotspot_breakpoints_non1p", "dynamic_hotspot_arms_used_all",
        "dynamic_hotspot_arms_used_non1p"
      ),
      names(features)
    )
  )
}

group_folds <- function(groups, k, seed) {
  set.seed(seed)
  ids <- sample(unique(as.character(groups)))
  map <- stats::setNames(rep(seq_len(k), length.out = length(ids)), ids)
  unname(map[as.character(groups)])
}

fit_enet_predict <- function(train, test, feature_names, y_col, group_col,
                             seed, alpha = 0.25,
                             type = c("response", "link")) {
  type <- match.arg(type)
  xtr <- as.matrix(train[, feature_names, drop = FALSE])
  xte <- as.matrix(test[, feature_names, drop = FALSE])
  storage.mode(xtr) <- storage.mode(xte) <- "double"
  med <- apply(xtr, 2, function(z) stats::median(z[is.finite(z)], na.rm = TRUE))
  med[!is.finite(med)] <- 0
  for (j in seq_along(feature_names)) {
    xtr[!is.finite(xtr[, j]), j] <- med[[j]]
    xte[!is.finite(xte[, j]), j] <- med[[j]]
  }
  inner_k <- min(5L, length(unique(train[[group_col]])))
  foldid <- group_folds(train[[group_col]], inner_k, seed)
  fit <- glmnet::cv.glmnet(
    xtr, train[[y_col]], family = "binomial", alpha = alpha,
    foldid = foldid, type.measure = "deviance", standardize = TRUE
  )
  drop(stats::predict(fit, xte, s = "lambda.1se", type = type))
}

evaluate_candidate <- function(data, scna_features, full_features,
                               repeats, folds, seed) {
  repeat_metrics <- vector("list", repeats)
  for (r in seq_len(repeats)) {
    outer <- group_folds(data$patient_id, folds, seed + r)
    ps <- pf <- rep(NA_real_, nrow(data))
    for (fold in seq_len(folds)) {
      tr <- data[outer != fold, , drop = FALSE]
      te <- data[outer == fold, , drop = FALSE]
      ps[outer == fold] <- fit_enet_predict(
        tr, te, scna_features, "endpoint", "patient_id",
        seed + 1000L * r + fold
      )
      pf[outer == fold] <- fit_enet_predict(
        tr, te, full_features, "endpoint", "patient_id",
        seed + 2000L * r + fold
      )
    }
    auc_scna <- auc_safe(data$endpoint, ps)
    auc_full <- auc_safe(data$endpoint, pf)
    repeat_metrics[[r]] <- tibble::tibble(
      repetition = r,
      auc_scna = auc_scna,
      auc_full = auc_full,
      auc_delta = auc_full - auc_scna,
      brier_scna = mean((ps - data$endpoint)^2),
      brier_full = mean((pf - data$endpoint)^2)
    )
  }
  dplyr::bind_rows(repeat_metrics)
}

fit_final_model <- function(data, feature_names, seed, alpha = 0.25) {
  x <- as.matrix(data[, feature_names, drop = FALSE])
  storage.mode(x) <- "double"
  imputation <- apply(x, 2, function(z) stats::median(z[is.finite(z)], na.rm = TRUE))
  imputation[!is.finite(imputation)] <- 0
  for (j in seq_along(feature_names)) x[!is.finite(x[, j]), j] <- imputation[[j]]
  foldid <- group_folds(data$patient_id, min(5L, length(unique(data$patient_id))), seed)
  fit <- glmnet::cv.glmnet(
    x, data$endpoint, family = "binomial", alpha = alpha,
    foldid = foldid, type.measure = "deviance", standardize = TRUE
  )
  co <- as.matrix(stats::coef(fit, s = "lambda.1se"))[, 1]
  raw <- drop(co[[1]] + x %*% co[-1])
  cal <- stats::glm(data$endpoint ~ raw, family = stats::binomial())
  list(
    intercept = unname(co[[1]]),
    coefficients = co[-1],
    calibration_intercept = unname(stats::coef(cal)[[1]]),
    calibration_slope = unname(stats::coef(cal)[[2]]),
    reference_scores = sort(raw),
    imputation = imputation,
    lambda = fit$lambda.1se,
    alpha = alpha
  )
}

run_bakeoff <- function(segments, outcomes, genome_build = "GRCh37",
                        v6_segments = NULL,
                        dgv_path = NULL, output_dir,
                        repeats = 10, folds = 5, seed = 20260703,
                        legacy_reproduction_verified = FALSE,
                        lock_model = FALSE) {
  assert_genome_build(genome_build)
  ensure_dir(output_dir)
  normalized <- normalize_segments(
    segments, format = "auto", genome_build = genome_build,
    allow_incompatible = TRUE
  )
  normalized_v6 <- if (!is.null(v6_segments)) {
    normalize_segments(
      v6_segments, format = "auto", genome_build = genome_build,
      allow_incompatible = TRUE
    )
  } else {
    NULL
  }
  out <- if (is.character(outcomes)) read_table_auto(outcomes) else
    tibble::as_tibble(outcomes)
  required <- c("sample_id", "patient_id", "recurrence_time", "recurrence_event")
  if (!all(required %in% names(out))) {
    mbpr_stop("Bakeoff outcomes require: ", paste(required, collapse = ", "))
  }
  candidates <- tidyr::crossing(
    preprocess = c("legacy", "v6"),
    arm_threshold = c(30, 50)
  )
  all_metrics <- list()
  feature_sets <- list()
  for (i in seq_len(nrow(candidates))) {
    prep <- candidates$preprocess[[i]]
    threshold <- candidates$arm_threshold[[i]]
    id <- paste(prep, threshold, sep = "_")
    feat <- if (prep == "v6" && !is.null(normalized_v6)) {
      extract_features(
        normalized_v6, genome_build, preprocess = "legacy",
        arm_threshold = threshold
      )
    } else {
      extract_features(
        normalized, genome_build, prep, threshold,
        dgv_path = if (prep == "v6") dgv_path else NULL
      )
    }
    feature_sets[[id]] <- feat
    merged <- dplyr::inner_join(out, feat, by = "sample_id")
    merged$endpoint <- fixed_horizon_endpoint(
      merged$recurrence_time, merged$recurrence_event, 5
    )
    merged <- merged[!is.na(merged$endpoint), , drop = FALSE]
    metrics <- evaluate_candidate(
      merged, scna_feature_names(feat), breakpoint_feature_names(feat),
      repeats, folds, seed + i * 10000L
    )
    metrics$candidate <- id
    metrics$preprocess <- prep
    metrics$arm_threshold <- threshold
    all_metrics[[id]] <- metrics
  }
  repeat_metrics <- dplyr::bind_rows(all_metrics)
  summary <- repeat_metrics %>%
    dplyr::group_by(.data$candidate, .data$preprocess, .data$arm_threshold) %>%
    dplyr::summarise(
      mean_auc_5y = mean(.data$auc_full),
      mean_scna_auc_5y = mean(.data$auc_scna),
      mean_auc_delta = mean(.data$auc_delta),
      repeats_auc_improved = mean(.data$auc_delta > 0),
      mean_brier_full = mean(.data$brier_full),
      mean_brier_scna = mean(.data$brier_scna),
      passes_breakpoint_rule = .data$repeats_auc_improved >= 0.80 &
        .data$mean_brier_full <= .data$mean_brier_scna,
      .groups = "drop"
    )
  best_auc <- max(summary$mean_auc_5y)
  preferred <- summary %>%
    dplyr::filter(.data$preprocess == "v6",
                  .data$mean_auc_5y >= best_auc - 0.01) %>%
    dplyr::arrange(.data$arm_threshold)
  selected <- if (nrow(preferred)) preferred[1, ] else
    summary[which.max(summary$mean_auc_5y), ]
  selected$protocol_deviation <- selected$preprocess != "v6"
  accepted <- isTRUE(selected$passes_breakpoint_rule) &&
    isTRUE(legacy_reproduction_verified)
  readr::write_csv(repeat_metrics, file.path(output_dir, "repeat_metrics.csv"))
  readr::write_csv(summary, file.path(output_dir, "candidate_summary.csv"))
  write_json(
    list(selected = as.list(selected), accepted = accepted,
         legacy_reproduction_verified = legacy_reproduction_verified,
         v6_reference_segments = !is.null(v6_segments),
         repeats = repeats, folds = folds, seed = seed),
    file.path(output_dir, "selection.json")
  )
  if (lock_model && !accepted) {
    mbpr_stop("Model locking blocked: acceptance rule or legacy reproduction failed.")
  }
  if (lock_model) {
    id <- selected$candidate[[1]]
    feat <- feature_sets[[id]]
    merged <- dplyr::inner_join(out, feat, by = "sample_id")
    full_names <- breakpoint_feature_names(feat)
    scna_names <- scna_feature_names(feat)
    artifact <- list(
      model_version = "mbpr-derivation-v1",
      preprocessing_lineage = list(
        preprocess = selected$preprocess[[1]],
        arm_threshold = selected$arm_threshold[[1]],
        dgv_reference_sha256 = if (!is.null(dgv_path)) sha256_file(dgv_path) else NULL
      ),
      training_population = list(samples = nrow(merged),
                                 patients = dplyr::n_distinct(merged$patient_id)),
      feature_order = full_names,
      imputation = stats::setNames(rep(0, length(full_names)), full_names),
      center = stats::setNames(rep(0, length(full_names)), full_names),
      scale = stats::setNames(rep(1, length(full_names)), full_names),
      horizons = list(), scna_only_horizons = list(),
      seeds = list(master = seed)
    )
    for (h in c(3, 5, 7)) {
      dh <- merged
      dh$endpoint <- fixed_horizon_endpoint(
        dh$recurrence_time, dh$recurrence_event, h
      )
      dh <- dh[!is.na(dh$endpoint), , drop = FALSE]
      ff <- fit_final_model(dh, full_names, seed + h)
      sf <- fit_final_model(dh, scna_names, seed + 100 + h)
      artifact$horizons[[as.character(h)]] <- ff[setdiff(names(ff), "imputation")]
      sc <- stats::setNames(rep(0, length(full_names)), full_names)
      sc[names(sf$coefficients)] <- sf$coefficients
      sf$coefficients <- sc
      artifact$scna_only_horizons[[as.character(h)]] <-
        sf[setdiff(names(sf), "imputation")]
      artifact$imputation[names(ff$imputation)] <- ff$imputation
    }
    artifact$checksums <- list(
      hotspot_windows = sha256_file(asset_path("hotspot_windows_v1.csv")),
      arm_reference = sha256_file(asset_path(
        paste0("arm_reference_", genome_build, ".csv")
      ))
    )
    artifact_path <- file.path(output_dir, "mbpr-derivation-v1.rds")
    saveRDS(artifact, artifact_path, version = 3)
  }
  list(summary = summary, selected = selected, accepted = accepted)
}
