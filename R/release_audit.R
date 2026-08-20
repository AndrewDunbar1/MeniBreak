derive_hotspot_windows <- function(boundaries, arm_calls, min_loss_samples = 20L,
                                   n_bins = 40L, window_bins = 4L,
                                   excluded_edge_bins = 2L) {
  loss_samples <- arm_calls %>%
    dplyr::filter(.data$has_loss) %>%
    dplyr::count(.data$arm_key, name = "loss_samples")
  counts <- boundaries %>%
    dplyr::mutate(
      bin = pmin(floor(.data$rel_tel_to_cent * n_bins), n_bins - 1L)
    ) %>%
    dplyr::inner_join(
      arm_calls %>%
        dplyr::filter(.data$has_loss) %>%
        dplyr::select("sample_id", "arm_key"),
      by = c("sample_id", "arm_key")
    ) %>%
    dplyr::count(.data$arm_key, .data$bin, name = "breakpoints")
  arms <- loss_samples$arm_key[loss_samples$loss_samples >= min_loss_samples]
  dplyr::bind_rows(lapply(arms, function(arm) {
    arm_counts <- counts[counts$arm_key == arm, , drop = FALSE]
    starts <- excluded_edge_bins:(n_bins - window_bins - excluded_edge_bins)
    totals <- vapply(starts, function(s) {
      sum(arm_counts$breakpoints[
        arm_counts$bin >= s & arm_counts$bin < s + window_bins
      ], na.rm = TRUE)
    }, numeric(1))
    best <- starts[which.max(totals)]
    tibble::tibble(
      arm_key = arm,
      loss_samples = loss_samples$loss_samples[loss_samples$arm_key == arm],
      hotspot_start_norm = best / n_bins,
      hotspot_end_norm = (best + window_bins) / n_bins,
      hotspot_breakpoints = max(totals)
    )
  }))
}

hotspot_stability_audit <- function(normalized, outcomes, repeats = 10,
                                    folds = 5, seed = 20260703,
                                    arm_threshold = 30) {
  arm_ref <- load_arm_reference(unique(normalized$genome_build))
  breakpoint_model <- normalized %>%
    dplyr::filter(.data$event %in% c("loss", "gain", "allelic_imbalance"))
  seg_model <- normalized %>%
    dplyr::filter(.data$event %in% c("loss", "gain"))
  frozen <- readr::read_csv(
    asset_path("hotspot_windows_v1.csv"), show_col_types = FALSE
  )
  audits <- list()
  for (r in seq_len(repeats)) {
    patient_fold <- group_folds(outcomes$patient_id, folds, seed + r)
    for (fold in seq_len(folds)) {
      training_samples <- outcomes$sample_id[patient_fold != fold]
      train_seg <- seg_model[seg_model$sample_id %in% training_samples, , drop = FALSE]
      train_bp <- breakpoint_model[
        breakpoint_model$sample_id %in% training_samples, , drop = FALSE
      ]
      calls <- make_arm_calls(train_seg, arm_ref, arm_threshold)
      boundaries <- make_boundaries(train_bp, arm_ref)
      derived <- derive_hotspot_windows(boundaries, calls)
      comparison <- frozen %>%
        dplyr::left_join(derived, by = "arm_key", suffix = c("_frozen", "_fold")) %>%
        dplyr::mutate(
          repetition = r,
          fold = fold,
          derived = !is.na(.data$hotspot_start_norm_fold),
          exact_match = .data$derived &
            abs(.data$hotspot_start_norm_frozen - .data$hotspot_start_norm_fold) < 1e-12 &
            abs(.data$hotspot_end_norm_frozen - .data$hotspot_end_norm_fold) < 1e-12,
          overlap = pmax(
            0,
            pmin(.data$hotspot_end_norm_frozen, .data$hotspot_end_norm_fold) -
              pmax(.data$hotspot_start_norm_frozen, .data$hotspot_start_norm_fold)
          ),
          union = pmax(.data$hotspot_end_norm_frozen, .data$hotspot_end_norm_fold) -
            pmin(.data$hotspot_start_norm_frozen, .data$hotspot_start_norm_fold),
          window_jaccard = dplyr::if_else(.data$derived, .data$overlap / .data$union, NA_real_)
        )
      audits[[length(audits) + 1L]] <- comparison
    }
  }
  dplyr::bind_rows(audits)
}

selected_oof_audit <- function(features, outcomes, source_manifest = NULL,
                               repeats = 10, folds = 5, seed = 20260703) {
  full_names <- breakpoint_feature_names(features)
  scna_names <- scna_feature_names(features)
  base <- dplyr::inner_join(outcomes, features, by = "sample_id")
  metrics <- predictions <- source_metrics <- list()
  for (h in c(3, 5, 7)) {
    d <- base
    d$endpoint <- fixed_horizon_endpoint(
      d$recurrence_time, d$recurrence_event, h
    )
    d <- d[!is.na(d$endpoint), , drop = FALSE]
    for (r in seq_len(repeats)) {
      outer <- group_folds(d$patient_id, folds, seed + h * 10000L + r)
      raw_scna <- raw_full <- rep(NA_real_, nrow(d))
      for (fold in seq_len(folds)) {
        tr <- d[outer != fold, , drop = FALSE]
        te <- d[outer == fold, , drop = FALSE]
        raw_scna[outer == fold] <- fit_enet_predict(
          tr, te, scna_names, "endpoint", "patient_id",
          seed + h * 100000L + r * 1000L + fold, type = "link"
        )
        raw_full[outer == fold] <- fit_enet_predict(
          tr, te, full_names, "endpoint", "patient_id",
          seed + h * 200000L + r * 1000L + fold, type = "link"
        )
      }
      prob_scna <- stats::plogis(raw_scna)
      prob_full <- stats::plogis(raw_full)
      cal <- calibration_stats(d$endpoint, prob_full)
      auc_scna <- auc_safe(d$endpoint, raw_scna)
      auc_full <- auc_safe(d$endpoint, raw_full)
      metrics[[length(metrics) + 1L]] <- tibble::tibble(
        horizon_years = h, repetition = r, n = nrow(d),
        events = sum(d$endpoint),
        auc_scna = auc_scna,
        auc_full = auc_full,
        auc_delta = auc_full - auc_scna,
        average_precision_scna = average_precision(d$endpoint, raw_scna),
        average_precision_full = average_precision(d$endpoint, raw_full),
        brier_scna = mean((prob_scna - d$endpoint)^2),
        brier_full = mean((prob_full - d$endpoint)^2),
        calibration_intercept = cal[["intercept"]],
        calibration_slope = cal[["slope"]],
        event_capture_top_10pct = event_capture_curve(
          d$endpoint, raw_full, 0.10
        )$event_capture[[1]],
        event_capture_top_25pct = event_capture_curve(
          d$endpoint, raw_full, 0.25
        )$event_capture[[1]]
      )
      pred <- tibble::tibble(
        sample_id = d$sample_id, patient_id = d$patient_id,
        horizon_years = h, repetition = r, endpoint = d$endpoint,
        raw_scna = raw_scna, raw_full = raw_full,
        probability_scna = prob_scna, probability_full = prob_full
      )
      predictions[[length(predictions) + 1L]] <- pred
      if (!is.null(source_manifest)) {
        ps <- dplyr::left_join(pred, source_manifest, by = "sample_id") %>%
          dplyr::group_by(.data$source) %>%
          dplyr::summarise(
            n = dplyr::n(), events = sum(.data$endpoint),
            auc_full = if (dplyr::n_distinct(.data$endpoint) == 2L)
              auc_safe(.data$endpoint, .data$raw_full) else NA_real_,
            auc_scna = if (dplyr::n_distinct(.data$endpoint) == 2L)
              auc_safe(.data$endpoint, .data$raw_scna) else NA_real_,
            .groups = "drop"
          ) %>%
          dplyr::mutate(horizon_years = h, repetition = r)
        source_metrics[[length(source_metrics) + 1L]] <- ps
      }
    }
  }
  list(
    metrics = dplyr::bind_rows(metrics),
    predictions = dplyr::bind_rows(predictions),
    source_metrics = dplyr::bind_rows(source_metrics)
  )
}

audit_selected_model <- function(segments, outcomes, output_dir,
                                 genome_build = "GRCh37",
                                 repeats = 10, folds = 5,
                                 seed = 20260703) {
  ensure_dir(output_dir)
  raw <- if (is.character(segments)) read_table_auto(segments) else
    tibble::as_tibble(segments)
  normalized <- normalize_segments(
    raw, genome_build = genome_build, allow_incompatible = TRUE
  )
  features <- extract_features(
    normalized, genome_build, preprocess = "legacy", arm_threshold = 30
  )
  out <- if (is.character(outcomes)) read_table_auto(outcomes) else
    tibble::as_tibble(outcomes)
  source_manifest <- NULL
  if ("source_file" %in% names(raw)) {
    id <- if ("cg" %in% names(raw)) trimws(as.character(raw$cg)) else
      trimws(as.character(raw$sample_id))
    if ("sample" %in% names(raw)) {
      fallback <- trimws(as.character(raw$sample))
      id[is.na(id) | !nzchar(id)] <- fallback[is.na(id) | !nzchar(id)]
    }
    source_manifest <- tibble::tibble(
      sample_id = id, source = as.character(raw$source_file)
    ) %>% dplyr::distinct(.data$sample_id, .keep_all = TRUE)
  }
  oof <- selected_oof_audit(
    features, out, source_manifest, repeats, folds, seed
  )
  hotspot <- hotspot_stability_audit(
    normalized, out, repeats, folds, seed, arm_threshold = 30
  )
  readr::write_csv(oof$metrics, file.path(output_dir, "secondary_metrics.csv"))
  readr::write_csv(oof$predictions, file.path(output_dir, "oof_predictions.csv"))
  readr::write_csv(oof$source_metrics, file.path(output_dir, "source_sensitivity.csv"))
  readr::write_csv(hotspot, file.path(output_dir, "hotspot_stability.csv"))
  summary <- list(
    candidate = "legacy_30",
    repeats = repeats, folds = folds, seed = seed,
    endpoint_summary = oof$metrics %>%
      dplyr::group_by(.data$horizon_years) %>%
      dplyr::summarise(
        mean_auc_full = mean(.data$auc_full),
        mean_auc_scna = mean(.data$auc_scna),
        mean_auc_delta = mean(.data$auc_delta),
        auc_improved_fraction = mean(.data$auc_delta > 0),
        mean_brier_full = mean(.data$brier_full),
        mean_average_precision_full = mean(.data$average_precision_full),
        .groups = "drop"
      ),
    hotspot_summary = hotspot %>%
      dplyr::summarise(
        fold_windows_derived_fraction = mean(.data$derived),
        exact_match_fraction = mean(.data$exact_match, na.rm = TRUE),
        median_window_jaccard = stats::median(.data$window_jaccard, na.rm = TRUE)
      )
  )
  write_json(summary, file.path(output_dir, "release_audit_summary.json"))
  invisible(list(oof = oof, hotspot = hotspot, summary = summary))
}
