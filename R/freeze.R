calibration_from_oof <- function(oof, outcomes, horizon, score_col) {
  averaged <- oof[oof$horizon_years == horizon, , drop = FALSE] %>%
    dplyr::group_by(.data$sample_id) %>%
    dplyr::summarise(raw = mean(.data[[score_col]]), .groups = "drop") %>%
    dplyr::inner_join(outcomes, by = "sample_id")
  averaged$endpoint <- fixed_horizon_endpoint(
    averaged$recurrence_time, averaged$recurrence_event, horizon
  )
  averaged <- averaged[!is.na(averaged$endpoint), , drop = FALSE]
  fit <- stats::glm(endpoint ~ raw, data = averaged, family = stats::binomial())
  list(
    intercept = unname(stats::coef(fit)[[1]]),
    slope = unname(stats::coef(fit)[[2]]),
    n = nrow(averaged),
    events = sum(averaged$endpoint)
  )
}

freeze_selected_model <- function(segments, outcomes, oof_predictions,
                                  legacy_audit_json, selection_json,
                                  release_audit_json, output_path,
                                  genome_build = "GRCh37",
                                  seed = 20260703) {
  legacy <- jsonlite::read_json(legacy_audit_json, simplifyVector = TRUE)
  selection <- jsonlite::read_json(selection_json, simplifyVector = TRUE)
  release <- jsonlite::read_json(release_audit_json, simplifyVector = TRUE)
  if (!isTRUE(legacy$passed)) mbpr_stop("Legacy reproduction gate failed.")
  if (!isTRUE(selection$accepted) ||
      !identical(selection$selected$candidate, "legacy_30")) {
    mbpr_stop("Prespecified bakeoff did not accept legacy_30.")
  }
  release_5y <- release$endpoint_summary[
    release$endpoint_summary$horizon_years == 5, , drop = FALSE
  ]
  if (!nrow(release_5y) ||
      release_5y$auc_improved_fraction[[1]] < 0.80 ||
      release_5y$mean_auc_delta[[1]] <= 0) {
    mbpr_stop("Release audit did not confirm the breakpoint acceptance rule.")
  }
  normalized <- normalize_segments(
    segments, genome_build = genome_build, allow_incompatible = TRUE
  )
  features <- extract_features(
    normalized, genome_build, preprocess = "legacy", arm_threshold = 30
  )
  out <- if (is.character(outcomes)) read_table_auto(outcomes) else
    tibble::as_tibble(outcomes)
  oof <- if (is.character(oof_predictions)) read_table_auto(oof_predictions) else
    tibble::as_tibble(oof_predictions)
  merged <- dplyr::inner_join(out, features, by = "sample_id")
  full_names <- breakpoint_feature_names(features)
  scna_names <- scna_feature_names(features)
  xall <- as.matrix(merged[, full_names, drop = FALSE])
  imputation <- apply(xall, 2, function(z) {
    value <- stats::median(z[is.finite(z)], na.rm = TRUE)
    if (is.finite(value)) value else 0
  })
  artifact <- list(
    model_version = "mbpr-derivation-v1",
    software_version = "0.1.0",
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    research_use_only = TRUE,
    primary_endpoint = "5-year recurrence",
    secondary_endpoints = c("3-year recurrence", "7-year recurrence"),
    preprocessing_lineage = list(
      candidate = "legacy_30",
      preprocess = "legacy",
      arm_threshold = 30,
      genome_build = genome_build,
      allelic_imbalance_breakpoints_included = TRUE,
      loh_excluded_from_model_features = TRUE,
      assay_nf2_excluded = TRUE,
      external_incompatible_overlaps = "sample-level scoring failure"
    ),
    training_population = list(
      segment_samples = nrow(features),
      outcome_linked_samples = nrow(merged),
      patients = dplyr::n_distinct(merged$patient_id)
    ),
    feature_order = full_names,
    scna_feature_order = scna_names,
    feature_reference = rbind(
      low = apply(xall, 2, stats::quantile, probs = 0.005, na.rm = TRUE),
      high = apply(xall, 2, stats::quantile, probs = 0.995, na.rm = TRUE)
    ),
    imputation = imputation,
    center = stats::setNames(rep(0, length(full_names)), full_names),
    scale = stats::setNames(rep(1, length(full_names)), full_names),
    horizons = list(),
    scna_only_horizons = list(),
    seeds = list(master = seed, final = c(`3` = seed + 3, `5` = seed + 5,
                                         `7` = seed + 7)),
    checksums = list()
  )
  for (h in c(3, 5, 7)) {
    d <- merged
    d$endpoint <- fixed_horizon_endpoint(
      d$recurrence_time, d$recurrence_event, h
    )
    d <- d[!is.na(d$endpoint), , drop = FALSE]
    ff <- fit_final_model(d, full_names, seed + h)
    sf <- fit_final_model(d, scna_names, seed + 100L + h)
    full_cal <- calibration_from_oof(oof, out, h, "raw_full")
    scna_cal <- calibration_from_oof(oof, out, h, "raw_scna")
    ff$calibration_intercept <- full_cal$intercept
    ff$calibration_slope <- full_cal$slope
    ff$calibration_n <- full_cal$n
    ff$calibration_events <- full_cal$events
    sf$calibration_intercept <- scna_cal$intercept
    sf$calibration_slope <- scna_cal$slope
    sf$calibration_n <- scna_cal$n
    sf$calibration_events <- scna_cal$events
    artifact$horizons[[as.character(h)]] <- ff[setdiff(names(ff), "imputation")]
    sc <- stats::setNames(rep(0, length(full_names)), full_names)
    sc[names(sf$coefficients)] <- sf$coefficients
    sf$coefficients <- sc
    artifact$scna_only_horizons[[as.character(h)]] <-
      sf[setdiff(names(sf), "imputation")]
    artifact$training_population[[paste0("endpoint_", h, "y")]] <- list(
      samples = nrow(d), patients = dplyr::n_distinct(d$patient_id),
      events = sum(d$endpoint)
    )
  }
  artifact$checksums <- list(
    input_segments_sha256 = if (is.character(segments)) sha256_file(segments) else NULL,
    outcomes_sha256 = if (is.character(outcomes)) sha256_file(outcomes) else NULL,
    legacy_audit_sha256 = sha256_file(legacy_audit_json),
    selection_sha256 = sha256_file(selection_json),
    release_audit_sha256 = sha256_file(release_audit_json),
    hotspot_windows_sha256 = sha256_file(asset_path("hotspot_windows_v1.csv")),
    arm_reference_GRCh37_sha256 = sha256_file(asset_path("arm_reference_GRCh37.csv")),
    arm_reference_GRCh38_sha256 = sha256_file(asset_path("arm_reference_GRCh38.csv"))
  )
  ensure_dir(dirname(output_path))
  saveRDS(artifact, output_path, version = 3)
  manifest <- list(
    model_version = artifact$model_version,
    artifact = basename(output_path),
    artifact_sha256 = sha256_file(output_path),
    feature_count = length(full_names),
    primary_endpoint = artifact$primary_endpoint,
    training_population = artifact$training_population,
    checksums = artifact$checksums
  )
  write_json(manifest, paste0(output_path, ".manifest.json"))
  invisible(manifest)
}
