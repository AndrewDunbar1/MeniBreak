model_registry <- function() {
  yaml::read_yaml(model_asset_path("registry.yml"))
}

model_info <- function(as_json = FALSE) {
  registry <- model_registry()
  artifact_path <- registry$artifact
  artifact <- NULL
  if (!is.null(artifact_path)) {
    full <- model_asset_path(artifact_path)
    if (file.exists(full)) artifact <- readRDS(full)
  }
  out <- list(
    software_version = registry$software_version,
    model_version = registry$model_version,
    status = registry$status,
    primary_endpoint = registry$primary_endpoint,
    secondary_endpoints = registry$secondary_endpoints,
    preprocessing_lineage = artifact$preprocessing_lineage %||% NULL,
    training_population = artifact$training_population %||% NULL,
    feature_order = artifact$feature_order %||% NULL,
    checksums = artifact$checksums %||% NULL,
    limitations = c(
      "Research use only; not for clinical decisions.",
      "Absolute probabilities are internally calibrated and provisional.",
      "External cohorts require an explicit genome build and validated preprocessing."
    )
  )
  if (as_json) jsonlite::toJSON(out, auto_unbox = TRUE, pretty = TRUE) else out
}

load_locked_model <- function() {
  registry <- model_registry()
  if (!identical(registry$status, "LOCKED") || is.null(registry$artifact)) {
    mbpr_stop(
      "Prediction model is not locked. Feature extraction and external validation ",
      "are available, but recurrence predictions are intentionally blocked. Run ",
      "the prespecified bakeoff and lock a model only if the acceptance rule passes."
    )
  }
  path <- model_asset_path(registry$artifact)
  if (!file.exists(path)) mbpr_stop("Locked model artifact is missing: ", path)
  observed <- sha256_file(path)
  if (!is.null(registry$artifact_sha256) &&
      !identical(observed, registry$artifact_sha256)) {
    mbpr_stop("Locked model checksum mismatch.")
  }
  readRDS(path)
}

matrix_for_model <- function(features, feature_order, imputation, center, scale) {
  missing <- setdiff(feature_order, names(features))
  if (length(missing)) mbpr_stop("Feature table is missing locked feature(s): ",
                                 paste(missing, collapse = ", "))
  x <- as.matrix(features[, feature_order, drop = FALSE])
  storage.mode(x) <- "double"
  for (j in seq_along(feature_order)) {
    x[!is.finite(x[, j]), j] <- imputation[[feature_order[[j]]]]
  }
  scale(x, center = center[feature_order], scale = scale[feature_order])
}

score_one_horizon <- function(x, spec) {
  raw <- drop(spec$intercept + x %*% spec$coefficients[colnames(x)])
  probability <- stats::plogis(spec$calibration_intercept +
                                 spec$calibration_slope * raw)
  percentile <- vapply(
    raw,
    function(z) mean(spec$reference_scores <= z, na.rm = TRUE) * 100,
    numeric(1)
  )
  list(raw = raw, probability = probability, percentile = percentile)
}

score_feature_table <- function(features, model = load_locked_model()) {
  x <- matrix_for_model(
    features, model$feature_order, model$imputation,
    model$center, model$scale
  )
  out <- tibble::tibble(sample_id = features$sample_id)
  for (h in c("3", "5", "7")) {
    s <- score_one_horizon(x, model$horizons[[h]])
    out[[paste0("raw_score_", h, "y")]] <- s$raw
    out[[paste0("reference_percentile_", h, "y")]] <- s$percentile
    out[[paste0("provisional_probability_", h, "y")]] <- s$probability
    if (!is.null(model$scna_only_horizons[[h]])) {
      ss <- score_one_horizon(x, model$scna_only_horizons[[h]])
      out[[paste0("scna_only_raw_score_", h, "y")]] <- ss$raw
      out[[paste0("scna_only_provisional_probability_", h, "y")]] <- ss$probability
    }
  }
  out$model_version <- model$model_version
  out
}

score_segments <- function(segments, format = "auto", genome_build,
                           output_dir = NULL, mapping_file = NULL,
                           caller_profile = NULL, features_format = "csv") {
  source_path <- if (is.character(segments) && length(segments) == 1L) segments else NULL
  normalized <- normalize_segments(
    segments, format = format, genome_build = genome_build,
    mapping_file = mapping_file, caller_profile = caller_profile,
    allow_incompatible = TRUE
  )
  incompatible_samples <- attr(normalized, "mbpr_incompatible_samples") %||%
    character()
  model <- load_locked_model()
  lineage <- model$preprocessing_lineage
  features <- extract_features(
    normalized, genome_build = genome_build,
    preprocess = lineage$preprocess,
    arm_threshold = lineage$arm_threshold,
    dgv_path = lineage$dgv_path %||% NULL
  )
  predictions <- score_feature_table(features, model)
  sample_qc <- attr(features, "mbpr_feature_qc")$per_sample
  predictions <- dplyr::left_join(
    predictions,
    sample_qc %>% dplyr::select("sample_id", "qc_status"),
    by = "sample_id"
  )
  if (length(incompatible_samples)) {
    bad <- predictions$sample_id %in% incompatible_samples
    prediction_columns <- grep(
      "^(raw_score|reference_percentile|provisional_probability|scna_only_raw_score|scna_only_provisional_probability)_",
      names(predictions), value = TRUE
    )
    predictions[bad, prediction_columns] <- NA_real_
    predictions$qc_status[bad] <- "FAIL_INCOMPATIBLE_OVERLAP"
  }
  failed_qc <- startsWith(predictions$qc_status, "FAIL_")
  if (any(failed_qc)) {
    prediction_columns <- grep(
      "^(raw_score|reference_percentile|provisional_probability|scna_only_raw_score|scna_only_provisional_probability)_",
      names(predictions), value = TRUE
    )
    predictions[failed_qc, prediction_columns] <- NA_real_
  }
  if (!is.null(model$feature_reference)) {
    refs <- model$feature_reference
    xref <- as.matrix(features[, model$feature_order, drop = FALSE])
    low <- refs["low", model$feature_order]
    high <- refs["high", model$feature_order]
    outside <- rowSums(
      sweep(xref, 2, low, "<") | sweep(xref, 2, high, ">"),
      na.rm = TRUE
    )
    predictions$ood_feature_count <- outside
    predictions$ood_flag <- outside > max(3L, ceiling(length(model$feature_order) * 0.10))
  } else {
    predictions$ood_feature_count <- NA_integer_
    predictions$ood_flag <- NA
  }
  if (!is.null(output_dir)) {
    ensure_dir(output_dir)
    readr::write_csv(predictions, file.path(output_dir, "predictions.csv"))
    if (identical(features_format, "parquet")) {
      nanoparquet::write_parquet(
        as.data.frame(features), file.path(output_dir, "features.parquet")
      )
    } else {
      readr::write_csv(features, file.path(output_dir, "features.csv"))
    }
    qc <- attr(features, "mbpr_feature_qc")
    qc$incompatible_samples <- incompatible_samples
    qc$out_of_distribution <- predictions %>%
      dplyr::select("sample_id", "ood_feature_count", "ood_flag")
    write_json(qc, file.path(output_dir, "qc.json"))
    manifest <- list(
      command = "mbpr score",
      arguments = commandArgs(trailingOnly = TRUE),
      timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      package_version = package_version_safe(),
      r_version = R.version.string,
      dependency_versions = as.list(vapply(
        c("data.table", "dplyr", "glmnet", "nanoparquet", "readr"),
        function(pkg) as.character(utils::packageVersion(pkg)),
        character(1)
      )),
      model_version = model$model_version,
      adapter = attr(normalized, "mbpr_adapter"),
      caller_profile = attr(normalized, "mbpr_caller_profile"),
      genome_build = genome_build,
      input_sha256 = if (!is.null(source_path)) sha256_file(source_path) else NULL,
      model_checksums = model$checksums
    )
    write_json(manifest, file.path(output_dir, "run_manifest.json"))
    write_score_report(predictions, attr(features, "mbpr_feature_qc"),
                       file.path(output_dir, "report.html"))
  }
  predictions
}
