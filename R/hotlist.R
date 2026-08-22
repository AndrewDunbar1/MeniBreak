hotlist_definition <- function() {
  readr::read_csv(asset_path("hotlist_v1.csv"), show_col_types = FALSE)
}

#' Extract the 20-bin breakpoint candidate matrix used by dynamic top-eight validation
#'
#' Each chromosome arm is divided into 20 equal arm-relative bins. A feature is
#' one when a sample has at least one eligible segment boundary in that bin and
#' zero otherwise. Candidate screening and selection are deliberately not done
#' here: they belong inside the outcome-linked validation workflow.
#'
#' @param segments Canonical output from `normalize_segments()`.
#' @param genome_build Explicit GRCh37 or GRCh38 build.
#' @param preprocess Legacy or V6 cleaning lineage.
#' @param dgv_path Optional DGV asset for V6 preprocessing.
#' @return One row per sample and one binary `bp_<arm>_binNN` feature per arm bin.
#' @export
extract_breakpoint_candidates <- function(segments, genome_build = NULL,
                                          preprocess = c("legacy", "v6"),
                                          dgv_path = NULL) {
  preprocess <- match.arg(preprocess)
  required <- c("sample_id", "chromosome", "start", "end", "event", "genome_build")
  if (!all(required %in% names(segments))) {
    mbpr_stop("extract_breakpoint_candidates requires output from normalize_segments().")
  }
  build <- genome_build %||% unique(segments$genome_build)
  assert_genome_build(build)
  clean <- preprocess_segments(segments, mode = preprocess, dgv_path = dgv_path)
  samples <- unique(clean$sample_id)
  eligible <- clean %>%
    dplyr::filter(.data$event %in% c("loss", "gain", "allelic_imbalance"))
  boundaries <- make_boundaries(eligible, load_arm_reference(build)) %>%
    dplyr::mutate(
      bin = pmin(20L, floor(.data$rel_tel_to_cent * 20) + 1L),
      feature = paste0("bp_", .data$arm_key, "_bin", sprintf("%02d", .data$bin))
    ) %>%
    dplyr::distinct(.data$sample_id, .data$feature) %>%
    dplyr::mutate(value = 1L)

  arms <- load_arm_reference(build) %>%
    dplyr::transmute(arm_key = paste0(sub("^chr", "", .data$chromosome), .data$arm)) %>%
    dplyr::pull(.data$arm_key) %>%
    unique()
  feature_order <- unlist(lapply(arms, function(arm) {
    paste0("bp_", arm, "_bin", sprintf("%02d", seq_len(20L)))
  }), use.names = FALSE)

  out <- tidyr::expand_grid(sample_id = samples, feature = feature_order) %>%
    dplyr::left_join(boundaries, by = c("sample_id", "feature")) %>%
    dplyr::mutate(value = dplyr::coalesce(.data$value, 0L)) %>%
    tidyr::pivot_wider(names_from = "feature", values_from = "value") %>%
    dplyr::select("sample_id", dplyr::all_of(feature_order))
  attr(out, "menibreak_candidate_matrix") <- list(
    version = "dynamic-top8-candidates-v1",
    genome_build = build,
    bins_per_arm = 20L,
    coordinate_system = "arm-relative telomere-to-centromere",
    encoding = "binary presence of one or more eligible boundaries",
    eligible_events = c("loss", "gain", "allelic_imbalance")
  )
  out
}

#' Extract the fixed MeniBreak nine-window hotlist
#'
#' Counts eligible segment boundaries in the nine arm-relative 5% windows
#' described in the MeniBreak collaboration document. This is feature
#' extraction, not application of a frozen prediction equation.
#'
#' @param segments Canonical output from `normalize_segments()`.
#' @param genome_build Explicit GRCh37 or GRCh38 build; inferred only from the
#'   already-normalized canonical table when omitted.
#' @param preprocess Legacy or V6 cleaning lineage.
#' @param dgv_path Optional DGV asset for V6 preprocessing.
#' @return One row per sample with a count and presence flag for each window.
#' @export
extract_hotlist_features <- function(segments, genome_build = NULL,
                                     preprocess = c("legacy", "v6"),
                                     dgv_path = NULL) {
  preprocess <- match.arg(preprocess)
  required <- c("sample_id", "chromosome", "start", "end", "event", "genome_build")
  if (!all(required %in% names(segments))) {
    mbpr_stop("extract_hotlist_features requires output from normalize_segments().")
  }
  build <- genome_build %||% unique(segments$genome_build)
  assert_genome_build(build)
  clean <- preprocess_segments(segments, mode = preprocess, dgv_path = dgv_path)
  samples <- unique(clean$sample_id)
  eligible <- clean %>%
    dplyr::filter(.data$event %in% c("loss", "gain", "allelic_imbalance"))
  boundaries <- make_boundaries(eligible, load_arm_reference(build))
  definition <- hotlist_definition()

  out <- tibble::tibble(sample_id = samples)
  for (i in seq_len(nrow(definition))) {
    d <- definition[i, ]
    counts <- boundaries %>%
      dplyr::filter(
        .data$arm_key == d$arm_key[[1]],
        .data$rel_tel_to_cent >= d$start_fraction_tel_to_cent[[1]],
        .data$rel_tel_to_cent < d$end_fraction_tel_to_cent[[1]]
      ) %>%
      dplyr::count(.data$sample_id, name = "value")
    count_name <- paste0(d$feature_id[[1]], "_count")
    present_name <- paste0(d$feature_id[[1]], "_present")
    names(counts)[names(counts) == "value"] <- count_name
    out <- dplyr::left_join(out, counts, by = "sample_id")
    out[[count_name]][is.na(out[[count_name]])] <- 0L
    out[[present_name]] <- as.integer(out[[count_name]] > 0)
  }
  attr(out, "menibreak_hotlist") <- list(
    version = "hotlist-v1",
    genome_build = build,
    coordinate_system = "arm-relative telomere-to-centromere",
    window_interval = "left-closed, right-open",
    eligible_events = c("loss", "gain", "allelic_imbalance"),
    excluded_events = c("loh", "neutral"),
    unit = "segment-boundary count"
  )
  out
}

write_hotlist_outputs <- function(segments, format, genome_build, output_dir,
                                  mapping_file = NULL, caller_profile = NULL,
                                  preprocess = "legacy", dgv_path = NULL) {
  ensure_dir(output_dir)
  normalized <- normalize_segments(
    segments, format = format, genome_build = genome_build,
    mapping_file = mapping_file, caller_profile = caller_profile,
    allow_incompatible = TRUE
  )
  incompatible_samples <- attr(normalized, "mbpr_incompatible_samples") %||% character()
  features <- extract_hotlist_features(
    normalized, genome_build = genome_build, preprocess = preprocess,
    dgv_path = dgv_path
  )
  hotlist_metadata <- attr(features, "menibreak_hotlist")
  features$qc_status <- ifelse(
    features$sample_id %in% incompatible_samples,
    "FAIL_INCOMPATIBLE_OVERLAP", "PASS"
  )
  readr::write_csv(features, file.path(output_dir, "hotlist_features.csv"))
  readr::write_csv(
    features %>% dplyr::select("sample_id", "qc_status"),
    file.path(output_dir, "hotlist_qc.csv")
  )
  file.copy(asset_path("hotlist_v1.csv"),
            file.path(output_dir, "hotlist_definition.csv"), overwrite = TRUE)
  write_json(
    c(hotlist_metadata, list(
      timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      input_sha256 = if (is.character(segments) && length(segments) == 1L)
        file_sha256(segments) else NA_character_,
      definition_sha256 = file_sha256(asset_path("hotlist_v1.csv")),
      incompatible_sample_count = length(incompatible_samples),
      scoring_rule = "samples with FAIL_INCOMPATIBLE_OVERLAP must be excluded from validation"
    )),
    file.path(output_dir, "hotlist_manifest.json")
  )
  invisible(features)
}

write_candidate_outputs <- function(segments, format, genome_build, output_dir,
                                    mapping_file = NULL, caller_profile = NULL,
                                    preprocess = "legacy", dgv_path = NULL) {
  ensure_dir(output_dir)
  normalized <- normalize_segments(
    segments, format = format, genome_build = genome_build,
    mapping_file = mapping_file, caller_profile = caller_profile,
    allow_incompatible = TRUE
  )
  incompatible_samples <- attr(normalized, "mbpr_incompatible_samples") %||% character()
  features <- extract_breakpoint_candidates(
    normalized, genome_build = genome_build, preprocess = preprocess,
    dgv_path = dgv_path
  )
  candidate_metadata <- attr(features, "menibreak_candidate_matrix")
  features$qc_status <- ifelse(
    features$sample_id %in% incompatible_samples,
    "FAIL_INCOMPATIBLE_OVERLAP", "PASS"
  )
  readr::write_csv(features, file.path(output_dir, "breakpoint_candidates.csv"))
  write_json(
    c(candidate_metadata, list(
      timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      input_sha256 = if (is.character(segments) && length(segments) == 1L)
        file_sha256(segments) else NA_character_,
      candidate_feature_count = sum(startsWith(names(features), "bp_")),
      incompatible_sample_count = length(incompatible_samples),
      outcome_use = "none; selection occurs only in the validation workflow"
    )),
    file.path(output_dir, "candidate_manifest.json")
  )
  invisible(features)
}
