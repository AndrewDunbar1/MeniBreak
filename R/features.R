arm_overlap_table <- function(seg_model, arm_ref) {
  seg_model %>%
    dplyr::inner_join(arm_ref, by = "chromosome", relationship = "many-to-many") %>%
    dplyr::mutate(
      overlap_start = pmax(.data$start.x, .data$start.y),
      overlap_end = pmin(.data$end.x, .data$end.y)
    ) %>%
    dplyr::filter(.data$overlap_end > .data$overlap_start)
}

make_arm_calls <- function(seg_model, arm_ref, threshold) {
  ov <- arm_overlap_table(seg_model, arm_ref)
  coverage <- ov %>%
    dplyr::group_by(.data$sample_id, .data$arm_key, .data$event, .data$arm_length) %>%
    dplyr::summarise(
      covered_bp = sum_union_length(.data$overlap_start, .data$overlap_end),
      .groups = "drop"
    ) %>%
    dplyr::mutate(fraction = .data$covered_bp / .data$arm_length)
  samples <- unique(seg_model$sample_id)
  grid <- tidyr::crossing(sample_id = samples, arm_key = allowed_arms)
  passing <- coverage %>%
    dplyr::filter(.data$fraction >= threshold / 100) %>%
    dplyr::group_by(.data$sample_id, .data$arm_key) %>%
    dplyr::summarise(
      has_loss = any(.data$event == "loss"),
      has_gain = any(.data$event == "gain"),
      .groups = "drop"
    )
  grid %>%
    dplyr::left_join(passing, by = c("sample_id", "arm_key")) %>%
    dplyr::mutate(
      has_loss = dplyr::coalesce(.data$has_loss, FALSE),
      has_gain = dplyr::coalesce(.data$has_gain, FALSE),
      mixed_call = .data$has_loss & .data$has_gain,
      cnv = dplyr::case_when(
        .data$has_loss ~ -1,
        .data$has_gain ~ 1,
        TRUE ~ 0
      )
    )
}

make_boundaries <- function(seg_model, arm_ref) {
  seg_model %>%
    dplyr::select("sample_id", "chromosome", "event",
                  segment_start = "start", segment_end = "end") %>%
    tidyr::pivot_longer(c("segment_start", "segment_end"),
                        names_to = "boundary", values_to = "position") %>%
    dplyr::inner_join(arm_ref, by = "chromosome", relationship = "many-to-many") %>%
    dplyr::filter(.data$position >= .data$start, .data$position <= .data$end) %>%
    dplyr::filter(!(.data$arm == "p" & .data$position == .data$start),
                  !(.data$arm == "q" & .data$position == .data$end)) %>%
    dplyr::mutate(
      rel_tel_to_cent = dplyr::if_else(
        .data$arm == "p",
        (.data$position - .data$start) / .data$arm_length,
        (.data$end - .data$position) / .data$arm_length
      ),
      rel_tel_to_cent = pmin(pmax(.data$rel_tel_to_cent, 0), 0.999999)
    )
}

make_breakpoint_aggregates <- function(boundaries, samples) {
  if (!nrow(boundaries)) {
    return(tibble::tibble(
      sample_id = samples,
      total_breakpoint_boundaries = 0, breakpoint_arms_touched = 0,
      breakpoint_route_entropy = 0, max_arm_breakpoint_fraction = 0,
      telomeric_breakpoint_fraction = 0, central_breakpoint_fraction = 0,
      centromeric_breakpoint_fraction = 0, mean_rel_tel_to_cent = 0,
      sd_rel_tel_to_cent = 0
    ))
  }
  arm_counts <- boundaries %>%
    dplyr::count(.data$sample_id, .data$arm_key, name = "n")
  counts <- arm_counts %>%
    dplyr::group_by(.data$sample_id) %>%
    dplyr::summarise(
      total_breakpoint_boundaries = sum(.data$n),
      breakpoint_arms_touched = sum(.data$n > 0),
      breakpoint_route_entropy = entropy_safe(.data$n),
      max_arm_breakpoint_fraction = max(.data$n / sum(.data$n)),
      .groups = "drop"
    )
  positions <- boundaries %>%
    dplyr::group_by(.data$sample_id) %>%
    dplyr::summarise(
      telomeric_breakpoint_fraction = mean(.data$rel_tel_to_cent < 0.25),
      central_breakpoint_fraction = mean(.data$rel_tel_to_cent >= 0.25 &
                                           .data$rel_tel_to_cent < 0.75),
      centromeric_breakpoint_fraction = mean(.data$rel_tel_to_cent >= 0.75),
      mean_rel_tel_to_cent = mean(.data$rel_tel_to_cent),
      sd_rel_tel_to_cent = dplyr::coalesce(stats::sd(.data$rel_tel_to_cent), 0),
      .groups = "drop"
    )
  tibble::tibble(sample_id = samples) %>%
    dplyr::left_join(counts, by = "sample_id") %>%
    dplyr::left_join(positions, by = "sample_id") %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ dplyr::coalesce(.x, 0)))
}

make_hotspot_features <- function(boundaries, arm_calls, samples) {
  hs <- readr::read_csv(asset_path("hotspot_windows_v1.csv"), show_col_types = FALSE)
  routed <- boundaries %>%
    dplyr::inner_join(hs, by = "arm_key") %>%
    dplyr::filter(.data$rel_tel_to_cent >= .data$hotspot_start_norm,
                  .data$rel_tel_to_cent < .data$hotspot_end_norm) %>%
    dplyr::inner_join(
      arm_calls %>% dplyr::filter(.data$has_loss) %>%
        dplyr::select("sample_id", "arm_key"),
      by = c("sample_id", "arm_key")
    )
  if (!nrow(routed)) {
    return(tibble::tibble(
      sample_id = samples,
      dynamic_hotspot_breakpoints_all = 0,
      dynamic_hotspot_breakpoints_non1p = 0,
      dynamic_hotspot_arms_used_all = 0,
      dynamic_hotspot_arms_used_non1p = 0
    ))
  }
  x <- routed %>%
    dplyr::group_by(.data$sample_id) %>%
    dplyr::summarise(
      dynamic_hotspot_breakpoints_all = dplyr::n(),
      dynamic_hotspot_breakpoints_non1p = sum(.data$arm_key != "1p"),
      dynamic_hotspot_arms_used_all = dplyr::n_distinct(.data$arm_key),
      dynamic_hotspot_arms_used_non1p = dplyr::n_distinct(.data$arm_key[.data$arm_key != "1p"]),
      .groups = "drop"
    )
  tibble::tibble(sample_id = samples) %>%
    dplyr::left_join(x, by = "sample_id") %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ dplyr::coalesce(.x, 0)))
}

extract_features <- function(segments, genome_build = NULL,
                             preprocess = c("legacy", "v6"),
                             arm_threshold = 50,
                             dgv_path = NULL) {
  preprocess <- match.arg(preprocess)
  if (!is.numeric(arm_threshold) || length(arm_threshold) != 1L ||
      !arm_threshold %in% c(30, 50)) {
    mbpr_stop("arm_threshold must be 30 or 50.")
  }
  if (!all(c("sample_id", "chromosome", "start", "end", "event") %in% names(segments))) {
    mbpr_stop("extract_features requires output from normalize_segments().")
  }
  build <- genome_build %||% unique(segments$genome_build)
  assert_genome_build(build)
  arm_ref <- load_arm_reference(build)
  clean <- preprocess_segments(segments, mode = preprocess, dgv_path = dgv_path)
  prep_qc <- attr(clean, "mbpr_preprocess_qc")
  samples <- unique(clean$sample_id)
  seg_model <- clean %>% dplyr::filter(.data$event %in% c("loss", "gain"))
  breakpoint_model <- clean %>%
    dplyr::filter(.data$event %in% c("loss", "gain", "allelic_imbalance"))
  if (!nrow(seg_model)) mbpr_stop("No loss or gain segments remain after preprocessing.")
  arm_calls <- make_arm_calls(seg_model, arm_ref, arm_threshold)
  boundaries <- make_boundaries(breakpoint_model, arm_ref)
  arm_wide <- arm_calls %>%
    dplyr::select("sample_id", "arm_key", "cnv") %>%
    tidyr::pivot_wider(names_from = "arm_key", values_from = "cnv",
                       names_prefix = "cnv_", values_fill = 0)
  aggregates <- make_breakpoint_aggregates(boundaries, samples)
  complexity <- breakpoint_model %>%
    dplyr::group_by(.data$sample_id) %>%
    dplyr::summarise(
      cna_segment_count = dplyr::n(),
      loss_segment_count = sum(.data$event == "loss"),
      gain_segment_count = sum(.data$event == "gain"),
      median_segment_length = stats::median(pmax(.data$end - .data$start, 0)),
      small_segment_count_lt5mb = sum((.data$end - .data$start) < 5e6),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      clean %>%
        dplyr::group_by(.data$sample_id) %>%
        dplyr::summarise(
          loh_excluded_segment_count = sum(.data$event %in% c("loh", "allelic_imbalance")),
          .groups = "drop"
        ),
      by = "sample_id"
    )
  hotspot <- make_hotspot_features(boundaries, arm_calls, samples)
  features <- tibble::tibble(sample_id = samples) %>%
    dplyr::left_join(arm_wide, by = "sample_id") %>%
    dplyr::left_join(aggregates, by = "sample_id") %>%
    dplyr::left_join(complexity, by = "sample_id") %>%
    dplyr::left_join(hotspot, by = "sample_id") %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ dplyr::coalesce(.x, 0))) %>%
    dplyr::mutate(
      scna_loss_count = rowSums(dplyr::across(dplyr::starts_with("cnv_"), ~ .x == -1)),
      scna_gain_count = rowSums(dplyr::across(dplyr::starts_with("cnv_"), ~ .x == 1)),
      scna_altered_arm_count = rowSums(dplyr::across(dplyr::starts_with("cnv_"), ~ .x != 0)),
      cnv_1p_loss = as.integer(.data$cnv_1p == -1),
      cnv_22q_loss = as.integer(.data$cnv_22q == -1)
    )
  per_sample_qc <- tibble::tibble(sample_id = samples) %>%
    dplyr::left_join(
      clean %>% dplyr::count(.data$sample_id, name = "segments_after_preprocessing"),
      by = "sample_id"
    ) %>%
    dplyr::left_join(
      seg_model %>%
        dplyr::filter(.data$chromosome %in% arm_ref$chromosome) %>%
        dplyr::count(.data$sample_id, name = "autosomal_model_segments"),
      by = "sample_id"
    ) %>%
    dplyr::left_join(
      clean %>%
        dplyr::group_by(.data$sample_id) %>%
        dplyr::summarise(
          chromosomes_with_events = dplyr::n_distinct(.data$chromosome),
          loss_segments = sum(.data$event == "loss"),
          gain_segments = sum(.data$event == "gain"),
          loh_segments = sum(.data$event == "loh"),
          allelic_imbalance_segments = sum(.data$event == "allelic_imbalance"),
          .groups = "drop"
        ),
      by = "sample_id"
    ) %>%
    dplyr::left_join(
      boundaries %>% dplyr::count(.data$sample_id, name = "modeled_breakpoint_boundaries"),
      by = "sample_id"
    ) %>%
    dplyr::left_join(
      arm_calls %>% dplyr::group_by(.data$sample_id) %>%
        dplyr::summarise(
          mixed_arm_calls = sum(.data$mixed_call),
          arms_called = sum(.data$cnv != 0),
          .groups = "drop"
        ),
      by = "sample_id"
    ) %>%
    dplyr::mutate(
      dplyr::across(dplyr::where(is.numeric), ~ dplyr::coalesce(.x, 0)),
      qc_status = dplyr::case_when(
        .data$autosomal_model_segments == 0 ~ "FAIL_NO_MODELED_AUTOSOMES",
        .data$modeled_breakpoint_boundaries == 0 ~ "WARN_NO_INTERNAL_BREAKPOINTS",
        TRUE ~ "PASS"
      )
    )
  attr(features, "mbpr_feature_qc") <- list(
    genome_build = build,
    preprocessing = prep_qc,
    arm_threshold = arm_threshold,
    samples = nrow(features),
    per_sample = per_sample_qc
  )
  features
}
