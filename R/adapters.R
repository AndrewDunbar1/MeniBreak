.event_map <- c(
  "LOSS" = "loss",
  "CN LOSS" = "loss",
  "COPY LOSS" = "loss",
  "HOMOZYGOUS COPY LOSS" = "loss",
  "HOMOZYGOUS DELETION" = "loss",
  "DEEP DELETION" = "loss",
  "-1" = "loss",
  "-2" = "loss",
  "GAIN" = "gain",
  "CN GAIN" = "gain",
  "COPY GAIN" = "gain",
  "HIGH COPY GAIN" = "gain",
  "AMPLIFICATION" = "gain",
  "AMP" = "gain",
  "1" = "gain",
  "2" = "gain",
  "LOH" = "loh",
  "CNLOH" = "loh",
  "COPY NEUTRAL LOH" = "loh",
  "ALLELIC IMBALANCE" = "allelic_imbalance",
  "AI" = "allelic_imbalance",
  "NEUTRAL" = "neutral",
  "DIPLOID" = "neutral",
  "0" = "neutral"
)

normalize_event <- function(x) {
  key <- toupper(gsub("[_\\-]+", " ", trimws(as.character(x))))
  unname(.event_map[key])
}

detect_segment_format <- function(nms) {
  low <- tolower(nms)
  has <- function(x) all(x %in% low)
  if (has(c("sample_id", "chromosome", "start", "end", "event"))) return("canonical_event")
  if (has(c("cg", "chr", "start", "end", "event_recoded"))) return("v6_processed")
  if (has(c("cg", "chr", "start", "end", "event"))) return("nexus_event")
  if (has(c("id", "chrom", "loc.start", "loc.end", "seg.mean"))) return("cbs_log2")
  if (has(c("sample", "chromosome", "start", "end", "segment_mean"))) return("gistic_seg")
  if (has(c("sample", "contig", "start", "end", "mean_log2_copy_ratio"))) return("gatk_modelsegments")
  "unknown"
}

mapping_for_format <- function(format, nms) {
  low <- stats::setNames(nms, tolower(nms))
  pick <- function(x) unname(low[[tolower(x)]])
  switch(
    format,
    canonical_event = c(
      sample_id = pick("sample_id"), chromosome = pick("chromosome"),
      start = pick("start"), end = pick("end"), event = pick("event")
    ),
    nexus_event = c(
      sample_id = pick("cg"), chromosome = pick("chr"),
      start = pick("start"), end = pick("end"), event = pick("event")
    ),
    v6_processed = c(
      sample_id = pick("cg"), chromosome = pick("chr"),
      start = pick("start"), end = pick("end"), event = pick("event_recoded")
    ),
    cbs_log2 = c(
      sample_id = pick("ID"), chromosome = pick("chrom"),
      start = pick("loc.start"), end = pick("loc.end"), segment_mean = pick("seg.mean")
    ),
    gistic_seg = c(
      sample_id = pick("sample"), chromosome = pick("chromosome"),
      start = pick("start"), end = pick("end"), segment_mean = pick("segment_mean")
    ),
    gatk_modelsegments = c(
      sample_id = pick("sample"), chromosome = pick("contig"),
      start = pick("start"), end = pick("end"),
      segment_mean = pick("mean_log2_copy_ratio")
    ),
    mbpr_stop("No built-in mapping for format: ", format)
  )
}

read_mapping_yaml <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    mbpr_stop("custom_mapping requires an existing mapping YAML.")
  }
  y <- yaml::read_yaml(path)
  mapping <- unlist(y$columns %||% y, use.names = TRUE)
  required <- c("sample_id", "chromosome", "start", "end")
  if (!all(required %in% names(mapping))) {
    mbpr_stop("Mapping YAML must define sample_id, chromosome, start, and end.")
  }
  if (!any(c("event", "segment_mean") %in% names(mapping))) {
    mbpr_stop("Mapping YAML must define either event or segment_mean.")
  }
  mapping
}

read_caller_profile <- function(profile) {
  if (is.null(profile)) {
    mbpr_stop(
      "Continuous segmentation requires --caller-profile with explicit ",
      "loss_threshold and gain_threshold. Thresholds are never guessed."
    )
  }
  path <- if (file.exists(profile)) profile else
    asset_path("caller_profiles", paste0(profile, ".yml"))
  y <- yaml::read_yaml(path)
  required <- c(
    "profile_id", "loss_threshold", "gain_threshold", "validation_status",
    "coordinate_system"
  )
  if (!all(required %in% names(y))) {
    mbpr_stop("Caller profile is missing required fields: ",
              paste(setdiff(required, names(y)), collapse = ", "))
  }
  if (!is.numeric(y$loss_threshold) || !is.numeric(y$gain_threshold) ||
      y$loss_threshold >= y$gain_threshold) {
    mbpr_stop("Caller profile thresholds are invalid.")
  }
  if (!y$coordinate_system %in% c("zero_based_half_open", "one_based_inclusive")) {
    mbpr_stop("Caller profile coordinate_system must be zero_based_half_open or one_based_inclusive.")
  }
  y
}

inspect_segments <- function(segments, format = "auto", genome_build,
                             mapping_file = NULL, caller_profile = NULL) {
  assert_genome_build(genome_build)
  d <- if (is.character(segments) && length(segments) == 1L) read_table_auto(segments) else
    tibble::as_tibble(segments)
  detected <- detect_segment_format(names(d))
  selected <- if (format == "auto") detected else format
  mapping <- if (selected == "custom_mapping") read_mapping_yaml(mapping_file) else
    if (selected == "unknown") NULL else mapping_for_format(selected, names(d))
  continuous <- !is.null(mapping) && "segment_mean" %in% names(mapping)
  profile <- if (continuous && !is.null(caller_profile)) read_caller_profile(caller_profile) else NULL
  result <- list(
    detected_format = detected,
    selected_format = selected,
    genome_build = genome_build,
    rows = nrow(d),
    columns = names(d),
    mapping = mapping,
    continuous = continuous,
    caller_profile = profile,
    score_ready = !is.null(mapping) && (!continuous || !is.null(profile)),
    issues = c(
      if (selected == "unknown") "No known schema detected; provide custom_mapping." else NULL,
      if (continuous && is.null(profile)) "Continuous values require an explicit caller profile." else NULL
    )
  )
  if (isTRUE(result$score_ready)) {
    normalized <- tryCatch(
      normalize_segments(
        d, format = selected, genome_build = genome_build,
        mapping_file = mapping_file, caller_profile = caller_profile,
        allow_incompatible = TRUE
      ),
      error = function(e) e
    )
    if (inherits(normalized, "error")) {
      result$score_ready <- FALSE
      result$issues <- c(result$issues, conditionMessage(normalized))
    } else {
      per_sample <- normalized %>%
        dplyr::group_by(.data$sample_id) %>%
        dplyr::summarise(
          segments = dplyr::n(),
          chromosomes_observed = dplyr::n_distinct(.data$chromosome),
          loss_segments = sum(.data$event == "loss"),
          gain_segments = sum(.data$event == "gain"),
          loh_segments = sum(.data$event == "loh"),
          allelic_imbalance_segments = sum(.data$event == "allelic_imbalance"),
          .groups = "drop"
        )
      result$qc <- list(
        samples = dplyr::n_distinct(normalized$sample_id),
        coordinate_system_input = attr(normalized, "mbpr_coordinate_system_input"),
        event_counts = as.list(table(normalized$event)),
        chromosome_counts = as.list(table(normalized$chromosome)),
        incompatible_overlap_rows = attr(normalized, "mbpr_incompatible_overlaps"),
        incompatible_samples = attr(normalized, "mbpr_incompatible_samples"),
        per_sample = per_sample
      )
      if (attr(normalized, "mbpr_incompatible_overlaps") > 0L) {
        result$issues <- c(
          result$issues,
          paste0(
            "Incompatible loss/gain overlaps affect ",
            length(attr(normalized, "mbpr_incompatible_samples")),
            " sample(s); scoring will emit sample-level QC failures."
          )
        )
      }
    }
  }
  result
}

normalize_segments <- function(segments, format = "auto", genome_build,
                               mapping_file = NULL, caller_profile = NULL,
                               allow_incompatible = FALSE) {
  assert_genome_build(genome_build)
  d <- if (is.character(segments) && length(segments) == 1L) read_table_auto(segments) else
    tibble::as_tibble(segments)
  selected <- if (format == "auto") detect_segment_format(names(d)) else format
  if (selected == "unknown") {
    mbpr_stop("Unknown segment schema. Run mbpr inspect and provide custom_mapping.")
  }
  mapping <- if (selected == "custom_mapping") read_mapping_yaml(mapping_file) else
    mapping_for_format(selected, names(d))
  missing_source <- setdiff(unname(mapping), names(d))
  if (length(missing_source)) {
    mbpr_stop("Input is missing mapped columns: ", paste(missing_source, collapse = ", "))
  }
  canonical <- tibble::tibble(
    sample_id = trimws(as.character(d[[mapping[["sample_id"]]]])),
    chromosome = normalize_chr(d[[mapping[["chromosome"]]]]),
    start = suppressWarnings(as.numeric(d[[mapping[["start"]]]])),
    end = suppressWarnings(as.numeric(d[[mapping[["end"]]]]))
  )
  if (selected %in% c("nexus_event", "v6_processed") &&
      "sample" %in% tolower(names(d))) {
    sample_col <- names(d)[match("sample", tolower(names(d)))]
    fallback <- trimws(as.character(d[[sample_col]]))
    use_fallback <- is.na(canonical$sample_id) | !nzchar(canonical$sample_id)
    canonical$sample_id[use_fallback] <- fallback[use_fallback]
  }
  profile <- NULL
  if ("event" %in% names(mapping)) {
    canonical$event_original <- as.character(d[[mapping[["event"]]]])
    canonical$event <- normalize_event(canonical$event_original)
  } else {
    profile <- read_caller_profile(caller_profile)
    canonical$segment_mean <- suppressWarnings(as.numeric(d[[mapping[["segment_mean"]]]]))
    canonical$event_original <- as.character(canonical$segment_mean)
    canonical$event <- ifelse(
      canonical$segment_mean <= profile$loss_threshold, "loss",
      ifelse(canonical$segment_mean >= profile$gain_threshold, "gain", "neutral")
    )
  }
  coordinate_system <- switch(
    selected,
    cbs_log2 = profile$coordinate_system,
    gistic_seg = profile$coordinate_system,
    gatk_modelsegments = profile$coordinate_system,
    custom_mapping = {
      y <- yaml::read_yaml(mapping_file)
      y$coordinate_system %||% mbpr_stop(
        "custom_mapping YAML must declare coordinate_system."
      )
    },
    canonical_event = "zero_based_half_open",
    nexus_event = "derivation_native_no_conversion",
    v6_processed = "derivation_native_no_conversion",
    "zero_based_half_open"
  )
  if (!coordinate_system %in% c(
    "zero_based_half_open", "one_based_inclusive",
    "derivation_native_no_conversion"
  )) {
    mbpr_stop("Unsupported coordinate_system: ", coordinate_system)
  }
  if (identical(coordinate_system, "one_based_inclusive")) {
    canonical$start <- canonical$start - 1
  }
  invalid <- !nzchar(canonical$sample_id) | is.na(canonical$chromosome) |
    !is.finite(canonical$start) | !is.finite(canonical$end) |
    canonical$start < 0 | canonical$end <= canonical$start
  if (any(invalid)) {
    mbpr_stop("Invalid sample IDs, chromosomes, or coordinates in ",
              sum(invalid), " row(s).")
  }
  noninteger <- abs(canonical$start - round(canonical$start)) > 1e-8 |
    abs(canonical$end - round(canonical$end)) > 1e-8
  if (any(noninteger)) {
    mbpr_stop("Non-integer genomic coordinates in ", sum(noninteger), " row(s).")
  }
  canonical$start <- as.numeric(round(canonical$start))
  canonical$end <- as.numeric(round(canonical$end))
  unknown_events <- is.na(canonical$event)
  if (any(unknown_events)) {
    vals <- unique(canonical$event_original[unknown_events])
    mbpr_stop("Unknown event label(s): ", paste(utils::head(vals, 12L), collapse = ", "))
  }
  arm_ref <- load_arm_reference(genome_build)
  chrom_end <- arm_ref %>%
    dplyr::group_by(.data$chromosome) %>%
    dplyr::summarise(chromosome_end = max(.data$end), .groups = "drop")
  checked <- dplyr::left_join(canonical, chrom_end, by = "chromosome")
  beyond <- !is.na(checked$chromosome_end) & checked$end > checked$chromosome_end
  if (any(beyond)) {
    mbpr_stop("Coordinate(s) exceed the declared ", genome_build,
              " chromosome reference in ", sum(beyond), " row(s).")
  }
  cna <- canonical[canonical$event %in% c("loss", "gain"), , drop = FALSE]
  groups <- split(seq_len(nrow(cna)), interaction(
    cna$sample_id, cna$chromosome, drop = TRUE
  ))
  incompatible <- 0L
  incompatible_samples <- character()
  for (idx in groups) {
    d0 <- cna[idx, , drop = FALSE]
    d0 <- d0[order(d0$start, d0$end), , drop = FALSE]
    if (nrow(d0) < 2L) next
    for (i in 2:nrow(d0)) {
      prior <- seq_len(i - 1L)
      overlap <- d0$end[prior] > d0$start[[i]]
      if (any(overlap & d0$event[prior] != d0$event[[i]])) {
        incompatible <- incompatible + 1L
        incompatible_samples <- c(incompatible_samples, d0$sample_id[[i]])
      }
    }
  }
  if (incompatible > 0L && !isTRUE(allow_incompatible)) {
    mbpr_stop("Overlapping incompatible loss/gain segments were detected.")
  }
  canonical$genome_build <- genome_build
  attr(canonical, "mbpr_adapter") <- selected
  attr(canonical, "mbpr_caller_profile") <- profile
  attr(canonical, "mbpr_coordinate_system_input") <- coordinate_system
  attr(canonical, "mbpr_incompatible_overlaps") <- incompatible
  attr(canonical, "mbpr_incompatible_samples") <- unique(incompatible_samples)
  canonical
}
