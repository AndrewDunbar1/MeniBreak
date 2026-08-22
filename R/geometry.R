safe_ratio <- function(numerator, denominator) {
  out <- numerator / denominator
  out[!is.finite(out)] <- 0
  out
}

row_entropy <- function(probability) {
  z <- probability
  z[z <= 0] <- 1
  -rowSums(ifelse(probability > 0, probability * log(z), 0))
}

geometry_route_block <- function(array, tag, arm_keys) {
  n <- dim(array)[1]
  bins <- seq_len(dim(array)[3])
  arm_counts <- apply(array, c(1, 2), sum)
  bin_counts <- apply(array, c(1, 3), sum)
  if (n == 1L) {
    arm_counts <- matrix(arm_counts, nrow = 1L,
                         dimnames = list(NULL, arm_keys))
    bin_counts <- matrix(bin_counts, nrow = 1L)
  }
  total <- rowSums(arm_counts)
  arm_probability <- safe_ratio(arm_counts, total)
  bin_probability <- safe_ratio(bin_counts, total)
  touched <- arm_counts > 0
  occupied <- array > 0
  arm_entropy <- row_entropy(arm_probability)
  radial_entropy <- row_entropy(bin_probability)
  position <- (bins - 0.5) / length(bins)
  radial_mean <- rowSums(sweep(bin_probability, 2, position, `*`))
  radial_variance <- rowSums(sweep(bin_probability, 2, position^2, `*`)) - radial_mean^2

  clusters <- spans <- isolated <- matrix(0, n, length(arm_keys))
  for (i in seq_len(n)) for (j in seq_along(arm_keys)) {
    z <- occupied[i, j, ]
    if (!any(z)) next
    runs <- rle(z)
    clusters[i, j] <- sum(runs$values)
    index <- which(z)
    spans[i, j] <- (max(index) - min(index) + 1) / length(bins)
    isolated[i, j] <- sum(runs$lengths[runs$values] == 1)
  }
  arm_n <- rowSums(touched)
  occupied_n <- apply(occupied, 1, sum)
  result <- tibble::tibble(
    arm_entropy_norm = safe_ratio(arm_entropy, log(pmax(arm_n, 2))),
    effective_arm_fraction = safe_ratio(exp(arm_entropy), pmax(arm_n, 1)),
    arm_hhi = rowSums(arm_probability^2),
    max_arm_fraction = apply(arm_probability, 1, max),
    p_arm_fraction = rowSums(arm_probability[, grepl("p$", arm_keys), drop = FALSE]),
    radial_mean = radial_mean,
    radial_sd = sqrt(pmax(radial_variance, 0)),
    radial_entropy_norm = radial_entropy / log(length(bins)),
    radial_hhi = rowSums(bin_probability^2),
    telomeric_fraction = rowSums(bin_probability[, 1:4, drop = FALSE]),
    outer_mid_fraction = rowSums(bin_probability[, 5:8, drop = FALSE]),
    middle_fraction = rowSums(bin_probability[, 9:12, drop = FALSE]),
    inner_mid_fraction = rowSums(bin_probability[, 13:16, drop = FALSE]),
    centromeric_fraction = rowSums(bin_probability[, 17:20, drop = FALSE]),
    cluster_per_occupied_bin = safe_ratio(rowSums(clusters), occupied_n),
    isolated_cluster_fraction = safe_ratio(rowSums(isolated), rowSums(clusters)),
    mean_within_arm_span = safe_ratio(rowSums(spans), arm_n)
  )
  names(result) <- paste0(tag, "__", names(result))
  list(features = result, arm_counts = arm_counts,
       bin_probability = bin_probability, total = total, touched = touched)
}

#' Extract the scale-free breakpoint-geometry feature set
#'
#' Produces the exact 56 normalized route-shape features used by the
#' `scale-free-geometry-v1` internal validation recipe. These are fixed feature
#' definitions, not a frozen prediction score; coefficients must be estimated
#' inside an outcome-linked validation workflow.
#'
#' @param segments Canonical output from `normalize_segments()`.
#' @param genome_build Explicit GRCh37 or GRCh38 build.
#' @param preprocess Legacy or V6 cleaning lineage.
#' @param dgv_path Optional DGV asset for V6 preprocessing.
#' @return One row per sample with 56 geometry features.
#' @export
extract_geometry_features <- function(segments, genome_build = NULL,
                                      preprocess = c("legacy", "v6"),
                                      dgv_path = NULL) {
  preprocess <- match.arg(preprocess)
  required <- c("sample_id", "chromosome", "start", "end", "event", "genome_build")
  if (!all(required %in% names(segments))) {
    mbpr_stop("extract_geometry_features requires output from normalize_segments().")
  }
  build <- genome_build %||% unique(segments$genome_build)
  assert_genome_build(build)
  clean <- preprocess_segments(segments, mode = preprocess, dgv_path = dgv_path)
  samples <- unique(clean$sample_id)
  arm_reference <- load_arm_reference(build)
  arm_keys <- arm_reference %>%
    dplyr::transmute(arm_key = paste0(sub("^chr", "", .data$chromosome), .data$arm)) %>%
    dplyr::pull(.data$arm_key) %>% unique()
  boundaries <- make_boundaries(
    clean %>% dplyr::filter(.data$event %in% c("loss", "gain", "allelic_imbalance")),
    arm_reference
  ) %>% dplyr::mutate(bin = pmin(20L, floor(.data$rel_tel_to_cent * 20) + 1L))

  make_array <- function(event) {
    selected <- if (event == "any") boundaries else
      boundaries[boundaries$event == event, , drop = FALSE]
    array_out <- array(0, dim = c(length(samples), length(arm_keys), 20L),
                       dimnames = list(samples, arm_keys, sprintf("%02d", 1:20)))
    if (nrow(selected)) {
      counts <- selected %>% dplyr::count(.data$sample_id, .data$arm_key, .data$bin)
      for (i in seq_len(nrow(counts))) {
        array_out[counts$sample_id[[i]], counts$arm_key[[i]], counts$bin[[i]]] <- counts$n[[i]]
      }
    }
    array_out
  }
  route <- lapply(c("any", "loss", "gain"), function(event) {
    geometry_route_block(make_array(event), event, arm_keys)
  })
  names(route) <- c("any", "loss", "gain")

  js_rows <- function(p, q) {
    midpoint <- (p + q) / 2
    kl <- function(a, b) rowSums(ifelse(a > 0, a * log(a / pmax(b, 1e-12)), 0))
    (kl(p, midpoint) + kl(q, midpoint)) / 2
  }
  loss_arm_probability <- safe_ratio(route$loss$arm_counts, route$loss$total)
  gain_arm_probability <- safe_ratio(route$gain$arm_counts, route$gain$total)
  cross_route <- tibble::tibble(
    loss_boundary_fraction = safe_ratio(route$loss$total, route$loss$total + route$gain$total),
    loss_gain_arm_js = js_rows(loss_arm_probability, gain_arm_probability),
    loss_gain_radial_js = js_rows(route$loss$bin_probability, route$gain$bin_probability),
    loss_gain_arm_jaccard = safe_ratio(
      rowSums(route$loss$touched & route$gain$touched),
      rowSums(route$loss$touched | route$gain$touched)
    ),
    loss_gain_arm_cosine = safe_ratio(
      rowSums(loss_arm_probability * gain_arm_probability),
      sqrt(rowSums(loss_arm_probability^2) * rowSums(gain_arm_probability^2))
    )
  )
  out <- dplyr::bind_cols(
    tibble::tibble(sample_id = samples), route$any$features,
    route$loss$features, route$gain$features, cross_route
  )
  stopifnot(ncol(out) == 57L)
  attr(out, "menibreak_geometry") <- list(
    version = "scale-free-geometry-v1", genome_build = build,
    feature_count = 56L, bins_per_arm = 20L,
    coordinate_system = "arm-relative telomere-to-centromere",
    eligible_events_any = c("loss", "gain", "allelic_imbalance"),
    event_specific_routes = c("loss", "gain")
  )
  out
}

write_geometry_outputs <- function(segments, format, genome_build, output_dir,
                                   mapping_file = NULL, caller_profile = NULL,
                                   preprocess = "legacy", dgv_path = NULL) {
  ensure_dir(output_dir)
  normalized <- normalize_segments(
    segments, format = format, genome_build = genome_build,
    mapping_file = mapping_file, caller_profile = caller_profile,
    allow_incompatible = TRUE
  )
  incompatible <- attr(normalized, "mbpr_incompatible_samples") %||% character()
  features <- extract_geometry_features(normalized, genome_build, preprocess, dgv_path)
  metadata <- attr(features, "menibreak_geometry")
  features$qc_status <- ifelse(features$sample_id %in% incompatible,
                               "FAIL_INCOMPATIBLE_OVERLAP", "PASS")
  readr::write_csv(features, file.path(output_dir, "geometry_features.csv"))
  write_json(c(metadata, list(
    timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    input_sha256 = if (is.character(segments) && length(segments) == 1L)
      file_sha256(segments) else NA_character_,
    incompatible_sample_count = length(incompatible),
    outcome_use = "none; fitting occurs only in the validation workflow"
  )), file.path(output_dir, "geometry_manifest.json"))
  invisible(features)
}
