parse_cli_options <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (!startsWith(token, "--")) mbpr_stop("Unexpected argument: ", token)
    key <- gsub("-", "_", substring(token, 3), fixed = TRUE)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

require_cli <- function(options, names) {
  missing <- names[!vapply(names, function(x) !is.null(options[[x]]), logical(1))]
  if (length(missing)) {
    mbpr_stop("Missing required option(s): --",
              paste(gsub("_", "-", missing), collapse = ", --"))
  }
}

write_schema_templates <- function(output_dir = ".") {
  ensure_dir(output_dir)
  assets <- c(
    "canonical_segments_template.csv", "outcomes_template.csv",
    "custom_mapping_template.yml"
  )
  for (a in assets) file.copy(asset_path(a), file.path(output_dir, a), overwrite = TRUE)
  file.copy(
    asset_path("caller_profiles", "example_custom.yml"),
    file.path(output_dir, "caller_profile_template.yml"),
    overwrite = TRUE
  )
  invisible(normalizePath(output_dir))
}

cli_usage <- function() {
  paste(
    "mbpr: deterministic SEG-to-recurrence research software",
    "",
    "Commands:",
    "  top8-candidates --segments FILE --format FORMAT --genome-build BUILD --output-dir DIR",
    "  hotlist-extract --segments FILE --format FORMAT --genome-build BUILD --output-dir DIR",
    "  inspect --segments FILE --format auto --genome-build BUILD [--mapping-file YAML] [--caller-profile YAML]",
    "  score --segments FILE --format FORMAT --genome-build BUILD --output-dir DIR [--mapping-file YAML] [--caller-profile YAML]",
    "  validate --predictions FILE --outcomes FILE --output-dir DIR [--bootstrap-reps N]",
    "  model-info",
    "  schema-template [--output-dir DIR]",
    sep = "\n"
  )
}

mbpr_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args) || args[[1]] %in% c("-h", "--help", "help")) {
    cat(cli_usage(), "\n")
    return(invisible(0L))
  }
  command <- args[[1]]
  o <- parse_cli_options(args[-1])
  if (command == "top8-candidates") {
    require_cli(o, c("segments", "genome_build", "output_dir"))
    write_candidate_outputs(
      o$segments, format = o$format %||% "auto",
      genome_build = o$genome_build, output_dir = o$output_dir,
      mapping_file = o$mapping_file, caller_profile = o$caller_profile,
      preprocess = o$preprocess %||% "legacy", dgv_path = o$dgv_path
    )
    cat("Wrote dynamic top-eight candidate features to", normalizePath(o$output_dir), "\n")
  } else if (command == "hotlist-extract") {
    require_cli(o, c("segments", "genome_build", "output_dir"))
    write_hotlist_outputs(
      o$segments, format = o$format %||% "auto",
      genome_build = o$genome_build, output_dir = o$output_dir,
      mapping_file = o$mapping_file, caller_profile = o$caller_profile,
      preprocess = o$preprocess %||% "legacy", dgv_path = o$dgv_path
    )
    cat("Wrote fixed nine-window hotlist features to", normalizePath(o$output_dir), "\n")
  } else if (command == "inspect") {
    require_cli(o, c("segments", "genome_build"))
    result <- inspect_segments(
      o$segments, format = o$format %||% "auto",
      genome_build = o$genome_build, mapping_file = o$mapping_file,
      caller_profile = o$caller_profile
    )
    if (identical(result$detected_format, "unknown")) {
      destination <- o$mapping_template %||% "custom_mapping_template.yml"
      file.copy(asset_path("custom_mapping_template.yml"), destination, overwrite = TRUE)
      result$mapping_template <- normalizePath(destination)
    }
    cat(jsonlite::toJSON(result, auto_unbox = TRUE, pretty = TRUE), "\n")
  } else if (command == "score") {
    require_cli(o, c("segments", "genome_build", "output_dir"))
    score_segments(
      o$segments, format = o$format %||% "auto",
      genome_build = o$genome_build, output_dir = o$output_dir,
      mapping_file = o$mapping_file, caller_profile = o$caller_profile,
      features_format = o$features_format %||% "csv"
    )
    cat("Wrote scoring outputs to", normalizePath(o$output_dir), "\n")
  } else if (command == "validate") {
    require_cli(o, c("predictions", "outcomes", "output_dir"))
    result <- validate_external(
      o$predictions, o$outcomes, o$output_dir,
      bootstrap_reps = as.integer(o$bootstrap_reps %||% 1000)
    )
    print(result$metrics)
  } else if (command == "model-info") {
    cat(model_info(as_json = TRUE), "\n")
  } else if (command == "schema-template") {
    path <- write_schema_templates(o$output_dir %||% ".")
    cat("Wrote schema templates to", path, "\n")
  } else {
    mbpr_stop("Unknown command: ", command, "\n\n", cli_usage())
  }
  invisible(0L)
}
