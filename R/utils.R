.mbpr_env <- new.env(parent = emptyenv())

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}

`%>%` <- dplyr::`%>%`

mbpr_stop <- function(..., call. = FALSE) {
  stop(paste0(...), call. = call.)
}

package_root <- function() {
  installed <- system.file(package = "mbpr")
  if (nzchar(installed)) return(installed)
  env_root <- Sys.getenv("MBPR_HOME", unset = "")
  if (nzchar(env_root)) return(normalizePath(env_root, mustWork = TRUE))
  candidates <- c(getwd(), file.path(getwd(), "software", "mbpr"))
  hit <- candidates[file.exists(file.path(candidates, "DESCRIPTION"))]
  if (!length(hit)) mbpr_stop("Cannot locate mbpr package root. Set MBPR_HOME.")
  normalizePath(hit[[1]], mustWork = TRUE)
}

asset_path <- function(...) {
  root <- package_root()
  installed_path <- file.path(root, "extdata", ...)
  source_path <- file.path(root, "inst", "extdata", ...)
  if (file.exists(installed_path)) return(installed_path)
  if (file.exists(source_path)) return(source_path)
  mbpr_stop("Missing packaged asset: ", paste(..., collapse = "/"))
}

model_path <- function(...) {
  root <- package_root()
  installed_path <- file.path(root, "models", ...)
  source_path <- file.path(root, "inst", "models", ...)
  if (file.exists(installed_path)) return(installed_path)
  if (file.exists(source_path)) return(source_path)
  mbpr_stop("Missing model asset: ", paste(..., collapse = "/"))
}

model_asset_path <- model_path

read_table_auto <- function(path) {
  if (!file.exists(path)) mbpr_stop("Input file does not exist: ", path)
  out <- tryCatch(
    data.table::fread(path, data.table = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(out) || ncol(out) < 2L) {
    mbpr_stop("Could not parse a delimited table from: ", path)
  }
  tibble::as_tibble(out)
}

normalize_chr <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("^CHR", "", x)
  x[x == "23"] <- "X"
  x[x == "24"] <- "Y"
  ok <- x %in% c(as.character(1:22), "X", "Y")
  x[!ok] <- NA_character_
  ifelse(is.na(x), NA_character_, paste0("chr", x))
}

assert_genome_build <- function(build) {
  allowed <- c("GRCh37", "GRCh38")
  if (length(build) != 1L || !build %in% allowed) {
    mbpr_stop("genome_build must be explicitly set to GRCh37 or GRCh38.")
  }
  build
}

file_sha256 <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

sha256_file <- file_sha256

safe_dir_create <- function(path) {
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    mbpr_stop("Unable to create output directory: ", path)
  }
  normalizePath(path, mustWork = TRUE)
}

ensure_dir <- safe_dir_create

write_json_pretty <- function(x, path) {
  jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, na = "null")
  invisible(path)
}

write_json <- write_json_pretty

package_version_safe <- function() {
  installed <- tryCatch(
    as.character(utils::packageVersion("mbpr")),
    error = function(e) NULL
  )
  if (!is.null(installed)) return(installed)
  description <- read.dcf(file.path(package_root(), "DESCRIPTION"))
  unname(description[1, "Version"])
}

entropy_safe <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (!length(x)) return(0)
  p <- x / sum(x)
  -sum(p * log(p))
}

sum_union_length <- function(starts, ends) {
  keep <- is.finite(starts) & is.finite(ends) & ends > starts
  starts <- starts[keep]
  ends <- ends[keep]
  if (!length(starts)) return(0)
  ord <- order(starts, ends)
  starts <- starts[ord]
  ends <- ends[ord]
  total <- 0
  cur_start <- starts[[1]]
  cur_end <- ends[[1]]
  if (length(starts) > 1L) {
    for (i in 2:length(starts)) {
      if (starts[[i]] <= cur_end) {
        cur_end <- max(cur_end, ends[[i]])
      } else {
        total <- total + cur_end - cur_start
        cur_start <- starts[[i]]
        cur_end <- ends[[i]]
      }
    }
  }
  total + cur_end - cur_start
}

average_precision <- function(y, score) {
  keep <- is.finite(y) & is.finite(score)
  y <- y[keep]
  score <- score[keep]
  if (!any(y == 1) || !any(y == 0)) return(NA_real_)
  ord <- order(score, decreasing = TRUE)
  y <- y[ord]
  precision <- cumsum(y == 1) / seq_along(y)
  sum(precision[y == 1]) / sum(y == 1)
}

auc_safe <- function(y, score) {
  keep <- is.finite(y) & is.finite(score)
  y <- y[keep]
  score <- score[keep]
  if (length(unique(y)) < 2L || length(unique(score)) < 2L) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE, direction = "<")))
}
