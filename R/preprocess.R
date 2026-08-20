allowed_arms <- c(
  "1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
  "7p","7q","8p","8q","9p","9q","10p","10q","11p","11q","12p","12q",
  "13q","14q","15q","16p","16q","17p","17q","18p","18q","19p","19q",
  "20p","20q","21q","22q"
)

load_arm_reference <- function(genome_build) {
  assert_genome_build(genome_build)
  path <- asset_path(paste0("arm_reference_", genome_build, ".csv"))
  ref <- readr::read_csv(path, show_col_types = FALSE)
  required <- c("chromosome", "arm", "arm_length", "start", "end")
  if (!all(required %in% names(ref))) mbpr_stop("Invalid arm reference: ", path)
  ref %>%
    dplyr::mutate(
      chromosome = normalize_chr(.data$chromosome),
      arm = tolower(as.character(.data$arm)),
      arm_key = paste0(sub("^chr", "", .data$chromosome), .data$arm)
    ) %>%
    dplyr::filter(.data$arm_key %in% allowed_arms)
}

load_dgv_reference <- function(path, genome_build) {
  if (is.null(path)) return(NULL)
  if (!file.exists(path)) mbpr_stop("DGV reference does not exist: ", path)
  ext <- tolower(tools::file_ext(path))
  d <- if (ext %in% c("xlsx", "xlsm", "xls")) {
    openxlsx::read.xlsx(path)
  } else {
    read_table_auto(path)
  }
  low <- stats::setNames(names(d), tolower(names(d)))
  required <- c("chr", "start", "end", "type")
  if (!all(required %in% names(low))) {
    mbpr_stop("DGV reference must contain chr, start, end, and type.")
  }
  out <- tibble::tibble(
    chromosome = normalize_chr(d[[low[["chr"]]]]),
    start = as.numeric(d[[low[["start"]]]]),
    end = as.numeric(d[[low[["end"]]]]),
    event = normalize_event(d[[low[["type"]]]]),
    genome_build = genome_build
  ) %>%
    dplyr::filter(.data$event %in% c("loss", "gain"))
  out
}

filter_dgv_containment <- function(seg, dgv) {
  if (is.null(dgv) || !nrow(dgv)) {
    attr(seg, "dgv_removed") <- 0L
    return(seg)
  }
  x <- data.table::as.data.table(seg)
  x[, mbpr_row_id := .I]
  y <- data.table::as.data.table(dgv)
  x <- x[event %in% c("loss", "gain")]
  y <- y[event %in% c("loss", "gain")]
  data.table::setkeyv(y, c("chromosome", "event", "start", "end"))
  hits <- data.table::foverlaps(
    x, y,
    by.x = c("chromosome", "event", "start", "end"),
    by.y = c("chromosome", "event", "start", "end"),
    type = "within", nomatch = 0L
  )
  remove <- seq_len(nrow(seg)) %in% unique(hits$mbpr_row_id)
  out <- seg[!remove, , drop = FALSE]
  attr(out, "dgv_removed") <- sum(remove)
  out
}

merge_same_event_intervals <- function(seg) {
  split_key <- interaction(seg$sample_id, seg$chromosome, seg$event, drop = TRUE)
  pieces <- split(seg, split_key)
  merged <- lapply(pieces, function(d) {
    d <- d[order(d$start, d$end), , drop = FALSE]
    if (!nrow(d)) return(d)
    out <- list()
    cur <- d[1, , drop = FALSE]
    if (nrow(d) > 1L) {
      for (i in 2:nrow(d)) {
        if (d$start[[i]] <= cur$end[[1]] + 1) {
          cur$end[[1]] <- max(cur$end[[1]], d$end[[i]])
        } else {
          out[[length(out) + 1L]] <- cur
          cur <- d[i, , drop = FALSE]
        }
      }
    }
    out[[length(out) + 1L]] <- cur
    dplyr::bind_rows(out)
  })
  dplyr::bind_rows(merged) %>%
    dplyr::arrange(.data$sample_id, .data$chromosome, .data$start, .data$end)
}

preprocess_segments <- function(seg, mode = c("legacy", "v6"),
                                dgv_path = NULL) {
  mode <- match.arg(mode)
  if (!all(c("sample_id", "chromosome", "start", "end", "event", "genome_build") %in%
           names(seg))) {
    mbpr_stop("preprocess_segments requires canonical normalized segments.")
  }
  input_rows <- nrow(seg)
  ai_n <- sum(seg$event == "allelic_imbalance")
  loh_n <- sum(seg$event == "loh")
  neutral_n <- sum(seg$event == "neutral")
  out <- seg
  dgv_removed <- 0L
  if (mode == "v6") {
    out <- out %>%
      dplyr::filter(.data$event != "allelic_imbalance") %>%
      dplyr::filter(!(.data$chromosome == "chr21" & .data$end <= 13200000)) %>%
      dplyr::distinct(dplyr::across(dplyr::all_of(
        c("sample_id", "chromosome", "start", "end", "event")
      )), .keep_all = TRUE)
    dgv <- load_dgv_reference(dgv_path, unique(out$genome_build))
    out <- filter_dgv_containment(out, dgv)
    dgv_removed <- attr(out, "dgv_removed") %||% 0L
    out <- merge_same_event_intervals(out)
  }
  attr(out, "mbpr_preprocess_qc") <- list(
    mode = mode,
    input_rows = input_rows,
    output_rows = nrow(out),
    allelic_imbalance_removed = if (mode == "v6") ai_n else 0L,
    loh_rows_retained_for_qc = loh_n,
    neutral_rows_retained_for_qc = neutral_n,
    dgv_removed = dgv_removed
  )
  out
}
