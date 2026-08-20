html_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

write_score_report <- function(predictions, qc, path) {
  p5 <- predictions$provisional_probability_5y
  finite_p5 <- p5[is.finite(p5)]
  breaks <- seq(0, 1, by = 0.1)
  histogram <- if (length(finite_p5)) {
    graphics::hist(finite_p5, breaks = breaks, plot = FALSE, include.lowest = TRUE)
  } else {
    list(counts = rep(0L, length(breaks) - 1L))
  }
  max_count <- max(histogram$counts, 1)
  bars <- paste0(
    "<div style='display:flex;align-items:center;gap:8px;margin:3px 0'>",
    "<span style='width:80px'>", sprintf("%.1f-%.1f", breaks[-length(breaks)], breaks[-1]),
    "</span><span style='display:inline-block;background:#3C5488;height:14px;width:",
    round(360 * histogram$counts / max_count), "px'></span><span>",
    histogram$counts, "</span></div>",
    collapse = "\n"
  )
  rows <- paste0(
    "<tr><td>", html_escape(predictions$sample_id), "</td><td>",
    html_escape(predictions$qc_status), "</td><td>",
    sprintf("%.3f", predictions$raw_score_5y), "</td><td>",
    sprintf("%.1f", predictions$reference_percentile_5y), "</td></tr>",
    collapse = "\n"
  )
  html <- paste0(
    "<!doctype html><html><head><meta charset='utf-8'><title>mbpr report</title>",
    "<style>body{font-family:Arial,sans-serif;max-width:1000px;margin:40px auto;",
    "color:#222}table{border-collapse:collapse;width:100%}th,td{padding:8px;",
    "border-bottom:1px solid #ddd;text-align:left}.note{background:#f5f5f5;",
    "padding:12px}</style></head><body><h1>mbpr scoring report</h1>",
    "<p class='note'>Research use only. Probabilities are internally calibrated ",
    "and provisional; no treatment labels are emitted.</p>",
    "<p>Samples: ", nrow(predictions), " &middot; median provisional 5-year ",
    "probability: ", sprintf("%.3f", stats::median(p5, na.rm = TRUE)), "</p>",
    "<h2>Five-year prediction distribution</h2>", bars,
    "<table><thead><tr><th>Sample</th><th>QC</th><th>5-year raw score</th>",
    "<th>Reference percentile</th></tr></thead><tbody>", rows,
    "</tbody></table></body></html>"
  )
  writeLines(html, path, useBytes = TRUE)
  invisible(path)
}

write_validation_report <- function(metrics, path) {
  rows <- paste0(
    "<tr><td>", metrics$horizon_years, "</td><td>", metrics$n,
    "</td><td>", metrics$events, "</td><td>", sprintf("%.3f", metrics$auc),
    "</td><td>", sprintf("%.3f-%.3f", metrics$auc_ci_low, metrics$auc_ci_high),
    "</td><td>", sprintf("%.3f", metrics$average_precision),
    "</td><td>", sprintf("%.3f", metrics$brier), "</td><td>",
    sprintf("%.3f", metrics$calibration_slope), "</td></tr>",
    collapse = "\n"
  )
  html <- paste0(
    "<!doctype html><html><head><meta charset='utf-8'><title>mbpr validation</title>",
    "<style>body{font-family:Arial,sans-serif;max-width:1000px;margin:40px auto;",
    "color:#222}table{border-collapse:collapse;width:100%}th,td{padding:8px;",
    "border-bottom:1px solid #ddd;text-align:left}.note{background:#f5f5f5;",
    "padding:12px}</style></head><body><h1>External-cohort validation</h1>",
    "<p class='note'>Research use only. Fixed-horizon analyses exclude earlier ",
    "censoring and use patient-clustered bootstrap intervals.</p>",
    "<table><thead><tr><th>Horizon</th><th>N</th><th>Events</th><th>AUC</th>",
    "<th>95% CI</th><th>Average precision</th><th>Brier</th>",
    "<th>Calibration slope</th></tr></thead><tbody>", rows,
    "</tbody></table></body></html>"
  )
  writeLines(html, path, useBytes = TRUE)
  invisible(path)
}
