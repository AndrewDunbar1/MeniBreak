fixed_horizon_endpoint <- function(time, event, horizon) {
  ifelse(event == 1 & time <= horizon, 1,
         ifelse(time >= horizon, 0, NA_real_))
}

calibration_stats <- function(y, p) {
  eps <- 1e-6
  lp <- stats::qlogis(pmin(pmax(p, eps), 1 - eps))
  fit <- try(stats::glm(y ~ lp, family = stats::binomial()), silent = TRUE)
  if (inherits(fit, "try-error")) return(c(intercept = NA_real_, slope = NA_real_))
  unname(stats::coef(fit))[1:2] |>
    stats::setNames(c("intercept", "slope"))
}

decision_curve <- function(y, p, thresholds = seq(0.05, 0.50, by = 0.05)) {
  n <- length(y)
  model <- tibble::tibble(
    threshold = thresholds,
    strategy = "model",
    net_benefit = vapply(
      thresholds,
      function(t) {
        pred <- p >= t
        (sum(pred & y == 1) / n) - (sum(pred & y == 0) / n) * t / (1 - t)
      },
      numeric(1)
    )
  )
  treat_all <- tibble::tibble(
    threshold = thresholds,
    strategy = "treat_all",
    net_benefit = mean(y == 1) - mean(y == 0) * thresholds / (1 - thresholds)
  )
  treat_none <- tibble::tibble(
    threshold = thresholds, strategy = "treat_none", net_benefit = 0
  )
  dplyr::bind_rows(model, treat_all, treat_none)
}

event_capture_curve <- function(y, score, fractions = seq(0.05, 1, by = 0.05)) {
  ord <- order(score, decreasing = TRUE)
  y <- y[ord]
  total <- sum(y == 1)
  tibble::tibble(
    highest_risk_fraction = fractions,
    event_capture = vapply(
      fractions,
    function(f) sum(utils::head(y, max(1L, ceiling(length(y) * f))) == 1) / total,
      numeric(1)
    )
  )
}

bootstrap_auc <- function(d, score_col, reps = 1000, seed = 20260703) {
  set.seed(seed)
  group <- d$patient_id %||% d$sample_id
  ids <- unique(group)
  vals <- replicate(reps, {
    draw <- sample(ids, length(ids), replace = TRUE)
    idx <- unlist(lapply(seq_along(draw), function(i) which(group == draw[[i]])))
    dd <- d[idx, , drop = FALSE]
    if (length(unique(dd$endpoint)) < 2L) return(NA_real_)
    auc_safe(dd$endpoint, dd[[score_col]])
  })
  stats::quantile(vals, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}

validate_external <- function(predictions, outcomes, output_dir = NULL,
                              bootstrap_reps = 1000, seed = 20260703) {
  pred <- if (is.character(predictions)) read_table_auto(predictions) else
    tibble::as_tibble(predictions)
  out <- if (is.character(outcomes)) read_table_auto(outcomes) else
    tibble::as_tibble(outcomes)
  required <- c("sample_id", "recurrence_time", "recurrence_event")
  if (!all(required %in% names(out))) {
    mbpr_stop("Outcomes require: ", paste(required, collapse = ", "))
  }
  if (!"patient_id" %in% names(out)) out$patient_id <- out$sample_id
  d <- dplyr::inner_join(pred, out, by = "sample_id")
  if (!nrow(d)) mbpr_stop("Predictions and outcomes have no matching sample_id values.")
  summaries <- list()
  curves <- list()
  for (h in c(3, 5, 7)) {
    score_col <- paste0("raw_score_", h, "y")
    prob_col <- paste0("provisional_probability_", h, "y")
    if (!all(c(score_col, prob_col) %in% names(d))) next
    dh <- d
    dh$endpoint <- fixed_horizon_endpoint(
      dh$recurrence_time, dh$recurrence_event, h
    )
    dh <- dh[!is.na(dh$endpoint), , drop = FALSE]
    if (length(unique(dh$endpoint)) < 2L) next
    cal <- calibration_stats(dh$endpoint, dh[[prob_col]])
    ci <- bootstrap_auc(dh, score_col, bootstrap_reps, seed + h)
    scna_col <- paste0("scna_only_raw_score_", h, "y")
    scna_prob_col <- paste0("scna_only_provisional_probability_", h, "y")
    summaries[[as.character(h)]] <- tibble::tibble(
      horizon_years = h,
      n = nrow(dh),
      events = sum(dh$endpoint),
      auc = auc_safe(dh$endpoint, dh[[score_col]]),
      auc_ci_low = ci[[1]], auc_ci_high = ci[[2]],
      average_precision = average_precision(dh$endpoint, dh[[score_col]]),
      brier = mean((dh[[prob_col]] - dh$endpoint)^2),
      calibration_intercept = cal[["intercept"]],
      calibration_slope = cal[["slope"]],
      scna_only_auc = if (scna_col %in% names(dh))
        auc_safe(dh$endpoint, dh[[scna_col]]) else NA_real_,
      scna_only_brier = if (scna_prob_col %in% names(dh))
        mean((dh[[scna_prob_col]] - dh$endpoint)^2) else NA_real_
    )
    curves[[paste0("decision_", h)]] <- decision_curve(
      dh$endpoint, dh[[prob_col]]
    ) %>% dplyr::mutate(horizon_years = h)
    if (scna_prob_col %in% names(dh)) {
      curves[[paste0("decision_scna_", h)]] <- decision_curve(
        dh$endpoint, dh[[scna_prob_col]]
      ) %>%
        dplyr::filter(.data$strategy == "model") %>%
        dplyr::mutate(strategy = "scna_only", horizon_years = h)
    }
    curves[[paste0("capture_", h)]] <- event_capture_curve(
      dh$endpoint, dh[[score_col]]
    ) %>% dplyr::mutate(model = "breakpoint_augmented", horizon_years = h)
    if (scna_col %in% names(dh)) {
      curves[[paste0("capture_scna_", h)]] <- event_capture_curve(
        dh$endpoint, dh[[scna_col]]
      ) %>% dplyr::mutate(model = "scna_only", horizon_years = h)
    }
  }
  metrics <- dplyr::bind_rows(summaries)
  if (!nrow(metrics)) mbpr_stop("No horizon had both events and non-events.")
  result <- list(
    metrics = metrics,
    decision_curves = dplyr::bind_rows(curves[grepl("^decision", names(curves))]),
    event_capture_curves = dplyr::bind_rows(curves[grepl("^capture", names(curves))])
  )
  if (!is.null(output_dir)) {
    ensure_dir(output_dir)
    readr::write_csv(metrics, file.path(output_dir, "validation_metrics.csv"))
    readr::write_csv(result$decision_curves,
                     file.path(output_dir, "decision_curves.csv"))
    readr::write_csv(result$event_capture_curves,
                     file.path(output_dir, "event_capture_curves.csv"))
    write_json(
      list(timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
           bootstrap_reps = bootstrap_reps, seed = seed),
      file.path(output_dir, "validation_manifest.json")
    )
    write_validation_report(metrics, file.path(output_dir, "report.html"))
  }
  result
}
