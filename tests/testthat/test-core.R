test_that("canonical event vocabulary and features are deterministic", {
  x <- normalize_segments(
    fixture("canonical_events.csv"), genome_build = "GRCh37"
  )
  expect_equal(
    x$event,
    c("loss", "gain", "loss", "loh", "loss", "gain", "allelic_imbalance")
  )
  a <- extract_features(x, preprocess = "legacy", arm_threshold = 50)
  b <- extract_features(x, preprocess = "legacy", arm_threshold = 50)
  expect_identical(a, b)
  expect_true(all(c("cnv_1p", "cnv_22q", "scna_altered_arm_count",
                    "dynamic_hotspot_breakpoints_all") %in% names(a)))
})

test_that("continuous schemas are blocked without an explicit profile", {
  x <- data.frame(
    ID = "A", chrom = 1, loc.start = 1, loc.end = 100000,
    seg.mean = -0.3
  )
  expect_error(
    normalize_segments(x, format = "cbs_log2", genome_build = "GRCh37"),
    "caller-profile|caller profile|Continuous"
  )
})

test_that("unknown events and incompatible overlaps fail explicitly", {
  unknown <- data.frame(
    sample_id = "A", chromosome = 1, start = 1, end = 100,
    event = "mystery"
  )
  expect_error(normalize_segments(unknown, genome_build = "GRCh37"),
               "Unknown event")
  overlap <- data.frame(
    sample_id = c("A", "A"), chromosome = c(1, 1),
    start = c(1, 50), end = c(100, 150), event = c("loss", "gain")
  )
  expect_error(normalize_segments(overlap, genome_build = "GRCh37"),
               "incompatible")
})

test_that("locked model checksum is enforced and scoring is deterministic", {
  expect_equal(model_info()$status, "LOCKED")
  a <- score_segments(fixture("canonical_events.csv"), genome_build = "GRCh37")
  b <- score_segments(fixture("canonical_events.csv"), genome_build = "GRCh37")
  expect_identical(a, b)
  expect_true(all(is.finite(a$raw_score_5y)))
  expect_true(all(a$provisional_probability_5y > 0 &
                    a$provisional_probability_5y < 1))
  expected <- read.csv(fixture("golden_predictions.csv"), check.names = FALSE)
  observed <- a[, names(expected)]
  expect_equal(as.data.frame(observed), expected, tolerance = 1e-12)
})

test_that("incompatible samples fail without blocking clean samples", {
  d <- data.frame(
    sample_id = c("BAD", "BAD", "GOOD"),
    chromosome = c(1, 1, 22),
    start = c(1, 50, 16000000),
    end = c(100, 150, 50000000),
    event = c("loss", "gain", "loss")
  )
  p <- score_segments(d, genome_build = "GRCh37")
  expect_equal(p$qc_status[p$sample_id == "BAD"], "FAIL_INCOMPATIBLE_OVERLAP")
  expect_true(is.na(p$raw_score_5y[p$sample_id == "BAD"]))
  expect_true(is.finite(p$raw_score_5y[p$sample_id == "GOOD"]))
})
