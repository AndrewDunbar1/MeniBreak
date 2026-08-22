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

test_that("fixed hotlist uses 5-percent arm-relative boundary counts", {
  d <- data.frame(
    sample_id = c("A", "A", "B"), chromosome = 1,
    start = c(112500000, 118750000, 1000000),
    end = c(118749999, 124999999, 2000000), event = "loss"
  )
  x <- normalize_segments(d, genome_build = "GRCh37")
  h <- extract_hotlist_features(x)
  expect_true(all(paste0(c("hotlist_1p_090_095", "hotlist_1p_095_100"), "_count") %in% names(h)))
  expect_gt(h$hotlist_1p_090_095_count[h$sample_id == "A"], 0)
  expect_equal(h$hotlist_1p_090_095_present[h$sample_id == "B"], 0)
})

test_that("dynamic top-eight candidates are binary 20-bin arm features", {
  x <- normalize_segments(fixture("canonical_events.csv"),
                          genome_build = "GRCh37")
  z <- extract_breakpoint_candidates(x)
  candidate_names <- grep("^bp_.+_bin[0-9]{2}$", names(z), value = TRUE)
  expect_gt(length(candidate_names), 700)
  expect_true(all(unlist(z[candidate_names], use.names = FALSE) %in% c(0L, 1L)))
  expect_true(all(c("bp_1p_bin19", "bp_1p_bin20") %in% candidate_names))
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
