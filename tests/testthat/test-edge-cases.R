test_that("GRCh37 and GRCh38 references produce deterministic arm calls", {
  d37 <- data.frame(
    sample_id = "A", chromosome = 1, start = 1, end = 70000000,
    event = "loss"
  )
  d38 <- data.frame(
    sample_id = "A", chromosome = 1, start = 1, end = 70000000,
    event = "loss"
  )
  f37 <- extract_features(
    normalize_segments(d37, genome_build = "GRCh37"),
    preprocess = "legacy", arm_threshold = 50
  )
  f38 <- extract_features(
    normalize_segments(d38, genome_build = "GRCh38"),
    preprocess = "legacy", arm_threshold = 50
  )
  expect_equal(f37$cnv_1p, -1)
  expect_equal(f38$cnv_1p, -1)
})

test_that("explicit continuous caller profiles are honored and labeled", {
  d <- data.frame(
    ID = c("A", "A", "A"), chrom = c(1, 2, 3),
    loc.start = c(1, 1, 1), loc.end = c(1000000, 1000000, 1000000),
    seg.mean = c(-0.3, 0, 0.3)
  )
  x <- normalize_segments(
    d, format = "cbs_log2", genome_build = "GRCh37",
    caller_profile = mbpr:::asset_path("caller_profiles", "example_custom.yml")
  )
  expect_equal(x$event, c("loss", "neutral", "gain"))
  expect_equal(x$start, c(0, 0, 0))
  expect_equal(attr(x, "mbpr_caller_profile")$validation_status, "unvalidated")
})

test_that("V6 removes 21p and contained DGV calls", {
  d <- data.frame(
    sample_id = c("A", "A", "A"), chromosome = c(21, 1, 1),
    start = c(1, 200000, 3000000), end = c(10000000, 1000000, 4000000),
    event = c("loss", "loss", "gain")
  )
  x <- normalize_segments(d, genome_build = "GRCh37")
  p <- mbpr:::preprocess_segments(
    x, mode = "v6", dgv_path = fixture("dgv_small.csv")
  )
  expect_false(any(p$chromosome == "chr21" & p$end <= 13200000))
  expect_false(any(p$chromosome == "chr1" & p$start == 200000))
})

test_that("arm-spanning events are apportioned by overlap", {
  d <- data.frame(
    sample_id = "A", chromosome = 1, start = 100000000,
    end = 150000000, event = "loss"
  )
  f <- extract_features(
    normalize_segments(d, genome_build = "GRCh37"),
    preprocess = "legacy", arm_threshold = 30
  )
  expect_equal(f$cnv_1p, 0)
  expect_equal(f$cnv_1q, 0)
})

test_that("patient-clustered validation accepts repeated tumors", {
  set.seed(7)
  pred <- data.frame(
    sample_id = paste0("S", 1:30),
    raw_score_3y = rnorm(30), raw_score_5y = rnorm(30),
    raw_score_7y = rnorm(30),
    provisional_probability_3y = runif(30, 0.1, 0.8),
    provisional_probability_5y = runif(30, 0.1, 0.8),
    provisional_probability_7y = runif(30, 0.1, 0.8)
  )
  out <- data.frame(
    sample_id = paste0("S", 1:30),
    patient_id = rep(paste0("P", 1:15), each = 2),
    recurrence_time = runif(30, 1, 9),
    recurrence_event = rep(c(0, 1), 15)
  )
  result <- validate_external(pred, out, bootstrap_reps = 10)
  expect_true(nrow(result$metrics) >= 1)
})
