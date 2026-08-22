if (requireNamespace("mbpr", quietly = TRUE)) {
  library(mbpr)
} else {
  Sys.setenv(MBPR_HOME = normalizePath("."))
  for (f in list.files("R", full.names = TRUE)) source(f)
}

fixture_path <- system.file(
  "extdata", "canonical_segments_template.csv", package = "mbpr"
)
if (!nzchar(fixture_path)) fixture_path <- "inst/extdata/canonical_segments_template.csv"
x <- normalize_segments(fixture_path, genome_build = "GRCh37")
stopifnot(nrow(x) == 2L)
f1 <- extract_features(x, preprocess = "legacy", arm_threshold = 50)
f2 <- extract_features(x, preprocess = "legacy", arm_threshold = 50)
stopifnot(identical(f1, f2), nrow(f1) == 1L)
g1 <- extract_geometry_features(x, preprocess = "legacy")
g2 <- extract_geometry_features(x, preprocess = "legacy")
stopifnot(identical(g1, g2), ncol(g1) == 57L)
stopifnot(identical(model_info()$status, "LOCKED"))
predictions <- score_segments(fixture_path, genome_build = "GRCh37")
stopifnot(nrow(predictions) == 1L, is.finite(predictions$raw_score_5y))
cat("mbpr smoke tests passed\n")
