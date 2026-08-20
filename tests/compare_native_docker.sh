#!/bin/sh
set -eu

if ! command -v docker >/dev/null 2>&1; then
  echo "SKIP: docker is not installed"
  exit 0
fi

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

"$ROOT/mbpr" score \
  --segments "$ROOT/tests/fixtures/canonical_events.csv" \
  --format canonical_event --genome-build GRCh37 \
  --output-dir "$TMP/native"

docker build -t mbpr:test "$ROOT"
docker run --rm \
  -v "$ROOT/tests/fixtures:/input:ro" \
  -v "$TMP/docker:/output" \
  mbpr:test score --segments /input/canonical_events.csv \
  --format canonical_event --genome-build GRCh37 --output-dir /output

Rscript --vanilla - "$TMP/native/predictions.csv" "$TMP/docker/predictions.csv" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
a <- read.csv(args[[1]], check.names = FALSE)
b <- read.csv(args[[2]], check.names = FALSE)
stopifnot(identical(names(a), names(b)), identical(a$sample_id, b$sample_id))
numeric <- names(a)[vapply(a, is.numeric, logical(1))]
delta <- max(abs(as.matrix(a[numeric]) - as.matrix(b[numeric])), na.rm = TRUE)
stopifnot(is.finite(delta), delta <= 1e-12)
cat("Native and Docker predictions match; maximum absolute delta:", delta, "\n")
RS
