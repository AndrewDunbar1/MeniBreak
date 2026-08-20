#!/usr/bin/env sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)

# tests/fixtures contains documented synthetic golden data. It is the only
# allow-listed location for a prediction-named file.
bad_files=$(find "$ROOT" -type f ! -path "$ROOT/tests/fixtures/*" \( -iname '*.xlsx' -o -iname '*.xls' -o -iname '*.idat' -o -iname '*.seg' -o -iname '*predictions*.csv' -o -iname '*oof*.csv' -o -iname '*fold_assignments*.csv' \) -print)
if [ -n "$bad_files" ]; then
  echo "Blocked potentially sensitive/sample-level files:"
  echo "$bad_files"
  exit 1
fi

if grep -RInE '/Users/[^/]+|/home/[^/]+|MRN|medical.record|date.of.birth|DOB|patient.name|@partners.org|@mgh.harvard.edu' "$ROOT" --exclude-dir=.git --exclude='audit_no_phi.sh'; then
  echo "Blocked possible local path, identifier, or institutional contact string."
  exit 1
fi

echo "Automated PHI/path audit passed. Manual review is still required."
