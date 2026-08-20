#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (requireNamespace("mbpr", quietly = TRUE)) {
  mbpr:::mbpr_cli(args)
} else {
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "scripts/mbpr.R"
  root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
  Sys.setenv(MBPR_HOME = root)
  for (f in list.files(file.path(root, "R"), full.names = TRUE)) source(f)
  mbpr_cli(args)
}
