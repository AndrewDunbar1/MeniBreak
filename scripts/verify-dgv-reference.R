#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !file.exists(args[[1]])) {
  stop("Usage: verify-dgv-reference.R PATH_TO_nat_rev_gen_human_normal_cnv.xlsx")
}
expected <- "f06e3b600b77373bd2fc72be0c8dd2c49eb299824d7148a5c11761af51ff6ad0"
observed <- digest::digest(file = args[[1]], algo = "sha256", serialize = FALSE)
if (!identical(observed, expected)) {
  stop("DGV checksum mismatch. Expected ", expected, "; observed ", observed)
}
cat("DGV reference checksum verified:", observed, "\n")
