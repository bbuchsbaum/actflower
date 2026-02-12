#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1L]]
}

threshold_file <- arg_value("--threshold-file", "tools/coverage_thresholds.json")

suppressPackageStartupMessages(library(covr))
suppressPackageStartupMessages(library(jsonlite))

fail <- function(msg) {
  cat(sprintf("[FAIL] %s\n", msg))
  quit(status = 1)
}

ok <- function(msg) {
  cat(sprintf("[OK] %s\n", msg))
}

th <- fromJSON(threshold_file)
min_percent <- as.numeric(th$minimum_percent)
file_min <- th$minimum_file_percent
if (is.null(file_min)) file_min <- list()

cov <- covr::package_coverage()
overall <- as.numeric(covr::percent_coverage(cov))
if (!is.finite(overall) || overall < min_percent) {
  fail(sprintf("overall coverage %.2f%% < %.2f%%", overall, min_percent))
}
ok(sprintf("overall coverage %.2f%% >= %.2f%%", overall, min_percent))

d <- as.data.frame(cov)
by_file <- aggregate(value ~ filename, data = d, FUN = function(v) mean(v > 0) * 100)
file_pct <- setNames(by_file$value, by_file$filename)

if (length(file_min)) {
  for (nm in names(file_min)) {
    need <- as.numeric(file_min[[nm]])
    got <- as.numeric(file_pct[[nm]])
    if (!is.finite(got)) {
      fail(sprintf("coverage file not found in report: %s", nm))
    }
    if (got < need) {
      fail(sprintf("coverage %s %.2f%% < %.2f%%", nm, got, need))
    }
    ok(sprintf("coverage %s %.2f%% >= %.2f%%", nm, got, need))
  }
}

cat("[PASS] coverage thresholds satisfied.\n")
