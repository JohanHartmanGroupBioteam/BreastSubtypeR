library(BreastSubtypeR)
library(testthat)
script <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1L])
source(file.path(dirname(script), "evidence-utils.R"))
mode <- commandArgs(TRUE)[1L]
stopifnot(normalizePath(find.package("BreastSubtypeR")) == normalizePath(file.path(Sys.getenv("R_LIBS"), "BreastSubtypeR")))
options(testthat.progress.max_fails = Inf)
filter <- switch(mode, primary = "^issue199-posterior-links$", ror = "^issue199-ror-bounds$",
    diagnostic = "^(aims-single-sample|expanded-examples|sspbc-ror|zero-evaluable-rules|ror-boundaries|ror-groups|ror-proliferation-coverage|ror-standardization)$", stop("unknown suite"))
result <- test_package("BreastSubtypeR", filter = filter, reporter = "summary", stop_on_failure = FALSE, stop_on_warning = FALSE)
d <- save_test_results(result, file.path(Sys.getenv("EVIDENCE_DIR"), "tests", mode))
if (mode == "primary") {
    stopifnot(nrow(d) == 9L, sum(d$passed) == 320L,
        !any(d$failed | d$error | d$warning | d$skipped))
} else if (mode == "ror") stopifnot(nrow(d) == 8L) else stopifnot(nrow(d) == 73L, length(unique(d$file)) == 8L)
# Diagnostic failures remain recorded failures; the separate provenance step classifies them.
