# CI-only: run the candidate's focused tests against its R CMD check installation.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
src <- normalizePath(args[1L], mustWork = TRUE)
rcheck <- normalizePath(args[2L], mustWork = TRUE)
out <- normalizePath(args[3L], mustWork = TRUE)
.libPaths(c(rcheck, .libPaths()))
stopifnot(as.character(packageVersion("BreastSubtypeR")) == "1.4.3",
          normalizePath(find.package("BreastSubtypeR")) ==
              normalizePath(file.path(rcheck, "BreastSubtypeR")))
setwd(out)
r <- testthat::test_dir(file.path(src, "tests/testthat"), filter = "release-",
                       package = "BreastSubtypeR", load_package = "installed",
                       reporter = "summary", stop_on_failure = FALSE,
                       stop_on_warning = FALSE)
saveRDS(r, "focused-results.rds")
d <- as.data.frame(r)
write.csv(d[, !vapply(d, is.list, logical(1))], "focused-results.csv", row.names = FALSE)
print(d[, c("test", "failed", "error", "warning", "skipped", "passed")])
stopifnot(nrow(d) == 7L, sum(d$passed) == 181L,
          !any(d$failed > 0 | d$error | d$warning > 0 | d$skipped))
cat("PASS: 181 installed-package API/Shiny/export assertions; no failures, warnings or skips.\n")
