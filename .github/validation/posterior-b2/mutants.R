library(BreastSubtypeR); library(testthat)
root <- system.file("tests/testthat", package = "BreastSubtypeR")
helper <- readLines(file.path(root, "helper-issue199-portability.R"), warn = FALSE)
target <- file.path(root, "test-issue199-posterior-links.R")
# test_file reloads nearby helpers; use an identical test copy without helpers.
dir.create("scratch/mutants", recursive = TRUE, showWarnings = FALSE)
scratch <- "scratch/mutants/test-issue199-posterior-links.R"
stopifnot(file.copy(target, scratch))
stopifnot(identical(readBin(target, "raw", file.info(target)$size),
    readBin(scratch, "raw", file.info(scratch)$size)))
rows <- list()
for (kind in c("control", "old_decision_gate", "old_walker_gate",
               "default_ties", "numeric_equal_selection")) {
    lines <- helper
    if (kind == "old_decision_gate") {
        line <- grep('testthat::expect_identical(a$cl, b$cl)', lines, fixed = TRUE)
        stopifnot(length(line) == 1L)
        lines <- append(lines, '            testthat::expect_identical(a$prob, b$prob)',
            after = line)
    }
    if (kind == "old_walker_gate") {
        line <- grep('} else if (p %in% derived) {', lines, fixed = TRUE)
        stopifnot(length(line) == 1L)
        lines <- append(lines, '            if (!identical(a, b, num.eq = FALSE)) fail(p, "restored exact prob gate")',
            after = line)
    }
    if (kind == "default_ties") {
        old <- 'max.col(p[!masked, , drop = FALSE], ties.method = "first")'
        stopifnot(sum(grepl(old, lines, fixed = TRUE)) == 1L)
        lines <- sub(old, 'max.col(p[!masked, , drop = FALSE])', lines, fixed = TRUE)
    }
    if (kind == "numeric_equal_selection") {
        old <- 'identical(unname(o$prob[, k]), probability, num.eq = FALSE)'
        stopifnot(sum(grepl(old, lines, fixed = TRUE)) == 1L)
        lines <- sub(old, 'identical(unname(o$prob[, k]), probability, num.eq = TRUE)',
            lines, fixed = TRUE)
    }
    e <- new.env(parent = globalenv())
    eval(parse(text = lines), envir = e)
    set.seed(199)
    result <- test_file(scratch, env = e, reporter = "silent",
        stop_on_failure = FALSE, stop_on_warning = FALSE)
    d <- as.data.frame(result)
    write.table(d[c("file", "test", "nb", "passed", "failed", "error", "warning", "skipped")],
        paste0("evidence/mutant-", kind, "-blocks.tsv"), sep = "\t", row.names = FALSE)
    totals <- colSums(d[c("passed", "failed", "error", "warning", "skipped")])
    rows[[kind]] <- data.frame(mutant = kind, t(totals), check.names = FALSE)
    if (kind == "control") stopifnot(!any(d$failed | d$error | d$warning | d$skipped))
    else stopifnot(any(d$failed | d$error), !any(d$warning | d$skipped))
}
summary <- do.call(rbind, rows)
write.table(summary, "evidence/mutation-summary.tsv", sep = "\t", row.names = FALSE)
print(summary)
