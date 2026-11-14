script <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1L])
source(file.path(dirname(script), "evidence-utils.R"))
failure <- tryCatch({source(file.path(dirname(script), "b1-native.R"));NULL}, error = identity)
if (exists("res")) save_test_results(res, file.path(Sys.getenv("EVIDENCE_DIR"), "b1", "secondary"))
if (!is.null(failure)) stop(failure)
stopifnot(identical(profile, "native-full-table"), nrow(d) == 14L, sum(d$passed) == 234L,
    sum(d$skipped) == 1L, !any(d$failed | d$error | d$warning))
