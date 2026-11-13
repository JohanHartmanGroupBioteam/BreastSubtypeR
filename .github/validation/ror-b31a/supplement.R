library(BreastSubtypeR)
library(testthat)
script <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1L])
here <- dirname(script)
source(file.path(here, "evidence-utils.R"))
root <- system.file("tests/testthat", package = "BreastSubtypeR")
sys.source(file.path(root, "helper-issue199-portability.R"), globalenv())
for (expr in parse(file.path(root, "test-issue199-ror-bounds.R"))) {
    if (!(is.call(expr) && identical(expr[[1L]], as.name("test_that")))) eval(expr, globalenv())
}
controls <- list()
refuse <- function(name, code, reason) {
    message <- tryCatch({ force(code); "ACCEPTED" }, error = conditionMessage)
    controls[[length(controls)+1L]] <<- data.frame(control = name, message,
        intended_reason = reason, detected = grepl(reason, message, fixed = TRUE))
}
set.seed(3112026); seed <- .Random.seed; kind <- RNGkind()
b <- issue199_ror_approved("TCGABRCAobj")$res$results$ssBC
for (j in seq_len(5L)) {
    a <- b; a$outList$distances[1L, j] <- ror31_next(a$outList$distances[1L, j])
    refuse(paste0("distance_column_", j), ror31_method(a, b), "five-distance identity")
}
f <- "ROR-P (Subtype + Prolif)"; cut <- 100*(-.25+.35)/.85
a <- b; z <- b
a$score.ROR[[f]][1L] <- cut + 2^-45
z$score.ROR[[f]][1L] <- cut - 2^-45
refuse("exported_cutoff_side_full_method", ror31_method(a, z), "ROR cutoff sides")
# A scratch function fault injection exercises the raw-side branch directly.
# It is not a scientific result and is never used in the primary native gate.
scratch <- new.env(parent = globalenv())
scratch$issue199_ror_domain <- local({
    original <- issue199_ror_domain; calls <- 0L
    function(...) { result <- original(...); calls <<- calls+1L
        if (calls == 1L) result$raw[1L] <- if (result$raw[1L] > .1) -.3 else .2
        result }
})
fault <- issue199_ror_method; environment(fault) <- scratch
refuse("raw_cutoff_guard_scratch_fault_injection", fault(b, b, b, "TCGABRCAobj", "ssBC"), "raw cutoff sides")
stopifnot(identical(seed, .Random.seed), identical(kind, RNGkind()))
controls <- do.call(rbind, controls)
write_tsv(controls, file.path(Sys.getenv("EVIDENCE_DIR"), "tests", "supplement-controls.tsv"))
print(controls)
stopifnot(nrow(controls) == 7L, all(controls$detected))
