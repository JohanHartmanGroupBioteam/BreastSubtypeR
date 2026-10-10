# Generates tests/testthat/fixtures/cihc_itr_1.5.5.rds: retained outputs of
# BreastSubtypeR 1.5.5 (commit 0dcc0c8a846290fcca3f266bd79248bff84fb58b,
# before target-aware cIHC.itr subsampling). Not run by the tests.
#
# Run with an installation of that commit only, e.g.
#   R_LIBS=<library with 1.5.5 from 0dcc0c8> Rscript make-cihc-itr-1.5.5.R <out.rds>
stopifnot(packageVersion("BreastSubtypeR") == "1.5.5")
suppressPackageStartupMessages({
    library(BreastSubtypeR)
    library(SummarizedExperiment)
})
# guard: must be the pre-fix code (no er_subsample argument)
stopifnot(!"er_subsample" %in% names(formals(BS_cIHC.itr)))
options(BreastSubtypeR.verbose = FALSE)

data("OSLO2EMIT0obj", package = "BreastSubtypeR")
se <- OSLO2EMIT0obj$data_input$se_NC
er <- colData(se)$ER
pick <- function(n_pos, n_neg) {
    c(colnames(se)[er == "ER+"][seq_len(n_pos)],
        colnames(se)[er == "ER-"][seq_len(n_neg)])
}
keep_parts <- function(res) {
    res[c("BS.all", "score.ROR", "BS.itr.keep", "BS.itr.keep.Subtype")]
}
run <- function(ids, ...) {
    res <- BS_cIHC.itr(se[, ids], Subtype = TRUE, hasClinical = TRUE, ...)
    c(keep_parts(res),
        list(testData = res$outList$testData, distances = res$outList$distances))
}
cases <- list(
    full_84_18 = list(ids = colnames(se)),
    erneg_7_18 = list(ids = pick(7, 18)),
    target_15_18 = list(ids = pick(15, 18)),
    near_16_18 = list(ids = pick(16, 18)),
    equal_18_18 = list(ids = pick(18, 18)),
    erpos_45_18 = list(ids = pick(45, 18)),
    erneg_7_18_seed7 = list(ids = pick(7, 18), seed = 7),
    erpos_45_18_seed7 = list(ids = pick(45, 18), seed = 7),
    erneg_7_18_ratio1 = list(ids = pick(7, 18), ratio = 1),
    erpos_45_18_ratio1 = list(ids = pick(45, 18), ratio = 1),
    erpos_45_18_ratio2 = list(ids = pick(45, 18), ratio = 2)
)
out <- list(cihc_itr = list(), multi = list())
for (nm in names(cases)) {
    cs <- cases[[nm]]
    args <- c(list(ids = cs$ids), cs[setdiff(names(cs), "ids")])
    out$cihc_itr[[nm]] <- c(list(ids = cs$ids), do.call(run, args))
}
methods_all <- c("parker.original", "genefu.scale", "genefu.robust", "ssBC",
    "ssBC.v2", "cIHC", "cIHC.itr", "PCAPAM50", "AIMS", "sspbc")
din <- OSLO2EMIT0obj$data_input
sub_input <- function(ids) list(se_NC = din$se_NC[, ids], se_SSP = din$se_SSP[, ids])
multi <- function(input, methods) {
    res <- suppressWarnings(BS_Multi(input, methods = methods, Subtype = TRUE))
    res[c("res_subtypes", "res_subtypes.Subtype")]
}
out$multi$full_manual <- multi(din, methods_all)
out$multi$full_auto <- multi(din, "AUTO")
out$multi$erneg_7_18_manual <- c(list(ids = pick(7, 18)),
    multi(sub_input(pick(7, 18)), methods_all))
attr(out, "provenance") <- list(
    package_version = "1.5.5",
    commit = "0dcc0c8a846290fcca3f266bd79248bff84fb58b",
    created = format(Sys.time(), tz = "UTC", usetz = TRUE),
    R = R.version.string,
    platform = R.version$platform,
    rng = RNGkind()
)
args <- commandArgs(TRUE)
saveRDS(out, args[1], compress = "xz")
cat("saved", args[1], ":", file.size(args[1]), "bytes\n")
