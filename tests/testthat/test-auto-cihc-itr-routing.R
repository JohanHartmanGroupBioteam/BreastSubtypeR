## AUTO selection of cIHC.itr by get_methods() and the cIHC.itr balancing plan
## for ER-negative-majority cohorts. Synthetic phenotype tables only; no
## expression data is fitted. The selection logic and its thresholds are
## unchanged; these tests record which AUTO cohorts can be affected by the
## target-aware ER subsampling of cIHC.itr (NEWS, 1.5.6).

make_pheno <- function(n_pos, n_neg, her2 = "HER2-") {
    er <- c(rep("ER+", n_pos), rep("ER-", n_neg))
    pheno <- data.frame(
        PatientID = sprintf("S%03d", seq_along(er)),
        ER = er,
        HER2 = rep(her2, length(er)),
        stringsAsFactors = FALSE
    )
    rownames(pheno) <- pheno$PatientID
    pheno
}

auto <- function(pheno) {
    old <- options(BreastSubtypeR.verbose = FALSE)
    on.exit(options(old), add = TRUE)
    BreastSubtypeR:::get_methods(pheno)
}

ten <- c(
    "parker.original", "genefu.scale", "genefu.robust", "ssBC", "ssBC.v2",
    "cIHC", "cIHC.itr", "PCAPAM50", "AIMS", "sspbc"
)
eight <- c(
    "genefu.robust", "ssBC", "ssBC.v2", "cIHC", "cIHC.itr", "PCAPAM50",
    "AIMS", "sspbc"
)

test_that("AUTO selects cIHC.itr in both mixed-cohort branches", {
    ## ER+ : ER- counts (all HER2-, no TN column) and the expected branch:
    ## ten methods when 0.39 < ER+ fraction < 0.69, eight methods otherwise.
    cases <- list(
        list(20, 80, eight), # fraction 0.20
        list(39, 61, eight), # fraction 0.39: not above the lower threshold
        list(40, 60, ten), # fraction 0.40
        list(45, 55, ten), # fraction 0.45
        list(50, 50, ten), # fraction 0.50
        list(68, 32, ten), # fraction 0.68
        list(69, 31, eight), # fraction 0.69: not below the upper threshold
        list(80, 20, eight) # fraction 0.80
    )
    for (cs in cases) {
        out <- auto(make_pheno(cs[[1]], cs[[2]]))
        expect_identical(out$cohort.select, "mixed")
        expect_identical(out$methods, cs[[3]])
        expect_true("cIHC.itr" %in% out$methods)
    }
})

test_that("ER composition alone does not make a cohort eligible for cIHC.itr", {
    ## Below the ER+ minimum of 15: the ER-/HER2- route, without cIHC.itr,
    ## although ER- samples outnumber ER+ samples.
    few_pos <- auto(make_pheno(10, 90))
    expect_false(identical(few_pos$cohort.select, "mixed"))
    expect_identical(few_pos$methods, c("ssBC.v2", "AIMS", "sspbc"))
    expect_false("cIHC.itr" %in% few_pos$methods)

    ## A HER2+ cohort with the same 20 : 80 ER counts takes the HER2+ route.
    her2_pos <- auto(make_pheno(20, 80, her2 = "HER2+"))
    expect_false(identical(her2_pos$cohort.select, "mixed"))
    expect_false("cIHC.itr" %in% her2_pos$methods)
})

test_that("the default-ratio plans differ for an ER-negative-majority cohort", {
    ## 20 ER+ / 80 ER- (ER+ fraction 0.20), default ratio 54 / 64. This is a
    ## plan check: the two modes centre on different subsets. Whether a
    ## subtype call changes depends on the data.
    plan <- BreastSubtypeR:::.cihc_itr_plan
    expect_identical(plan(20, 80, er_subsample = "larger"),
        list(group = "ER-", size = 17L)) # 20 ER+ : 17 ER-
    expect_identical(plan(20, 80, er_subsample = "target"),
        list(group = "ER-", size = 24L)) # 20 ER+ : 24 ER-
})
