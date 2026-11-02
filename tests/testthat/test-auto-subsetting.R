# Historical example input is pinned independently of public data; see fixtures/README.md.
## Regression tests for the AUTO-mode ssBC / ssBC.v2 sample subsets returned by
## get_methods(). Samples whose ER (or HER2) status is missing must not appear
## in the subsets: a missing value used to produce an NA sample name, which
## either broke the SummarizedExperiment subsetting in BS_Multi() or, when the
## padded vector was as long as the cohort, silently disabled the subsetting.

make_pheno <- function(er, her2) {
    pheno <- data.frame(
        PatientID = sprintf("S%03d", seq_along(er)),
        ER = er,
        HER2 = her2,
        stringsAsFactors = FALSE
    )
    rownames(pheno) <- pheno$PatientID
    pheno
}

expect_valid_subset <- function(subset, pheno, expected_ids) {
    expect_false(anyNA(subset))
    expect_true(all(subset %in% rownames(pheno)))
    expect_setequal(subset, expected_ids)
}

test_that("missing ER values do not leak NA sample names into the ssBC subsets", {
    old <- options(BreastSubtypeR.verbose = FALSE)
    on.exit(options(old), add = TRUE)

    pheno <- make_pheno(
        er = c(rep("ER+", 50), rep("ER-", 45), rep(NA_character_, 5)),
        her2 = rep(c("HER2-", "HER2+"), 50)
    )
    out <- BreastSubtypeR:::get_methods(pheno)

    expect_true(all(c("ssBC", "ssBC.v2") %in% out$methods))
    known_er <- rownames(pheno)[!is.na(pheno$ER)]
    expect_valid_subset(out$samples_ER.icd, pheno, known_er)
    expect_valid_subset(out$samples_ERHER2.icd, pheno, known_er)
    expect_length(out$samples_ER.icd, 95L)
})

test_that("missing HER2 values do not leak NA sample names into the ssBC.v2 subset", {
    old <- options(BreastSubtypeR.verbose = FALSE)
    on.exit(options(old), add = TRUE)

    pheno <- make_pheno(
        er = c(rep("ER+", 10), rep("ER-", 80), rep(NA_character_, 10)),
        her2 = c(rep("HER2-", 45), rep("HER2+", 45), rep(NA_character_, 10))
    )
    out <- BreastSubtypeR:::get_methods(pheno)

    expect_identical(out$cohort.select, "ERneg")
    expect_identical(out$methods, c("ssBC", "ssBC.v2", "AIMS", "sspbc"))
    ## ssBC: only the ER- group reaches its threshold (ER+ has 10 < 15).
    er_neg <- rownames(pheno)[which(pheno$ER == "ER-")]
    expect_valid_subset(out$samples_ER.icd, pheno, er_neg)
    ## ssBC.v2: every ER/HER2 subgroup with known status is large enough
    ## (ER+/HER2- 10 >= 8, ER-/HER2+ 45 >= 9, ER-/HER2- 35 >= 9).
    known_both <- rownames(pheno)[!is.na(pheno$ER) & !is.na(pheno$HER2)]
    expect_valid_subset(out$samples_ERHER2.icd, pheno, known_both)
    expect_length(out$samples_ERHER2.icd, 90L)
    ## The subset is shorter than the cohort, so BS_Multi() will subset the
    ## SummarizedExperiment with it; every name must be a real sample.
    expect_lt(length(out$samples_ER.icd), nrow(pheno))
})

test_that("the AUTO fallback to AIMS and sspbc is announced", {
    old <- options(BreastSubtypeR.verbose = TRUE)
    on.exit(options(old), add = TRUE)

    ## ER+ only, 17 samples, HER2+ 5 and HER2- 6 (both below 8), 6 missing:
    ## no ER+ sub-rule matches, so get_methods() falls back to AIMS and sspbc.
    pheno <- make_pheno(
        er = rep("ER+", 17),
        her2 = c(rep("HER2+", 5), rep("HER2-", 6), rep(NA_character_, 6))
    )
    expect_message(
        out <- BreastSubtypeR:::get_methods(pheno),
        "running the single-sample predictors AIMS and sspbc only",
        fixed = TRUE
    )
    expect_identical(out$methods, c("AIMS", "sspbc"))
    expect_null(out$samples_ER.icd)
    expect_null(out$samples_ERHER2.icd)
})

test_that("subgroups exactly at their minimum are included in the ssBC subsets", {
    old <- options(BreastSubtypeR.verbose = FALSE)
    on.exit(options(old), add = TRUE)

    ## 15 ER+ (minimum 15) and 40 ER-: both ER groups qualify for ssBC.
    pheno <- make_pheno(
        er = c(rep("ER+", 15), rep("ER-", 40)),
        her2 = c(rep(c("HER2-", "HER2+"), 7), "HER2-", rep(c("HER2-", "HER2+"), 20))
    )
    out <- BreastSubtypeR:::get_methods(pheno)
    expect_identical(out$cohort.select, "mixed")
    expect_true("ssBC" %in% out$methods)
    expect_setequal(out$samples_ER.icd, rownames(pheno))

    ## ER+ only with HER2+ 9 and HER2- 8 (minimum 8): both HER2 subgroups
    ## qualify for ssBC.v2.
    pheno2 <- make_pheno(
        er = rep("ER+", 17),
        her2 = c(rep("HER2+", 9), rep("HER2-", 8))
    )
    out2 <- BreastSubtypeR:::get_methods(pheno2)
    expect_identical(out2$methods, c("ssBC", "ssBC.v2", "AIMS", "sspbc"))
    expect_setequal(out2$samples_ERHER2.icd, rownames(pheno2))
    expect_setequal(out2$samples_ER.icd, rownames(pheno2))
})

test_that("the packaged OSLO2-EMIT0 cohort runs ssBC on all tumors", {
    old <- options(BreastSubtypeR.verbose = FALSE)
    on.exit(options(old), add = TRUE)

    load_historical_example("OSLO2EMIT0obj")
    pheno <- as.data.frame(SummarizedExperiment::colData(OSLO2EMIT0obj$data_input$se_NC))
    ## 84 ER+ and 18 ER- tumors; the ER- count equals the ER- minimum of 18.
    expect_identical(sum(pheno$ER == "ER-"), 18L)
    out <- BreastSubtypeR:::get_methods(pheno)
    expect_identical(out$cohort.select, "mixed")
    expect_setequal(out$samples_ER.icd, pheno$PatientID)
    ## ssBC.v2 subsets: ER+/HER2- (80) and ER-/HER2- (14); the two HER2+
    ## subgroups (4 each) stay below their minimums of 8 and 9.
    expect_length(out$samples_ERHER2.icd, 94L)
    expect_true(all(pheno$HER2[match(out$samples_ERHER2.icd, pheno$PatientID)] == "HER2-"))
})
