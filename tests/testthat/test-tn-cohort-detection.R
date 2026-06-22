make_pheno <- function(er, her2, tn = NULL) {
    pheno <- data.frame(
        PatientID = paste0("S", seq_along(er)),
        ER = er,
        HER2 = her2,
        stringsAsFactors = FALSE
    )
    if (!is.null(tn)) {
        pheno$TN <- tn
    }
    rownames(pheno) <- pheno$PatientID
    pheno
}

expect_auto_panel <- function(pheno, methods, cohort) {
    old <- options(BreastSubtypeR.verbose = FALSE)
    on.exit(options(old), add = TRUE)
    out <- BreastSubtypeR:::get_methods(pheno)
    expect_identical(out$methods, methods)
    expect_identical(out$cohort.select, cohort)
}

tn_panel <- c("ssBC", "ssBC.v2", "AIMS", "sspbc")
full_mixed_panel <- c(
    "parker.original", "genefu.scale", "genefu.robust",
    "ssBC", "ssBC.v2", "cIHC", "cIHC.itr", "PCAPAM50", "AIMS", "sspbc"
)

test_that("TN-only cohorts with a TN column select the TNBC panel", {
    pheno <- make_pheno(
        er = rep("ER-", 20),
        her2 = rep("HER2-", 20),
        tn = rep("TN", 20)
    )

    expect_auto_panel(pheno, tn_panel, "TNBC")
})

test_that("mixed cohorts with some TN samples are not treated as TNBC", {
    pheno <- make_pheno(
        er = c(rep("ER+", 50), rep("ER-", 30)),
        her2 = c(rep("HER2-", 30), rep("HER2+", 20), rep("HER2-", 20), rep("HER2+", 10)),
        tn = c(rep("nonTN", 60), rep("TN", 20))
    )

    expect_auto_panel(pheno, full_mixed_panel, "mixed")
})

test_that("UNC-like mixed cohorts with TN annotations are not treated as TNBC", {
    pheno <- make_pheno(
        er = c(rep("ER+", 62), rep("ER-", 38)),
        her2 = c(rep("HER2+", 9), rep("HER2-", 53), rep("HER2+", 6), rep("HER2-", 32)),
        tn = c(rep("nonTN", 70), rep("TN", 30))
    )

    expect_auto_panel(pheno, full_mixed_panel, "mixed")
})

test_that("missing or empty TN values do not trigger TNBC handling", {
    pheno <- make_pheno(
        er = c(rep("ER+", 25), rep("ER-", 25)),
        her2 = c(rep("HER2-", 20), rep("HER2+", 5), rep("HER2-", 20), rep("HER2+", 5)),
        tn = rep(c(NA_character_, "", " "), length.out = 50)
    )

    expect_auto_panel(pheno, full_mixed_panel, "mixed")
})

test_that("partly missing TN values with TN and nonTN remain non-TNBC", {
    pheno <- make_pheno(
        er = c(rep("ER+", 25), rep("ER-", 25)),
        her2 = c(rep("HER2-", 20), rep("HER2+", 5), rep("HER2-", 20), rep("HER2+", 5)),
        tn = c(rep("TN", 10), rep("nonTN", 10), rep(NA_character_, 10), rep("", 10), rep(" ", 10))
    )

    expect_auto_panel(pheno, full_mixed_panel, "mixed")
})

test_that("partly missing TN values select TNBC when all evaluable values are TN", {
    pheno <- make_pheno(
        er = rep("ER-", 20),
        her2 = rep("HER2-", 20),
        tn = c(rep("TN", 18), NA_character_, "")
    )

    expect_auto_panel(pheno, tn_panel, "TNBC")
})
