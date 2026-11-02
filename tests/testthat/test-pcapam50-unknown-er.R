# Historical example input is pinned independently of public data; see fixtures/README.md.
## PCAPAM50 and samples whose ER status is missing or non-canonical.
##
## BS_PCAPAM50() derives the IHC label from the ER column, so unknown ER
## becomes NA. The reference implementation (CRAN package PCAPAM50) tests
## ER-negativity as !grepl("^L", IHC), which is TRUE for NA, so such samples
## used to be counted as ER-negative in the PC1 axis check, the cutoff search
## and the ER-balanced centering set. They are now excluded from those steps
## and still classified.

test_that(".ihc_status() reports missing and empty labels as not evaluable", {
    status <- BreastSubtypeR:::.ihc_status(
        c("Luminal", "luminal", "LA", "Non-Luminal", "TN", "Her2+", NA, "", "   ")
    )
    expect_identical(
        status,
        c("luminal", "luminal", "luminal", "nonluminal", "nonluminal",
          "nonluminal", NA, NA, NA)
    )
    expect_identical(BreastSubtypeR:::.ihc_status(character(0)), character(0))
    expect_identical(BreastSubtypeR:::.ihc_status(factor("LB1")), "luminal")
})

test_that("a complete ER column reproduces the packaged PCAPAM50 calls", {
    load_historical_example("OSLO2EMIT0obj")
    se <- OSLO2EMIT0obj$data_input$se_NC
    expect_identical(
        sum(is.na(SummarizedExperiment::colData(se)$ER)),
        0L
    )
    res <- suppressMessages(BS_PCAPAM50(se, Subtype = FALSE, hasClinical = FALSE))
    fresh <- res$BS.all$BS[match(
        rownames(OSLO2EMIT0obj$res$res_subtypes),
        res$BS.all$PatientID
    )]
    expect_identical(
        unname(fresh),
        OSLO2EMIT0obj$res$res_subtypes$PCAPAM50
    )
})

test_that("samples with unknown ER are reported, excluded and still classified", {
    load_historical_example("OSLO2EMIT0obj")
    se <- OSLO2EMIT0obj$data_input$se_NC
    unknown <- colnames(se)[1:6]
    SummarizedExperiment::colData(se)$ER[1:6] <- NA_character_

    expect_message(
        res <- suppressWarnings(BS_PCAPAM50(se, Subtype = FALSE, hasClinical = FALSE)),
        "6 sample(s) without an evaluable ER/IHC status are excluded",
        fixed = TRUE
    )
    ## every sample keeps a call, including those with unknown ER
    expect_identical(nrow(res$BS.all), ncol(se))
    expect_false(anyNA(res$BS.all$BS))
    expect_true(all(unknown %in% res$BS.all$PatientID))

    ## a non-canonical label is treated like a missing one
    se2 <- OSLO2EMIT0obj$data_input$se_NC
    SummarizedExperiment::colData(se2)$ER[1:6] <- "unknown"
    res2 <- suppressMessages(suppressWarnings(
        BS_PCAPAM50(se2, Subtype = FALSE, hasClinical = FALSE)
    ))
    expect_identical(res2$BS.all$BS, res$BS.all$BS)
})
