# Historical example input is pinned independently of public data; see fixtures/README.md.
## The ER-balancing methods need both ER groups. Without one of them they
## used to fail deep inside the median computation with messages such as
## "undefined columns selected" or "arguments imply differing number of rows".

single_group_input <- function(er_value) {
    load_historical_example("OSLO2EMIT0obj")
    se <- OSLO2EMIT0obj$data_input$se_NC
    SummarizedExperiment::colData(se)$ER <- er_value
    se
}

test_that("cIHC, cIHC.itr and PCAPAM50 report a missing ER group clearly", {
    se_pos <- single_group_input("ER+")
    expect_error(
        suppressMessages(BS_cIHC(se_pos, Subtype = FALSE, hasClinical = FALSE)),
        "cIHC requires both ER+ and ER- samples for ER balancing; found 102 ER+ and 0 ER- samples",
        fixed = TRUE
    )
    expect_error(
        suppressMessages(BS_cIHC.itr(se_pos, iteration = 2, Subtype = FALSE, hasClinical = FALSE)),
        "cIHC.itr requires both ER+ and ER- samples for ER balancing; found 102 ER+ and 0 ER- samples",
        fixed = TRUE
    )
    expect_error(
        suppressMessages(BS_PCAPAM50(se_pos, Subtype = FALSE, hasClinical = FALSE)),
        "PCAPAM50 requires both luminal (ER+) and non-luminal (ER-) IHC classes",
        fixed = TRUE
    )

    se_neg <- single_group_input("ER-")
    expect_error(
        suppressMessages(BS_cIHC(se_neg, Subtype = FALSE, hasClinical = FALSE)),
        "found 0 ER+ and 102 ER- samples",
        fixed = TRUE
    )
})

test_that("cIHC still runs on a cohort with both ER groups", {
    load_historical_example("OSLO2EMIT0obj")
    se <- OSLO2EMIT0obj$data_input$se_NC
    res <- suppressMessages(BS_cIHC(se, Subtype = FALSE, hasClinical = FALSE))
    expect_identical(
        unname(res$BS.all$BS[match(rownames(OSLO2EMIT0obj$res$res_subtypes), res$BS.all$PatientID)]),
        OSLO2EMIT0obj$res$res_subtypes$cIHC
    )
})
