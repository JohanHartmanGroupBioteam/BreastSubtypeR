# Historical example input is pinned independently of public data; see fixtures/README.md.
test_that("BS_parker describes the existing numeric NODE requirement", {
    load_historical_example("OSLO2EMIT0obj", envir = environment())
    se <- OSLO2EMIT0obj$data_input$se_NC
    SummarizedExperiment::colData(se)$NODE <- as.character(
        SummarizedExperiment::colData(se)$NODE)
    expect_error(BS_parker(se, hasClinical = TRUE),
        "colData(se_obj)$NODE must be numeric.", fixed = TRUE)
})
