test_that("BS_parker describes the existing numeric NODE requirement", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = environment())
    se <- OSLO2EMIT0obj$data_input$se_NC
    SummarizedExperiment::colData(se)$NODE <- as.character(
        SummarizedExperiment::colData(se)$NODE)
    expect_error(BS_parker(se, hasClinical = TRUE),
        "colData(se_obj)$NODE must be numeric.", fixed = TRUE)
})
