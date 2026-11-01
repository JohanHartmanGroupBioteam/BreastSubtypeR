# Historical example input is pinned independently of public data; see fixtures/README.md.
test_that("BS_Multi performs subtyping correctly", {
    # Load required dataset
    load_historical_example("OSLO2EMIT0obj")

    # Define methods to use for subtyping
    methods <- c("genefu.robust", "PCAPAM50", "ssBC", "AIMS")

    # Perform subtyping
    res.test <- BS_Multi(
        data_input = OSLO2EMIT0obj$data_input,
        methods = methods,
        Subtype = FALSE,
        hasClinical = FALSE
    )

    # Validate results
    expect_true(identical(res.test$res_subtypes[, methods], OSLO2EMIT0obj$res$res_subtypes[, methods])) # Validate results
})
