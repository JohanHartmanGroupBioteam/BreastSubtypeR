# Historical example input is pinned independently of public data; see fixtures/README.md.
## Raw Shannon entropy of the per-method calls.
##
## table() drops missing values, so a row in which every executed method
## returned NA used to yield -sum(numeric(0)) = 0, the same value as unanimous
## agreement. Such a row has no call distribution and is now reported as NA.
## Rows with at least one call are unchanged, and the definition of the
## statistic (unnormalized, in bits) is untouched.

test_that("entropy is NA only when no method returned a call", {
    expect_identical(BreastSubtypeR:::get_entropy(c(NA, NA, NA)), NA_real_)
    expect_identical(BreastSubtypeR:::get_entropy(character(0)), NA_real_)
    expect_identical(BreastSubtypeR:::get_entropy(NA_character_), NA_real_)

    ## one call is a degenerate but defined distribution: entropy 0
    expect_equal(BreastSubtypeR:::get_entropy(c("LumA", NA, NA)), 0)
    expect_equal(BreastSubtypeR:::get_entropy(rep("LumA", 8)), 0)
})

test_that("entropy values for rows with calls are the raw Shannon entropy in bits", {
    expect_equal(BreastSubtypeR:::get_entropy(c("LumA", "Basal")), 1)
    expect_equal(BreastSubtypeR:::get_entropy(c("LumA", "LumA", "Basal", "Her2")), 1.5)
    expect_equal(
        BreastSubtypeR:::get_entropy(c("LumA", "LumA", "LumB", NA)),
        -(2 / 3 * log2(2 / 3) + 1 / 3 * log2(1 / 3))
    )
    ## missing calls are dropped, not counted as a category
    expect_equal(
        BreastSubtypeR:::get_entropy(c("LumA", "Basal", NA, NA)),
        BreastSubtypeR:::get_entropy(c("LumA", "Basal"))
    )
})

test_that("the packaged example keeps its stored entropy column", {
    load_historical_example("OSLO2EMIT0obj")
    rs <- OSLO2EMIT0obj$res$res_subtypes
    methods <- setdiff(colnames(rs), "entropy")
    recomputed <- apply(rs[, methods, drop = FALSE], 1, BreastSubtypeR:::get_entropy)
    expect_equal(unname(recomputed), rs$entropy)
    ## every tumor in the packaged cohort has calls, so none is NA
    expect_false(anyNA(recomputed))
})

test_that("a manual run that calls nothing for a sample reports NA entropy", {
    load_historical_example("OSLO2EMIT0obj")
    input <- OSLO2EMIT0obj$data_input
    ## ssBC and ssBC.v2 center within ER (and ER/HER2) subgroups, so a sample
    ## whose status is unknown gets no call from either; with only these two
    ## methods executed, such a sample has no calls at all. AUTO never produces
    ## this situation because every AUTO panel also runs AIMS and sspbc.
    SummarizedExperiment::colData(input$se_NC)$ER[1:4] <- NA_character_
    SummarizedExperiment::colData(input$se_NC)$HER2[1:4] <- NA_character_

    res <- suppressMessages(suppressWarnings(BS_Multi(
        data_input = input,
        methods = c("ssBC", "ssBC.v2"),
        Subtype = FALSE,
        hasClinical = FALSE
    )))
    unknown <- rownames(res$res_subtypes)[1:4]
    expect_true(all(is.na(res$res_subtypes[unknown, "ssBC"])))
    expect_true(all(is.na(res$res_subtypes[unknown, "ssBC.v2"])))
    expect_true(all(is.na(res$res_subtypes[unknown, "entropy"])))
    ## samples that were called keep a numeric entropy
    called <- setdiff(rownames(res$res_subtypes), unknown)
    expect_false(anyNA(res$res_subtypes[called, "entropy"]))
})
