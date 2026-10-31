## Public NC wrappers keep rejecting nonnumeric NODE (test-node-message.R)
## and pass numeric invalid counts through to the RORgroup() warning and
## missing additional group (#187).
test_that("BS_parker warns on invalid numeric NODE and leaves the group missing", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = environment())
    se <- OSLO2EMIT0obj$data_input$se_NC
    node <- as.numeric(SummarizedExperiment::colData(se)$NODE)
    node[1] <- -1
    node[2] <- 2.5
    SummarizedExperiment::colData(se)$NODE <- node
    expect_warning(
        res <- suppressMessages(expect_experimental_ror_pc(BS_parker(se, hasClinical = TRUE))),
        "finite nonnegative integer"
    )
    group <- res$score.ROR[["ROR-PC Group (Subtype + Clinic + Prolif.Subtype)"]]
    expect_type(group, "character")
    expect_true(all(is.na(group[1:2])))
    expect_false(any(group %in% c("NA", "NaN")))
    expect_false(any(is.na(group[-(1:2)])))
    expect_true(all(group[-(1:2)] %in% c("low", "med", "high")))
    expect_false(anyNA(res$score.ROR[["ROR-PC (Subtype + Clinic + Prolif.Subtype)"]]))
})
