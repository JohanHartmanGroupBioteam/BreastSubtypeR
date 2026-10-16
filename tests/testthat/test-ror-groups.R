## Regression tests for the ROR risk-group columns produced by RORgroup().
## The synthetic 'out' object isolates the combined (ROR-C) and combined +
## proliferation (ROR-PC) models: only the LumA correlation and the
## proliferation genes vary, so the expected risk groups follow directly from
## the published coefficients and the thresholds hard-coded in RORgroup()
## (ROR-C: low < -0.1, high > 0.2; ROR-PC: low < -0.2, high > 0.2).

make_ror_input <- function() {
    ids <- c("S1", "S2", "S3")
    genes <- c("ESR1", "ERBB2", "CCNB1", "UBE2C", "MKI67")
    prolif <- c(-1, 1, 0)
    testData <- rbind(
        ESR1 = c(0.5, 0.4, 0.3),
        ERBB2 = c(-0.5, -0.4, -0.3),
        CCNB1 = prolif,
        UBE2C = prolif,
        MKI67 = prolif
    )
    colnames(testData) <- ids
    centroids <- matrix(
        c(
            -1.0, 0.2, 1.2, 0.8, 0.1,
            0.1, 1.1, -0.4, 0.2, -0.2,
            0.9, 0.4, -0.8, 0.6, -0.1,
            0.8, 0.3, -0.6, 0.9, 0.0,
            0.7, 0.5, -0.7, 0.7, 0.2
        ),
        nrow = length(genes),
        byrow = TRUE,
        dimnames = list(genes, c("Basal", "Her2", "LumA", "LumB", "Normal"))
    )
    ## Basal, Her2 and LumB correlations are zero, so the combined score is
    ## -0.2608388 * LumA and the combined + proliferation score is
    ## -0.090436516 * LumA + 0.327259375 * prolif.
    distances <- matrix(
        0,
        nrow = length(ids),
        ncol = 5,
        dimnames = list(ids, colnames(centroids))
    )
    distances[, "LumA"] <- c(-0.9, 0.9, 0)
    distances[, "Normal"] <- c(0.1, 0.2, 0.3)
    predictions <- setNames(c("Basal", "LumA", "Normal"), ids)
    out <- list(
        predictions = predictions,
        testData = testData,
        distances = distances,
        dist.RORSubtype = distances[, c("Basal", "Her2", "LumA", "LumB")],
        centroids = centroids
    )
    df.cln <- data.frame(
        PatientID = ids,
        TSIZE = c(0, 0, 0),
        NODE = c(0, 0, 0),
        stringsAsFactors = FALSE
    )
    list(out = out, df.cln = df.cln)
}

test_that("ROR-C group follows the ROR-C score, not the ROR-PC score", {
    input <- make_ror_input()
    res <- suppressMessages(BreastSubtypeR:::RORgroup(
        input$out,
        input$df.cln,
        Subtype = FALSE,
        hasClinical = TRUE
    ))

    ## S1: combined = 0.2348 (high), combined + prolif = -0.2459 (low)
    ## S2: combined = -0.2348 (low), combined + prolif = 0.2459 (high)
    ## S3: both scores are 0 (med)
    expect_identical(
        unname(as.character(res[["ROR-C Group (Subtype + Clinic)"]])),
        c("high", "low", "med")
    )
    expect_identical(
        unname(as.character(res[["ROR-PC Group (Subtype + Clinic + Prolif)"]])),
        c("low", "high", "med")
    )

    ## The group must be reproducible from the reported ROR-C score. The score
    ## is reported on the 100 * (x + 0.35) / 0.85 scale used for all ROR columns.
    combined <- res[["ROR-C (Subtype + Clinic)"]] * 0.85 / 100 - 0.35
    expected <- ifelse(combined > 0.2, "high", ifelse(combined < -0.1, "low", "med"))
    expect_identical(
        unname(as.character(res[["ROR-C Group (Subtype + Clinic)"]])),
        unname(expected)
    )
})

test_that("ROR-C and ROR-PC groups are unchanged without clinical data", {
    input <- make_ror_input()
    res <- suppressMessages(BreastSubtypeR:::RORgroup(
        input$out,
        input$df.cln,
        Subtype = FALSE,
        hasClinical = FALSE
    ))
    expect_false("ROR-C Group (Subtype + Clinic)" %in% colnames(res))
    expect_identical(
        unname(as.character(res[["ROR-S Group (Subtype Only)"]])),
        c("high", "low", "med")
    )
})

make_node_ror_input <- function() {
    input <- make_ror_input()
    ids <- paste0("N", seq_len(4))
    input$out$testData <- input$out$testData[, c(1, 2, 3, 1)]
    colnames(input$out$testData) <- ids
    input$out$testData[c("CCNB1", "UBE2C", "MKI67"), ] <-
        rep(c(-2, -0.5, 0, 1), each = 3)
    input$out$distances <- input$out$distances[c(1, 2, 3, 1), ]
    rownames(input$out$distances) <- ids
    input$out$distances[, "LumA"] <- c(0.9, 0.5, -0.2, -0.9)
    input$out$dist.RORSubtype <- input$out$distances[, seq_len(4)]
    input$out$predictions <- setNames(c("LumA", "LumA", "Normal", "Basal"), ids)
    input$df.cln <- data.frame(PatientID = ids, TSIZE = 0, NODE = 0)
    input
}

test_that("NODE changes only the subtype-specific ROR-PC risk group", {
    input <- make_node_ror_input()
    baseline <- suppressMessages(BreastSubtypeR:::RORgroup(
        input$out, input$df.cln, hasClinical = TRUE
    ))
    score_cols <- c(
        "ROR-S (Subtype Only)", "ROR-P (Subtype + Prolif)",
        "ROR-C (Subtype + Clinic)", "ROR-PC (Subtype + Clinic + Prolif)",
        "ROR-PC (Subtype + Clinic + Prolif.Subtype)"
    )
    group_cols <- c(
        "ROR-S Group (Subtype Only)", "ROR-P Group (Subtype + Prolif)",
        "ROR-C Group (Subtype + Clinic)",
        "ROR-PC Group (Subtype + Clinic + Prolif)"
    )
    node_group <- "ROR-PC Group (Subtype + Clinic + Prolif.Subtype)"

    ## Fixed numerical expectations also protect the NODE-independent scores.
    expect_equal(unname(as.matrix(baseline[, score_cols])), rbind(
        c(3.86822953, -48.54809056, 13.55824471, -45.40136640, 11.00418748),
        c(20.44967000, 15.66693021, 25.83301176, 16.60612406, 34.22788655),
        c(49.46719082, 43.41530702, 47.31385412, 43.30438861, 51.81366476),
        c(78.48471165, 91.07613314, 68.79469647, 89.25320464, 74.82431242)
    ), tolerance = 1e-7)
    expect_identical(unname(as.matrix(baseline[, group_cols])), rbind(
        c("low", "low", "low", "low"),
        c("low", "med", "low", "low"),
        c("med", "med", "med", "med"),
        c("high", "high", "high", "high")
    ))

    ## Scores are away from all cutoffs: exact-cutoff handling is tracked
    ## separately in issue #154 and is not an accepted result in these tests.
    nodes <- c(0, 1, 3, 4, NA)
    expected_groups <- list(
        c("low", "low", "med", "high"),
        c("low", "med", "high", "high"),
        c("low", "med", "high", "high"),
        rep("high", 4),
        rep(NA_character_, 4)
    )
    for (i in seq_along(nodes)) {
        input$df.cln$NODE <- nodes[i]
        res <- suppressMessages(BreastSubtypeR:::RORgroup(
            input$out, input$df.cln, hasClinical = TRUE
        ))
        info <- paste("NODE =", nodes[i])
        expect_identical(res[, score_cols], baseline[, score_cols], info = info)
        expect_identical(res[, group_cols], baseline[, group_cols], info = info)
        expect_identical(unname(as.character(res[[node_group]])),
            expected_groups[[i]], info = info
        )
    }
})

test_that("missing TSIZE removes clinical ROR results regardless of NODE", {
    input <- make_node_ror_input()
    genomic_cols <- c(
        "ROR-S (Subtype Only)", "ROR-S Group (Subtype Only)",
        "ROR-P (Subtype + Prolif)", "ROR-P Group (Subtype + Prolif)"
    )
    clinical_cols <- c(
        "ROR-C (Subtype + Clinic)", "ROR-C Group (Subtype + Clinic)",
        "ROR-PC (Subtype + Clinic + Prolif)",
        "ROR-PC Group (Subtype + Clinic + Prolif)",
        "ROR-PC (Subtype + Clinic + Prolif.Subtype)",
        "ROR-PC Group (Subtype + Clinic + Prolif.Subtype)"
    )
    for (node in c(0, 1, 3, 4, NA)) {
        input$df.cln$NODE <- node
        input$df.cln$TSIZE <- 0
        baseline <- suppressMessages(BreastSubtypeR:::RORgroup(
            input$out, input$df.cln, hasClinical = TRUE
        ))
        input$df.cln$TSIZE <- c(NA, 0, NA, 0)
        res <- suppressMessages(BreastSubtypeR:::RORgroup(
            input$out, input$df.cln, hasClinical = TRUE
        ))
        info <- paste("NODE =", node)
        expect_identical(res[, genomic_cols], baseline[, genomic_cols], info = info)
        expect_true(all(is.na(res[c(1, 3), clinical_cols])), info = info)
        expect_identical(res[c(2, 4), clinical_cols],
            baseline[c(2, 4), clinical_cols], info = info
        )
    }
})
