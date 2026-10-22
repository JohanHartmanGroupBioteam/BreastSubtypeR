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

    ## These scores are away from cutoffs. Resolved NODE boundaries (#154,
    ## implemented in #157) are tested below. The original Parker equality
    ## limitation remains separately deferred in #172.
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

test_that("subtype-specific ROR-PC endpoints include the upper cutoff", {
    ## Adjacent IEEE double values make equality distinct from either side.
    ## Test categorization directly: the unchanged final score multiplication
    ## skips exact doubles 15 and 60, whereas exact 40 is attainable below.
    score <- c(15 - 2^-49, 15, 15 + 2^-49,
        40 - 2^-47, 40, 40 + 2^-47,
        60 - 2^-47, 60, 60 + 2^-47, NA_real_)
    names(score) <- paste0("B", seq_along(score))
    expected_n0 <- c("low", "low", "low", "low", "low", "med",
        "med", "med", "high", NA_character_)
    expected_n1_3 <- c("low", "low", "med", "med", "med", "high",
        "high", "high", "high", NA_character_)
    expect_identical(BreastSubtypeR:::.ror_pc_subtype_labels(score, 40, 60),
        setNames(expected_n0, names(score)))
    expect_identical(BreastSubtypeR:::.ror_pc_subtype_labels(score, 15, 40),
        setNames(expected_n1_3, names(score)))
    expect_identical(BreastSubtypeR:::.ror_pc_subtype_labels(
        c(NA_real_, NaN, -Inf, Inf), 15, 40),
        c(NA_character_, NA_character_, "low", "high"))
})

test_that("undefined additional ROR-PC scores have missing categories", {
    input <- make_node_ror_input()
    score <- "ROR-PC (Subtype + Clinic + Prolif.Subtype)"
    group <- "ROR-PC Group (Subtype + Clinic + Prolif.Subtype)"
    for (node in c(0, 1, 2, 3, 4, NA_real_)) {
        input$df.cln$NODE <- node
        for (tsize in c(0, NA_real_)) {
            input$df.cln$TSIZE <- tsize
            reference <- suppressMessages(BreastSubtypeR:::RORgroup(
                input$out, input$df.cln, hasClinical = TRUE))
            undefined <- input$out
            undefined$dist.RORSubtype[1, 1] <- NaN
            undefined$dist.RORSubtype[2, 1] <- NA_real_
            actual <- suppressMessages(BreastSubtypeR:::RORgroup(
                undefined, input$df.cln, hasClinical = TRUE))
            info <- paste("NODE", node, "TSIZE", tsize)
            other <- setdiff(names(actual), c(score, group))
            expect_identical(actual[, other], reference[, other], info = info)
            expect_identical(actual[3:4, c(score, group)],
                reference[3:4, c(score, group)], info = info)
            if (!is.na(tsize)) {
                expect_true(is.nan(actual[[score]][1]), info = info)
                expect_true(is.na(actual[[score]][2]), info = info)
                expect_false(is.nan(actual[[score]][2]), info = info)
            }
            if (!is.na(node) && node >= 4 && !is.na(tsize)) {
                ## The independent NODE >= 4 rule is deliberately unchanged.
                expect_identical(unname(actual[[group]][1:2]),
                    c("high", "high"), info = info)
            } else {
                expect_true(all(is.na(actual[[group]][1:2])), info = info)
                expect_false(any(actual[[group]][1:2] %in% c("NA", "NaN")),
                    info = info)
            }
            path <- tempfile(fileext = ".txt")
            utils::write.table(actual, path, sep = "\t", quote = FALSE,
                row.names = FALSE)
            restored <- utils::read.delim(path, check.names = FALSE)
            unlink(path)
            expect_identical(is.na(restored[[group]]), is.na(actual[[group]]),
                info = info)
            expect_false(any(restored[[group]] %in% c("NA", "NaN")), info = info)
        }
    }
})

make_boundary_ror_input <- function() {
    input <- make_ror_input()
    ## Fixed synthetic proliferation values yield scores immediately around
    ## 15/40/60, exactly 40, and four nonboundary controls. No score is rounded.
    prolif <- c(-3.0728038441551813, -3.0728038441551808,
        -0.76860092918918199, -0.76860092918918088, -0.76860092918918044,
        1.0747614027836181, 1.074761402783619,
        -3.5336444271483809, -2.1511226781687807,
        0.15308023679721919, 2.4572831517632188)
    ids <- paste0("B", seq_along(prolif))
    input$out$testData <- input$out$testData[, rep(1L, length(ids))]
    colnames(input$out$testData) <- ids
    input$out$testData[c("CCNB1", "UBE2C", "MKI67"), ] <- rep(prolif, each = 3)
    input$out$distances <- matrix(0, nrow = length(ids), ncol = 5,
        dimnames = list(ids, colnames(input$out$centroids)))
    input$out$dist.RORSubtype <- input$out$distances[, seq_len(4)]
    input$out$predictions <- setNames(rep("LumA", length(ids)), ids)
    input$df.cln <- data.frame(PatientID = ids, TSIZE = 0, NODE = 0)
    input
}

test_that("NODE boundary labels leave scores and other outputs unchanged", {
    input <- make_boundary_ror_input()
    score_col <- "ROR-PC (Subtype + Clinic + Prolif.Subtype)"
    group_col <- "ROR-PC Group (Subtype + Clinic + Prolif.Subtype)"
    baseline <- suppressMessages(BreastSubtypeR:::RORgroup(
        input$out, input$df.cln, hasClinical = TRUE
    ))
    ## Compare exact doubles, including the genuinely exact score 40.
    expect_identical(unname(baseline[[score_col]]),
        c(15 - 2^-48, 15 + 2^-49, 40 - 2^-47, 40, 40 + 2^-47,
            60 - 2^-47, 60 + 2^-47, 10 - 2^-49, 25, 50, 75))
    expect_identical(unname(as.character(baseline$Call)), rep("LumA", 11))
    standard_groups <- c("ROR-S Group (Subtype Only)",
        "ROR-P Group (Subtype + Prolif)", "ROR-C Group (Subtype + Clinic)",
        "ROR-PC Group (Subtype + Clinic + Prolif)")
    expected_prolif <- c(rep("low", 5), "high", "high",
        "low", "low", "med", "high")
    expect_identical(unname(as.matrix(baseline[, standard_groups])),
        unname(cbind(rep("med", 11), expected_prolif, rep("med", 11), expected_prolif)))
    expected_n0 <- c("low", "low", "low", "low", "med", "med", "high",
        "low", "low", "med", "high")
    expected_n1_3 <- c("low", "med", "med", "med", "high", "high", "high",
        "low", "med", "high", "high")
    nodes <- c(0, 1, 2, 3, 4, NA_real_)
    expected <- list(expected_n0, expected_n1_3, expected_n1_3,
        expected_n1_3, rep("high", 11), rep(NA_character_, 11))
    other_cols <- setdiff(names(baseline), group_col)
    for (i in seq_along(nodes)) {
        input$df.cln$NODE <- nodes[i]
        res <- suppressMessages(BreastSubtypeR:::RORgroup(
            input$out, input$df.cln, hasClinical = TRUE
        ))
        info <- paste("NODE =", nodes[i])
        expect_identical(res[, other_cols], baseline[, other_cols], info = info)
        expect_identical(unname(as.character(res[[group_col]])), expected[[i]],
            info = info)
        expect_true(all(is.na(res[[group_col]]) |
            res[[group_col]] %in% c("low", "med", "high")), info = info)
    }
})

test_that("boundary correction preserves missing clinical output", {
    input <- make_boundary_ror_input()
    group_col <- "ROR-PC Group (Subtype + Clinic + Prolif.Subtype)"
    clinical_cols <- c("ROR-C (Subtype + Clinic)", "ROR-C Group (Subtype + Clinic)",
        "ROR-PC (Subtype + Clinic + Prolif)",
        "ROR-PC Group (Subtype + Clinic + Prolif)",
        "ROR-PC (Subtype + Clinic + Prolif.Subtype)", group_col)
    for (node in c(0, 1, 2, 3, 4, NA_real_)) {
        input$df.cln$NODE <- node
        input$df.cln$TSIZE <- 0
        complete <- suppressMessages(BreastSubtypeR:::RORgroup(
            input$out, input$df.cln, hasClinical = TRUE
        ))
        input$df.cln$TSIZE <- NA_real_
        missing <- suppressMessages(BreastSubtypeR:::RORgroup(
            input$out, input$df.cln, hasClinical = TRUE
        ))
        info <- paste("NODE =", node)
        expect_true(all(is.na(missing[, clinical_cols])), info = info)
        other_cols <- setdiff(names(complete), clinical_cols)
        expect_identical(missing[, other_cols], complete[, other_cols], info = info)
        no_clinical <- suppressMessages(BreastSubtypeR:::RORgroup(
            input$out, input$df.cln, hasClinical = FALSE
        ))
        expect_false(any(clinical_cols %in% names(no_clinical)), info = info)
    }
})

## Invalid NODE counts (#187): negative, fractional or infinite values must
## give a genuine missing additional group with one warning per call, while
## NA/NaN stay silently missing and every other output is untouched.
node_group_col <- "ROR-PC Group (Subtype + Clinic + Prolif.Subtype)"
invalid_node_pattern <- "finite nonnegative integer"

run_ror_capturing <- function(input, hasClinical = TRUE) {
    warnings <- character()
    res <- withCallingHandlers(
        suppressMessages(BreastSubtypeR:::RORgroup(
            input$out, input$df.cln, hasClinical = hasClinical
        )),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    list(res = res, warnings = warnings)
}

expect_missing_group <- function(x, info = NULL) {
    expect_type(x, "character")
    expect_true(all(is.na(x)), info = info)
    expect_false(any(x %in% c("NA", "NaN")), info = info)
    expect_false(any(grepl("[0-9]", x)), info = info)
}

test_that("valid NODE counts are accepted as integers and integer-valued doubles", {
    input <- make_node_ror_input()
    baseline <- run_ror_capturing(input)
    expect_length(baseline$warnings, 0)
    expected <- list(
        `0` = c("low", "low", "med", "high"),
        `1` = c("low", "med", "high", "high"),
        `2` = c("low", "med", "high", "high"),
        `3` = c("low", "med", "high", "high"),
        `4` = rep("high", 4),
        `12` = rep("high", 4)
    )
    other <- setdiff(names(baseline$res), node_group_col)
    for (node in names(expected)) {
        for (storage in c("double", "integer")) {
            value <- as.numeric(node)
            if (storage == "integer") value <- as.integer(value)
            input$df.cln$NODE <- value
            out <- run_ror_capturing(input)
            info <- paste("NODE =", node, storage)
            expect_length(out$warnings, 0)
            expect_type(out$res[[node_group_col]], "character")
            expect_identical(unname(out$res[[node_group_col]]),
                expected[[node]], info = info)
            expect_identical(out$res[, other], baseline$res[, other], info = info)
        }
    }
})

test_that("missing NODE gives a silent missing character group", {
    input <- make_node_ror_input()
    reference <- run_ror_capturing(input)
    other <- setdiff(names(reference$res), node_group_col)
    for (node in list(NA_real_, NaN, c(NA, 0, NaN, 4))) {
        input$df.cln$NODE <- node
        out <- run_ror_capturing(input)
        info <- paste("NODE =", paste(node, collapse = ","))
        expect_length(out$warnings, 0)
        expect_type(out$res[[node_group_col]], "character")
        expect_identical(out$res[, other], reference$res[, other], info = info)
    }
    input$df.cln$NODE <- NA_real_
    expect_missing_group(run_ror_capturing(input)$res[[node_group_col]])
    input$df.cln$NODE <- NaN
    expect_missing_group(run_ror_capturing(input)$res[[node_group_col]])
    input$df.cln$NODE <- c(NA, 0, NaN, 4)
    mixed <- run_ror_capturing(input)$res[[node_group_col]]
    expect_identical(unname(mixed), c(NA, "low", NA, "high"))
})

test_that("invalid NODE values give one warning and missing groups", {
    input <- make_node_ror_input()
    reference <- run_ror_capturing(input)
    other <- setdiff(names(reference$res), node_group_col)
    invalid <- c(-1, -0.5, 0.5, 2.5, 4.5, Inf, -Inf)
    for (value in invalid) {
        input$df.cln$NODE <- value
        out <- run_ror_capturing(input)
        info <- paste("NODE =", value)
        expect_length(out$warnings, 1)
        expect_match(out$warnings, invalid_node_pattern, info = info)
        expect_match(out$warnings, "4 sample", info = info)
        expect_missing_group(out$res[[node_group_col]], info = info)
        expect_identical(out$res[, other], reference$res[, other], info = info)
        expect_identical(out$res$patientID, reference$res$patientID, info = info)
        expect_identical(dim(out$res), dim(reference$res), info = info)
        expect_identical(names(out$res), names(reference$res), info = info)
    }
})

test_that("mixed valid and invalid NODE values are handled per sample", {
    input <- make_node_ror_input()
    reference <- run_ror_capturing(input)
    other <- setdiff(names(reference$res), node_group_col)
    cases <- list(
        list(node = c(-1, 0.5, 4.5, Inf), groups = rep(NA_character_, 4), n = 4),
        list(node = c(-0.5, 2.5, -Inf, 2), groups = c(NA, NA, NA, "high"), n = 3),
        list(node = c(0, -1, 3, 4), groups = c("low", NA, "high", "high"), n = 1),
        list(node = c(2.5, 1, NA, 4), groups = c(NA, "med", NA, "high"), n = 1)
    )
    for (case in cases) {
        input$df.cln$NODE <- case$node
        out <- run_ror_capturing(input)
        info <- paste("NODE =", paste(case$node, collapse = ","))
        expect_length(out$warnings, 1)
        expect_match(out$warnings, paste(case$n, "sample"), info = info)
        group <- out$res[[node_group_col]]
        expect_type(group, "character")
        expect_identical(unname(group), case$groups, info = info)
        expect_false(any(group %in% c("NA", "NaN")), info = info)
        expect_false(any(grepl("[0-9]", group)), info = info)
        expect_identical(out$res[, other], reference$res[, other], info = info)
    }
})

test_that("invalid NODE groups stay missing after text round-trip", {
    input <- make_node_ror_input()
    input$df.cln$NODE <- c(-1, 2.5, Inf, 2)
    out <- run_ror_capturing(input)$res
    path <- tempfile(fileext = ".txt")
    utils::write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)
    restored <- utils::read.delim(path, check.names = FALSE)
    unlink(path)
    expect_identical(is.na(restored[[node_group_col]]), c(TRUE, TRUE, TRUE, FALSE))
    expect_false(any(restored[[node_group_col]] %in% c("NA", "NaN")))
})

test_that("shuffled clinical rows keep NODE validity aligned by PatientID", {
    input <- make_node_ror_input()
    input$df.cln$NODE <- c(0, -1, 3, 4)
    ordered <- run_ror_capturing(input)
    shuffled <- input
    shuffled$df.cln <- shuffled$df.cln[c(3, 1, 4, 2), ]
    rownames(shuffled$df.cln) <- NULL
    out <- run_ror_capturing(shuffled)
    expect_length(out$warnings, 1)
    expect_identical(out$res, ordered$res)
    expect_identical(unname(out$res[[node_group_col]]),
        c("low", NA, "high", "high"))
    expect_identical(out$res$patientID, paste0("N", 1:4))
})

test_that("missing TSIZE overrides valid and invalid NODE alike", {
    input <- make_node_ror_input()
    input$df.cln$NODE <- c(4, -1, 2.5, 0)
    input$df.cln$TSIZE <- c(NA, NA, 0, 0)
    out <- run_ror_capturing(input)
    expect_length(out$warnings, 1)
    group <- out$res[[node_group_col]]
    expect_type(group, "character")
    expect_identical(unname(group), c(NA, NA, NA, "high"))
    expect_true(is.na(out$res[["ROR-C (Subtype + Clinic)"]][1]))
    expect_true(is.na(out$res[["ROR-PC (Subtype + Clinic + Prolif.Subtype)"]][2]))
})

test_that("valid NODE >= 4 keeps high risk with an undefined additional score", {
    input <- make_node_ror_input()
    input$out$dist.RORSubtype[1, 1] <- NaN
    input$out$dist.RORSubtype[2, 1] <- NA_real_
    score <- "ROR-PC (Subtype + Clinic + Prolif.Subtype)"
    input$df.cln$NODE <- c(4, 5L, -1, 4.5)
    out <- run_ror_capturing(input)
    expect_length(out$warnings, 1)
    expect_true(is.nan(out$res[[score]][1]))
    expect_true(is.na(out$res[[score]][2]))
    expect_identical(unname(out$res[[node_group_col]]),
        c("high", "high", NA, NA))
})

test_that("invalid NODE is irrelevant when hasClinical = FALSE", {
    input <- make_node_ror_input()
    reference <- run_ror_capturing(input, hasClinical = FALSE)
    input$df.cln$NODE <- c(-1, 2.5, Inf, -Inf)
    out <- run_ror_capturing(input, hasClinical = FALSE)
    expect_length(out$warnings, 0)
    expect_identical(out$res, reference$res)
    expect_false(node_group_col %in% names(out$res))
})

test_that("all-invalid and all-missing batches return a character group column", {
    input <- make_node_ror_input()
    reference <- run_ror_capturing(input)
    other <- setdiff(names(reference$res), node_group_col)
    for (node in list(rep(-1, 4), rep(Inf, 4), rep(0.5, 4), rep(NA_real_, 4))) {
        input$df.cln$NODE <- node
        out <- run_ror_capturing(input)
        info <- paste("NODE =", paste(node, collapse = ","))
        expect_missing_group(out$res[[node_group_col]], info = info)
        expect_identical(out$res[, other], reference$res[, other], info = info)
        expect_identical(nrow(out$res), 4L, info = info)
    }
})
