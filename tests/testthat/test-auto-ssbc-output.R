auto233_input <- function() {
    e <- new.env(parent = emptyenv())
    utils::data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = e)
    e$OSLO2EMIT0obj$data_input
}
auto233_subset <- function(x, i) lapply(x, function(z) z[, i, drop = FALSE])
auto233_annotate <- function(x, er, her2, tn = "nonTN") {
    for (n in names(x)) {
        SummarizedExperiment::colData(x[[n]])$ER <- er
        SummarizedExperiment::colData(x[[n]])$HER2 <- her2
        SummarizedExperiment::colData(x[[n]])$TN <- rep(tn, length.out = ncol(x[[n]]))
    }
    x
}
auto233_run <- function(input, subtype = TRUE, clinical = TRUE, methods = "AUTO") {
    suppressMessages(expect_experimental_ror_pc(BS_Multi(input, methods = methods,
        Subtype = subtype, hasClinical = clinical), n = as.integer(clinical &&
            any(c("ssBC", "ssBC.v2") %in% if (identical(methods, "AUTO"))
                BreastSubtypeR:::get_methods(as.data.frame(SummarizedExperiment::colData(input$se_NC)))$methods
            else methods))))
}

test_that("AUTO restores each component's own identities without changing payloads", {
    # Independent artificial outputs deliberately have different coverage/order
    # per component. These are assembly tests, not invented scientific scores.
    ids <- c("a", "B", "c", "D")
    b <- data.frame(PatientID = c("c", "a"), BS = c("Basal", NA_character_), row.names = c("c", "a"))
    scores <- data.frame(patientID = c("B", NA_character_), value = c(2, NA_real_),
        count = c(1L, NA_integer_), group = factor(c("med", NA), levels = c("low", "med", "high")),
        flag = c(TRUE, NA), row.names = c("B", "c"))
    mat <- matrix(c(NA_real_, 2, 3, NaN), 2L, dimnames = list(c("B", "a"), c("g2", "g1")))
    native <- list(BS.all = b, score.ROR = scores, mdns = "gene metadata",
        outList = list(predictions = c(c = "Basal", a = NA_character_),
            predictions.FourSubtype = c(B = "Her2", a = NA_character_),
            distances = mat, dist.RORSubtype = mat[2:1, , drop = FALSE],
            testData = t(mat), centroids = mat))
    before <- serialize(native, NULL)
    actual <- BreastSubtypeR:::.restore_auto_ssbc(native, ids)
    expect_auto_ssbc_axes(actual, native, ids)
    expect_identical(actual$score.ROR$patientID, c("a", "B", NA_character_, "D"))
    expect_identical(actual$score.ROR$value, c(NA_real_, 2, NA_real_, NA_real_))
    expect_identical(actual$score.ROR$count, c(NA_integer_, 1L, NA_integer_, NA_integer_))
    expect_identical(actual$score.ROR$flag, c(NA, TRUE, NA, NA))
    expect_true(is.na(actual$BS.all["a", "BS"]))
    expect_true(any(!is.na(actual$outList$testData[, "a"])))
    expect_identical(serialize(native, NULL), before)
    native$outList$predictions.FourSubtype <- NULL
    expect_false("predictions.FourSubtype" %in%
        names(BreastSubtypeR:::.restore_auto_ssbc(native, ids)$outList))
    for (idx in list(integer(), 1L)) {
        tiny <- native
        tiny$BS.all <- b[idx, , drop = FALSE]
        tiny$score.ROR <- scores[idx, , drop = FALSE]
        tiny$outList$predictions <- native$outList$predictions[idx]
        tiny$outList$distances <- mat[idx, , drop = FALSE]
        tiny$outList$dist.RORSubtype <- mat[idx, , drop = FALSE]
        tiny$outList$testData <- t(mat)[, idx, drop = FALSE]
        expect_auto_ssbc_axes(BreastSubtypeR:::.restore_auto_ssbc(tiny, ids), tiny, ids)
    }
    for (bad in list(c("a", "a"), c("a", NA_character_), c("a", "")))
        expect_error(BreastSubtypeR:::.restore_auto_ssbc(native, bad), "sample IDs")
    for (field in c("predictions", "distances", "testData")) {
        for (bad in list(c("a", "a"), c("a", "unknown"), c("a", NA_character_), c("a", ""), NULL)) {
            broken <- native
            if (field == "predictions") names(broken$outList[[field]]) <- bad else
                if (field == "testData") colnames(broken$outList[[field]]) <- bad else
                    rownames(broken$outList[[field]]) <- bad
            expect_error(BreastSubtypeR:::.restore_auto_ssbc(broken, ids), "sample IDs")
        }
    }
    broken <- native; broken$score.ROR$patientID[1] <- "a"
    expect_error(BreastSubtypeR:::.restore_auto_ssbc(broken, ids), "conflicting IDs")
    broken <- native; broken$BS.all$PatientID <- c("a", "a")
    expect_error(BreastSubtypeR:::.restore_auto_ssbc(broken, ids), "sample IDs")
    broken <- native; broken$outList$testData <- mat[1, ]
    expect_error(BreastSubtypeR:::.restore_auto_ssbc(broken, ids), "native matrix dimensions")
})

test_that("OSLO AUTO preserves native evaluated values for all flags and ID orders", {
    withr::local_seed(233)
    original <- auto233_input()
    orders <- list(seq_len(102), 102:1, c(seq(2, 102, 2), seq(1, 101, 2)))
    for (i in orders) for (factor_ids in c(FALSE, TRUE))
        for (subtype in c(FALSE, TRUE)) for (clinical in c(FALSE, TRUE)) {
            input <- auto233_subset(original, i)
            if (factor_ids) for (n in names(input))
                SummarizedExperiment::colData(input[[n]])$PatientID <-
                    factor(colnames(input[[n]]), levels = c("unused", rev(colnames(input[[n]]))))
            before <- input; seed <- .Random.seed; locale <- Sys.getlocale(); kind <- RNGkind()
            actual <- auto233_run(input, subtype, clinical)
            # Disable only the new post-calculation assembly to compare unchanged
            # BS.all, aggregates, entropy and every unaffected method as objects.
            baseline <- testthat::with_mocked_bindings(auto233_run(input, subtype, clinical),
                .restore_auto_ssbc = function(result, ids) result, .package = "BreastSubtypeR")
            expect_identical(actual[names(actual) != "results"], baseline[names(baseline) != "results"])
            expect_identical(names(actual$results), names(baseline$results))
            for (method in names(actual$results)) {
                if (method %in% c("ssBC", "ssBC.v2")) {
                    expect_auto_ssbc_axes(actual$results[[method]], baseline$results[[method]], colnames(input$se_NC))
                    expect_identical(actual$results[[method]]$BS.all, baseline$results[[method]]$BS.all)
                } else expect_identical(actual$results[[method]], baseline$results[[method]])
            }
            expect_identical(input, before)
            expect_identical(.Random.seed, seed)
            expect_identical(RNGkind(), kind)
            expect_identical(Sys.getlocale(), locale)
        }
})

test_that("AUTO thresholds, unavailable annotations and absent methods keep their contracts", {
    original <- auto233_input()
    missing <- original
    for (n in names(missing)) {
        SummarizedExperiment::colData(missing[[n]])$ER[c(1, 3)] <- NA_character_
        SummarizedExperiment::colData(missing[[n]])$HER2[c(4, 6)] <- NA_character_
    }
    cases <- list(missing = missing,
        both = auto233_annotate(auto233_subset(original, 1:30), rep(c("ER+", "ER-"), c(20, 10)),
            rep(c("HER2+", "HER2-", "HER2+", "HER2-"), c(10, 10, 5, 5))),
        below = auto233_annotate(auto233_subset(original, 1:14), rep("ER+", 14), rep(c("HER2+", "HER2-"), c(7, 7))),
        at = auto233_annotate(auto233_subset(original, 1:15), rep("ER+", 15), rep(c("HER2+", "HER2-"), c(8, 7))),
        neg_below = auto233_annotate(auto233_subset(original, 1:17), rep("ER-", 17), rep("HER2-", 17)),
        neg_at = auto233_annotate(auto233_subset(original, 1:18), rep("ER-", 18), rep("HER2-", 18)),
        tn_unknown = auto233_annotate(auto233_subset(original, 1:19), rep("ER-", 19), rep("HER2-", 19), c(rep("TN", 18), NA)),
        her2 = auto233_annotate(auto233_subset(original, 1:17), rep(c("ER+", "ER-"), c(8, 9)), rep("HER2+", 17)))
    expected_methods <- list(missing = c("ssBC", "ssBC.v2"), both = c("ssBC", "ssBC.v2"),
        below = character(), at = "ssBC.v2", neg_below = character(), neg_at = "ssBC.v2",
        tn_unknown = c("ssBC", "ssBC.v2"), her2 = "ssBC.v2")
    for (case in names(cases)) for (factor_ids in c(FALSE, TRUE)) {
        input <- cases[[case]]
        if (factor_ids) for (n in names(input))
            SummarizedExperiment::colData(input[[n]])$PatientID <- factor(colnames(input[[n]]),
                levels = c("unused", rev(colnames(input[[n]]))))
        actual <- suppressWarnings(auto233_run(input))
        baseline <- suppressWarnings(testthat::with_mocked_bindings(auto233_run(input),
            .restore_auto_ssbc = function(result, ids) result, .package = "BreastSubtypeR"))
        expect_identical(intersect(c("ssBC", "ssBC.v2"), names(actual$results)), expected_methods[[case]])
        expect_identical(actual[names(actual) != "results"], baseline[names(baseline) != "results"])
        for (method in expected_methods[[case]])
            expect_auto_ssbc_axes(actual$results[[method]], baseline$results[[method]], colnames(input$se_NC))
        for (method in setdiff(names(actual$results), c("ssBC", "ssBC.v2")))
            expect_identical(actual$results[[method]], baseline$results[[method]])
        if (case == "tn_unknown") {
            for (method in expected_methods[[case]]) {
                a <- actual$results[[method]]
                expect_true(is.na(a$BS.all$BS[19]))
                expect_true(is.na(a$score.ROR$patientID[19]))
            }
        }
    }
    # Invalid input labels must fail before any classifier can run.
    bad <- original
    SummarizedExperiment::colData(bad$se_NC)$PatientID <- rev(colnames(bad$se_NC))
    expect_error(BS_Multi(bad), "PatientID labels to match")
})

test_that("plot safeguards remove only wholly unavailable measurements from copies", {
    result <- auto233_run(auto233_input())$results$ssBC.v2
    before <- serialize(result, NULL)
    x <- result$outList$testData
    correlations <- result$outList$distances
    out <- data.frame(PatientID = result$BS.all$PatientID, Subtype = result$BS.all$BS)
    keep <- colSums(!is.na(x)) > 0L
    expect_identical(sum(!keep), 8L)
    # A missing call with available expression/correlations must be retained.
    out$Subtype[which(keep)[1]] <- NA_character_
    shuffled <- out[nrow(out):1, , drop = FALSE]
    shuffled$PatientID <- factor(shuffled$PatientID, levels = c("unused", out$PatientID))
    expect_message(p <- Vis_PCA(x, shuffled), "excluded 8 sample")
    reference <- stats::prcomp(t(x[, keep, drop = FALSE]), center = TRUE, scale. = TRUE)
    expect_identical(p$data[, colnames(reference$x)], as.data.frame(reference$x))
    expect_identical(p$data$Subtype, out$Subtype[keep])
    expect_message(h <- Vis_heatmap(x, shuffled), "excluded 8 sample")
    selected <- as.character(shuffled$PatientID[shuffled$PatientID %in% colnames(x)[keep]])
    expect_identical(h@matrix, t(scale(t(x[, selected, drop = FALSE]))))
    expect_identical(as.character(get("value", h@top_annotation@anno_list$Subtype@fun@var_env)),
        shuffled$Subtype[shuffled$PatientID %in% selected])
    pdf <- tempfile(fileext = ".pdf")
    grDevices::pdf(pdf)
    tryCatch(expect_s4_class(ComplexHeatmap::draw(h), "HeatmapList"),
             finally = { grDevices::dev.off(); unlink(pdf) })
    expect_message(b <- Vis_boxplot(shuffled, correlations), "excluded 8 sample")
    expect_identical(b$data$predictions, out$Subtype[keep])
    expect_identical(b$data$cor, unname(apply(correlations[keep, , drop = FALSE], 1, max)))
    expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
    expect_s3_class(suppressWarnings(ggplot2::ggplot_build(b)), "ggplot_built")
    expect_identical(serialize(result, NULL), before)
    # Partial missingness is retained: preserve existing PCA/heatmap errors and
    # row-maximum NA semantics, rather than quietly doing complete-case analysis.
    partial <- x[, keep, drop = FALSE]; partial[1, 1] <- NA_real_
    expect_error(Vis_PCA(partial, out[keep, ]))
    expect_error(Vis_heatmap(partial, out[keep, ]))
    partial_cor <- correlations[keep, , drop = FALSE]; partial_cor[1, 1] <- NA_real_
    expect_true(is.na(Vis_boxplot(out[keep, ], partial_cor)$data$cor[1]))
    for (remaining in c(0L, 1L)) {
        too_small <- x; too_small[,] <- NA_real_
        if (remaining) too_small[, which(keep)[1]] <- x[, which(keep)[1]]
        expect_error(suppressMessages(Vis_PCA(too_small, out)), "at least two")
        expect_error(suppressMessages(Vis_heatmap(too_small, out)), "at least two")
    }
    correlations[,] <- NA_real_
    expect_error(suppressMessages(Vis_boxplot(out, correlations)), "at least one")
})
