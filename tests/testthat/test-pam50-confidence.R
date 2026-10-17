# These synthetic inputs distinguish the original Spearman p-value default
# from exact = FALSE without copying the author-supplied example files.
confidence_test_input <- function() {
    genes <- paste0("gene", seq_len(5))
    ids <- c("sample 1/A", "sample-2", "sample.3")
    test_data <- cbind(1:5, c(4, 1, 5, 2, 3), c(2, 5, 3, 1, 4))
    dimnames(test_data) <- list(genes, ids)
    centroids <- cbind(LumA = c(1, 2, 3, 5, 4), Basal = c(4, 1, 5, 3, 2))
    rownames(centroids) <- genes
    list(
        predictions = setNames(c("LumA", "Basal", "LumA"), ids),
        testData = test_data,
        centroids = centroids
    )
}

confidence_direct <- function(x, y) {
    suppressWarnings(round(1 - cor.test(x, y, method = "spearman")$p.value, 2))
}

test_that("Confidence supports matrix and data.frame centroids", {
    out <- confidence_test_input()
    expected <- c(
        confidence_direct(out$testData[, 1], out$centroids[, "LumA"]),
        confidence_direct(out$testData[, 2], out$centroids[, "Basal"]),
        confidence_direct(out$testData[, 3], out$centroids[, "LumA"])
    )
    for (centroids in list(out$centroids, as.data.frame(out$centroids))) {
        out$centroids <- centroids
        actual <- BreastSubtypeR:::.pam50_confidence(out)
        expect_identical(actual, expected)
        expect_type(actual, "double")
        expect_null(names(actual))
        expect_length(actual, length(out$predictions))
    }
    expect_null(names(as.data.frame(out$centroids)[, "LumA"]))
})

test_that("Confidence preserves default p-values and two-decimal rounding", {
    out <- confidence_test_input()
    actual <- BreastSubtypeR:::.pam50_confidence(out)
    approximate <- round(1 - cor.test(1:5, c(1, 2, 3, 5, 4),
        method = "spearman", exact = FALSE)$p.value, 2)
    expect_identical(actual[1], 0.92)
    expect_identical(approximate, 0.96)
    expect_false(identical(actual[1], approximate))
    expect_identical(actual, round(actual, 2))
})

test_that("Confidence retains single-sample dimensions and values", {
    out <- confidence_test_input()
    expected <- BreastSubtypeR:::.pam50_confidence(out)
    for (i in seq_along(out$predictions)) {
        single <- out
        single$predictions <- out$predictions[i]
        single$testData <- out$testData[, i, drop = FALSE]
        for (centroids in list(out$centroids, as.data.frame(out$centroids))) {
            single$centroids <- centroids
            expect_identical(BreastSubtypeR:::.pam50_confidence(single),
                expected[i])
        }
    }
})

test_that("Confidence aligns reordered genes and centroid columns", {
    out <- confidence_test_input()
    expected <- BreastSubtypeR:::.pam50_confidence(out)
    out$testData <- out$testData[c(4, 2, 5, 1, 3), , drop = FALSE]
    out$centroids <- out$centroids[c(3, 1, 4, 5, 2), 2:1, drop = FALSE]
    expect_identical(BreastSubtypeR:::.pam50_confidence(out), expected)
    out$centroids <- as.data.frame(out$centroids)
    expect_identical(BreastSubtypeR:::.pam50_confidence(out), expected)
})

test_that("Confidence uses only verified common genes with partial overlap", {
    out <- confidence_test_input()
    out$testData <- rbind(out$testData, test_only = c(100, 200, 300))
    out$centroids <- rbind(out$centroids[5:2, , drop = FALSE],
        centroid_only = c(10, 20))
    expected <- vapply(seq_along(out$predictions), function(j) {
        confidence_direct(
            out$testData[paste0("gene", 2:5), j],
            out$centroids[paste0("gene", 2:5), out$predictions[j]]
        )
    }, numeric(1))
    expect_identical(BreastSubtypeR:::.pam50_confidence(out), expected)
    out$centroids <- as.data.frame(out$centroids)
    expect_identical(BreastSubtypeR:::.pam50_confidence(out), expected)
})

test_that("Confidence preserves missing-pair and tied-rank semantics", {
    out <- confidence_test_input()
    out$testData[2, 1] <- NA_real_
    out$testData[4, 2] <- NaN
    out$testData[, 3] <- c(1, 1, 2, 3, 4)
    out$centroids[5, "LumA"] <- NA_real_
    expected <- vapply(seq_along(out$predictions), function(j) {
        confidence_direct(out$testData[, j],
            out$centroids[, out$predictions[j]])
    }, numeric(1))
    expect_silent(actual <- BreastSubtypeR:::.pam50_confidence(out))
    expect_identical(actual, expected)
    expect_true(all(is.finite(actual)))
})

test_that("Confidence ranks infinities as the original calculation does", {
    out <- confidence_test_input()
    out$testData[3, 1] <- Inf
    out$testData[1, 2] <- -Inf
    out$centroids[2, "LumA"] <- Inf
    expected <- vapply(seq_along(out$predictions), function(j) {
        confidence_direct(out$testData[, j],
            out$centroids[, out$predictions[j]])
    }, numeric(1))
    expect_identical(BreastSubtypeR:::.pam50_confidence(out), expected)
    expect_true(all(is.finite(expected)))
})

test_that("Undefined Confidence is NA while a computed zero remains zero", {
    out <- confidence_test_input()
    out$testData[, 1] <- NA_real_
    out$testData[, 2] <- c(1, rep(NA_real_, 4))
    out$testData[, 3] <- 1
    expect_silent(actual <- BreastSubtypeR:::.pam50_confidence(out))
    expect_identical(actual, rep(NA_real_, 3))

    out <- confidence_test_input()
    out$centroids[, "LumA"] <- 1
    actual <- BreastSubtypeR:::.pam50_confidence(out)
    expect_identical(actual[c(1, 3)], rep(NA_real_, 2))
    expect_true(is.finite(actual[2]))

    out <- confidence_test_input()
    out$testData <- out$testData[1:2, , drop = FALSE]
    expected <- vapply(seq_along(out$predictions), function(j) {
        confidence_direct(out$testData[, j],
            out$centroids[1:2, out$predictions[j]])
    }, numeric(1))
    expect_identical(expected, rep(0, 3))
    expect_identical(BreastSubtypeR:::.pam50_confidence(out), expected)
})

test_that("Confidence rejects absent and ambiguous gene identifiers", {
    for (component in c("testData", "centroids")) {
        for (bad_names in list(NULL, c("", paste0("gene", 2:5)),
            c(NA_character_, paste0("gene", 2:5)),
            c("gene2", paste0("gene", 2:5)))) {
            out <- confidence_test_input()
            rownames(out[[component]]) <- bad_names
            expect_identical(BreastSubtypeR:::.pam50_confidence(out),
                rep(NA_real_, 3))
        }
    }
    for (overlap in 0:1) {
        out <- confidence_test_input()
        rownames(out$centroids) <- paste0("other", seq_len(5))
        if (overlap == 1L) rownames(out$centroids)[1] <- "gene1"
        expect_identical(BreastSubtypeR:::.pam50_confidence(out),
            rep(NA_real_, 3))
    }
})

test_that("Confidence rejects missing and ambiguous predicted centroids", {
    out <- confidence_test_input()
    out$predictions[] <- c("unknown", NA_character_, "")
    expect_identical(BreastSubtypeR:::.pam50_confidence(out),
        rep(NA_real_, 3))
    out <- confidence_test_input()
    out$centroids <- cbind(out$centroids, LumA = out$centroids[, "LumA"])
    actual <- BreastSubtypeR:::.pam50_confidence(out)
    expect_identical(actual[c(1, 3)], rep(NA_real_, 2))
    expect_true(is.finite(actual[2]))
})

test_that("Confidence does not mutate the prediction or centroid inputs", {
    out <- confidence_test_input()
    original <- out
    BreastSubtypeR:::.pam50_confidence(out)
    expect_identical(out, original)
})

test_that("Bundled PAM50 results have original Confidence and unchanged ROR", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = environment())
    se <- OSLO2EMIT0obj$data_input$se_NC
    result <- suppressMessages(BS_parker(
        se, Subtype = TRUE, hasClinical = TRUE
    ))
    out <- result$outList
    expected <- vapply(seq_along(out$predictions), function(j) {
        confidence_direct(out$testData[, j],
            out$centroids[rownames(out$testData), out$predictions[j]])
    }, numeric(1))
    expect_identical(result$score.ROR$Confidence, expected)
    expect_identical(rownames(result$score.ROR), names(out$predictions))
    expect_true(any(expected > 0))

    matrix_out <- out
    matrix_out$centroids <- as.matrix(out$centroids)
    matrix_out$centroids <- matrix_out$centroids[
        rev(seq_len(nrow(matrix_out$centroids))), , drop = FALSE]
    matrix_scores <- suppressMessages(BreastSubtypeR:::RORgroup(
        matrix_out, as.data.frame(SummarizedExperiment::colData(se)),
        Subtype = TRUE, hasClinical = TRUE
    ))
    expect_identical(matrix_scores$Confidence, expected)
    other_columns <- setdiff(colnames(result$score.ROR), "Confidence")
    expect_identical(matrix_scores[, other_columns, drop = FALSE],
        result$score.ROR[, other_columns, drop = FALSE])
    expect_identical(matrix_out$predictions, out$predictions)
    expect_identical(matrix_out$distances, out$distances)
    expect_identical(matrix_out$dist.RORSubtype, out$dist.RORSubtype)
})
