# Frozen outputs were captured from unmodified public commit
# 6fbb1fb0c47c899d9159141136c481ac283ba644 under R 4.6.0 / Bioc 3.24.
# The fixture contains complete AIMS and SSPBC results for deterministic
# OSLO2-EMIT0 subsets (2, 12 and 102 samples), plus SSPBC for one sample.
# Tests only read this fixture; they must not regenerate expected results.
aims_regression_fixture <- function() {
    readRDS(test_path("fixtures", "aims_1.5.6_multisample.rds"))
}

aims_regression_input <- function() {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = environment())
    suppressMessages(Mapping(OSLO2EMIT0obj$se_obj, RawCounts = FALSE))$se_SSP
}

test_that("AIMS accepts one mapped sample with complete result dimensions", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR")
    input <- suppressMessages(Mapping(
        OSLO2EMIT0obj$se_obj[, 1], RawCounts = FALSE
    ))$se_SSP
    result <- suppressMessages(BS_AIMS(input))
    frozen <- aims_regression_fixture()$aims$two

    expect_identical(names(result), names(frozen))
    expect_type(result$cl, "character")
    expect_identical(dim(result$cl), c(1L, 1L))
    expect_identical(rownames(result$cl), colnames(input))
    expect_identical(colnames(result$cl), "20")
    expect_identical(result$cl, frozen$cl[1, , drop = FALSE])
    expect_identical(dim(result$prob), c(1L, 1L))
    expect_identical(dimnames(result$prob), dimnames(result$cl))
    expect_equal(result$prob, frozen$prob[1, , drop = FALSE], tolerance = 1e-12)
    expect_true(all(is.finite(result$prob)))
    expect_true(all(result$prob >= 0 & result$prob <= 1))
    expect_type(result$all.probs, "list")
    expect_identical(names(result$all.probs), "20")
    expect_identical(dim(result$all.probs[["20"]]), c(1L, 5L))
    expect_equal(
        result$all.probs[["20"]],
        frozen$all.probs[["20"]][1, , drop = FALSE], tolerance = 1e-12
    )
    expect_equal(unname(rowSums(result$all.probs[["20"]])), 1)
    expect_equal(unname(result$prob[1, 1]), max(result$all.probs[["20"]]))
    expect_identical(dim(result$rules.matrix), c(100L, 1L))
    expect_identical(result$rules.matrix, frozen$rules.matrix[, 1, drop = FALSE])
    expect_true(all(result$rules.matrix %in% c(0, 1)))
    expect_identical(dim(result$data.used), c(151L, 1L))
    expect_equal(result$data.used, frozen$data.used[, 1, drop = FALSE])
    expect_identical(colnames(result$data.used), colnames(input))
    expect_identical(length(result$EntrezID.used), nrow(result$data.used))
    expect_identical(result$EntrezID.used, frozen$EntrezID.used)
})

test_that("AIMS preserves an arbitrary single-sample identifier exactly", {
    input <- aims_regression_input()[, 1, drop = FALSE]
    sample_id <- "patient 01/A-B.001"
    colnames(input) <- sample_id
    result <- suppressMessages(BS_AIMS(input))
    expect_identical(rownames(result$cl), sample_id)
    expect_identical(rownames(result$prob), sample_id)
    expect_identical(colnames(result$data.used), sample_id)
})

for (cohort in c("two", "multiple", "full")) {
    test_that(paste("AIMS preserves frozen", cohort, "sample predictions"), {
        fixture <- aims_regression_fixture()
        input <- aims_regression_input()[, fixture$indices[[cohort]], drop = FALSE]
        result <- suppressMessages(BS_AIMS(input))
        frozen <- fixture$aims[[cohort]]
        expect_identical(result$cl, frozen$cl)
        expect_identical(rownames(result$cl), colnames(input))
        expect_identical(result$rules.matrix, frozen$rules.matrix)
        expect_identical(result$EntrezID.used, frozen$EntrezID.used)
        # Tolerate only floating-point variation across supported platforms.
        expect_equal(result, frozen, tolerance = 1e-12)
    })
}

test_that("AIMS retains missing and invalid input checks", {
    expect_error(BS_AIMS(), "se_obj")
    expect_error(BS_AIMS(NULL), "Input must be a SummarizedExperiment object")
    expect_error(BS_AIMS(matrix(1, 2, 2)), "SummarizedExperiment")
    input <- aims_regression_input()[, 1:2, drop = FALSE]
    negative <- input
    SummarizedExperiment::assay(negative) <- -abs(SummarizedExperiment::assay(input))
    expect_error(BS_AIMS(negative), "AIMS needs to be applied")
    nonnumeric <- input
    SummarizedExperiment::assay(nonnumeric) <- matrix(
        "bad", nrow(input), ncol(input), dimnames = dimnames(input)
    )
    expect_error(BS_AIMS(nonnumeric), "Verify D is a numeric matrix")
})

test_that("SSPBC predictions are unchanged for single and multiple samples", {
    fixture <- aims_regression_fixture()
    input <- aims_regression_input()
    indices <- c(fixture$indices, list(single = 1L))
    for (cohort in names(indices)) {
        result <- suppressMessages(BS_sspbc(input[, indices[[cohort]], drop = FALSE]))
        expect_identical(result$cl, fixture$sspbc[[cohort]]$cl, info = cohort)
        expect_equal(result, fixture$sspbc[[cohort]], tolerance = 1e-12, info = cohort)
    }
})
