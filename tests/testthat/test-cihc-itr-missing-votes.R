# Synthetic controls from the independent NCP-03 investigation. With ratio = 1
# the two ER+ and two ER- patients all enter each iteration's calibration set.
ncp03_input <- function() {
    data("BreastSubtypeRobj", envir = environment(), package = "BreastSubtypeR")
    genes <- rownames(BreastSubtypeRobj$centroid)
    ids <- c("T-9", "P-2", "S-7", "A-3")
    x <- outer(seq_along(genes), seq_along(ids),
        function(g, s) 5 + 0.01 * g + sin(g * s / 7) + cos(g / 5 + s))
    dimnames(x) <- list(genes, ids)
    pheno <- S4Vectors::DataFrame(PatientID = factor(ids, levels = c("unused", rev(ids))),
        ER = c("ER+", "ER-", "ER-", "ER+"), TSIZE = c(1, 0, 1, 0),
        NODE = c(1, 0, 2, 4), row.names = ids)
    SummarizedExperiment::SummarizedExperiment(assays = list(expression = x), colData = pheno)
}

ncp03_result <- function(input, four, clinical, ratio = 1) {
    before <- unserialize(serialize(input, NULL))
    rng <- .Random.seed
    result <- expect_experimental_ror_pc(BS_cIHC.itr(input, iteration = 2L, ratio = ratio,
        Subtype = four, hasClinical = clinical), n = as.integer(clinical))
    expect_identical(.Random.seed, rng)
    expect_identical(input, before)
    result
}

ncp03_expect_calls <- function(result, input, calls, four) {
    ids <- colnames(input)
    expect_identical(result$BS.all$PatientID, ids)
    expect_identical(rownames(result$BS.all), ids)
    expect_identical(result$score.ROR$patientID, ids)
    expect_identical(rownames(result$score.ROR), ids)
    expect_identical(result$outList$predictions, setNames(calls, ids))
    expect_identical(result$BS.all$BS, calls)
    expect_identical(result$score.ROR$Call, calls)
    votes <- matrix(calls, nrow = length(ids), ncol = 2L,
        dimnames = list(ids, c("itr.1", "itr.2")))
    expect_identical(result$BS.itr.keep, votes)
    if (four) {
        expect_identical(result$BS.itr.keep.Subtype, votes)
        expect_identical(result$outList$predictions.FourSubtype, setNames(calls, ids))
        expect_identical(result$BS.all$BS.Subtype, calls)
        expect_identical(result$score.ROR$Call.Subtype, calls)
    } else {
        expect_null(result$BS.itr.keep.Subtype)
        expect_null(result$outList$predictions.FourSubtype)
        expect_null(result$BS.all$BS.Subtype)
        expect_null(result$score.ROR$Call.Subtype)
    }
    # Each iteration uses all patients here, so the oracle centers each raw
    # gene on its observed cohort median before matching by gene name.
    x <- SummarizedExperiment::assay(input)
    expected <- sweep(x, 1L, apply(x, 1L, median, na.rm = TRUE))
    expected <- expected[sort(rownames(x)), , drop = FALSE]
    expect_identical(result$outList$testData, expected)
    expect_identical(result$score.ROR$ER, unname(expected["ESR1", ]))
    expect_identical(result$score.ROR$HER2, unname(expected["ERBB2", ]))

    classes <- c("Basal", "Her2", "LumA", "LumB", "Normal")
    for (field in c("distances", "dist.RORSubtype")) {
        mat <- result$outList[[field]]
        labels <- if (field == "distances") classes else classes[1:4]
        expect_identical(dimnames(mat), list(ids, labels))
        expect_type(mat, "double")
        missing <- unname(as.vector(mat[is.na(calls), , drop = FALSE]))
        expect_identical(missing, rep(NA_real_, length(missing)))
    }
}

test_that("iterative cIHC preserves each sample and each evaluable vote at sparse coverage", {
    withr::local_seed(178)
    withr::local_options(list(BreastSubtypeR.verbose = FALSE))
    input <- ncp03_input()
    # Ordered per-sample votes observed in the unmodified iteration outputs.
    expected <- list("50" = c("Basal", "Basal", "LumB", "Basal"),
        "2" = c("Basal", "LumA", "LumA", "Her2"),
        "1" = c("Basal", NA_character_, "LumA", "Her2"),
        "0" = c("Basal", NA_character_, "LumA", "LumB"))
    for (k in c(50L, 2L, 1L, 0L)) {
        se <- input
        if (k < 50L) {
            SummarizedExperiment::assay(se)[, "P-2"] <- NA_real_
            idx <- seq_len(k)
            SummarizedExperiment::assay(se)[c("ESR1", "ERBB2")[idx], "P-2"] <- c(12, 1)[idx]
        }
        for (four in c(FALSE, TRUE)) for (clinical in c(FALSE, TRUE)) {
            result <- expect_silent(ncp03_result(se, four, clinical))
            ncp03_expect_calls(result, se, expected[[as.character(k)]], four)
            rr <- result$score.ROR
            pheno <- if (clinical) as.data.frame(SummarizedExperiment::colData(se))
            ror198_expect_reference(result, pheno)
            expect_identical(ncol(rr), 14L + as.integer(four) + 6L * as.integer(clinical))
            if (k < 2L) {
                fields <- c("Confidence", grep("^ROR-", names(rr), value = TRUE))
                expect_true(all(is.na(rr["P-2", fields])))
            }
        }
    }
})

test_that("iterative cIHC returns typed no-calls for all patients and retains the nodal rule", {
    withr::local_seed(178)
    withr::local_options(list(BreastSubtypeR.verbose = FALSE))
    input <- ncp03_input()
    SummarizedExperiment::assay(input)[] <- NA_real_
    missing <- rep(NA_character_, ncol(input))
    for (four in c(FALSE, TRUE)) for (clinical in c(FALSE, TRUE)) {
        result <- expect_silent(ncp03_result(input, four, clinical))
        ncp03_expect_calls(result, input, missing, four)
        rr <- result$score.ROR
        risk <- grep("^ROR-", names(rr), value = TRUE)
        groups <- grep(" Group ", risk, value = TRUE)
        for (field in c("Confidence", setdiff(risk, groups)))
            expect_identical(rr[[field]], rep(NA_real_, ncol(input)))
        for (field in setdiff(groups, ror221_fields[2]))
            expect_identical(rr[[field]], missing)
        if (clinical) expect_identical(rr[[ror221_fields[2]]], c(NA, NA, NA, "high"))
    }
})

test_that("a fully measured patient at the calibration medians returns a missing vote", {
    withr::local_seed(178)
    withr::local_envvar(c(LANGUAGE = "C"))
    withr::local_options(list(BreastSubtypeR.verbose = FALSE))
    se <- ncp03_input()[, 1:3]
    SummarizedExperiment::colData(se)$ER <- c("ER+", "ER+", "ER-")
    g <- seq_len(nrow(se))
    offset <- 1 + sin(g / 4)
    SummarizedExperiment::assay(se)[] <- cbind(g + offset, g, g - offset)
    # Spearman reports zero variance for each affected centroid comparison.
    warnings <- character()
    result <- withCallingHandlers(ncp03_result(se, TRUE, FALSE, ratio = 2),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    expect_identical(warnings, rep("the standard deviation is zero", 26L))
    ncp03_expect_calls(result, se, c("Basal", NA_character_, "LumA"), TRUE)
})

test_that("the existing rank behavior for explicit infinities is preserved", {
    withr::local_seed(178)
    withr::local_options(list(BreastSubtypeR.verbose = FALSE))
    input <- ncp03_input()
    x <- outer(seq_len(nrow(input)), seq_len(ncol(input)),
        function(g, j) cos(g * j / 7) + sin(g + j))
    SummarizedExperiment::colData(input)$ER <- c("ER+", "ER-", "ER+", "ER-")
    for (k in 0:1) {
        x[, 2L] <- rep(c(Inf, -Inf), length.out = nrow(x))
        if (k == 1L) x[1L, 2L] <- 0
        SummarizedExperiment::assay(input)[] <- x
        rr <- expect_silent(ncp03_result(input, TRUE, FALSE))$score.ROR
        expect_identical(rr$Call[2], "LumA")
        expect_identical(rr$Call.Subtype[2], "LumA")
        expect_identical(rr$Confidence[2], c(0.77, 0.76)[k + 1L])
        expect_equal(rr[["ROR-S (Subtype Only)"]][2],
            c(29.9227845478982, 30.4728533978102)[k + 1L], tolerance = 1e-12)
    }
})

test_that("iterative vote ties use the established table ordering", {
    votes <- c("LumA", NA_character_, "Basal", "LumA", "Basal")
    expect_identical(BreastSubtypeR:::get_consensus_subtype(votes), "Basal")
    expect_identical(BreastSubtypeR:::get_consensus_subtype(rev(votes)), "Basal")
    expect_identical(BreastSubtypeR:::get_consensus_subtype(c(NA_character_, NA_character_)), NA_character_)
})
