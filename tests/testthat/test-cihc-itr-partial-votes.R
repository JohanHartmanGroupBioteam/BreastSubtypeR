# Reproduce the archived 3 ER+ / 1 ER- NCP-03 cohort. The reference replays
# its original draws, using the native per-iteration classifier before either
# aggregation helper. Patient T remains in every calibration pair.
check_partial_vote_averages <- function(measured) {
    withr::local_seed(7, .rng_kind = "Mersenne-Twister", .rng_sample_kind = "Rejection")
    withr::local_envvar(c(LANGUAGE = "C"))
    data("BreastSubtypeRobj", package = "BreastSubtypeR", envir = environment())
    centers <- BreastSubtypeRobj$centroid
    ids <- c("T", "Z", "B", "K")
    x <- outer(seq_len(nrow(centers)), seq_along(ids),
        function(i, j) sin(i * j / 4) + cos(i / (j + 1)))
    dimnames(x) <- list(rownames(centers), ids)
    x[c("ESR1", "ERBB2"), ] <- rbind(c(10, 5, 13, 3), c(20, 15, 12, 23))
    x[!rownames(x) %in% measured, "T"] <- NA_real_
    clinical <- S4Vectors::DataFrame(PatientID = ids, ER = c("ER-", "ER+", "ER+", "ER+"),
        row.names = ids)
    input <- SummarizedExperiment::SummarizedExperiment(assays = list(expression = x),
        colData = clinical)
    draws <- c("K", "K", "Z", "Z", "B", "K", "Z", "B", "K", "B", "Z", "K")
    # Only three distinct sets of medians are possible with ratio 1.
    by_peer <- lapply(setNames(ids[-1], ids[-1]), function(peer) {
        medians <- apply(x[, c(peer, "T")], 1, stats::median, na.rm = TRUE)
        calibrated <- sweep(x, 1, medians)
        suppressWarnings(BreastSubtypeR:::sspPredict(centers, calibrated, Subtype = TRUE))
    })
    iterations <- setNames(by_peer[draws], paste0("itr.", seq_along(draws)))
    calls <- vapply(iterations, "[[", setNames(character(length(ids)), ids), "predictions")
    calls4 <- vapply(iterations, "[[", setNames(character(length(ids)), ids), "predictions.FourSubtype")
    consensus <- if (length(measured) == 2L) {
        c(T = "LumA", Z = "Her2", B = "LumB", K = "Her2")
    } else {
        c(T = NA_character_, Z = "Her2", B = "LumB", K = "Basal")
    }

    warnings <- character()
    before <- .Random.seed
    actual <- withCallingHandlers(BS_cIHC.itr(input, iteration = 12L, ratio = 1,
        seed = 118L, Subtype = TRUE), warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    expect_identical(.Random.seed, before)
    if (length(measured) == 2L) {
        expect_identical(unique(warnings), "the standard deviation is zero")
        expect_identical(unname(calls["T", ]),
            c("LumA", "LumA", NA, NA, "Basal", "LumA", NA, "Basal", "LumA", "Basal", NA, "LumA"))
    } else {
        expect_length(warnings, 0L)
        expect_true(all(is.na(calls["T", ])))
    }
    expect_identical(actual$BS.itr.keep, calls)
    expect_identical(actual$BS.itr.keep.Subtype, calls4)
    expect_identical(actual$outList$predictions, consensus)
    average_expression <- apply(simplify2array(lapply(iterations, "[[", "testData")), 1:2, mean)
    expect_equal(actual$outList$testData, average_expression, tolerance = 1e-12)

    for (field in c("distances", "dist.RORSubtype")) {
        expected <- iterations[[1]][[field]]
        for (id in ids) {
            # Both similarity matrices use the five-class vote, including in
            # four-class mode. Base-R which() drops the missing comparisons.
            support <- which(calls[id, ] == consensus[id])
            if (length(support)) {
                rows <- lapply(iterations[support], function(it) it[[field]][id, ])
                expected[id, ] <- -colMeans(do.call(rbind, rows), na.rm = TRUE)
            } else {
                expected[id, ] <- NA_real_
                expect_identical(actual$outList[[field]][id, ], expected[id, ])
            }
        }
        expect_equal(actual$outList[[field]], expected, tolerance = 1e-12)
    }
}

test_that("a partially evaluable patient retains only consensus-supporting iterations", {
    check_partial_vote_averages(c("ESR1", "ERBB2"))
})

test_that("an uncalled patient remains in the peer calibration cohort", {
    check_partial_vote_averages("ESR1")
})

test_that("undefined similarity cells retain NaN in an otherwise evaluable patient", {
    ids <- c("observed", "uncalled")
    classes <- c("Basal", "Her2", "LumA", "LumB", "Normal")
    d <- matrix(c(0.2, NA, NA, NA, 0.1, NA, 0.3, NA, 0.4, NA), nrow = 2,
        dimnames = list(ids, classes))
    calls <- c(observed = "Basal", uncalled = NA_character_)
    it <- list(predictions = calls, distances = d, dist.RORSubtype = d[, 1:4],
        testData = matrix(c(1, NA_real_), nrow = 1, dimnames = list("ESR1", ids)))
    averages <- BreastSubtypeR:::get_average_subtype(list(it), calls)
    for (field in c("mean_distance", "mean_distance.Subtype")) {
        values <- as.matrix(averages[[field]])
        expect_true(is.nan(values["observed", "Her2"]))
        expect_identical(values["observed", "Basal"], 0.2)
        expect_identical(unname(values["uncalled", ]), rep(NA_real_, ncol(values)))
    }
})
