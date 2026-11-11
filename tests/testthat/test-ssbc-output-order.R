# Public expression data, with explicitly synthetic TN/receptor annotations
# for coverage of the four supported strategies and subgroup edge cases.
# Expectations concern output assembly, not equality between differently
# calibrated cohorts. Historical scientific fixtures remain independent.
ssbc234_inputs <- function() {
    e <- new.env(parent = emptyenv())
    utils::data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = e)
    x <- e$OSLO2EMIT0obj$data_input$se_NC
    SummarizedExperiment::colData(x)$TN <- rep(c("TN", "nonTN", "nonTN"), length.out = ncol(x))
    annotate <- function(i, er, her2, tn) {
        z <- x[, i, drop = FALSE]
        SummarizedExperiment::colData(z)$ER <- er
        SummarizedExperiment::colData(z)$HER2 <- her2
        SummarizedExperiment::colData(z)$TN <- tn
        z
    }
    list(original = x, reversed = x[, 102:1, drop = FALSE],
        permuted = x[, c(seq(2, 102, 2), seq(1, 101, 2)), drop = FALSE],
        singleton = annotate(1L, "ER+", "HER2-", "TN"),
        singleton_subgroups = annotate(1:6, c(rep("ER+", 5), "ER-"),
            c(rep("HER2-", 4), "HER2+", "HER2+"), c(rep("TN", 5), "nonTN")),
        empty_subgroups = annotate(1:6, rep("ER+", 6), rep("HER2-", 6), rep("TN", 6)),
        missing = annotate(1:6, c(NA, "ER+", "ER-", "ER+", "ER-", "ER+"),
            c("HER2-", "HER2-", NA, "HER2+", "HER2-", "HER2-"),
            c("TN", NA, "TN", "nonTN", "nonTN", "TN")),
        all_missing = annotate(1:4, rep(NA_character_, 4), rep(NA_character_, 4),
            rep(NA_character_, 4)))
}

ssbc234_singletons <- function(cohort, s) {
    # Independently enumerated singleton strata, confirmed on installed main
    # 930ade96 before the output-only correction. No warning is ignored merely
    # because a small cohort might warn: exact subgroup, ID and count must match.
    strata <- switch(cohort,
        singleton = switch(s, ER = c(ER_pos = 1L), ER.v2 = c(ERpos_HER2neg = 1L),
            TN = c(TN = 1L), TN.v2 = c(TNBC = 1L)),
        singleton_subgroups = switch(s, ER = c(ER_neg = 6L),
            ER.v2 = c(HER2pos_ERneg = 6L, HER2pos_ERpos = 5L),
            TN = c(nonTN = 6L), TN.v2 = integer()),
        missing = switch(s, ER.v2 = c(ERneg_HER2neg = 5L, HER2pos_ERpos = 4L), integer()),
        integer())
    if (!length(strata)) return(character())
    sprintf(paste0("ssBC/ssBC.v2: subgroup '%s' has n=1 ",
        "(sample OSLO2EMIT0.%03d); no subtype or risk estimate."), names(strata), strata)
}

ssbc234_with_warnings <- function(expr, expected, clinical) {
    seen <- character()
    result <- expect_experimental_ror_pc(withCallingHandlers(force(expr),
        warning = function(w) {
            msg <- conditionMessage(w)
            if (msg %in% expected) {
                expect_identical(class(w), c("simpleWarning", "warning", "condition"))
                expect_null(conditionCall(w))
                seen <<- c(seen, msg)
                invokeRestart("muffleWarning")
            }
            # Unexpected warnings remain visible to testthat.
        }), n = as.integer(clinical))
    expect_identical(seen, expected)
    result
}

test_that("ssBC preserves ordered sample axes and matrix dimensions for every strategy", {
    inputs <- ssbc234_inputs()
    e <- new.env(parent = emptyenv())
    utils::data("BreastSubtypeRobj", package = "BreastSubtypeR", envir = e)
    for (cohort in names(inputs)) for (factor_ids in c(FALSE, TRUE))
        for (s in c("ER", "ER.v2", "TN", "TN.v2"))
            for (subtype in c(FALSE, TRUE)) for (clinical in c(FALSE, TRUE)) {
                x <- inputs[[cohort]]
                ids <- colnames(x)
                SummarizedExperiment::colData(x)$PatientID <- if (factor_ids)
                    factor(ids, levels = c("unused", rev(ids))) else ids
                before <- x
                expected_warnings <- ssbc234_singletons(cohort, s)
                run <- function() BS_ssBC(x, s = s, Subtype = subtype, hasClinical = clinical)
                a <- ssbc234_with_warnings(run(), expected_warnings, clinical)
                if (cohort == "all_missing" && subtype && s != "TN.v2") {
                    # NCP-02: no sample has a selected receptor stratum. The
                    # four-class return is the five-class empty return plus
                    # only the four-class fields, in their canonical places.
                    five <- ssbc234_with_warnings(BS_ssBC(x, s = s, Subtype = FALSE,
                        hasClinical = clinical), expected_warnings, clinical)
                    expect_identical(a, ncp02_four_from_five(five))
                    expect_identical(a, ncp02_explicit_empty(ids, clinical))
                    expect_false("patientID" %in% names(a$score.ROR))
                    expect_identical(ncol(a$score.ROR), if (clinical) 20L else 14L)
                    expect_identical(x, before)
                }
                out <- a$outList
                expect_identical(names(a), c("BS.all", "score.ROR", "mdns", "outList"))
                expect_identical(names(out), c("predictions", if (subtype) "predictions.FourSubtype",
                    "testData", "distances", "dist.RORSubtype", "centroids"))
                for (field in c("BS.all", "score.ROR")) {
                    expect_s3_class(a[[field]], "data.frame")
                    expect_identical(rownames(a[[field]]), ids)
                }
                expect_identical(as.character(a$BS.all$PatientID), ids)
                expect_identical(names(out$predictions), ids)
                # The native all-missing table uses logical NA, while its
                # named vector uses character NA. Preserve both payload types.
                if (cohort == "all_missing") {
                    expect_identical(typeof(a$BS.all$BS), "logical")
                    expect_identical(unname(out$predictions), rep(NA_character_, length(ids)))
                } else expect_identical(unname(out$predictions), a$BS.all$BS)
                expect_identical(a$score.ROR$Call, a$BS.all$BS)
                expect_identical("predictions.FourSubtype" %in% names(out), subtype)
                expect_identical("BS.Subtype" %in% names(a$BS.all), subtype)
                if (subtype) {
                    expect_identical(names(out$predictions.FourSubtype), ids)
                    if (cohort == "all_missing") {
                        expect_identical(typeof(a$BS.all$BS.Subtype), "logical")
                        expect_identical(unname(out$predictions.FourSubtype), rep(NA_character_, length(ids)))
                    } else expect_identical(unname(out$predictions.FourSubtype), a$BS.all$BS.Subtype)
                    expect_identical(a$score.ROR$Call.Subtype, a$BS.all$BS.Subtype)
                }
                for (field in c("distances", "dist.RORSubtype")) {
                    k <- if (field == "distances") 5L else 4L
                    expect_true(is.matrix(out[[field]]))
                    expect_identical(typeof(out[[field]]), "double")
                    expect_identical(dim(out[[field]]), c(length(ids), k))
                    expect_identical(rownames(out[[field]]), ids)
                    expect_identical(colnames(out[[field]]), c("Basal", "Her2", "LumA", "LumB", "Normal")[seq_len(k)])
                }
                expect_true(is.matrix(out$testData))
                expect_identical(typeof(out$testData), "double")
                expect_identical(colnames(out$testData), ids)
                # All-missing ER/ER.v2/TN already return zero gene rows with
                # Subtype=FALSE. Preserve that distinct native empty result.
                genes <- if (cohort == "all_missing" && s != "TN.v2") NULL else rownames(out$centroids)
                expect_identical(rownames(out$testData), genes)
                expect_identical(dim(out$testData), c(length(genes), length(ids)))
                expect_identical(a$mdns, e$BreastSubtypeRobj$ssBC.subgroupQuantile)
                expect_identical(out$centroids, e$BreastSubtypeRobj$centroid)
                unavailable <- switch(s,
                    ER = is.na(SummarizedExperiment::colData(x)$ER),
                    ER.v2 = is.na(SummarizedExperiment::colData(x)$ER) | is.na(SummarizedExperiment::colData(x)$HER2),
                    TN = is.na(SummarizedExperiment::colData(x)$TN),
                    TN.v2 = is.na(SummarizedExperiment::colData(x)$TN) | SummarizedExperiment::colData(x)$TN != "TN")
                expect_true(all(is.na(a$BS.all$BS[unavailable])))
                expect_true(all(is.na(out$testData[, unavailable, drop = FALSE])))
                if ("patientID" %in% names(a$score.ROR)) {
                    expect_true(all(is.na(a$score.ROR$patientID[unavailable])))
                    expect_identical(as.character(a$score.ROR$patientID[!unavailable]), ids[!unavailable])
                }
                for (msg in expected_warnings) {
                    id <- sub(".*sample (OSLO2EMIT0[.][0-9]+).*", "\\1", msg)
                    expect_identical(out$predictions[[id]], NA_character_)
                    if (subtype) expect_identical(out$predictions.FourSubtype[[id]], NA_character_)
                    expect_true(all(is.na(out$distances[id, ])))
                    expect_true(all(is.na(out$dist.RORSubtype[id, ])))
                    expect_true(all(is.na(out$testData[, id])))
                    scores <- a$score.ROR[id, , drop = FALSE]
                    expect_identical(as.character(scores$patientID), id)
                    for (key in setdiff(names(scores), "patientID")) {
                        expected <- if (grepl("Group|^Call", key)) NA_character_ else NA_real_
                        expect_identical(scores[[key]], expected)
                    }
                }
                expect_identical(x, before)
            }
})
