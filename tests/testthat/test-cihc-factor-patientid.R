cihc_id_membership <- function(actual, expected) {
    complete <- is.character(actual) && length(actual) == length(expected) &&
        !anyNA(actual) && !anyDuplicated(actual) && all(nzchar(actual)) &&
        !anyNA(expected) && !anyDuplicated(expected) && setequal(actual, expected)
    expect_true(complete)
    if (!complete) stop("Incomplete sample identities in cIHC regression comparison")
}

cihc_id_outputs <- function(result, ids, subtype) {
    axes <- list(result$BS.all$PatientID, as.character(result$score.ROR$patientID),
        names(result$outList$predictions), rownames(result$outList$distances),
        rownames(result$outList$dist.RORSubtype), colnames(result$outList$testData))
    if (subtype) axes <- c(axes, list(names(result$outList$predictions.FourSubtype)))
    for (axis in axes) {
        cihc_id_membership(axis, ids)
        expect_identical(axis, ids)
    }
    expect_false(anyNA(result$BS.all$BS))
    expect_identical(unname(result$outList$predictions), result$BS.all$BS)
    if (subtype) {
        expected <- result$outList$predictions
        normal <- expected == "Normal"
        scores <- result$outList$distances[, c("Basal", "Her2", "LumA", "LumB"), drop = FALSE]
        expected[normal] <- colnames(scores)[max.col(scores[normal, , drop = FALSE], ties.method = "first")]
        expect_identical(result$outList$predictions.FourSubtype, expected)
        expect_identical(unname(expected), result$BS.all$BS.Subtype)
    }
}

test_that("native cIHC centering selects expression values by biological ID", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR")
    original <- OSLO2EMIT0obj$data_input$se_NC
    ids <- colnames(original)
    assay <- SummarizedExperiment::assay(original)
    reference <- BS_cIHC(original, Subtype = TRUE, hasClinical = TRUE)
    captures <- new.env(parent = emptyenv())
    native <- getFromNamespace("makeCalls_ihc", "BreastSubtypeR")
    record <- function(e) {
        captures$native <- list(input = e$.input_expression, clinical = e$df.cln,
            larger = as.character(e$ERP.ihc$PatientID),
            smaller = as.character(e$ERN.ihc$PatientID), draw = e$i,
            intended = c(as.character(e$ERP.ihc$PatientID[e$i]), as.character(e$ERN.ihc$PatientID)),
            selected = e$mbal.ihc, medians = e$medians.all)
    }
    environment(native) <- list2env(list(.record = record), parent = environment(native))
    body(native) <- bquote({
        .input_expression <- mat
        on.exit(.record(environment()), add = TRUE)
        .(body(native))
    })
    local_mocked_bindings(makeCalls_ihc = native, .package = "BreastSubtypeR")
    reference_capture <- NULL
    for (factor_ids in c(FALSE, TRUE)) {
        input <- original
        if (factor_ids) SummarizedExperiment::colData(input)$PatientID <- factor(ids, levels = rev(ids))
        result <- BS_cIHC(input, Subtype = TRUE, hasClinical = TRUE)
        cihc_id_outputs(result, ids, TRUE)
        expect_identical(result, reference)
        z <- captures$native
        cihc_id_membership(colnames(z$input), ids)
        cihc_id_membership(as.character(z$clinical$PatientID), ids)
        expect_identical(z$input, assay)
        expect_false(anyNA(z$intended))
        expect_identical(anyDuplicated(z$intended), 0L)
        expect_true(all(z$intended %in% ids))
        er <- SummarizedExperiment::colData(original)$ER
        expect_length(z$intended, 2L * min(table(er)))
        expect_identical(z$larger, ids[er == "ER+"])
        expect_identical(z$smaller, ids[er == "ER-"])
        expect_identical(colnames(z$selected), z$intended)
        expect_identical(z$selected, assay[, z$intended, drop = FALSE])
        matched <- match(ids, as.character(z$clinical$PatientID))
        expect_false(anyNA(matched))
        for (column in c("ER", "TSIZE", "NODE")) {
            expect_identical(z$clinical[[column]][matched], SummarizedExperiment::colData(original)[[column]])
        }
        if (factor_ids) {
            expect_identical(z$intended, reference_capture$intended)
            expect_identical(z$draw, reference_capture$draw)
            expect_identical(z$medians, reference_capture$medians)
        } else reference_capture <- z
    }
})

test_that("cIHC factor IDs preserve complete character results in both column orders", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR")
    original <- OSLO2EMIT0obj$data_input$se_NC
    original_ids <- colnames(original)
    for (indices in list(seq_len(ncol(original)), rev(seq_len(ncol(original))))) {
        reference_input <- original[, indices, drop = FALSE]
        ids <- colnames(reference_input)
        for (subtype in c(FALSE, TRUE)) for (clinical in c(FALSE, TRUE)) {
            reference <- BS_cIHC(reference_input, Subtype = subtype, hasClinical = clinical)
            cihc_id_outputs(reference, ids, subtype)
            encodings <- list(factor(ids), factor(ids, levels = rev(ids)),
                factor(ids, levels = c("unused-before", rev(ids), "unused-after")),
                factor(original_ids, levels = original_ids)[indices])
            for (encoding in encodings) {
                input <- reference_input
                SummarizedExperiment::colData(input)$PatientID <- encoding
                saved <- input
                result <- BS_cIHC(input, Subtype = subtype, hasClinical = clinical)
                expect_identical(input, saved)
                cihc_id_outputs(result, ids, subtype)
                expect_identical(is.na(result$score.ROR), is.na(reference$score.ROR))
                expect_identical(result$BS.all, reference$BS.all)
                expect_identical(result$mdns, reference$mdns)
                expect_identical(result$score.ROR, reference$score.ROR)
                expect_identical(result$outList, reference$outList)
            }
        }
    }
})

test_that("cIHC names remain aligned with matrix assays and smaller ER-positive groups", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR")
    original <- OSLO2EMIT0obj$data_input$se_NC
    ids <- paste0(rep(c("a-", "A_", "z.", "2"), length.out = ncol(original)), seq_len(ncol(original)))
    colnames(original) <- ids
    SummarizedExperiment::colData(original)$PatientID <- ids
    SummarizedExperiment::assay(original) <- as.matrix(SummarizedExperiment::assay(original))
    er <- SummarizedExperiment::colData(original)$ER
    keep <- rev(c(which(er == "ER+")[seq_len(10)], which(er == "ER-")))
    original <- original[, keep, drop = FALSE]
    reference <- BS_cIHC(original, Subtype = TRUE, hasClinical = TRUE)
    factor_input <- original
    SummarizedExperiment::colData(factor_input)$PatientID <- factor(ids, levels = c("unused", rev(ids)))[keep]
    result <- BS_cIHC(factor_input, Subtype = TRUE, hasClinical = TRUE)
    cihc_id_outputs(result, ids[keep], TRUE)
    expect_identical(result, reference)
    numeric_ids <- seq_len(ncol(original)) + 1000L
    colnames(original) <- as.character(numeric_ids)
    SummarizedExperiment::colData(original)$PatientID <- as.character(numeric_ids)
    reference <- BS_cIHC(original)
    SummarizedExperiment::colData(original)$PatientID <- numeric_ids
    expect_identical(BS_cIHC(original), reference)
})

test_that("cIHC rejects invalid identifiers before scientific calculation", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR")
    original <- OSLO2EMIT0obj$data_input$se_NC
    ids <- colnames(original)
    reached <- FALSE
    local_mocked_bindings(makeCalls_ihc = function(...) {
        reached <<- TRUE
        stop("unexpected scientific calculation")
    }, .package = "BreastSubtypeR")
    cases <- list(missing = replace(ids, 5, NA_character_), empty = replace(ids, 5, ""),
        duplicated = replace(ids, 5, ids[4]), unmatched = replace(ids, 5, "NOT_A_SAMPLE"))
    errors <- c(missing = "must not contain missing or empty values",
        empty = "must not contain missing or empty values", duplicated = "must be unique",
        unmatched = "not found among assay column names")
    for (name in names(cases)) for (factor_ids in c(FALSE, TRUE)) {
        input <- original
        SummarizedExperiment::colData(input)$PatientID <- if (factor_ids) factor(cases[[name]]) else cases[[name]]
        expect_error(BS_cIHC(input), errors[[name]], fixed = TRUE)
        expect_false(reached)
    }
})

test_that("cIHC factor IDs preserve caller RNG and locale", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR")
    input <- OSLO2EMIT0obj$data_input$se_NC
    ids <- colnames(input)
    SummarizedExperiment::colData(input)$PatientID <- factor(ids, levels = rev(ids))
    withr::local_preserve_seed()
    original_kind <- RNGkind()
    withr::defer(do.call(RNGkind, as.list(original_kind)))
    locale <- Sys.getlocale()
    for (kind in c("Mersenne-Twister", "L'Ecuyer-CMRG")) {
        RNGkind(kind, "Inversion", "Rejection")
        for (present in c(FALSE, TRUE)) {
            if (present) set.seed(200) else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
            before <- if (present) get(".Random.seed", .GlobalEnv) else NULL
            before_kind <- RNGkind()
            result <- BS_cIHC(input, Subtype = TRUE, hasClinical = TRUE)
            expect_identical(exists(".Random.seed", .GlobalEnv, inherits = FALSE), present)
            if (present) expect_identical(get(".Random.seed", .GlobalEnv), before)
            expect_identical(RNGkind(), before_kind)
            expect_identical(Sys.getlocale(), locale)
            cihc_id_outputs(result, ids, TRUE)
        }
    }
})
