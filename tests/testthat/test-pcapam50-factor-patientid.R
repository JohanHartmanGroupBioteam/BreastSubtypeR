pcapam50_id_membership <- function(actual, expected) {
    complete <- is.character(actual) && length(actual) == length(expected) &&
        !anyNA(actual) && !anyDuplicated(actual) && all(nzchar(actual)) &&
        !anyNA(expected) && !anyDuplicated(expected) && setequal(actual, expected)
    expect_true(complete)
    if (!complete) stop("Incomplete sample identities in regression comparison")
}

pcapam50_id_outputs <- function(result, ids, subtype) {
    axes <- list(result$BS.all$PatientID, as.character(result$score.ROR$patientID),
        names(result$outList$predictions), rownames(result$outList$distances),
        rownames(result$outList$dist.RORSubtype), colnames(result$outList$testData))
    if (subtype) axes <- c(axes, list(names(result$outList$predictions.FourSubtype)))
    for (axis in axes) {
        pcapam50_id_membership(axis, ids)
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

test_that("factor PatientID preserves complete character-ID results", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR")
    original <- OSLO2EMIT0obj$data_input$se_NC
    ids <- colnames(original)
    for (subtype in c(FALSE, TRUE)) for (clinical in c(FALSE, TRUE)) {
        reference <- expect_experimental_ror_pc(BS_PCAPAM50(original,
            Subtype = subtype, hasClinical = clinical), n = as.integer(clinical))
        pcapam50_id_outputs(reference, ids, subtype)
        encodings <- list(factor(ids), factor(ids, levels = rev(ids)),
            factor(ids, levels = c("unused-before", rev(ids), "unused-after")))
        for (encoding in encodings) {
            input <- original
            SummarizedExperiment::colData(input)$PatientID <- encoding
            saved <- input
            result <- expect_experimental_ror_pc(BS_PCAPAM50(input,
                Subtype = subtype, hasClinical = clinical), n = as.integer(clinical))
            expect_identical(input, saved)
            pcapam50_id_outputs(result, ids, subtype)
            expect_identical(is.na(result$score.ROR), is.na(reference$score.ROR))
            expect_identical(result$BS.all, reference$BS.all)
            expect_identical(result$mdns.fl, reference$mdns.fl)
            expect_identical(result$score.ROR, reference$score.ROR)
            expect_identical(result$outList, reference$outList)
        }
    }
})

test_that("factor IDs retain names through reordering, subsets and unused levels", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR")
    original <- OSLO2EMIT0obj$data_input$se_NC
    names <- paste0(rep(c("a-", "A_", "z.", "2"), length.out = ncol(original)), seq_len(ncol(original)))
    colnames(original) <- names
    SummarizedExperiment::colData(original)$PatientID <- names
    er <- SummarizedExperiment::colData(original)$ER
    keep <- c(which(er == "ER+")[seq_len(40)], which(er == "ER-")[seq_len(15)])
    for (indices in list(rev(seq_len(ncol(original))), rev(keep))) {
        input <- original
        SummarizedExperiment::colData(input)$PatientID <- factor(names,
            levels = c("unused-before", rev(names), "unused-after"))
        input <- input[, indices, drop = FALSE]
        reference_input <- original[, indices, drop = FALSE]
        reference <- expect_experimental_ror_pc(BS_PCAPAM50(reference_input, Subtype = TRUE, hasClinical = TRUE))
        result <- expect_experimental_ror_pc(BS_PCAPAM50(input, Subtype = TRUE, hasClinical = TRUE))
        pcapam50_id_outputs(result, names[indices], TRUE)
        expect_identical(result, reference)
    }
    numeric_ids <- seq_len(ncol(original)) + 1000L
    colnames(original) <- as.character(numeric_ids)
    SummarizedExperiment::colData(original)$PatientID <- as.character(numeric_ids)
    reference <- BS_PCAPAM50(original)
    SummarizedExperiment::colData(original)$PatientID <- numeric_ids
    expect_identical(BS_PCAPAM50(original), reference)
})

test_that("native PCA-PAM50 centering selects expression values by biological ID", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR")
    original <- OSLO2EMIT0obj$data_input$se_NC
    ids <- colnames(original)
    assay <- as.matrix(SummarizedExperiment::assay(original))
    reference <- expect_experimental_ror_pc(BS_PCAPAM50(original, Subtype = TRUE, hasClinical = TRUE))
    captures <- new.env(parent = emptyenv())
    instrument <- function(name, stage) {
        fn <- getFromNamespace(name, "BreastSubtypeR")
        record <- function(e) {
            first <- stage == "first"
            captures[[stage]] <- list(input = e$.input_matrix,
                clinical = if (first) e$df.cln else e$df.pam,
                intended = if (first) c(as.character(e$ERP.pc1ihc$PatientID[e$i]),
                    as.character(e$ERN.pc1ihc$PatientID)) else
                    c(as.character(e$ERP.pam_sampled$PatientID), as.character(e$ERN.pam_sampled$PatientID)),
                selected = if (first) e$mbal.pc1ihc else e$mbal.pam)
        }
        environment(fn) <- list2env(list(.record = record), parent = environment(fn))
        body(fn) <- bquote({
            .input_matrix <- mat
            on.exit(.record(environment()), add = TRUE)
            .(body(fn))
        })
        fn
    }
    local_mocked_bindings(
        makeCalls.PC1ihc = instrument("makeCalls.PC1ihc", "first"),
        makeCalls.v1PAM = instrument("makeCalls.v1PAM", "second"),
        .package = "BreastSubtypeR")
    reference_capture <- NULL
    for (factor_ids in c(FALSE, TRUE)) {
        input <- original
        if (factor_ids) SummarizedExperiment::colData(input)$PatientID <- factor(ids, levels = rev(ids))
        result <- expect_experimental_ror_pc(BS_PCAPAM50(input, Subtype = TRUE, hasClinical = TRUE))
        expect_identical(result, reference)
        expect_true(all(c("first", "second") %in% ls(captures)))
        for (stage in c("first", "second")) {
            z <- captures[[stage]]
            pcapam50_id_membership(colnames(z$input), ids)
            pcapam50_id_membership(as.character(z$clinical$PatientID), ids)
            expect_false(anyNA(z$intended))
            expect_identical(anyDuplicated(z$intended), 0L)
            expect_true(all(z$intended %in% ids))
            expect_identical(rownames(z$selected), rownames(assay))
            expect_identical(colnames(z$selected), z$intended)
            expect_identical(z$selected, assay[, z$intended, drop = FALSE])
            matched <- match(ids, as.character(z$clinical$PatientID))
            expect_false(anyNA(matched))
            for (column in c("TSIZE", "NODE")) {
                expect_identical(z$clinical[[column]][matched], SummarizedExperiment::colData(original)[[column]])
            }
        }
        er <- SummarizedExperiment::colData(original)$ER
        er_status <- ifelse(er == "ER+", "pos", "neg")
        expect_identical(colnames(captures$first$input), ids[order(er_status, decreasing = TRUE)])
        expect_identical(colnames(captures$first$input), as.character(captures$first$clinical$PatientID))
        expect_identical(captures$second$clinical$PatientID, sort(ids))
        if (factor_ids) expect_identical(as.list(captures), reference_capture)
        else reference_capture <- as.list(captures)
    }
})

test_that("PCA-PAM50 rejects invalid identifiers before scientific calculation", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR")
    original <- OSLO2EMIT0obj$data_input$se_NC
    ids <- colnames(original)
    reached <- FALSE
    local_mocked_bindings(makeCalls.PC1ihc = function(...) {
        reached <<- TRUE
        stop("unexpected scientific calculation")
    }, .package = "BreastSubtypeR")
    cases <- list(missing = replace(ids, 5, NA_character_), empty = replace(ids, 5, ""),
        duplicated = replace(ids, 5, ids[4]), unmatched = replace(ids, 5, "NOT_A_SAMPLE"))
    errors <- c(missing = "must not contain missing or empty values", empty = "must not contain missing or empty values",
        duplicated = "must be unique", unmatched = "not found among assay column names")
    for (name in names(cases)) for (factor_ids in c(FALSE, TRUE)) {
        input <- original
        SummarizedExperiment::colData(input)$PatientID <- if (factor_ids) factor(cases[[name]]) else cases[[name]]
        expect_error(BS_PCAPAM50(input), errors[[name]], fixed = TRUE)
        expect_false(reached)
    }
})

test_that("factor-ID correction preserves caller RNG and locale", {
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
            if (present) set.seed(197) else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
            before <- if (present) get(".Random.seed", .GlobalEnv) else NULL
            before_kind <- RNGkind()
            result <- expect_experimental_ror_pc(BS_PCAPAM50(input, Subtype = TRUE, hasClinical = TRUE))
            expect_identical(exists(".Random.seed", .GlobalEnv, inherits = FALSE), present)
            if (present) expect_identical(get(".Random.seed", .GlobalEnv), before)
            expect_identical(RNGkind(), before_kind)
            expect_identical(Sys.getlocale(), locale)
            pcapam50_id_outputs(result, ids, TRUE)
        }
    }
})
