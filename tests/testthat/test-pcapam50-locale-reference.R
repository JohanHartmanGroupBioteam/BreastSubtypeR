# Original-reference oracle and provenance: fixtures/pcapam50-locale/README.md.
# These tests preserve per-locale compatibility, not cross-locale invariance.
pcap232_fixture <- function() readRDS(test_path("fixtures", "pcapam50-locale", "reference.rds"))

pcap232_input <- function(fixture, idkind = "original", indices = seq_len(102L), factor_ids = FALSE) {
    env <- new.env(parent = emptyenv())
    data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = env)
    se <- env$OSLO2EMIT0obj$data_input$se_NC
    expect_identical(colnames(se), fixture$ids$original)
    expect_identical(rownames(se), fixture$genes)
    cd <- as.data.frame(SummarizedExperiment::colData(se))
    rownames(cd) <- NULL
    expect_identical(cd[, names(fixture$clinical)], fixture$clinical)
    ids <- if (idkind == "original") colnames(se) else
        ifelse(seq_len(ncol(se)) %% 2L == 0L, toupper(colnames(se)), tolower(colnames(se)))
    expect_identical(ids, fixture$ids[[idkind]])
    colnames(se) <- ids
    SummarizedExperiment::colData(se)$PatientID <- if (factor_ids)
        factor(ids, levels = c("unused-before", rev(ids), "unused-after")) else ids
    se[, indices, drop = FALSE]
}

pcap232_locale <- function(profile) {
    original <- Sys.getlocale("LC_COLLATE")
    on.exit(Sys.setlocale("LC_COLLATE", original))
    candidates <- if (profile == "C") "C" else
        c("en_US.UTF-8", "en_US.utf8", "English_United States.utf8")
    for (name in candidates) {
        effective <- suppressWarnings(Sys.setlocale("LC_COLLATE", name))
        if (!is.na(effective) && nzchar(effective)) {
            expect_identical(Sys.getlocale("LC_COLLATE"), effective)
            return(list(requested = name, effective = effective))
        }
    }
    if (profile == "C") stop("Required C collation is unavailable")
    skip("No alternate en_US UTF-8 collation is available; C tests run separately")
}

pcap232_sentinel <- function(fixture, profile, locale) {
    expect_identical(Sys.getlocale("LC_COLLATE"), locale$effective)
    # Full character-key sentinel, not the locale name alone: OS collation
    # implementations need not agree even when their locale names match.
    expected <- fixture$cases[[paste("mixed", profile, sep = "__")]]$first$ordered_ids
    actual <- sort(fixture$ids$mixed)
    expect_identical(actual, expected, info = "Effective sample-ID collation sentinel")
    expect_identical(sort(fixture$genes), fixture$genes)
    expect_identical(sort(c("neg", "pos"), decreasing = TRUE), fixture$er_order)
    if (!identical(actual, expected)) stop("Available locale has an unverified collation profile")
}

pcap232_capture <- function(.env = parent.frame()) {
    captures <- new.env(parent = emptyenv())
    instrument <- function(name, stage) {
        fn <- getFromNamespace(name, "BreastSubtypeR")
        record <- function(e) {
            first <- stage == "first"
            captures[[stage]] <- list(
                input = e$.pcap232_input,
                clinical = if (first) e$df.cln else e$df.pam,
                ordered_ids = if (first) e$df.pca1$PatientID else e$df.pam$PatientID,
                positive = if (first) e$ERP.pc1ihc$PatientID else e$ERP.pam$PatientID,
                negative = if (first) e$ERN.pc1ihc$PatientID else e$ERN.pam$PatientID,
                draw = e$i,
                selected_matrix = if (first) e$mbal.pc1ihc else e$mbal.pam,
                medians = e$mdns, calls = e$out$predictions,
                correlations = e$out$distances)
        }
        environment(fn) <- list2env(list(.pcap232_record = record), parent = environment(fn))
        # Observe the real function, with its original body/formals/environment
        # parent. No sampling, median, PCA or scoring expression is replaced.
        body(fn) <- bquote({
            .pcap232_input <- mat
            on.exit(.pcap232_record(environment()), add = TRUE)
            .(body(fn))
        })
        fn
    }
    local_mocked_bindings(
        makeCalls.PC1ihc = instrument("makeCalls.PC1ihc", "first"),
        makeCalls.v1PAM = instrument("makeCalls.v1PAM", "second"),
        .package = "BreastSubtypeR", .env = .env)
    captures
}

pcap232_ids <- function(actual, expected) {
    expect_type(actual, "character")
    expect_length(actual, length(expected))
    expect_false(anyNA(actual))
    expect_identical(anyDuplicated(actual), 0L)
    expect_setequal(actual, expected)
}

pcap232_numeric <- function(actual, expected) {
    expect_identical(dim(actual), dim(expected))
    expect_identical(dimnames(actual), dimnames(expected))
    expect_identical(names(actual), names(expected))
    expect_true(all(is.finite(actual)))
    expect_true(all(is.finite(expected)))
    # Absolute tolerance; does not change any decision, rank or sample ID.
    # Justification and reference decision margins are in the fixture README.
    expect_lte(max(abs(actual - expected)), 1e-10)
}

pcap232_stage <- function(actual, expected, se, first) {
    ids <- colnames(se)
    pcap232_ids(actual$ordered_ids, ids)
    expect_identical(actual$ordered_ids, expected$ordered_ids)
    expect_identical(actual$positive, expected$positive)
    expect_identical(actual$negative, expected$negative)
    expect_identical(actual$draw, expected$draw)
    expect_identical(colnames(actual$selected_matrix), expected$selected)
    expect_identical(actual$selected_matrix,
        as.matrix(SummarizedExperiment::assay(se))[, expected$selected, drop = FALSE])
    pcap232_numeric(actual$medians, expected$medians)
    pcap232_ids(names(actual$calls), ids)
    expect_identical(actual$calls[names(expected$calls)], expected$calls)
    pcap232_ids(rownames(actual$correlations), ids)
    pcap232_numeric(actual$correlations[rownames(expected$correlations), , drop = FALSE],
        expected$correlations)
    cd <- as.data.frame(SummarizedExperiment::colData(se))
    pcap232_ids(as.character(actual$clinical$PatientID), ids)
    index <- match(ids, as.character(actual$clinical$PatientID))
    for (column in if (first) c("ER", "PR", "HER2", "TSIZE", "NODE") else c("TSIZE", "NODE"))
        expect_identical(actual$clinical[[column]][index], cd[[column]])
    # ER preparation preserves within-class input order; this differs from the
    # character-key ordering used by both calibration sampling stages.
    prepared_ids <- c(ids[cd$ER == "ER+"], ids[cd$ER == "ER-"])
    expect_identical(colnames(actual$input), prepared_ids)
    if (first) expect_identical(actual$clinical$PatientID, prepared_ids)
}

pcap232_output <- function(result, se, expected) {
    ids <- colnames(se)
    axes <- list(result$BS.all$PatientID, as.character(result$score.ROR$patientID),
        names(result$outList$predictions), names(result$outList$predictions.FourSubtype),
        rownames(result$outList$distances), rownames(result$outList$dist.RORSubtype),
        colnames(result$outList$testData))
    for (axis in axes) {
        pcap232_ids(axis, ids)
        expect_identical(axis, ids)
    }
    expect_identical(result$outList$predictions, expected$calls[ids])
    expect_identical(result$BS.all$BS, unname(expected$calls[ids]))
    pcap232_numeric(setNames(result$mdns.fl[names(expected$medians), "v1PAM.mdns"],
        names(expected$medians)), expected$medians)
    pcap232_numeric(result$outList$distances,
        expected$correlations[ids, , drop = FALSE])
    # Clinical scoring is checked by the existing independent Parker formula,
    # not by claiming that the original reference supports package extensions.
    ror198_expect_reference(result, as.data.frame(SummarizedExperiment::colData(se)))
}

pcap232_align <- function(result, ids) {
    result$BS.all <- result$BS.all[match(ids, result$BS.all$PatientID), , drop = FALSE]
    result$score.ROR <- result$score.ROR[match(ids, result$score.ROR$patientID), , drop = FALSE]
    rownames(result$BS.all) <- rownames(result$score.ROR) <- NULL
    for (field in c("predictions", "predictions.FourSubtype"))
        result$outList[[field]] <- result$outList[[field]][ids]
    for (field in c("distances", "dist.RORSubtype"))
        result$outList[[field]] <- result$outList[[field]][ids, , drop = FALSE]
    result$outList$testData <- result$outList$testData[, ids, drop = FALSE]
    result
}

for (profile in c("C", "en_US.UTF-8")) {
    test_that(paste("both PCAPAM50 stages match the original reference in", profile), {
        fixture <- pcap232_fixture()
        before_locale <- Sys.getlocale("LC_COLLATE")
        locale <- pcap232_locale(profile)
        expect_identical(Sys.getlocale("LC_COLLATE"), before_locale)
        withr::local_preserve_seed()
        kind <- RNGkind()
        withr::defer(do.call(RNGkind, as.list(kind)))
        do.call(RNGkind, as.list(fixture$rng_kind))
        captures <- pcap232_capture()
        withr::with_locale(c(LC_COLLATE = locale$requested), {
            pcap232_sentinel(fixture, profile, locale)
            for (idkind in c("original", "mixed")) {
                expected <- fixture$cases[[paste(idkind, profile, sep = "__")]]
                baseline <- NULL
                orders <- list(original = seq_len(102L), reverse = 102L:1L,
                    permuted = fixture$permutation)
                for (indices in orders) for (factor_ids in c(FALSE, TRUE)) {
                    se <- pcap232_input(fixture, idkind, indices, factor_ids)
                    saved <- se
                    set.seed(97232)
                    before_seed <- .Random.seed
                    result <- expect_experimental_ror_pc(BS_PCAPAM50(se,
                        Subtype = TRUE, hasClinical = TRUE, seed = fixture$seed))
                    expect_identical(.Random.seed, before_seed)
                    expect_identical(RNGkind(), fixture$rng_kind)
                    expect_identical(Sys.getlocale("LC_COLLATE"), locale$effective)
                    expect_identical(se, saved)
                    expect_setequal(ls(captures), c("first", "second"))
                    pcap232_stage(captures$first, expected$first, se, TRUE)
                    pcap232_stage(captures$second, expected$second, se, FALSE)
                    pcap232_output(result, se, expected$second)
                    aligned <- pcap232_align(result, fixture$ids[[idkind]])
                    if (is.null(baseline)) baseline <- aligned else {
                        expect_identical(aligned$BS.all, baseline$BS.all)
                        expect_equal(aligned, baseline, tolerance = 1e-10)
                    }
                }
            }
        })
        expect_identical(Sys.getlocale("LC_COLLATE"), before_locale)
    })

    test_that(paste("PCAPAM50 repeats and preserves caller state in", profile), {
        fixture <- pcap232_fixture()
        before_locale <- Sys.getlocale("LC_COLLATE")
        locale <- pcap232_locale(profile)
        se <- pcap232_input(fixture, "mixed", fixture$permutation, TRUE)
        withr::local_preserve_seed()
        kind <- RNGkind()
        withr::defer(do.call(RNGkind, as.list(kind)))
        # Only Mersenne-Twister is pinned by the reference fixture. The second
        # kind checks the caller-state contract, not unverified reference calls.
        for (generator in c("Mersenne-Twister", "L'Ecuyer-CMRG")) {
            RNGkind(generator, "Inversion", "Rejection")
            expected_kind <- RNGkind()
            for (present in c(FALSE, TRUE)) {
                if (present) set.seed(232) else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
                    rm(".Random.seed", envir = .GlobalEnv)
                before_seed <- if (present) get(".Random.seed", .GlobalEnv) else NULL
                run <- function() withr::with_locale(c(LC_COLLATE = locale$requested), {
                    pcap232_sentinel(fixture, profile, locale)
                    expect_experimental_ror_pc(BS_PCAPAM50(se,
                        Subtype = TRUE, hasClinical = TRUE, seed = fixture$seed))
                })
                first <- run()
                expect_identical(run(), first)
                expect_identical(exists(".Random.seed", .GlobalEnv, inherits = FALSE), present)
                if (present) expect_identical(get(".Random.seed", .GlobalEnv), before_seed)
                expect_identical(RNGkind(), expected_kind)
                expect_identical(Sys.getlocale("LC_COLLATE"), before_locale)
                # A real error after both scientific stages (the documented
                # experimental warning promoted to an error) also restores state.
                expect_error(withr::with_locale(c(LC_COLLATE = locale$requested),
                    withr::with_options(list(warn = 2), BS_PCAPAM50(se,
                        Subtype = TRUE, hasClinical = TRUE, seed = fixture$seed))),
                    "Experimental 18-gene ROR-PC", fixed = TRUE)
                expect_identical(exists(".Random.seed", .GlobalEnv, inherits = FALSE), present)
                if (present) expect_identical(get(".Random.seed", .GlobalEnv), before_seed)
                expect_identical(RNGkind(), expected_kind)
                expect_identical(Sys.getlocale("LC_COLLATE"), before_locale)
            }
        }
    })
}

test_that("the original-reference fixture records the observed locale dependence", {
    f <- pcap232_fixture()
    a <- f$cases[["mixed__C"]]
    b <- f$cases[["mixed__en_US.UTF-8"]]
    expect_identical(names(a$second$calls)[a$second$calls != b$second$calls],
        c("oslo2emit0.019", "OSLO2EMIT0.093", "OSLO2EMIT0.097"))
    expect_equal(length(intersect(a$first$selected, b$first$selected)), 17L)
    expect_equal(length(intersect(a$second$selected, b$second$selected)), 16L)
    expect_false(identical(a$first$medians, b$first$medians))
    expect_false(identical(a$second$medians, b$second$medians))
    # This describes this captured fixture; it is not a universal three-call
    # guarantee for other data, RNG settings, R versions or collation libraries.
})
