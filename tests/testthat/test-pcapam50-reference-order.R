test_that("second-stage candidates follow the reference character-key merge", {
    ids <- c("a2", "A10", "a-1", "a_1", "a.2", "a.10", "B2", "b1", "10", "2")
    labels <- setNames(c(rep("LumA", 6), rep("Basal", 4)), ids)
    original_locale <- Sys.getlocale("LC_COLLATE")
    captured <- NULL
    selected <- NULL
    draw <- function(tab, seed) {
        positive <- as.character(tab$PatientID[tab$PAM50 == "LumA"])
        negative <- as.character(tab$PatientID[tab$PAM50 == "Basal"])
        withr::with_seed(seed, c(positive[sample(length(positive), length(negative))], negative))
    }
    withr::local_seed(17)
    testthat::local_mocked_bindings(
        makeCalls.PC1ihc = function(mat, df.cln, ...) {
            expect_identical(colnames(mat), as.character(df.cln$PatientID))
            # The native helper returns IDs from names(out$predictions).
            list(BS.all = data.frame(PatientID = colnames(mat),
                 BS = unname(labels[colnames(mat)])))
        },
        makeCalls.v1PAM = function(mat, df.pam, seed, ...) {
            captured <<- df.pam
            selected <<- draw(df.pam, seed)
            stop("captured second-stage input", call. = FALSE)
        },
        .package = "BreastSubtypeR"
    )
    for (locale in unique(c(original_locale, "C"))) {
        withr::with_locale(c(LC_COLLATE = locale), {
            expect_identical(Sys.getlocale("LC_COLLATE"), locale)
            for (clinical in c(FALSE, TRUE)) {
                pheno <- data.frame(
                    PatientID = ids, ER = rep(c("ER+", "ER-"), 5),
                    TSIZE = rep(c(0, 1), 5), NODE = as.numeric(seq_along(ids) - 1L),
                    row.names = ids
                )
                se <- SummarizedExperiment::SummarizedExperiment(
                    assays = list(expr = matrix(seq_len(20), nrow = 2,
                        dimnames = list(c("gene1", "gene2"), ids))),
                    colData = pheno
                )
                captured <- selected <- NULL
                before <- get(".Random.seed", .GlobalEnv)
                kind <- RNGkind()
                expect_error(BS_PCAPAM50(se, hasClinical = clinical, seed = 118),
                             "captured second-stage input", fixed = TRUE)
                expect_identical(get(".Random.seed", .GlobalEnv), before)
                expect_identical(RNGkind(), kind)
                expect_true(is.data.frame(captured))
                if (!is.data.frame(captured)) stop("Second-stage capture did not occur")
                complete <- is.character(captured$PatientID) &&
                    length(captured$PatientID) == length(ids) &&
                    !anyNA(captured$PatientID) && !anyDuplicated(captured$PatientID) &&
                    setequal(captured$PatientID, ids)
                expect_true(complete)
                if (!complete) stop("Incomplete second-stage identities")
                matched <- match(ids, captured$PatientID)
                expect_false(anyNA(matched))
                expect_identical(captured$PAM50[matched], unname(labels))
                expect_identical(names(captured), c("PatientID", "PAM50",
                                                   if (clinical) c("TSIZE", "NODE")))
                if (clinical) {
                    expect_identical(captured$TSIZE[matched], pheno$TSIZE)
                    expect_identical(captured$NODE[matched], pheno$NODE)
                }
                expected <- merge(data.frame(PatientID = ids, PAM50 = unname(labels)),
                                  pheno, by = "PatientID")
                expect_identical(captured$PatientID, expected$PatientID,
                                 info = "intermediate candidate order")
                expect_identical(selected, draw(expected, 118),
                                 info = "second-stage selected IDs")
                expect_length(selected, 8L)
                expect_false(anyNA(selected))
                expect_identical(anyDuplicated(selected), 0L)
                expect_true(all(selected %in% ids))
            }
        })
        expect_identical(Sys.getlocale("LC_COLLATE"), original_locale)
    }
})

test_that("intermediate character and factor keys retain merge order and payloads", {
    ids <- c("a2", "A10", "a-1", "a_1", "a.2", "a.10", "B2", "b1", "10", "2")
    original_locale <- Sys.getlocale("LC_COLLATE")
    draw <- function(tab) {
        positive <- as.character(tab$PatientID[tab$PAM50 == "LumA"])
        negative <- as.character(tab$PatientID[tab$PAM50 == "Basal"])
        withr::with_seed(118, c(positive[sample(length(positive), length(negative))], negative))
    }
    for (locale in unique(c(original_locale, "C"))) {
        withr::with_locale(c(LC_COLLATE = locale), {
            expect_identical(Sys.getlocale("LC_COLLATE"), locale)
            for (factor_keys in c(FALSE, TRUE)) {
                tab <- data.frame(
                    PatientID = if (factor_keys) factor(ids, levels = rev(ids)) else ids,
                    PAM50 = c(rep("LumA", 6), rep("Basal", 4)),
                    TSIZE = rep(c(0, 1), 5), NODE = as.numeric(seq_along(ids) - 1L)
                )
                # The reference merge has a character left key. This checks
                # table semantics, not native factor-ID input to the wrapper.
                expected <- merge(data.frame(PatientID = ids), tab, by = "PatientID")
                actual <- tab[order(as.character(tab$PatientID)), , drop = FALSE]
                actual_ids <- as.character(actual$PatientID)
                complete <- length(actual_ids) == length(ids) && !anyNA(actual_ids) &&
                    !anyDuplicated(actual_ids) && setequal(actual_ids, ids) &&
                    length(expected$PatientID) == length(ids) &&
                    !anyNA(expected$PatientID) && !anyDuplicated(expected$PatientID) &&
                    setequal(expected$PatientID, ids)
                expect_true(complete)
                if (!complete) stop("Incomplete intermediate table identities")
                expect_identical(actual_ids, expected$PatientID)
                matched <- match(ids, actual_ids)
                expect_false(anyNA(matched))
                for (column in c("PAM50", "TSIZE", "NODE")) {
                    expect_identical(actual[[column]], expected[[column]])
                    expect_identical(actual[[column]][matched], tab[[column]])
                }
                selected <- draw(actual)
                expect_identical(selected, draw(expected))
                expect_length(selected, 8L)
                expect_false(anyNA(selected))
                expect_identical(anyDuplicated(selected), 0L)
                expect_true(all(selected %in% ids))
            }
        })
        expect_identical(Sys.getlocale("LC_COLLATE"), original_locale)
    }
})

test_that("reference ordering preserves caller RNG state and package four-class logic", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR")
    se <- OSLO2EMIT0obj$data_input$se_NC
    se <- se[, rev(seq_len(ncol(se))), drop = FALSE]
    withr::local_preserve_seed()
    original_kind <- RNGkind()
    withr::defer(do.call(RNGkind, as.list(original_kind)))
    RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
    kind <- RNGkind()
    locale <- Sys.getlocale()
    for (seed_exists in c(FALSE, TRUE)) {
        if (seed_exists) set.seed(17) else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv)
        before <- if (seed_exists) get(".Random.seed", .GlobalEnv) else NULL
        res <- expect_experimental_ror_pc(BS_PCAPAM50(se,
            Subtype = TRUE, hasClinical = TRUE, seed = 118))
        expect_identical(RNGkind(), kind)
        expect_identical(Sys.getlocale(), locale)
        expect_identical(exists(".Random.seed", .GlobalEnv, inherits = FALSE), seed_exists)
        if (seed_exists) expect_identical(get(".Random.seed", .GlobalEnv), before)
        expect_identical(res$BS.all$PatientID, colnames(se))
        expect_identical(names(res$outList$predictions.FourSubtype), colnames(se))
        expected <- res$outList$predictions
        normal <- expected == "Normal"
        distances <- res$outList$distances[, c("Basal", "Her2", "LumA", "LumB"), drop = FALSE]
        expected[normal] <- colnames(distances)[max.col(distances[normal, , drop = FALSE], ties.method = "first")]
        expect_identical(res$outList$predictions.FourSubtype, expected)
        expect_false(anyNA(res$BS.all$BS))
    }
})
