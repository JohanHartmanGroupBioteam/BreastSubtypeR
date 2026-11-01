mapping218_input <- function(x, ids = rep("596", nrow(x))) {
    storage.mode(x) <- "double"
    rownames(x) <- paste0("p", seq_len(nrow(x)))
    if (is.null(colnames(x)) && ncol(x)) {
        colnames(x) <- paste0("sample ", seq_len(ncol(x)))
    }
    y <- data.frame(probe = rownames(x), ENTREZID = ids)
    list(x = x, y = y, se = SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = x),
        rowData = S4Vectors::DataFrame(y, Length = rep(1000, nrow(x))),
        colData = S4Vectors::DataFrame(PatientID = colnames(x),
            row.names = colnames(x))))
}

mapping218_expected <- function(input, row) {
    out <- input$x[row, , drop = FALSE]
    rownames(out) <- as.character(input$y$ENTREZID[match(
        rownames(input$x)[row], input$y$probe)])
    out
}

mapping218_reference <- function() {
    env <- new.env(parent = emptyenv())
    data("BreastSubtypeRobj", package = "BreastSubtypeR", envir = env)
    ref <- env$BreastSubtypeRobj
    ids <- intersect(as.character(ref$genes.sig50$EntrezGene.ID),
        as.character(ref$genes.signature$EntrezGene.ID[
            ref$genes.signature$SSP_based == "Yes"]))
    list(id = ids[[1L]], second = ids[[2L]], ref = ref,
        symbol = ref$genes.sig50$Symbol[match(ids[[1L]],
            ref$genes.sig50$EntrezGene.ID)])
}

test_that("variability selection skips undefined statistics in either order", {
    undefined <- list(c(NA, NA, NA), c(NaN, NaN, NaN),
        c(Inf, Inf, Inf))
    for (method in c("stdev", "iqr")) {
        for (invalid in undefined) {
            for (valid in list(c(1, 2, 3), c(7, 7, 7))) {
                for (bad_first in c(TRUE, FALSE)) {
                    x <- if (bad_first) rbind(invalid, valid) else rbind(valid, invalid)
                    input <- mapping218_input(x)
                    expect_identical(BreastSubtypeR:::duplicate_genes(
                        input$x, input$y, method),
                        mapping218_expected(input, if (bad_first) 2L else 1L))
                }
            }
        }
        input <- mapping218_input(rbind(c(Inf, Inf, Inf), c(NA, NA, NA),
            c(1, 2, 3), c(11, 12, 13)))
        expect_identical(BreastSubtypeR:::duplicate_genes(
            input$x, input$y, method), mapping218_expected(input, 3L))
    }
})

test_that("all-undefined duplicate statistics raise an identifiable diagnostic", {
    for (method in c("stdev", "iqr")) {
        for (x in list(rbind(c(NA, NA), c(NaN, NaN)),
            rbind(c(Inf, Inf), c(-Inf, -Inf)),
            rbind(c(Inf, Inf), c(NA, NA)), matrix(NA_real_, 3L, 2L))) {
            input <- mapping218_input(x, rep("100", nrow(x)))
            error <- expect_error(BreastSubtypeR:::duplicate_genes(
                input$x, input$y, method),
                class = "BreastSubtypeR_duplicate_undefined")
            expect_match(conditionMessage(error), "Entrez ID 100", fixed = TRUE)
            expect_match(conditionMessage(error), paste0('method = "', method, '"'),
                fixed = TRUE)
            expect_match(conditionMessage(error), paste(nrow(x), "probes"), fixed = TRUE)
            expect_match(conditionMessage(error), "NA/NaN", fixed = TRUE)
            expect_identical(error$entrez_id, "100")
            expect_identical(error$method, method)
            expect_identical(error$n_probes, nrow(x))
            expect_null(error$call)
        }
    }
})

test_that("zero, one and two samples have explicit variability boundaries", {
    empty <- mapping218_input(matrix(numeric(), nrow = 2L, ncol = 0L))
    for (method in c("stdev", "iqr")) {
        expect_error(BreastSubtypeR:::duplicate_genes(empty$x, empty$y, method),
            class = "BreastSubtypeR_duplicate_undefined")
    }
    for (values in list(c(9, 7), c(NA, 7), c(7, NA))) {
        input <- mapping218_input(matrix(values, ncol = 1L))
        expect_error(BreastSubtypeR:::duplicate_genes(input$x, input$y, "stdev"),
            class = "BreastSubtypeR_duplicate_undefined")
        row <- if (is.na(values[[1L]])) 2L else 1L
        expect_identical(BreastSubtypeR:::duplicate_genes(input$x, input$y, "iqr"),
            mapping218_expected(input, row))
    }
    for (x in list(rbind(c(9, NA), c(2, 4)), rbind(c(2, 4), c(10, 16)))) {
        input <- mapping218_input(x)
        for (method in c("stdev", "iqr")) {
            expect_identical(BreastSubtypeR:::duplicate_genes(
                input$x, input$y, method), mapping218_expected(input, 2L))
        }
    }
})

test_that("singletons retain all observations and zero-column dimensions", {
    x <- rbind(c(9, NA, NaN), c(NA, NA, NA), c(Inf, -Inf, Inf), c(1, 2, 3))
    for (columns in list(1:3, 1L, integer())) {
        input <- mapping218_input(x[, columns, drop = FALSE],
            c("200", "100", "300", "400"))
        for (method in c("stdev", "iqr")) {
            expect_identical(BreastSubtypeR:::duplicate_genes(
                input$x, input$y, method),
                mapping218_expected(input, seq_len(nrow(x))))
        }
    }
})

test_that("ties and gene ordering follow expression rows, not annotation rows", {
    for (duplicate_rows in list(rbind(c(1, 2, 3), c(11, 12, 13)),
        rbind(c(8, 8, 8), c(3, 3, 3)))) {
        input <- mapping218_input(rbind(c(7, 8, 9), duplicate_rows),
            c("200", "100", "100"))
        for (method in c("stdev", "iqr")) {
            expected <- mapping218_expected(input, 1:2)
            expect_identical(BreastSubtypeR:::duplicate_genes(
                input$x, input$y[3:1, ], method), expected)
            expect_identical(BreastSubtypeR:::duplicate_genes(
                input$x[c(3, 2, 1), , drop = FALSE], input$y, method),
                mapping218_expected(input, c(3L, 1L)))
        }
    }
})

test_that("SD, type-7 IQR and available-observation calculations are retained", {
    input <- mapping218_input(rbind(c(0, 0, 0, 8), c(0, 3, 3, 6)))
    ## Type-7 IQRs are 2 and 1.5; the type-1 ordering would be reversed.
    expect_identical(BreastSubtypeR:::duplicate_genes(input$x, input$y, "iqr"),
        mapping218_expected(input, 1L))
    input <- mapping218_input(rbind(c(1, NA, 3), c(10, 11, 12)))
    ## SDs are sqrt(2) and 1; IQRs tie at 1. Original holes stay missing.
    for (method in c("stdev", "iqr")) {
        expect_identical(BreastSubtypeR:::duplicate_genes(
            input$x, input$y, method), mapping218_expected(input, 1L))
    }
})

test_that("infinite statistics rank normally, including finite-input SD overflow", {
    cases <- list(
        list(method = "stdev", x = rbind(c(-1e308, 0, 1e308), c(-1, 0, 1))),
        list(method = "iqr", x = rbind(c(1, 2, Inf), c(1, 2, 3))),
        list(method = "iqr", x = rbind(c(1, 2, 3, 4, Inf), c(1, 2, 3, 4, 5))))
    for (case in cases) {
        input <- mapping218_input(case$x)
        expect_identical(BreastSubtypeR:::duplicate_genes(
            input$x, input$y, case$method), mapping218_expected(input, 1L))
    }
    for (method in c("stdev", "iqr")) {
        defined <- if (method == "stdev") c(-1e308, 0, 1e308) else c(1, 2, Inf)
        input <- mapping218_input(rbind(c(NA, NA, NA), defined, defined))
        expect_identical(BreastSubtypeR:::duplicate_genes(
            input$x, input$y, method), mapping218_expected(input, 2L))
    }
})

test_that("public Mapping selects the defined probe before optional KNN", {
    ref <- mapping218_reference()
    local_mocked_bindings(
        .mapping_donor_matrix = function(...) stop("unexpected donor preparation"),
        knn_vectorized_reference = function(...) stop("unexpected estimator"),
        .package = "BreastSubtypeR")
    for (method in c("stdev", "iqr")) {
        for (invalid in list(c(NA, NA, NA), c(Inf, Inf, Inf), c(9, NA, NA))) {
            input <- mapping218_input(rbind(invalid, c(1, 2, 3)),
                rep(ref$id, 2))
            for (impute in c(FALSE, TRUE)) {
                result <- Mapping(input$se, method = method,
                    impute = impute, verbose = FALSE)
                expect_identical(names(result), c("se_NC", "se_SSP"))
                nc <- SummarizedExperiment::assay(result$se_NC)
                ssp <- SummarizedExperiment::assay(result$se_SSP)
                expect_identical(as.numeric(nc[ref$symbol, ]), c(1, 2, 3))
                expect_identical(as.numeric(ssp[ref$id, ]), c(2, 4, 8))
                for (out in result) {
                    expect_identical(dim(out), c(1L, 3L))
                    expect_identical(colnames(out), colnames(input$x))
                    expect_identical(SummarizedExperiment::colData(out),
                        SummarizedExperiment::colData(input$se))
                }
            }
        }
    }
})

test_that("public Mapping propagates undefined protected groups and boundaries", {
    ref <- mapping218_reference()
    for (method in c("stdev", "iqr")) {
        for (x in list(matrix(NA_real_, 2L, 3L), matrix(numeric(), 2L, 0L))) {
            input <- mapping218_input(x, rep(ref$id, 2))
            for (impute in c(FALSE, TRUE)) {
                error <- expect_error(Mapping(input$se, method = method,
                    impute = impute, verbose = FALSE),
                    class = "BreastSubtypeR_duplicate_undefined")
                expect_identical(error$entrez_id, ref$id)
            }
        }
        singleton <- mapping218_input(matrix(7, 1L, 1L), ref$id)
        result <- Mapping(singleton$se, method = method, verbose = FALSE)
        expect_identical(as.numeric(SummarizedExperiment::assay(result$se_NC)), 7)
        empty <- mapping218_input(matrix(numeric(), 1L, 0L), ref$id)
        result <- Mapping(empty$se, method = method, verbose = FALSE)
        expect_identical(dim(result$se_NC), c(1L, 0L))
    }
    input <- mapping218_input(matrix(c(9, 7), 2L, 1L), rep(ref$id, 2))
    for (raw in c(FALSE, TRUE)) {
        expect_error(Mapping(input$se, RawCounts = raw, method = "stdev",
            verbose = FALSE), class = "BreastSubtypeR_duplicate_undefined")
        result <- Mapping(input$se, RawCounts = raw, method = "iqr", verbose = FALSE)
        expect_identical(dim(result$se_NC), c(1L, 1L))
        if (!raw) expect_identical(as.numeric(
            SummarizedExperiment::assay(result$se_NC)), 9)
    }
})

test_that("KNN excludes recognized background failures and protects required groups", {
    input <- mapping218_input(rbind(c(1, NA, 3), c(NA, NA, NA), c(NA, NA, NA)),
        c("596", "900001", "900001"))
    for (method in c("max", "stdev", "iqr")) {
        actual <- NULL
        warning <- expect_warning(actual <- BreastSubtypeR:::.mapping_donor_matrix(
            input$x, input$y, method, "596"),
            class = "BreastSubtypeR_knn_donor_exclusion")
        expect_identical(actual$matrix, mapping218_expected(input, 1L))
        expect_identical(actual$exclusions,
            data.frame(gene = "900001", reason = "duplicate_resolution_failed"))
        expect_identical(warning$exclusions, actual$exclusions)
        if (method == "max") {
            expect_error(BreastSubtypeR:::.mapping_donor_matrix(
                input$x, input$y, method, c("596", "900001")),
                'Cannot resolve Entrez ID 900001 with method = "max"', fixed = TRUE)
        } else {
            expect_error(BreastSubtypeR:::.mapping_donor_matrix(
                input$x, input$y, method, c("596", "900001")),
                class = "BreastSubtypeR_duplicate_undefined")
        }
    }
})

test_that("KNN propagates unrelated errors rather than matching variability messages", {
    input <- mapping218_input(rbind(c(NA, NA), c(NA, NA)), c("900001", "900001"))
    specific <- tryCatch(BreastSubtypeR:::duplicate_genes(input$x, input$y, "stdev"),
        error = identity)
    unrelated <- simpleError(conditionMessage(specific))
    local_mocked_bindings(duplicate_genes = function(...) stop(unrelated),
        .package = "BreastSubtypeR")
    for (method in c("max", "stdev", "iqr")) {
        for (required in list(character(), "900001")) {
            error <- tryCatch(BreastSubtypeR:::.mapping_donor_matrix(
                input$x, input$y, method, required), error = identity)
            expect_identical(error, unrelated)
        }
    }
})

test_that("public optional KNN excludes background groups without changing estimation or RNG", {
    ref <- mapping218_reference()
    donors <- matrix(rep(c(9, 20, 30), 10L), nrow = 10L, byrow = TRUE)
    input <- mapping218_input(rbind(c(1, NA, 3), donors,
        c(NA, NA, NA), c(NA, NA, NA)),
        c(ref$id, as.character(900001:900010), "900011", "900011"))
    clean <- input$se[seq_len(11L), ]
    withr::local_seed(218)
    for (method in c("max", "stdev", "iqr")) {
        seed <- .Random.seed
        expected <- Mapping(clean, method = method, impute = TRUE, verbose = FALSE)
        actual <- NULL
        expect_warning(actual <- Mapping(input$se, method = method,
            impute = TRUE, verbose = FALSE),
            class = "BreastSubtypeR_knn_donor_exclusion")
        expect_identical(actual, expected)
        expect_identical(as.numeric(SummarizedExperiment::assay(actual$se_NC)), c(1, 20, 3))
        expect_identical(as.numeric(SummarizedExperiment::assay(actual$se_SSP)), 2^c(1, 20, 3))
        expect_identical(.Random.seed, seed)
    }
    ## A protected nontarget group must fail during donor preparation too.
    required <- BreastSubtypeR:::.mapping_required_ids(ref$ref)
    local_mocked_bindings(.mapping_required_ids = function(...) c(required, "900011"),
        .package = "BreastSubtypeR")
    for (method in c("stdev", "iqr")) {
        error <- expect_error(Mapping(input$se, method = method,
            impute = TRUE, verbose = FALSE),
            class = "BreastSubtypeR_duplicate_undefined")
        expect_identical(error$entrez_id, "900011")
    }
})

test_that("the max default, its legacy error and RNG behavior are unchanged", {
    ref <- mapping218_reference()
    input <- mapping218_input(rbind(c(9, 9, 9), c(1, 2, 3)), rep(ref$id, 2))
    withr::local_seed(218)
    seed <- .Random.seed
    expect_identical(Mapping(input$se, verbose = FALSE),
        Mapping(input$se, method = "max", impute = FALSE, verbose = FALSE))
    expect_identical(as.numeric(SummarizedExperiment::assay(
        Mapping(input$se, verbose = FALSE)$se_NC)), c(9, 9, 9))
    input <- mapping218_input(matrix(NA_real_, 2L, 3L))
    error <- tryCatch(BreastSubtypeR:::duplicate_genes(input$x, input$y, "max"),
        error = identity)
    expect_s3_class(error, "simpleError")
    expect_false(inherits(error, "BreastSubtypeR_duplicate_undefined"))
    expect_identical(conditionMessage(error), paste0(
        'Cannot resolve Entrez ID 596 with method = "max": ',
        'none of its 2 probes has a defined total across samples ',
        '(all row sums are NA/NaN). Review missing or nonfinite ',
        'measurements and the preprocessing strategy. ',
        'No probe was selected for this Entrez ID.'))
    for (method in c("stdev", "iqr")) {
        tryCatch(BreastSubtypeR:::duplicate_genes(input$x, input$y, method),
            error = function(e) NULL)
    }
    expect_identical(.Random.seed, seed)
})
