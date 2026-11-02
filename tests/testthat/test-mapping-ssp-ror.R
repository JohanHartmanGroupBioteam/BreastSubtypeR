mapping131_models <- function() {
    env <- new.env(parent = emptyenv())
    data(list = c("BreastSubtypeRobj", "sspbc.models"),
        package = "BreastSubtypeR", envir = env)
    pairs <- strsplit(env$sspbc.models$ssp.ror$all.pairs, "<", fixed = TRUE)
    ids <- unique(sub("^e", "", unlist(pairs, use.names = FALSE)))
    ref <- env$BreastSubtypeRobj
    legacy <- as.character(ref$genes.signature$EntrezGene.ID)
    shared <- intersect(as.character(ref$genes.sig50$EntrezGene.ID),
        as.character(ref$genes.signature$EntrezGene.ID[
            ref$genes.signature$SSP_based == "Yes"]))[[1L]]
    list(model = env$sspbc.models$ssp.ror, ids = ids, legacy = legacy,
        ror_only = setdiff(ids, legacy)[[1L]], shared = shared, ref = ref)
}

mapping131_se <- function(x, ids = rownames(x)) {
    ids <- as.character(ids)
    storage.mode(x) <- "double"
    rownames(x) <- paste0("probe", seq_len(nrow(x)))
    if (is.null(colnames(x))) colnames(x) <- paste0("sample ", seq_len(ncol(x)))
    SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = x),
        rowData = S4Vectors::DataFrame(probe = rownames(x), ENTREZID = ids,
            Length = 1000 + 100 * seq_len(nrow(x))),
        colData = S4Vectors::DataFrame(PatientID = colnames(x),
            row.names = colnames(x)))
}

mapping131_complete <- function() {
    model <- mapping131_models()
    ids <- unique(c(model$legacy, model$ids))
    x <- outer(seq_along(ids), 1:3, function(i, j) 4 + sin(i * j) + j / 10)
    dimnames(x) <- list(ids, c("sample B", "sample/A", "s.3"))
    mapping131_se(x)
}

test_that("SSP-ROR retention is opt-in and validates its logical option", {
    expect_identical(formals(Mapping)$include_ssp_ror, FALSE)
    for (bad in list(NULL, NA, logical(), c(TRUE, FALSE), 1, "TRUE")) {
        expect_error(Mapping(stop("input evaluated"), include_ssp_ror = bad),
            "include_ssp_ror must be one nonmissing logical value.", fixed = TRUE)
    }
    se <- mapping131_complete()
    expected <- Mapping(se, verbose = FALSE)
    expect_identical(names(expected), c("se_NC", "se_SSP"))
    expect_identical(Mapping(se, verbose = FALSE, include_ssp_ror = FALSE), expected)
    # Appending the option preserves the existing positional arguments.
    expect_identical(Mapping(se, FALSE, "max", FALSE, FALSE), expected)
})

test_that("the separate ROR output retains all bundled model genes and rules", {
    model <- mapping131_models()
    se <- mapping131_complete()
    before <- se
    expect_length(model$ids, 296L)
    expect_length(model$model$all.pairs, 361L)
    for (impute in c(FALSE, TRUE)) {
        legacy <- Mapping(se, impute = impute, verbose = FALSE)
        actual <- Mapping(se, impute = impute, verbose = FALSE, include_ssp_ror = TRUE)
        expect_identical(names(actual), c("se_NC", "se_SSP", "se_SSP_ROR"))
        expect_identical(actual[c("se_NC", "se_SSP")], legacy)
        ror <- actual$se_SSP_ROR
        expect_s4_class(ror, "SummarizedExperiment")
        expect_identical(SummarizedExperiment::assayNames(ror), "counts")
        expect_identical(rownames(ror), model$ids)
        expect_identical(colnames(ror), colnames(se))
        expect_identical(SummarizedExperiment::colData(ror), SummarizedExperiment::colData(se))
        expected <- 2^SummarizedExperiment::assay(se)[
            match(model$ids, SummarizedExperiment::rowData(se)$ENTREZID), , drop = FALSE]
        rownames(expected) <- model$ids
        expect_identical(SummarizedExperiment::assay(ror), expected)
        native <- suppressMessages(BS_sspbc(ror, ssp.name = "ssp.ror"))
        expect_identical(rownames(native$rules.matrix), model$model$all.pairs)
        expect_identical(dim(native$rules.matrix), c(361L, ncol(se)))
        expect_false(anyNA(native$rules.matrix))
        direct <- suppressMessages(BreastSubtypeR:::applySSP(gex = expected,
            id = rownames(expected), id.type = "EntrezGene", ssp.name = "ssp.ror"))
        expect_identical(native, direct)
    }
    expect_identical(se, before)
})

test_that("existing subtype results remain identical with a separate ROR output", {
    se <- mapping131_complete()
    legacy <- Mapping(se, verbose = FALSE)
    actual <- Mapping(se, verbose = FALSE, include_ssp_ror = TRUE)
    expect_identical(BS_parker(actual$se_NC), BS_parker(legacy$se_NC))
    expect_identical(BS_AIMS(actual$se_SSP), BS_AIMS(legacy$se_SSP))
    for (key in c("ssp.pam50", "ssp.subtype")) {
        expect_identical(BS_sspbc(actual$se_SSP, key), BS_sspbc(legacy$se_SSP, key))
    }
    run_multi <- function(input) suppressMessages(BS_Multi(input,
        methods = c("parker.original", "AIMS", "sspbc"), Subtype = TRUE))
    expect_identical(run_multi(actual), run_multi(legacy))
    expect_identical(InputCoverage(actual), InputCoverage(legacy))
})

test_that("ROR output keeps sample dimensions and does not invent absent genes", {
    model <- mapping131_models()
    singleton <- mapping131_complete()[, 2L, drop = FALSE]
    for (method in c("max", "mean", "median", "iqr", "stdev")) {
        for (impute in c(FALSE, TRUE)) {
            actual <- Mapping(singleton, method = method, impute = impute,
                verbose = FALSE, include_ssp_ror = TRUE)
            expect_identical(dim(actual$se_SSP_ROR), c(296L, 1L))
            expect_identical(colnames(actual$se_SSP_ROR), colnames(singleton))
            expect_identical(actual[c("se_NC", "se_SSP")],
                Mapping(singleton, method = method, impute = impute, verbose = FALSE))
        }
    }
    ssp <- as.character(model$ref$genes.signature$EntrezGene.ID[
        model$ref$genes.signature$SSP_based == "Yes"])
    absent <- mapping131_se(matrix(c(1, 2, 3), nrow = 1L), setdiff(ssp, model$ids)[[1L]])
    for (impute in c(FALSE, TRUE)) {
        actual <- Mapping(absent, impute = impute, verbose = FALSE, include_ssp_ror = TRUE)
        expect_identical(dim(actual$se_SSP_ROR), c(0L, 3L))
        expect_identical(colnames(actual$se_SSP_ROR), colnames(absent))
        expect_identical(actual[c("se_NC", "se_SSP")],
            Mapping(absent, impute = impute, verbose = FALSE))
    }
})

test_that("ROR uses the unchanged duplicate rules on both established scales", {
    model <- mapping131_models()
    x <- rbind(c(2, 8, 3), c(6, 3, 7), c(4, 9, 2), c(7, 3, 8), c(3, 4, 5))
    ids <- c(model$shared, model$shared, model$ror_only, model$ror_only, "900001")
    se <- mapping131_se(x, ids)
    annotation <- as.data.frame(SummarizedExperiment::rowData(se))
    for (raw in c(FALSE, TRUE)) {
        source <- SummarizedExperiment::assay(se)
        if (raw) source <- edgeR::rpkm(edgeR::DGEList(counts = round(source)),
            gene.length = annotation$Length, normalized.lib.sizes = FALSE, log = FALSE)
        for (method in c("max", "mean", "median", "iqr", "stdev")) {
            collapsed <- BreastSubtypeR:::duplicate_genes(source, annotation, method)
            expected <- collapsed[intersect(model$ids, rownames(collapsed)), , drop = FALSE]
            if (!raw) expected <- 2^expected
            for (impute in c(FALSE, TRUE)) {
                legacy <- Mapping(se, RawCounts = raw, method = method,
                    impute = impute, verbose = FALSE)
                actual <- Mapping(se, RawCounts = raw, method = method,
                    impute = impute, verbose = FALSE, include_ssp_ror = TRUE)
                expect_identical(actual[c("se_NC", "se_SSP")], legacy)
                expect_identical(SummarizedExperiment::assay(actual$se_SSP_ROR), expected)
            }
        }
    }
})

test_that("ROR targets use original full-input donors without changing existing estimates", {
    model <- mapping131_models()
    x <- rbind(c(0, NA, 2, 3), c(1, 2, NA, 4),
        t(vapply(1:12, function(i) c(i / 10, 2 + i / 10, 3 + i / 10, 4), numeric(4))))
    rownames(x) <- c(model$shared, model$ror_only, as.character(900001:900012))
    colnames(x) <- c("sample B", "sample/A", "s.3", "four")
    se <- mapping131_se(x)
    before <- se
    frozen <- BreastSubtypeR:::knn_vectorized_reference
    observed <- NULL
    local_mocked_bindings(knn_vectorized_reference = function(x, targets) {
        result <- frozen(x, targets)
        observed <<- list(x = x, targets = targets, result = result)
        result
    }, .package = "BreastSubtypeR")
    legacy <- Mapping(se, impute = TRUE, verbose = FALSE)
    actual <- Mapping(se, impute = TRUE, verbose = FALSE, include_ssp_ror = TRUE)
    expect_identical(actual[c("se_NC", "se_SSP")], legacy)
    expect_identical(observed$x, x)
    expect_identical(observed$targets, c(model$shared, model$ror_only))
    expect_identical(observed$result, frozen(x, observed$targets))
    expected <- 2^observed$result$matrix[
        intersect(model$ids, rownames(x)), , drop = FALSE]
    expect_identical(SummarizedExperiment::assay(actual$se_SSP_ROR), expected)
    disabled <- Mapping(se, impute = FALSE, verbose = FALSE, include_ssp_ror = TRUE)
    expect_true(is.na(SummarizedExperiment::assay(disabled$se_SSP_ROR)[model$ror_only, 3]))
    expect_setequal(rownames(actual$se_SSP_ROR), intersect(model$ids, rownames(x)))
    expect_identical(se, before)
})

test_that("only opted-in represented ROR holes trigger KNN and its support checks", {
    model <- mapping131_models()
    x <- rbind(c(1, 2, 3), c(0, NA, 2),
        matrix(rep(c(0, 5, 2), 10L), 10L, 3L, byrow = TRUE))
    ids <- c(model$shared, model$ror_only, as.character(900001:900010))
    se <- mapping131_se(x, ids)
    legacy <- Mapping(se, impute = TRUE, verbose = FALSE)
    actual <- Mapping(se, impute = TRUE, verbose = FALSE, include_ssp_ror = TRUE)
    expect_identical(actual[c("se_NC", "se_SSP")], legacy)
    expect_identical(as.numeric(SummarizedExperiment::assay(actual$se_SSP_ROR)[
        model$ror_only, ]), 2^c(0, 5, 2))
    SummarizedExperiment::assay(se)[2, ] <- NA_real_
    error <- expect_error(Mapping(se, impute = TRUE, verbose = FALSE,
        include_ssp_ror = TRUE), class = "BreastSubtypeR_knn_unsupported")
    expect_true(model$ror_only %in% error$targets)
    expect_identical(Mapping(se, impute = TRUE, verbose = FALSE), legacy)
    # Estimates must remain finite after the same SSP back-transformation.
    SummarizedExperiment::assay(se)[2, ] <- c(0, NA, 2)
    SummarizedExperiment::assay(se)[3:12, 2] <- 2048
    expect_error(Mapping(se, impute = TRUE, verbose = FALSE, include_ssp_ror = TRUE),
        class = "BreastSubtypeR_knn_unsupported")
})

test_that("ROR required duplicate errors propagate only when requested", {
    model <- mapping131_models()
    se <- mapping131_se(rbind(c(1, 2, 3), rep(NA_real_, 3), rep(NA_real_, 3)),
        c(model$shared, model$ror_only, model$ror_only))
    for (method in c("max", "iqr", "stdev")) {
        for (impute in c(FALSE, TRUE)) {
            expect_no_error(Mapping(se, method = method, impute = impute, verbose = FALSE))
            error <- tryCatch(Mapping(se, method = method, impute = impute,
                verbose = FALSE, include_ssp_ror = TRUE), error = identity)
            expect_s3_class(error, "error")
            expect_match(conditionMessage(error), paste("Entrez ID", model$ror_only), fixed = TRUE)
            if (method != "max") expect_s3_class(error, "BreastSubtypeR_duplicate_undefined")
        }
    }
})

test_that("raw ROR gaps are never imputed and missing genes are never synthesized", {
    model <- mapping131_models()
    se <- mapping131_se(rbind(c(100, 200, 300), c(200, 300, 100)),
        c(model$shared, model$ror_only))
    local_mocked_bindings(knn_vectorized_reference = function(...) stop("raw KNN called"),
        .package = "BreastSubtypeR")
    SummarizedExperiment::rowData(se)$Length[2] <- NA_real_
    actual <- Mapping(se, RawCounts = TRUE, impute = FALSE, verbose = FALSE,
        include_ssp_ror = TRUE)
    expect_true(all(is.na(SummarizedExperiment::assay(actual$se_SSP_ROR)[model$ror_only, ])))
    expect_error(Mapping(se, RawCounts = TRUE, impute = TRUE, verbose = FALSE,
        include_ssp_ror = TRUE), "neither counts nor linear FPKM will be imputed", fixed = TRUE)
    expect_no_error(Mapping(se, RawCounts = TRUE, impute = TRUE, verbose = FALSE))
    SummarizedExperiment::assay(se)[2, 1] <- NA_real_
    for (impute in c(FALSE, TRUE)) {
        expect_error(Mapping(se, RawCounts = TRUE, impute = impute, verbose = FALSE,
            include_ssp_ror = TRUE), "raw-count imputation is not supported", fixed = TRUE)
    }
})

test_that("ROR retention and imputation preserve caller RNG state", {
    model <- mapping131_models()
    x <- rbind(c(1, 2, 3), c(0, NA, 2),
        matrix(rep(c(0, 5, 2), 10L), 10L, 3L, byrow = TRUE))
    se <- mapping131_se(x, c(model$shared, model$ror_only, as.character(900001:900010)))
    withr::local_preserve_seed()
    kind <- RNGkind()
    for (seed_present in c(TRUE, FALSE)) {
        if (seed_present) set.seed(131) else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv)
        seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        for (impute in c(FALSE, TRUE)) {
            Mapping(se, impute = impute, verbose = FALSE, include_ssp_ror = TRUE)
            expect_identical(RNGkind(), kind)
            expect_identical(exists(".Random.seed", .GlobalEnv, inherits = FALSE), seed_present)
            expect_identical(get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE), seed)
        }
        SummarizedExperiment::assay(se)[2, ] <- NA_real_
        expect_error(Mapping(se, impute = TRUE, verbose = FALSE, include_ssp_ror = TRUE),
            class = "BreastSubtypeR_knn_unsupported")
        expect_identical(RNGkind(), kind)
        expect_identical(get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE), seed)
        SummarizedExperiment::assay(se)[2, ] <- c(0, NA, 2)
    }
})
