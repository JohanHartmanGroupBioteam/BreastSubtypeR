coverage_test_models <- function() {
    env <- new.env(parent = emptyenv())
    data(list = c("BreastSubtypeRobj", "AIMSmodel", "sspbc.models", "Gene.ID.ann"),
        package = "BreastSubtypeR", envir = env)
    list(
        nc = env$BreastSubtypeRobj$genes.sig50,
        aims_signature = as.character(env$BreastSubtypeRobj$genes.signature$EntrezGene.ID[
            env$BreastSubtypeRobj$genes.signature$AIMS == "Yes"]),
        models = list(AIMS = env$AIMSmodel,
            ssp.subtype = env$sspbc.models$ssp.subtype,
            ssp.pam50 = env$sspbc.models$ssp.pam50),
        annotation = env$Gene.ID.ann
    )
}

coverage_test_genes <- function(model) {
    unique(unlist(strsplit(model$all.pairs, "<", fixed = TRUE), use.names = FALSE))
}

coverage_test_input <- function(n = 3L) {
    models <- coverage_test_models()
    genes <- unique(unlist(lapply(models$models, function(model) {
        sub("^e", "", coverage_test_genes(model))
    }), use.names = FALSE))
    make_se <- function(ids) {
        values <- matrix(seq_len(length(ids) * n), nrow = length(ids), ncol = n,
            dimnames = list(ids, paste0("sample ", seq_len(n))))
        SummarizedExperiment::SummarizedExperiment(assays = list(counts = values))
    }
    list(se_NC = make_se(as.character(models$nc$Symbol)), se_SSP = make_se(genes))
}

# This reference follows the existing classifier path, without calling any
# reporting helper. Native rules.matrix values define evaluability.
coverage_native_reference <- function(input, key, models = coverage_test_models()) {
    x <- as.matrix(SummarizedExperiment::assay(input$se_SSP))
    ids <- rownames(x)
    model <- models$models[[key]]
    if (key == "AIMS") {
        keep <- ids %in% models$aims_signature
        x <- x[keep, , drop = FALSE]
        ids <- ids[keep]
        compare <- BreastSubtypeR:::.comp.sel.pairs_AIMS
    } else {
        ids <- BreastSubtypeR:::translate_id2entrez(
            ids, models$annotation, id.type = "EntrezGene", e = TRUE)
        compare <- BreastSubtypeR:::.comp.sel.pairs
    }
    genes <- coverage_test_genes(model)
    keep <- ids %in% genes
    collapsed <- BreastSubtypeR:::.removeDuplicatedEntrezPerPatients(
        x[keep, , drop = FALSE], ids[keep])
    rules <- suppressMessages(compare(
        list(exprs = collapsed$dataset, GeneName = collapsed$EntrezID),
        model$all.pairs))
    list(genes_present = as.integer(length(intersect(genes, ids))),
        rules = rules, evaluable = colSums(!is.na(rules)))
}

expect_coverage_native <- function(input) {
    models <- coverage_test_models()
    actual <- InputCoverage(input)
    for (i in seq_along(models$models)) {
        key <- names(models$models)[i]
        native <- coverage_native_reference(input, key, models)
        expect_identical(actual$genes_present[i + 1L], native$genes_present, info = key)
        expect_equal(actual$rules_evaluable_min[i + 1L], min(native$evaluable), info = key)
        expect_equal(actual$rules_evaluable_max[i + 1L], max(native$evaluable), info = key)
    }
    invisible(actual)
}

test_that("coverage has four stable rows and derives independent model denominators", {
    models <- coverage_test_models()
    actual <- InputCoverage(coverage_test_input())
    expect_s3_class(actual, "data.frame")
    expect_identical(names(actual), c("model", "genes_present", "genes_total",
        "genes_absent", "rules_evaluable_min", "rules_evaluable_max", "rules_total"))
    expect_identical(actual$model, c("NC-based (PAM50)", "AIMS",
        "SSPBC (ssp.subtype)", "SSPBC (ssp.pam50)"))
    expect_identical(actual$genes_total, unname(c(as.integer(length(unique(models$nc$Symbol))),
        vapply(models$models, function(model) length(coverage_test_genes(model)), integer(1)))))
    expect_identical(actual$rules_total, unname(c(NA_integer_,
        vapply(models$models, function(model) length(model$all.pairs), integer(1)))))
    # Audited bundled-model regression anchors; production derives these
    # counts from the packaged objects, as independently checked above.
    expect_identical(actual$genes_total, c(50L, 151L, 153L, 216L))
    expect_identical(actual$rules_total, c(NA_integer_, 100L, 84L, 120L))
    expect_identical(actual$genes_present, actual$genes_total)
    expect_identical(actual$genes_absent, rep(0L, 4L))
    expect_identical(actual$rules_evaluable_min, actual$rules_total)
    expect_identical(actual$rules_evaluable_max, actual$rules_total)
    expect_true(all(is.na(actual[1L, c("rules_evaluable_min", "rules_evaluable_max", "rules_total")])))
    expect_true(all(vapply(actual[-1L], is.integer, logical(1))))
    for (model in models$models) {
        expect_gt(length(model$all.pairs), max(model$k))
        for (k in model$k) {
            expect_setequal(names(model$one.vs.all.tsp[[k]]$tables), make.names(model$all.pairs))
        }
    }
    expect_false(identical(coverage_test_genes(models$models$ssp.subtype),
        coverage_test_genes(models$models$ssp.pam50)))
})

test_that("structural gene coverage counts IDs, including present missing-valued rows", {
    input <- coverage_test_input()
    input$se_NC <- input$se_NC[-seq_len(7L), , drop = FALSE]
    SummarizedExperiment::assay(input$se_NC)[] <- NA_real_
    actual <- InputCoverage(input)
    expect_identical(actual$genes_present[1L], nrow(input$se_NC))
    expect_identical(actual$genes_absent[1L], 7L)
    expect_true(all(is.na(actual[1L, 5:7])))
    SummarizedExperiment::assay(input$se_SSP)[] <- NA_real_
    actual <- expect_coverage_native(input)
    expect_identical(actual$genes_present[-1L], actual$genes_total[-1L])
    expect_identical(actual$rules_evaluable_min[-1L], rep(0L, 3L))
    expect_identical(actual$rules_evaluable_max[-1L], rep(0L, 3L))
})

test_that("missing model-specific genes remove their actual incident rules", {
    models <- coverage_test_models()
    plain_genes <- lapply(models$models, function(model) sub("^e", "", coverage_test_genes(model)))
    for (i in seq_along(models$models)) {
        model <- models$models[[i]]
        other_genes <- unique(unlist(plain_genes[-i], use.names = FALSE))
        unique_gene <- setdiff(plain_genes[[i]], other_genes)[1L]
        expect_false(is.na(unique_gene), info = names(models$models)[i])
        input <- coverage_test_input()
        input$se_SSP <- input$se_SSP[rownames(input$se_SSP) != unique_gene, , drop = FALSE]
        actual <- expect_coverage_native(input)
        expect_identical(actual$genes_absent[-1L], as.integer(seq_along(models$models) == i))
        rule_parts <- strsplit(model$all.pairs, "<", fixed = TRUE)
        affected <- sum(vapply(rule_parts, function(pair) {
            unique_gene %in% sub("^e", "", pair)
        }, logical(1)))
        expect_identical(actual$rules_evaluable_min[i + 1L], length(model$all.pairs) - affected)

        multiplicities <- table(unlist(rule_parts, use.names = FALSE))
        selected <- c(names(which.min(multiplicities)), names(which.max(multiplicities)))
        expect_gt(max(multiplicities), min(multiplicities))
        for (gene in selected) {
            input <- coverage_test_input(1L)
            input$se_SSP <- input$se_SSP[rownames(input$se_SSP) != sub("^e", "", gene), , drop = FALSE]
            actual <- expect_coverage_native(input)
            affected <- sum(vapply(rule_parts, function(pair) gene %in% pair, logical(1)))
            expect_identical(actual$rules_evaluable_min[i + 1L], length(model$all.pairs) - affected)
        }
    }
})

test_that("coverage respects native ties, missing comparisons and infinities per sample", {
    input <- coverage_test_input(8L)
    x <- SummarizedExperiment::assay(input$se_SSP)
    x[1L, 2L] <- NA_real_
    x[2L, 3L] <- NaN
    x[, 4L] <- NA_real_
    x[, 5L] <- 0
    x[, 6L] <- rep(c(Inf, -Inf), length.out = nrow(x))
    x[, 7L] <- Inf
    x[, 8L] <- -Inf
    SummarizedExperiment::assay(input$se_SSP) <- x
    actual <- expect_coverage_native(input)
    expect_identical(actual$rules_evaluable_min[-1L], rep(0L, 3L))
    expect_identical(actual$rules_evaluable_max[-1L], actual$rules_total[-1L])
    for (j in seq_len(ncol(input$se_SSP))) {
        single <- lapply(input, function(se) se[, j, drop = FALSE])
        observed <- expect_coverage_native(single)
        expect_identical(observed$rules_evaluable_min, observed$rules_evaluable_max)
        if (j %in% c(1L, 5L, 6L, 7L, 8L)) {
            expect_identical(observed$rules_evaluable_min[-1L], observed$rules_total[-1L])
        }
    }
    for (key in names(coverage_test_models()$models)) {
        native <- coverage_native_reference(input, key)$rules
        expect_true(all(native[, 5L] == 0))
        expect_true(all(native[, 7L] == 0))
        expect_true(all(native[, 8L] == 0))
    }
})

test_that("per-sample ranges retain differing nonzero counts", {
    input <- coverage_test_input(2L)
    models <- coverage_test_models()
    for (model in models$models) {
        gene <- sub("^e", "", coverage_test_genes(model)[1L])
        SummarizedExperiment::assay(input$se_SSP)[gene, 2L] <- NA_real_
    }
    actual <- expect_coverage_native(input)
    expect_true(all(actual$rules_evaluable_min[-1L] > 0L))
    expect_true(all(actual$rules_evaluable_min[-1L] < actual$rules_evaluable_max[-1L]))
    expect_identical(actual$rules_evaluable_max[-1L], actual$rules_total[-1L])
})

test_that("rules can be entirely absent or only partly represented", {
    input <- coverage_test_input(1L)
    for (i in c(0L, 1L, 2L)) {
        subset <- input
        subset$se_SSP <- subset$se_SSP[seq_len(i), , drop = FALSE]
        actual <- expect_coverage_native(subset)
        if (i < 2L) expect_identical(actual$rules_evaluable_max[-1L], rep(0L, 3L))
    }
    empty <- lapply(input, function(se) se[integer(), , drop = FALSE])
    actual <- expect_coverage_native(empty)
    expect_identical(actual$genes_present, rep(0L, 4L))
    expect_identical(actual$genes_absent, actual$genes_total)
    expect_identical(actual$rules_evaluable_min[-1L], rep(0L, 3L))
    expect_true(all(is.na(actual[1L, 5:7])))
})

test_that("coverage aligns identifiers and preserves native duplicate handling", {
    input <- coverage_test_input(3L)
    expected <- InputCoverage(input)
    unnamed <- lapply(input, function(se) {
        colnames(se) <- NULL
        se
    })
    expect_identical(InputCoverage(unnamed), expected)
    permuted <- lapply(input, function(se) se[rev(seq_len(nrow(se))), , drop = FALSE])
    expect_identical(InputCoverage(permuted), expected)
    input$se_NC <- input$se_NC[c(seq_len(nrow(input$se_NC)), 1L), , drop = FALSE]
    input$se_SSP <- input$se_SSP[c(seq_len(nrow(input$se_SSP)), 1L), , drop = FALSE]
    x <- SummarizedExperiment::assay(input$se_SSP)
    x[1L, ] <- c(NA, 0, -Inf)
    x[nrow(x), ] <- c(1, NA, Inf)
    SummarizedExperiment::assay(input$se_SSP) <- x
    actual <- expect_coverage_native(input)
    expect_identical(actual$genes_present, expected$genes_present)
    expect_identical(actual$rules_evaluable_min[-1L], actual$rules_total[-1L])
    expect_identical(InputCoverage(lapply(input, function(se) {
        se[rev(seq_len(nrow(se))), , drop = FALSE]
    })), actual)
})

test_that("unknown and ambiguous identifiers are not silently normalized", {
    input <- coverage_test_input(1L)
    rownames(input$se_NC) <- paste0("unknown_", seq_len(nrow(input$se_NC)))
    rownames(input$se_SSP) <- paste0("e", rownames(input$se_SSP))
    actual <- expect_coverage_native(input)
    expect_identical(actual$genes_present, rep(0L, 4L))
    input <- coverage_test_input(1L)
    ids <- rownames(input$se_SSP)
    ids[1L] <- paste(ids[1L], ids[2L], sep = ";")
    rownames(input$se_SSP) <- ids
    expect_coverage_native(input)
})

test_that("malformed or unavailable inputs fail explicitly", {
    input <- coverage_test_input()
    for (bad in list(NULL, matrix(1), list(), input["se_NC"],
        list(se_NC = NULL, se_SSP = input$se_SSP),
        list(se_NC = SummarizedExperiment::SummarizedExperiment(), se_SSP = input$se_SSP),
        list(se_NC = matrix(1), se_SSP = input$se_SSP))) {
        expect_error(InputCoverage(bad), "InputCoverage|se_NC|se_SSP")
    }
    duplicated_name <- input
    names(duplicated_name) <- c("se_NC", "se_NC")
    expect_error(InputCoverage(duplicated_name), "InputCoverage|se_NC|se_SSP")
    for (key in names(input)) {
        bad <- input
        SummarizedExperiment::assay(bad[[key]]) <- matrix("bad", nrow(bad[[key]]),
            ncol(bad[[key]]), dimnames = dimnames(bad[[key]]))
        expect_error(InputCoverage(bad), "numeric")
        for (id in c("", NA_character_)) {
            bad <- input
            rownames(bad[[key]])[1L] <- id
            expect_error(InputCoverage(bad), "gene identifiers")
        }
        bad <- input
        rownames(bad[[key]]) <- NULL
        expect_error(InputCoverage(bad), "gene identifiers")
    }
    bad <- lapply(input, function(se) se[, integer(), drop = FALSE])
    expect_error(InputCoverage(bad), "sample")
    bad <- input
    bad$se_SSP <- bad$se_SSP[, -1L, drop = FALSE]
    expect_error(InputCoverage(bad), "sample")
    bad <- input
    bad$se_SSP <- bad$se_SSP[, rev(seq_len(ncol(bad$se_SSP))), drop = FALSE]
    expect_error(InputCoverage(bad), "sample")
})

test_that("coverage does not call classifiers or repeat preprocessing", {
    input <- coverage_test_input()
    expected <- InputCoverage(input)
    forbidden <- function(...) stop("reporting called a scientific processing step")
    testthat::local_mocked_bindings(
        Mapping = forbidden, domapping = forbidden, impute_missing = forbidden,
        BS_AIMS = forbidden, BS_sspbc = forbidden, applyAIMS_AIMS = forbidden,
        applySSP = forbidden, .apply.nbc_AIMS = forbidden, .apply.nbc = forbidden,
        .package = "BreastSubtypeR")
    expect_identical(InputCoverage(input), expected)
})

test_that("reporting does not modify inputs, models, RNG or native predictions", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR")
    input <- suppressMessages(Mapping(OSLO2EMIT0obj$se_obj[, 1:3], verbose = FALSE))
    models <- coverage_test_models()
    input_before <- serialize(input, NULL)
    models_before <- serialize(models, NULL)
    predict_all <- function() suppressMessages(list(
        AIMS = BS_AIMS(input$se_SSP),
        ssp.subtype = BS_sspbc(input$se_SSP, ssp.name = "ssp.subtype"),
        ssp.pam50 = BS_sspbc(input$se_SSP, ssp.name = "ssp.pam50")))
    before <- predict_all()
    withr::local_preserve_seed()
    set.seed(161L)
    seed <- .Random.seed
    actual <- InputCoverage(input)
    expect_identical(.Random.seed, seed)
    expect_identical(serialize(input, NULL), input_before)
    expect_identical(serialize(models, NULL), models_before)
    expect_identical(predict_all(), before)
    for (i in seq_along(before)) {
        counts <- colSums(!is.na(before[[i]]$rules.matrix))
        expect_equal(actual$rules_evaluable_min[i + 1L], min(counts))
        expect_equal(actual$rules_evaluable_max[i + 1L], max(counts))
    }
    rm(".Random.seed", envir = .GlobalEnv)
    InputCoverage(input)
    expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("Mapping uses the same descriptive report without changing its return", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR")
    source <- OSLO2EMIT0obj$se_obj[, 1:2]
    quiet_messages <- testthat::capture_messages(quiet <- Mapping(source, verbose = FALSE))
    messages <- testthat::capture_messages(loud <- Mapping(source, verbose = TRUE))
    expect_identical(loud, quiet)
    expect_identical(names(quiet), c("se_NC", "se_SSP"))
    actual <- InputCoverage(quiet)
    expect_false(any(grepl("reference genes represented|model genes represented", quiet_messages)))
    for (i in seq_len(nrow(actual))) {
        expect_true(any(grepl(actual$model[i], messages, fixed = TRUE)), info = actual$model[i])
        expect_true(any(grepl(paste0(actual$genes_present[i], "/", actual$genes_total[i]),
            messages, fixed = TRUE)), info = actual$model[i])
    }
    expect_identical(sum(grepl("NC-based (PAM50)", messages, fixed = TRUE)), 1L)
    expect_false(any(grepl("adequate|inadequate|compatible|threshold", messages, ignore.case = TRUE)))
})
