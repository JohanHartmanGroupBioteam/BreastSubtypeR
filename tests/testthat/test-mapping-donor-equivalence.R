# Exact issue-219 regression evidence against the pre-optimization implementation.
# Conditions are compared intact, including their classes, calls and payloads.
mapping219_capture <- function(expr) {
    warnings <- messages <- list()
    value <- tryCatch(withCallingHandlers(force(expr),
        warning = function(w) {
            warnings[[length(warnings) + 1L]] <<- w
            invokeRestart("muffleWarning")
        }, message = function(m) {
            messages[[length(messages) + 1L]] <<- m
            invokeRestart("muffleMessage")
        }), error = identity)
    list(value = value, warnings = warnings, messages = messages)
}

mapping219_input <- function(x, ids) {
    rownames(x) <- if (nrow(x)) paste0("p", seq_len(nrow(x))) else NULL
    if (ncol(x)) colnames(x) <- paste0("sample ", rev(seq_len(ncol(x))))
    list(x = x, y = data.frame(probe = rownames(x), ENTREZID = ids))
}

mapping219_compare <- function(x, y, method, required) {
    before <- list(x = x, y = y)
    expected <- mapping219_capture(mapping219_baseline(x, y, method, required))
    actual <- mapping219_capture(BreastSubtypeR:::.mapping_donor_matrix(
        x, y, method, required))
    expect_identical(actual, expected)
    expect_identical(list(x = x, y = y), before)
    actual
}

test_that("donor preparation is identical for mixed groups and all five methods", {
    input <- mapping219_input(rbind(
        c(1, 2, 3, 4), c(1, NA, NaN, 4), rep(NA_real_, 4), rep(NaN, 4),
        c(Inf, 2, 3, 4), c(-Inf, 2, 3, 4),
        c(1, 2, 3, 4), c(4, 3, 2, 1), # Ties in total, SD and IQR.
        c(NA, NaN, NA, NaN), rep(NA_real_, 4),
        c(1, NA, 3, 4), c(11, 12, 13, 14),
        c(Inf, -Inf, 3, 4), c(0, 1, 2, 3),
        c(5, 6, 7, 8)),
        c("20", "3", "40", "5", "60", "7", "80", "80", "9", "9",
            "100", "100", "11", "11", NA_character_))
    for (method in c("max", "mean", "median", "stdev", "iqr")) {
        for (required in list(character(), c("3", "40", "5", "60", "7"),
            "9", unique(na.omit(input$y$ENTREZID)))) {
            for (ordering in 1:3) {
                x <- input$x
                y <- input$y
                if (ordering == 2L) {
                    x <- x[rev(seq_len(nrow(x))), c(4, 2, 1, 3), drop = FALSE]
                }
                if (ordering == 3L) y <- y[rev(seq_len(nrow(y))), , drop = FALSE]
                mapping219_compare(x, y, method, required)
            }
        }
    }
})

test_that("singleton and duplicate boundaries preserve storage and dimensions", {
    for (n in c(0L, 1L, 3L)) {
        for (duplicate in c(FALSE, TRUE)) {
            ids <- if (duplicate) c("20", "3", "3", "40") else c("20", "3", "4", "40")
            input <- mapping219_input(matrix(c(1, NA, 2, 3, 4, NaN, 5, 6,
                Inf, -Inf, 7, 8), 4L)[, seq_len(n), drop = FALSE], ids)
            for (method in c("max", "mean", "median", "stdev", "iqr")) {
                for (required in list(character(), ids)) {
                    mapping219_compare(input$x, input$y, method, required)
                }
            }
        }
    }
    for (storage in c("integer", "logical", "double")) {
        input <- mapping219_input(matrix(c(1, NA, 3, 4, 5, 6), 3L), c("10", "2", "30"))
        storage.mode(input$x) <- storage
        for (method in c("max", "mean", "median", "stdev", "iqr")) {
            mapping219_compare(input$x, input$y, method, "2")
        }
    }
    input <- mapping219_input(matrix(numeric(), 0L, 3L), character())
    for (method in c("max", "mean", "median", "stdev", "iqr")) {
        mapping219_compare(input$x, input$y, method, character())
    }
})

test_that("annotation and expression exceptions retain the original behavior", {
    base <- mapping219_input(matrix(seq_len(12), 4L), c("20", "3", "3", "40"))
    cases <- list(base)
    factor_ids <- numeric_ids <- empty_id <- duplicate_probes <- duplicate_rows <- absent <- base
    factor_ids$y$ENTREZID <- factor(factor_ids$y$ENTREZID, levels = c("40", "3", "20"))
    numeric_ids$y$ENTREZID <- as.integer(numeric_ids$y$ENTREZID)
    empty_id$y$ENTREZID[1] <- ""
    duplicate_probes$y$probe[2] <- duplicate_probes$y$probe[1]
    rownames(duplicate_rows$x)[2] <- rownames(duplicate_rows$x)[1]
    absent$y$probe[1] <- "absent probe"
    factor_probes <- integer_probes <- empty_probe <- classed <- named_axes <-
        duplicate_samples <- no_samples <- missing_sample <- empty_sample <- base
    factor_probes$y$probe <- factor(factor_probes$y$probe,
        levels = rev(factor_probes$y$probe))
    integer_probes$y$probe <- seq_len(nrow(integer_probes$y))
    empty_probe$y$probe[1] <- ""
    class(classed$x) <- c("mapping219_matrix", "matrix", "array")
    names(dimnames(named_axes$x)) <- c("features", "patients")
    colnames(duplicate_samples$x) <- c("same", "same", "other")
    colnames(no_samples$x) <- NULL
    colnames(missing_sample$x) <- c("first", NA_character_, "last")
    colnames(empty_sample$x) <- c("first", "", "last")
    cases <- c(cases, list(factor_ids, numeric_ids, empty_id, duplicate_probes, duplicate_rows,
        absent, factor_probes, integer_probes, empty_probe, classed, named_axes,
        duplicate_samples, no_samples, missing_sample, empty_sample))
    for (case in cases) {
        for (method in c("max", "mean", "median", "stdev", "iqr")) {
            for (required in list(character(), c("20", "3", "40"))) {
                mapping219_compare(case$x, case$y, method, required)
            }
        }
    }
    # Real subsetting failures must propagate for either protection policy.
    for (required in list(character(), "20")) {
        result <- mapping219_compare(absent$x, absent$y, "max", required)
        expect_s3_class(result$value, "error")
        expect_length(result$warnings, 0L)
    }
    for (method in list("unknown", NA_character_, character(), c("max", "mean"), 1)) {
        mapping219_compare(base$x, base$y, method, character())
    }
})

test_that("empty donor sets and extreme finite singletons are identical", {
    excluded <- mapping219_input(rbind(c(NA, NaN), c(Inf, -Inf)), c("20", "3"))
    missing_ids <- excluded
    missing_ids$y$ENTREZID[] <- NA_character_
    extremes <- mapping219_input(rbind(c(-1e308, 1e308),
        c(.Machine$double.xmin / 2, -.Machine$double.xmin / 2), c(-0, 0)),
        c("20", "3", "40"))
    for (method in c("max", "mean", "median", "stdev", "iqr")) {
        for (input in list(excluded, missing_ids)) {
            actual <- mapping219_compare(input$x, input$y, method, character())
            expect_identical(dim(actual$value$matrix), c(0L, 2L))
            expect_identical(typeof(actual$value$matrix), "double")
        }
        mapping219_compare(extremes$x, extremes$y, method, character())
    }
})

test_that("classed annotation metadata retains per-group warnings and errors", {
    input <- mapping219_input(matrix(seq_len(12L), 4L), c("20", "3", "3", "40"))
    input$y$metadata <- structure(seq_len(nrow(input$y)), class = "mapping219_metadata")
    previous <- utils::getS3method("[", "mapping219_metadata", optional = TRUE)
    check <- function(kind, method, required) {
        # Registration is local to this helper and restored even on test failure.
        registry <- get(".__S3MethodsTable__.", envir = asNamespace("base"))
        key <- "[.mapping219_metadata"
        registered <- exists(key, envir = registry, inherits = FALSE)
        if (registered) old <- get(key, envir = registry, inherits = FALSE)
        withr::defer(if (registered) assign(key, old, envir = registry) else
            rm(list = key, envir = registry))
        registerS3method("[", "mapping219_metadata", function(x, i, ...) {
            # The initial full-annotation filter is uneventful; only singleton
            # group subsetting signals the condition whose behavior we freeze.
            if (!missing(i) && length(i) == 1L) {
                condition <- structure(list(message = paste("metadata subset", kind),
                    call = sys.call(), rows = i),
                    class = c(paste0("mapping219_annotation_", kind), kind, "condition"))
                if (kind == "warning") warning(condition) else stop(condition)
            }
            structure(NextMethod("["), class = class(x))
        }, envir = baseenv())
        actual <- mapping219_compare(input$x, input$y, method, required)
        if (kind == "warning") {
            expect_length(actual$warnings, 2L)
            expect_s3_class(actual$warnings[[1L]], "mapping219_annotation_warning")
            expect_identical(actual$warnings[[1L]]$rows, 1L)
            expect_identical(actual$warnings[[2L]]$rows, 4L)
            expect_false(inherits(actual$value, "error"))
        } else {
            expect_length(actual$warnings, 0L)
            expect_s3_class(actual$value, "mapping219_annotation_error")
            expect_identical(actual$value$rows, 1L)
        }
    }
    for (kind in c("warning", "error")) {
        for (method in c("max", "mean", "median", "stdev", "iqr")) {
            for (required in list(character(), "20")) check(kind, method, required)
        }
    }
    expect_identical(utils::getS3method("[", "mapping219_metadata", optional = TRUE),
        previous)
})

test_that("prepared donors produce identical frozen estimates and diagnostics", {
    x <- rbind(c(0, NA, 2, NaN),
        t(vapply(seq_len(12L), function(i) c(i / 10, i, 2 + i / 10, i + 10), numeric(4))),
        c(1, NA, NA, NA), c(NA, NA, NA, NA), c(Inf, 1, 2, 3),
        c(2, 3, 4, 5), c(3, 2, 5, 4))
    input <- mapping219_input(x, c("target", as.character(900001:900015), "dup", "dup"))
    for (method in c("max", "mean", "median", "stdev", "iqr")) {
        expected <- mapping219_capture(mapping219_baseline(input$x, input$y, method, "target"))
        actual <- mapping219_compare(input$x, input$y, method, "target")
        frozen <- BreastSubtypeR:::knn_vectorized_reference(expected$value$matrix, "target")
        observed <- BreastSubtypeR:::knn_vectorized_reference(actual$value$matrix, "target")
        expect_identical(observed, frozen)
        accepted <- BreastSubtypeR:::.mapping_impute_log(actual$value$matrix, "target",
            exclusions = actual$value$exclusions)
        expect_identical(accepted$diagnostics, frozen)
        expect_identical(accepted$matrix, frozen$matrix)
        expect_identical(accepted$exclusions, expected$value$exclusions)
        for (target in list(rep(NA_real_, 4), c(Inf, NA, 2, 3))) {
            changed <- input$x
            changed[1, ] <- target
            reference <- mapping219_capture(mapping219_baseline(
                changed, input$y, method, "target"))$value
            prepared <- mapping219_compare(changed, input$y, method, "target")$value
            expect_identical(BreastSubtypeR:::knn_vectorized_reference(
                prepared$matrix, "target"), BreastSubtypeR:::knn_vectorized_reference(
                reference$matrix, "target"))
            expected_error <- mapping219_capture(BreastSubtypeR:::.mapping_impute_log(
                reference$matrix, "target", exclusions = reference$exclusions))
            actual_error <- mapping219_capture(BreastSubtypeR:::.mapping_impute_log(
                prepared$matrix, "target", exclusions = prepared$exclusions))
            expect_identical(actual_error, expected_error)
            expect_s3_class(actual_error$value, "BreastSubtypeR_knn_unsupported")
        }
    }
})

mapping219_se <- function(x, ids) {
    input <- mapping219_input(x, ids)
    SummarizedExperiment::SummarizedExperiment(assays = list(counts = input$x),
        rowData = S4Vectors::DataFrame(input$y, Length = rep(1000, nrow(x))),
        colData = S4Vectors::DataFrame(PatientID = colnames(input$x),
            row.names = colnames(input$x)))
}

mapping219_mapping <- function(donor, se, ...) {
    local_mocked_bindings(.mapping_donor_matrix = donor, .package = "BreastSubtypeR")
    mapping219_capture(Mapping(se, verbose = FALSE, ...))
}

test_that("public Mapping including opt-in SSP-ROR is exactly unchanged", {
    models <- new.env(parent = emptyenv())
    data("BreastSubtypeRobj", package = "BreastSubtypeR", envir = models)
    ref <- models$BreastSubtypeRobj
    shared <- intersect(as.character(ref$genes.sig50$EntrezGene.ID),
        as.character(ref$genes.signature$EntrezGene.ID[
            ref$genes.signature$SSP_based == "Yes"]))[[1L]]
    ror <- setdiff(BreastSubtypeR:::.mapping_ssp_ror_ids(),
        as.character(ref$genes.signature$EntrezGene.ID))[[1L]]
    x <- rbind(c(0, NA, 2, 3), c(1, 2, NaN, 4),
        t(vapply(seq_len(12L), function(i) c(i / 10, i, 2 + i / 10, 4), numeric(4))),
        rep(NA_real_, 4), c(Inf, 2, 3, 4), c(1, 2, 3, 4), c(4, 3, 2, 1))
    se <- mapping219_se(x, c(shared, ror, as.character(900001:900014), "900015", "900015"))
    candidate <- BreastSubtypeR:::.mapping_donor_matrix
    original <- se
    for (method in c("max", "mean", "median", "stdev", "iqr")) {
        for (include in c(FALSE, TRUE)) {
            for (impute in c(FALSE, TRUE)) {
                expected <- mapping219_mapping(mapping219_baseline, se,
                    method = method, impute = impute, include_ssp_ror = include)
                observed <- mapping219_mapping(candidate, se,
                    method = method, impute = impute, include_ssp_ror = include)
                expect_identical(observed, expected)
                expect_false(inherits(observed$value, "error"))
            }
        }
    }
    expect_identical(se, original)
    reordered <- se[rev(seq_len(nrow(se))), c(4L, 2L, 1L, 3L)]
    expect_identical(mapping219_mapping(candidate, reordered, impute = TRUE,
        include_ssp_ror = TRUE), mapping219_mapping(mapping219_baseline, reordered,
        impute = TRUE, include_ssp_ror = TRUE))
    # Complete targets still bypass unusable background groups.
    complete <- se
    SummarizedExperiment::assay(complete)[1:2, ] <- rbind(1:4, 2:5)
    for (include in c(FALSE, TRUE)) {
        expect_identical(mapping219_mapping(candidate, complete, impute = TRUE,
            include_ssp_ror = include), mapping219_mapping(mapping219_baseline,
            complete, impute = TRUE, include_ssp_ror = include))
        expect_identical(mapping219_mapping(candidate, complete, impute = TRUE,
            include_ssp_ror = include), mapping219_mapping(candidate, complete,
            include_ssp_ror = include))
    }
    expect_identical(mapping219_mapping(candidate, se),
        mapping219_mapping(candidate, se, impute = FALSE))
    # Genuine unsupported KNN: fewer than ten donors, with structured diagnostics.
    small <- se[1:4, ]
    expected <- mapping219_mapping(mapping219_baseline, small, impute = TRUE,
        include_ssp_ror = TRUE)
    actual <- mapping219_mapping(candidate, small, impute = TRUE, include_ssp_ror = TRUE)
    expect_identical(actual, expected)
    expect_s3_class(actual$value, "BreastSubtypeR_knn_unsupported")
})

test_that("donor preparation and public imputation preserve RNG state and absence", {
    withr::local_preserve_seed()
    input <- mapping219_input(rbind(c(0, NA, 2),
        t(vapply(seq_len(12L), function(i) c(i / 10, i, 2 + i / 10), numeric(3)))),
        c("596", as.character(900001:900012)))
    se <- mapping219_se(input$x, input$y$ENTREZID)
    candidate <- BreastSubtypeR:::.mapping_donor_matrix
    for (donor in list(mapping219_baseline, candidate)) {
        for (method in c("max", "mean", "median", "stdev", "iqr")) {
            set.seed(219)
            seed <- .Random.seed
            mapping219_capture(donor(input$x, input$y, method, "596"))
            expect_identical(.Random.seed, seed)
            mapped <- mapping219_mapping(donor, se, method = method,
                impute = TRUE, include_ssp_ror = TRUE)
            expect_false(inherits(mapped$value, "error"))
            expect_identical(.Random.seed, seed)
            rm(".Random.seed", envir = .GlobalEnv)
            mapping219_capture(donor(input$x, input$y, method, "596"))
            expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            mapped <- mapping219_mapping(donor, se, method = method,
                impute = TRUE, include_ssp_ror = TRUE)
            expect_false(inherits(mapped$value, "error"))
            expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        }
    }
})
