# New B6 integration evidence. These are not the missing historical baseline
# or replacements for any frozen B2/B3/B4 estimator expectations.

b6_reference <- function() {
    env <- new.env(parent = emptyenv())
    data("BreastSubtypeRobj", package = "BreastSubtypeR", envir = env)
    env$BreastSubtypeRobj
}

b6_target <- function() {
    ref <- b6_reference()
    intersect(as.character(ref$genes.sig50$EntrezGene.ID),
        as.character(ref$genes.signature$EntrezGene.ID[
            ref$genes.signature$SSP_based == "Yes"]))[[1L]]
}

b6_se <- function(x, ids = rownames(x)) {
    stopifnot(is.matrix(x), length(ids) == nrow(x))
    rownames(x) <- paste0("b6_probe_", seq_len(nrow(x)))
    if (is.null(colnames(x))) colnames(x) <- paste0("synthetic sample ", seq_len(ncol(x)))
    SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = x),
        rowData = S4Vectors::DataFrame(probe = rownames(x), ENTREZID = ids,
            Length = rep(1000, nrow(x))),
        colData = S4Vectors::DataFrame(PatientID = colnames(x)))
}

b6_matrix <- function(target = c(0, NA_real_, 2, NA_real_)) {
    x <- rbind(target, t(vapply(seq_len(10L), function(i)
        c(i / 10, i, 2 + i / 10, i + 10), numeric(4))),
        background = c(100, NaN, 100, NA_real_))
    rownames(x) <- c("target", paste0("d", seq_len(10L)), "background")
    colnames(x) <- c("sample 3", "s/1", "S.2", "sample 4")
    x
}

b6_save_interface <- function(case, x, targets, observed, expected,
    exclusions = NULL) {
    directory <- Sys.getenv("B6_EVIDENCE_DIR", unset = "")
    if (!nzchar(directory)) return(invisible(NULL))
    stopifnot(dir.exists(directory))
    record <- list(label = "NEW B6 synthetic integration interface evidence",
        base = "f65b1bff312c5ccf453f248ab8d1f4dc336c5057",
        estimator_sha256 = "18293960a4fdd625c473f9bed587830f5f1027b08d445e93cfa502f7263a160c",
        case = case, donor_matrix = x, target_ids = targets,
        integrated_diagnostics = observed, frozen_B3_diagnostics = expected,
        diagnostics_identical = identical(observed, expected),
        exclusions = exclusions)
    # Every execution gets a distinct file; prior validation attempts survive.
    saveRDS(record, tempfile(paste0(case, "-"), tmpdir = directory, fileext = ".rds"))
    invisible(NULL)
}

b6_adapter_equal <- function(x, targets, case, ssp = character()) {
    original <- x
    expected <- BreastSubtypeR:::knn_vectorized_reference(x, targets)
    actual <- BreastSubtypeR:::.mapping_impute_log(x, targets, ssp_targets = ssp)
    expect_identical(actual$diagnostics, expected)
    expect_identical(actual$matrix, expected$matrix)
    expect_identical(x, original)
    b6_save_interface(case, x, targets, actual$diagnostics, expected)
    actual
}

b6_adapter_error <- function(x, targets, case, status) {
    original <- x
    expected <- BreastSubtypeR:::knn_vectorized_reference(x, targets)
    error <- expect_error(BreastSubtypeR:::.mapping_impute_log(x, targets),
        class = "BreastSubtypeR_knn_unsupported")
    expect_identical(error$diagnostics, expected)
    expect_identical(error$diagnostics$targets[[targets[[1L]]]]$status, status)
    expect_true(targets[[1L]] %in% error$targets)
    expect_match(conditionMessage(error), targets[[1L]], fixed = TRUE)
    expect_identical(x, original)
    b6_save_interface(case, x, targets, error$diagnostics, expected, error$exclusions)
    error
}

test_that("B6 defaults and complete-output paths bypass both preparation and estimator", {
    id <- b6_target()
    x <- rbind(c(1, 2, 3), c(1, NA, 3), c(NA, 2, 3))
    se <- b6_se(x, c(id, "900001", "900001"))
    local_mocked_bindings(
        .mapping_donor_matrix = function(...) stop("unexpected donor preparation"),
        knn_vectorized_reference = function(...) stop("unexpected estimator"),
        .package = "BreastSubtypeR")
    expect_identical(formals(Mapping)$impute, FALSE)
    expect_identical(formals(BreastSubtypeR:::domapping)$impute, FALSE)
    expected <- Mapping(se, impute = FALSE, verbose = FALSE)
    expect_identical(Mapping(se, impute = TRUE, verbose = FALSE), expected)
    expect_identical(Mapping(se, verbose = FALSE), expected)
    SummarizedExperiment::assay(se)[1, 2] <- NA_real_
    observed <- Mapping(se, verbose = FALSE)
    expect_identical(observed, Mapping(se, impute = FALSE, verbose = FALSE))
    expect_true(anyNA(SummarizedExperiment::assay(observed$se_NC)))
    expect_false("knn_vectorized_reference" %in% getNamespaceExports("BreastSubtypeR"))
})

test_that("B6 duplicate collapse resolves holes before any optional preparation", {
    id <- b6_target()
    local_mocked_bindings(
        .mapping_donor_matrix = function(...) stop("unexpected donor preparation"),
        knn_vectorized_reference = function(...) stop("unexpected estimator"),
        .package = "BreastSubtypeR")
    for (method in c("mean", "median")) {
        se <- b6_se(rbind(c(1, NA, 3), c(NA, 2, NA)), rep(id, 2))
        actual <- Mapping(se, method = method, impute = TRUE, verbose = FALSE)
        expect_identical(unname(as.numeric(as.matrix(SummarizedExperiment::assay(actual$se_NC)))),
            c(1, 2, 3))
    }
    for (method in c("max", "iqr", "stdev")) {
        se <- b6_se(rbind(c(1, NA, 1), c(2, 4, 6)), rep(id, 2))
        actual <- Mapping(se, method = method, impute = TRUE, verbose = FALSE)
        expect_identical(unname(as.numeric(as.matrix(SummarizedExperiment::assay(actual$se_NC)))),
            c(2, 4, 6))
    }
    se <- b6_se(rbind(c(1, NA, 3), c(NA, 2, 3)), rep(id, 2))
    expect_error(Mapping(se, method = "max", impute = TRUE, verbose = FALSE),
        paste0("Cannot resolve Entrez ID ", id), fixed = TRUE)
})

test_that("B6 actual output targets are ordered represented IDs rather than all protected IDs", {
    reference <- list(genes.signature = data.frame(EntrezGene.ID = c("1", "2", "3"),
        SSP_based = c("Yes", "No", "No")),
        genes.sig50 = data.frame(EntrezGene.ID = "3"))
    x <- matrix(1, 4, 2, dimnames = list(c("2", "3", "4", "1"), c("s2", "s1")))
    expect_identical(BreastSubtypeR:::.mapping_output_targets(x, reference), c("3", "1"))
    expect_identical(BreastSubtypeR:::.mapping_output_targets(x[1:2, , drop = FALSE], reference), "3")
    actual_reference <- b6_reference()
    protected <- BreastSubtypeR:::.mapping_required_ids(actual_reference)
    env <- new.env(parent = emptyenv())
    data(list = c("AIMSmodel", "sspbc.models"), package = "BreastSubtypeR", envir = env)
    expected <- unique(c(as.character(actual_reference$genes.signature$EntrezGene.ID),
        as.character(actual_reference$genes.sig50$EntrezGene.ID),
        sub("^e", "", unlist(strsplit(c(env$AIMSmodel$all.pairs,
            env$sspbc.models$ssp.pam50$all.pairs,
            env$sspbc.models$ssp.subtype$all.pairs), "<", fixed = TRUE), use.names = FALSE))))
    expect_identical(protected, expected)
})

test_that("B6 Mapping supplies full-input donors and only represented output targets", {
    id <- b6_target()
    x <- b6_matrix()
    rownames(x) <- c(id, as.character(900001:900011))
    x[11, 2] <- NA_real_ # Ten selected neighbors; nine contributors at this hole.
    se <- b6_se(x)
    original <- se
    frozen <- BreastSubtypeR:::knn_vectorized_reference
    observed <- NULL
    local_mocked_bindings(knn_vectorized_reference = function(x, targets) {
        result <- frozen(x, targets)
        observed <<- list(x = x, targets = targets, result = result)
        result
    }, .package = "BreastSubtypeR")
    actual <- Mapping(se, method = "max", impute = TRUE, verbose = FALSE)
    expect_identical(observed$x, x)
    expect_identical(observed$targets, id)
    direct <- frozen(observed$x, observed$targets)
    expect_identical(observed$result, direct)
    expect_identical(observed$result$targets[[id]]$neighbors$gene, as.character(900001:900010))
    expect_identical(observed$result$targets[[id]]$cells$contributors, c(9L, 10L))
    expect_identical(unname(as.numeric(as.matrix(SummarizedExperiment::assay(actual$se_NC)))),
        unname(direct$matrix[id, ]))
    expect_identical(unname(as.numeric(as.matrix(SummarizedExperiment::assay(actual$se_SSP)))),
        unname(2^direct$matrix[id, ]))
    expect_identical(names(actual), c("se_NC", "se_SSP"))
    expect_identical(nrow(actual$se_NC), 1L)
    expect_identical(nrow(actual$se_SSP), 1L)
    expect_identical(se, original)
    expect_identical(colnames(actual$se_NC), colnames(x))
    expect_identical(colnames(actual$se_SSP), colnames(x))
    b6_save_interface("mapping-full-universe-nine-contributors", observed$x,
        observed$targets, observed$result, direct)
})

test_that("B6 retains useful high-missingness rows and removes all legacy support cutoffs", {
    id <- b6_target()
    x <- rbind(matrix(c(0, NA, NA, NA), 10, 4, byrow = TRUE), c(0, 1, 2, 3))
    rownames(x) <- c(id, as.character(900001:900010))
    colnames(x) <- c("s4", "s3", "s2", "s1")
    y <- data.frame(probe = rownames(x), ENTREZID = rownames(x))
    prepared <- BreastSubtypeR:::.mapping_donor_matrix(x, y, "max", id)
    expect_identical(prepared$matrix, x)
    expect_identical(nrow(prepared$exclusions), 0L)
    expect_true(mean(is.na(x[1, ])) > 0.5)
    expect_true(all(colMeans(is.na(x))[-1] > 0.8))
    actual <- b6_adapter_equal(prepared$matrix, id, "retired-cutoffs-one-contributor", id)
    expect_identical(actual$diagnostics$targets[[id]]$cells$contributors, rep(1L, 3))
    expect_identical(unname(actual$matrix[id, ]), c(0, 1, 2, 3))
    mapped <- Mapping(b6_se(x), impute = TRUE, verbose = FALSE)
    expect_identical(unname(as.numeric(as.matrix(SummarizedExperiment::assay(mapped$se_NC)))), c(0, 1, 2, 3))
})

test_that("B6 donor exclusions are explicit and cannot hide protected or unrelated errors", {
    x <- rbind(target = c(0, NA, 2), good = c(1, 2, 3),
        mostly = c(1, NA, NA), missing = rep(NA_real_, 3),
        infinite = c(Inf, 2, 3), a = c(1, NA, 3), b = c(NA, 2, 3))
    colnames(x) <- c("s3", "s1", "s2")
    y <- data.frame(probe = rownames(x), ENTREZID = c("1", "2", "3", "4", "5", "6", "6"))
    expect_warning(prepared <- BreastSubtypeR:::.mapping_donor_matrix(x, y, "max", "1"),
        class = "BreastSubtypeR_knn_donor_exclusion")
    expect_identical(rownames(prepared$matrix), c("1", "2", "3"))
    expect_identical(prepared$exclusions$gene, c("4", "5", "6"))
    expect_identical(prepared$exclusions$reason,
        c("entirely_missing", "input_infinity", "duplicate_resolution_failed"))
    expect_identical(unname(prepared$matrix["3", ]), c(1, NA_real_, NA_real_))
    expect_error(BreastSubtypeR:::.mapping_donor_matrix(x, y, "max", c("1", "6")),
        "Cannot resolve Entrez ID 6", fixed = TRUE)
    protected <- suppressWarnings(BreastSubtypeR:::.mapping_donor_matrix(x, y, "max", c("1", "4", "5")))
    expect_true(all(c("4", "5") %in% rownames(protected$matrix)))
    expect_true(all(is.na(protected$matrix["4", ])))
    expect_identical(unname(protected$matrix["5", 1]), Inf)
    local_mocked_bindings(duplicate_genes = function(...) stop("unexpected internal failure"),
        .package = "BreastSubtypeR")
    expect_error(BreastSubtypeR:::.mapping_donor_matrix(x, y, "max", character()),
        "unexpected internal failure", fixed = TRUE)
})

test_that("B6 protected nontarget missingness stays distinct from target failure", {
    x <- b6_matrix()
    x["background", ] <- c(NA, NaN, NA, NaN)
    y <- data.frame(probe = rownames(x), ENTREZID = rownames(x))
    prepared <- BreastSubtypeR:::.mapping_donor_matrix(x, y, "max", c("target", "background"))
    expect_identical(prepared$matrix, x)
    actual <- b6_adapter_equal(prepared$matrix, "target", "protected-nontarget-missing")
    expect_identical(actual$matrix["background", ], x["background", ])
    expect_identical(names(actual$diagnostics$targets), "target")
    expect_identical(actual$diagnostics$targets$target$candidates$reason[12], "no_overlap")
})

test_that("B6 accepts genuine zero estimates with observed selected contributors", {
    x <- b6_matrix()
    x[2:11, 2] <- c(-5:-1, 1:5)
    x[2:11, 4] <- 0
    actual <- b6_adapter_equal(x, "target", "genuine-zero-estimates", "target")
    expect_identical(unname(actual$matrix["target", c(2, 4)]), c(0, 0))
    expect_identical(actual$diagnostics$targets$target$cells$contributors, c(10L, 10L))
    expect_identical(actual$diagnostics$targets$target$cells$reason, rep("estimated", 2))
})

test_that("B6 unresolved supported, unsupported and partial targets are structured failures", {
    x <- b6_matrix()
    x["target", ] <- c(NA, NaN, NA, NaN)
    b6_adapter_error(x, "target", "entirely-missing-target", "entirely_missing")
    b6_adapter_error(b6_matrix()[1:10, , drop = FALSE], "target", "insufficient-neighbors", "insufficient_neighbors")
    x <- b6_matrix()
    x[2:11, c(2, 4)] <- NA_real_
    error <- b6_adapter_error(x, "target", "no-selected-contributors", "unresolved_coordinates")
    expect_identical(error$diagnostics$targets$target$cells$contributors, c(0L, 0L))
    expect_identical(error$diagnostics$targets$target$cells$reason, rep("no_contributors", 2))
    x <- b6_matrix()
    x[2:11, 4] <- NA_real_
    error <- b6_adapter_error(x, "target", "partially-resolved-target", "partial")
    expect_identical(error$diagnostics$targets$target$cells$resolved, c(TRUE, FALSE))
    expect_true(is.finite(error$diagnostics$matrix["target", 2]))
    expect_true(is.na(error$diagnostics$matrix["target", 4]))
    disjoint <- matrix(NA_real_, 11, 2,
        dimnames = list(c("target", paste0("d", 1:10)), c("s1", "s2")))
    disjoint[1, 1] <- 0
    disjoint[-1, 2] <- 1
    error <- b6_adapter_error(disjoint, "target", "disjoint-observations", "insufficient_neighbors")
    expect_true(all(error$diagnostics$targets$target$candidates$reason[-1] == "no_overlap"))
    x <- rbind(c(0, NA_real_), matrix(rep(c(0, .Machine$double.xmax), each = 10), 10, 2))
    dimnames(x) <- list(c("target", paste0("d", 1:10)), c("s1", "s2"))
    error <- b6_adapter_error(x, "target", "nonfinite-estimate", "unresolved_coordinates")
    expect_identical(error$diagnostics$targets$target$cells$reason, "nonfinite_estimate")
})

test_that("B6 public Mapping rejects partial or entirely missing represented targets", {
    id <- b6_target()
    for (entire in c(FALSE, TRUE)) {
        x <- b6_matrix()
        rownames(x) <- c(id, as.character(900001:900011))
        if (entire) x[1, ] <- c(NA, NaN, NA, NaN) else x[2:11, 4] <- NA_real_
        se <- b6_se(x)
        original <- se
        error <- expect_error(Mapping(se, impute = TRUE, verbose = FALSE),
            class = "BreastSubtypeR_knn_unsupported")
        expect_identical(error$targets, id)
        expect_identical(error$diagnostics$targets[[id]]$status,
            if (entire) "entirely_missing" else "partial")
        expect_identical(se, original)
        expect_null(error$se_NC)
        expect_null(error$se_SSP)
    }
})

test_that("B6 public Mapping excludes bad background but never protected infinities", {
    id <- b6_target()
    x <- b6_matrix()
    rownames(x) <- c(id, as.character(900001:900011))
    x[12, 1] <- Inf
    se <- b6_se(x)
    expect_warning(mapped <- Mapping(se, impute = TRUE, verbose = FALSE),
        class = "BreastSubtypeR_knn_donor_exclusion")
    expect_false(anyNA(SummarizedExperiment::assay(mapped$se_NC)))
    # A second real output gene is protected even if it has no missing cells.
    ref <- b6_reference()
    second <- setdiff(as.character(ref$genes.sig50$EntrezGene.ID), id)[[1L]]
    x[12, ] <- c(Inf, 2, 100, 4)
    ids <- c(id, as.character(900001:900010), second)
    se <- b6_se(x, ids)
    error <- expect_error(Mapping(se, impute = TRUE, verbose = FALSE),
        class = "BreastSubtypeR_knn_unsupported")
    expect_true(id %in% error$targets)
    expect_identical(error$diagnostics$input_counts[["posInf"]], 1L)
    expect_identical(error$diagnostics$targets[[id]]$status, "invalid_input_infinity")
})

test_that("B6 unsupported Mapping retains the exact prepared interface and exclusions", {
    id <- b6_target()
    x <- b6_matrix()
    rownames(x) <- c(id, as.character(900001:900011))
    x[2:11, 4] <- NA_real_
    x[12, 1] <- Inf
    observed <- NULL
    frozen <- BreastSubtypeR:::knn_vectorized_reference
    local_mocked_bindings(knn_vectorized_reference = function(x, targets) {
        observed <<- list(x = x, targets = targets)
        frozen(x, targets)
    }, .package = "BreastSubtypeR")
    expect_warning(error <- expect_error(Mapping(b6_se(x), impute = TRUE, verbose = FALSE),
        class = "BreastSubtypeR_knn_unsupported"), class = "BreastSubtypeR_knn_donor_exclusion")
    expect_identical(observed$x, x[-12, , drop = FALSE])
    expect_identical(observed$targets, id)
    expect_identical(error$exclusions$gene, "900011")
    expect_identical(error$exclusions$reason, "input_infinity")
    expected <- frozen(observed$x, observed$targets)
    expect_identical(error$diagnostics, expected)
    b6_save_interface("mapping-partial-with-exclusion", observed$x, observed$targets,
        error$diagnostics, expected, error$exclusions)
})

test_that("B6 NA and NaN stay separately diagnosed and absent rows are not created", {
    x <- b6_matrix()
    x[1, 2] <- NaN
    actual <- b6_adapter_equal(x, "target", "NA-NaN-preservation")
    expect_identical(actual$diagnostics$input_counts[c("NA_only", "NaN")], c(NA_only = 2L, `NaN` = 2L))
    expect_identical(actual$matrix["background", ], x["background", ])
    expect_identical(rownames(actual$matrix), rownames(x))
    expect_error(BreastSubtypeR:::.mapping_impute_log(x, c("target", "absent")),
        class = "BreastSubtypeR_knn_boundary")
    # Public Mapping already established this target set from represented IDs.
    id <- b6_target()
    rownames(x) <- c(id, as.character(900001:900011))
    mapped <- Mapping(b6_se(x), impute = TRUE, verbose = FALSE)
    expect_identical(nrow(mapped$se_NC), 1L)
    expect_identical(nrow(mapped$se_SSP), 1L)
})

test_that("B6 retained infinities block estimation and distance overflow remains candidate-local", {
    x <- b6_matrix()
    x["background", 1] <- Inf
    error <- b6_adapter_error(x, "target", "protected-positive-infinity", "invalid_input_infinity")
    expect_identical(error$diagnostics$input_counts[["posInf"]], 1L)
    x["background", 1] <- -Inf
    error <- b6_adapter_error(x, "target", "protected-negative-infinity", "invalid_input_infinity")
    expect_identical(error$diagnostics$input_counts[["negInf"]], 1L)
    x <- b6_matrix()
    x["background", c(1, 3)] <- .Machine$double.xmax
    actual <- b6_adapter_equal(x, "target", "nonfinite-distance-excluded")
    expect_identical(actual$diagnostics$targets$target$candidates$reason[12], "nonfinite_distance")
    expect_identical(actual$diagnostics$targets$target$status, "imputed")
})

test_that("B6 boundary rejects changed dimensions, names, observed and nontarget values", {
    x <- b6_matrix()
    original <- BreastSubtypeR:::knn_vectorized_reference
    corruptions <- list(
        dimensions = function(z) {z$matrix <- z$matrix[, -1, drop = FALSE]; z},
        names = function(z) {colnames(z$matrix) <- rev(colnames(z$matrix)); z},
        observed = function(z) {z$matrix[1, 1] <- 999; z},
        background_hole = function(z) {z$matrix[12, 2] <- 999; z},
        target_order = function(z) {names(z$targets) <- "other"; z})
    for (corrupt in corruptions) {
        local_mocked_bindings(knn_vectorized_reference = function(x, targets)
            corrupt(original(x, targets)), .package = "BreastSubtypeR")
        expect_error(BreastSubtypeR:::.mapping_impute_log(x, "target"),
            class = "BreastSubtypeR_knn_boundary")
    }
})

test_that("B6 target order and exact ties follow supplied row and column ordering", {
    x <- b6_matrix()
    x[2:11, c(1, 3)] <- matrix(rep(c(0, 2), each = 10), 10, 2)
    x <- x[c(1, 11:2, 12), c(3, 2, 1, 4)]
    actual <- b6_adapter_equal(x, c("d3", "target"), "reordered-inputs-exact-ties")
    expect_identical(names(actual$diagnostics$targets), c("d3", "target"))
    expect_identical(actual$diagnostics$targets$target$neighbors$gene, paste0("d", 10:1))
    expect_identical(dimnames(actual$matrix), dimnames(x))
})

test_that("B6 success and failure preserve RNG kind and existing or absent seed", {
    withr::local_preserve_seed()
    original_kind <- RNGkind()
    for (seed_present in c(TRUE, FALSE)) {
        if (seed_present) set.seed(62026) else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv)
        before <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        BreastSubtypeR:::.mapping_impute_log(b6_matrix(), "target")
        expect_error(BreastSubtypeR:::.mapping_impute_log(b6_matrix()[1:10, ], "target"),
            class = "BreastSubtypeR_knn_unsupported")
        expect_identical(RNGkind(), original_kind)
        expect_identical(exists(".Random.seed", .GlobalEnv, inherits = FALSE), seed_present)
        expect_identical(get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE), before)
    }
})

test_that("B6 raw counts and missing raw-derived expression can never enter KNN", {
    id <- b6_target()
    local_mocked_bindings(knn_vectorized_reference = function(...) stop("raw estimator invoked"),
        .mapping_donor_matrix = function(...) stop("raw donor preparation invoked"),
        .package = "BreastSubtypeR")
    for (hole in list(NA_real_, NaN)) {
        se <- b6_se(matrix(c(100, hole, 300), 1), id)
        for (setting in c(FALSE, TRUE))
            expect_error(Mapping(se, RawCounts = TRUE, impute = setting, verbose = FALSE),
                "raw-count imputation is not supported", fixed = TRUE)
    }
    se <- b6_se(matrix(c(100, 200, 300), 1), id)
    SummarizedExperiment::rowData(se)$Length <- NA_real_
    expect_error(Mapping(se, RawCounts = TRUE, impute = TRUE, verbose = FALSE),
        "neither counts nor linear FPKM will be imputed", fixed = TRUE)
    actual <- Mapping(se, RawCounts = TRUE, impute = FALSE, verbose = FALSE)
    expect_true(all(is.na(SummarizedExperiment::assay(actual$se_SSP))))
})

test_that("B6 new finite log estimates must also have finite existing SSP projection", {
    x <- rbind(c(0, NA), matrix(rep(c(0, 2048), each = 10), 10, 2))
    dimnames(x) <- list(c("target", paste0("d", 1:10)), c("observed", "missing"))
    expected <- BreastSubtypeR:::knn_vectorized_reference(x, "target")
    expect_identical(expected$matrix[1, 2], 2048)
    error <- expect_error(BreastSubtypeR:::.mapping_impute_log(x, "target", ssp_targets = "target"),
        class = "BreastSubtypeR_knn_unsupported")
    expect_identical(error$reasons, c(target = "nonfinite_ssp_projection"))
    expect_identical(error$diagnostics, expected)
    b6_save_interface("SSP-projection-overflow", x, "target", error$diagnostics, expected)
    # NC-only private targets keep their finite log scale; no clipping occurs.
    expect_identical(BreastSubtypeR:::.mapping_impute_log(x, "target")$matrix, expected$matrix)
    id <- b6_target()
    rownames(x) <- c(id, as.character(900001:900010))
    expect_error(Mapping(b6_se(x), impute = TRUE, verbose = FALSE),
        class = "BreastSubtypeR_knn_unsupported")
})

test_that("B6 integer capacity arithmetic is checked without large allocations", {
    limit <- as.double(.Machine$integer.max)
    zero <- c(NA_only = 0, `NaN` = 0, posInf = 0, negInf = 0)
    check <- BreastSubtypeR:::.mapping_knn_capacity
    exact <- check(limit, limit, 1, 1, zero)
    expect_identical(exact[["donor_pairs"]], limit - 1)
    expect_identical(check(46341, 10, 46341, 46341, zero)[["donor_pairs"]], 46341 * 46340)
    expect_error(check(46342, 10, 46342, 46342, zero), class = "BreastSubtypeR_knn_capacity")
    expect_error(check(limit + 1, 1, 1, 1, zero), class = "BreastSubtypeR_knn_capacity")
    expect_error(check(1, limit + 1, 1, 1, zero), class = "BreastSubtypeR_knn_capacity")
    expect_error(check(10, 10, 11, 1, zero), class = "BreastSubtypeR_knn_capacity")
    expect_error(check(10, 10, 1, 2, zero), class = "BreastSubtypeR_knn_capacity")
    counts <- zero; counts[["NA_only"]] <- limit + 1
    expect_error(check(10, 10, 1, 1, counts), class = "BreastSubtypeR_knn_capacity")
    counts <- zero; counts[["posInf"]] <- limit; counts[["negInf"]] <- 1
    expect_error(check(limit, 1, 1, 0, counts), class = "BreastSubtypeR_knn_capacity")
    counts[["negInf"]] <- 0
    expect_identical(check(limit, 1, 1, 0, counts)[["infinity_count"]], limit)
    for (bad in list(NA_real_, Inf, -1, 1.5))
        expect_error(check(bad, 1, 0, 0, zero), class = "BreastSubtypeR_knn_capacity")
    # Full/entirely-missing targets do not count as actual distance computations.
    x <- b6_matrix()
    x["background", ] <- NA_real_
    recorded <- NULL
    original <- BreastSubtypeR:::.mapping_knn_capacity
    local_mocked_bindings(.mapping_knn_capacity = function(n_genes, n_samples,
        n_targets, n_computed, input_counts) {
        recorded <<- c(n_targets = n_targets, n_computed = n_computed)
        original(n_genes, n_samples, n_targets, n_computed, input_counts)
    }, .package = "BreastSubtypeR")
    expect_error(BreastSubtypeR:::.mapping_impute_log(x, c("d1", "target", "background")),
        class = "BreastSubtypeR_knn_unsupported")
    expect_identical(recorded, c(n_targets = 3, n_computed = 1))
})

# Consolidated corrections: new evidence, preserving all original B6 assertions.
test_that("B6 corrections reject malformed diagnostics and finite cell disagreements", {
    x <- b6_matrix()
    frozen <- BreastSubtypeR:::knn_vectorized_reference
    corruptions <- list(
        no_targets = function(z) {z$targets <- NULL; z},
        no_work = function(z) {z$work <- NULL; z},
        bad_counts = function(z) {z$input_counts <- as.double(z$input_counts); z},
        bad_positions = function(z) {z$input_positions$NA_only <- NULL; z},
        no_record = function(z) {z$targets$target <- NULL; z},
        no_status = function(z) {z$targets$target$status <- NULL; z},
        status_length = function(z) {z$targets$target$status <- c("imputed", "partial"); z},
        status_na = function(z) {z$targets$target$status <- NA_character_; z},
        index_type = function(z) {z$targets$target$index <- 1; z},
        no_cells = function(z) {z$targets$target$cells <- NULL; z},
        cell_length = function(z) {z$targets$target$cells <- z$targets$target$cells[1, ]; z},
        duplicate_coordinate = function(z) {z$targets$target$cells$column[] <- 2L; z},
        reordered_coordinates = function(z) {z$targets$target$cells <- z$targets$target$cells[2:1, ]; z},
        missing_coordinate = function(z) {z$targets$target$missing_columns <- 2L; z},
        coordinate_na = function(z) {z$targets$target$cells$column[1] <- NA_integer_; z},
        no_estimate = function(z) {z$targets$target$cells$estimate <- NULL; z},
        estimate_type = function(z) {z$targets$target$cells$estimate <- as.character(z$targets$target$cells$estimate); z},
        resolved_type = function(z) {z$targets$target$cells$resolved <- c(1L, 1L); z},
        resolved_na = function(z) {z$targets$target$cells$resolved[1] <- NA; z},
        reason_na = function(z) {z$targets$target$cells$reason[1] <- NA_character_; z},
        contributors_type = function(z) {z$targets$target$cells$contributors <- c(10, 10); z},
        contributors_na = function(z) {z$targets$target$cells$contributors[1] <- NA_integer_; z},
        short_neighbors = function(z) {z$targets$target$neighbors <- z$targets$target$neighbors[-1, ]; z},
        short_candidates = function(z) {z$targets$target$candidates <- z$targets$target$candidates[-1, ]; z},
        no_neighbors = function(z) {z$targets$target$neighbors <- NULL; z},
        bad_candidates = function(z) {z$targets$target$candidates$distance <- NULL; z},
        matrix_above = function(z) {z$matrix[1, 2] <- z$matrix[1, 2] + 1; z},
        matrix_below = function(z) {z$matrix[1, 2] <- z$matrix[1, 2] - 1; z},
        estimate_above = function(z) {z$targets$target$cells$estimate[1] <- z$targets$target$cells$estimate[1] + 1; z},
        estimate_below = function(z) {z$targets$target$cells$estimate[1] <- z$targets$target$cells$estimate[1] - 1; z})
    for (name in names(corruptions)) {
        corrupt <- corruptions[[name]]
        local_mocked_bindings(knn_vectorized_reference = function(x, targets)
            corrupt(frozen(x, targets)), .package = "BreastSubtypeR")
        error <- expect_error(BreastSubtypeR:::.mapping_impute_log(x, "target"),
            class = "BreastSubtypeR_knn_boundary", info = name)
        expect_true(!is.null(error$diagnostics), info = name)
    }
})

b6_two_target_input <- function() {
    ref <- b6_reference()
    ids <- intersect(as.character(ref$genes.sig50$EntrezGene.ID),
        as.character(ref$genes.signature$EntrezGene.ID[ref$genes.signature$SSP_based == "Yes"]))[1:2]
    x <- rbind(c(0, NA, 0, 0), c(100, 100, 100, NA),
        t(vapply(1:12, function(i) c(i, 10 + i, i, 110 + i), numeric(4))))
    dimnames(x) <- list(c(ids, as.character(910001:910012)), c("sample B", "sample/A", "s.3", "four"))
    list(x = x, ids = ids)
}

test_that("B6 corrections map two incomplete output targets in each supplied ordering", {
    fixture <- b6_two_target_input()
    frozen <- BreastSubtypeR:::knn_vectorized_reference
    observed <- NULL
    local_mocked_bindings(knn_vectorized_reference = function(x, targets) {
        result <- frozen(x, targets)
        observed <<- list(x = x, targets = targets, result = result)
        result
    }, .package = "BreastSubtypeR")
    for (ordering in c("original", "reordered")) {
        x <- fixture$x
        if (ordering == "reordered") x <- x[c(14:3, 2, 1), c(4, 2, 1, 3)]
        se <- b6_se(x); before <- se
        mapped <- Mapping(se, impute = TRUE, verbose = FALSE)
        targets <- rownames(x)[rownames(x) %in% fixture$ids]
        expect_identical(observed$x, x)
        expect_identical(observed$targets, targets)
        expected <- frozen(x, targets)
        expect_identical(observed$result, expected)
        values <- vapply(targets, function(id) expected$targets[[id]]$cells$estimate[[1]], numeric(1))
        expect_false(values[[1]] == values[[2]])
        ref <- b6_reference()
        nc <- as.matrix(SummarizedExperiment::assay(mapped$se_NC))
        ssp <- SummarizedExperiment::assay(mapped$se_SSP)
        for (id in targets) {
            symbol <- ref$genes.sig50$Symbol[match(id, as.character(ref$genes.sig50$EntrezGene.ID))]
            expect_identical(unname(nc[symbol, ]), unname(expected$matrix[id, ]))
            expect_identical(unname(ssp[id, ]), unname(2^expected$matrix[id, ]))
            present <- !is.na(x[id, ])
            expect_identical(unname(nc[symbol, present]), unname(x[id, present]))
        }
        expect_identical(colnames(mapped$se_NC), colnames(x))
        expect_identical(colnames(mapped$se_SSP), colnames(x))
        expect_identical(se, before)
        b6_save_interface(paste0("correction-two-targets-", ordering), observed$x,
            observed$targets, observed$result, expected)
    }
})

test_that("B6 corrections reject mixed resolution across two public output targets", {
    fixture <- b6_two_target_input()
    frozen <- BreastSubtypeR:::knn_vectorized_reference
    observed <- NULL
    local_mocked_bindings(knn_vectorized_reference = function(x, targets) {
        observed <<- list(x = x, targets = targets)
        frozen(x, targets)
    }, .package = "BreastSubtypeR")
    for (ordering in c("original", "reordered")) {
        x <- fixture$x; x[fixture$ids[2], ] <- NA_real_
        if (ordering == "reordered") x <- x[c(14:3, 2, 1), c(4, 2, 1, 3)]
        se <- b6_se(x); before <- se
        error <- expect_error(Mapping(se, impute = TRUE, verbose = FALSE),
            class = "BreastSubtypeR_knn_unsupported")
        targets <- rownames(x)[rownames(x) %in% fixture$ids]
        expect_identical(observed$x, x)
        expect_identical(observed$targets, targets)
        expected <- frozen(x, targets)
        expect_identical(error$diagnostics, expected)
        expect_identical(expected$targets[[fixture$ids[1]]]$status, "imputed")
        expect_identical(expected$targets[[fixture$ids[2]]]$status, "entirely_missing")
        expect_identical(error$targets, fixture$ids[2])
        expect_null(error$se_NC); expect_null(error$se_SSP)
        expect_identical(se, before)
        b6_save_interface(paste0("correction-mixed-targets-", ordering), observed$x,
            observed$targets, error$diagnostics, expected)
    }
})
