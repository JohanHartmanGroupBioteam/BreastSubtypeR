# Expanded public examples have their own current-code oracle. The older tests
# retain immutable historical inputs and every original expectation.
issue199_rng <- function(seed, code) {
    withr::with_preserve_seed({
        kind <- RNGkind()
        tryCatch({
            RNGkind("Mersenne-Twister", "Inversion", "Rejection")
            set.seed(seed)
            force(code)
        }, finally = do.call(RNGkind, as.list(kind)))
    })
}

issue199_fresh <- local({
    cache <- new.env(parent = emptyenv())
    function(name, historical = FALSE) {
        key <- paste(name, historical)
        if (!exists(key, cache, inherits = FALSE)) {
            obj <- if (historical) historical_example(name) else current_example(name)
            mapped <- issue199_rng(1992026, suppressMessages(Mapping(
                obj$se_obj, RawCounts = name == "TCGABRCAobj",
                method = "max", impute = TRUE, verbose = FALSE)))
            result <- issue199_rng(20260831, suppressMessages(BS_Multi(
                data_input = mapped, methods = "AUTO", Subtype = FALSE,
                hasClinical = FALSE, er_subsample = "target")))
            assign(key, list(data_input = mapped, res = result), cache)
        }
        get(key, cache, inherits = FALSE)
    }
})

issue199_calls <- function() readRDS(test_path("fixtures",
    "issue199-cb7a7c70-calls.rds"))

test_that("the 19 portability paths enforce elementwise bounds and exact decisions", {
    policy <- issue199_portability_policy()
    expect_identical(nrow(policy), 19L)
    expect_identical(sum(policy$dataset == "OSLO2EMIT0obj"), 5L)
    expect_identical(sum(policy$dataset == "TCGABRCAobj" &
        policy$scope == "C_vs_approved_C"), 13L)
    expect_identical(sum(policy$scope == "B_mapping_vs_frozen_A"), 1L)
    expect_false(anyDuplicated(policy[c("scope", "dataset", "path")]) > 0L)
    for (i in seq_len(nrow(policy))) {
        rule <- policy[i, ]
        posterior <- grepl("/all.probs/", rule$path, fixed = TRUE)
        expected <- if (posterior) matrix(c(.4, .3, .2, .09, .01), nrow = 1L) else
            matrix(1, nrow = 1L)
        compare <- function(actual) issue199_portable_compare(actual, expected,
            rule$dataset, rule$scope, rule$path)
        actual <- expected
        bound <- rule$atol + rule$rtol * abs(expected[1L])
        actual[1L] <- expected[1L] + bound / 2
        expect_true(compare(actual), info = rule$path)
        actual[1L] <- expected[1L] + 4 * bound
        expect_error(compare(actual), "numeric bound exceeded", info = rule$path)
    }
    p <- "/res/results/AIMS/all.probs/20"
    expected <- matrix(c(.5, .5, 0, 0, 0), nrow = 1L)
    actual <- expected
    actual[1L, 1:2] <- c(.5 + 2^-53, .5 - 2^-53)
    expect_error(issue199_portable_compare(actual, expected, "OSLO2EMIT0obj",
        path = p), "posterior ranks/ties/argmax")
    expected[1L, 1:2] <- rev(actual[1L, 1:2])
    expect_error(issue199_portable_compare(actual, expected, "OSLO2EMIT0obj",
        path = p), "posterior ranks/ties/argmax")
    cut <- 100 * (-.25 + .35) / .85
    expect_error(issue199_portable_compare(cut + 2^-49, cut, "TCGABRCAobj",
        path = "/res/results/ssBC/score.ROR/ROR-P (Subtype + Prolif)"),
        "ROR cutoff sides")
})

test_that("portability never relaxes structure, metadata, masks or unlisted values", {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1, nrow = 1L)))
    posterior <- matrix(c(.4, .3, .2, .09, .01), nrow = 1L)
    expected <- list(data_input = list(se_NC = se, se_SSP = se),
        res = list(results = list(
            AIMS = list(all.probs = list(`20` = posterior), data.used = matrix(1:2 + 0, 1L)),
            sspbc = list(all.probs = list(`24` = posterior), data.used = matrix(1, 1L))),
            untouched = 1, call = "LumA"))
    compare <- function(actual, reference = expected) issue199_portable_compare(
        actual, reference, "OSLO2EMIT0obj")
    expect_true(compare(expected))
    actual <- expected; actual$res$untouched <- 1 + .Machine$double.eps
    expect_error(compare(actual), "unlisted value must be exact")
    actual <- expected; reference <- expected
    actual$res$untouched <- -0; reference$res$untouched <- 0
    expect_error(compare(actual, reference), "unlisted value must be exact")
    actual <- expected; actual$res$call <- "LumB"
    expect_error(compare(actual), "unlisted value must be exact")
    actual <- expected; dim(actual$res$results$AIMS$data.used) <- c(2L, 1L)
    expect_error(compare(actual), "attributes/names/dimensions")
    actual <- expected; S4Vectors::metadata(actual$data_input$se_NC) <- list(changed = TRUE)
    expect_error(compare(actual), "type/class/length")
    actual <- expected
    SummarizedExperiment::assay(actual$data_input$se_SSP, withDimnames = FALSE) <-
        as.data.frame(SummarizedExperiment::assay(actual$data_input$se_SSP, withDimnames = FALSE))
    expect_error(compare(actual), "type/class/length")
    actual <- expected; SummarizedExperiment::assayNames(actual$data_input$se_SSP) <- "other"
    expect_error(compare(actual), "attributes/names/dimensions")
    actual <- expected
    SummarizedExperiment::assays(actual$data_input$se_SSP)$extra <- matrix(1, 1L)
    expect_error(compare(actual), "type/class/length")
    actual <- expected; attr(actual$res$results$AIMS$data.used, "changed") <- TRUE
    expect_error(compare(actual), "attributes/names/dimensions")
    actual <- expected; class(actual$res$results$AIMS$data.used) <- "changed"
    expect_error(compare(actual), "type/class/length")
    actual <- expected; storage.mode(actual$res$results$AIMS$data.used) <- "integer"
    expect_error(compare(actual), "type/class/length")
    for (value in c(NA_real_, NaN, Inf)) {
        actual <- expected; actual$res$results$AIMS$data.used[1L] <- value
        expect_error(compare(actual), "nonfinite masks/signs")
    }
    # A misspelled named posterior child must leave its exact policy path unused.
    actual <- expected; names(actual$res$results$AIMS$all.probs) <- "2"
    expect_error(compare(actual, actual), "unused.*policy path")
    p <- "/data_input/se_NC/assay"
    expect_error(issue199_portable_compare(1e-320, 0, "TCGABRCAobj", path = p),
        "zero support/signs")
    expect_error(issue199_portable_compare(-0, 0, "TCGABRCAobj", path = p), "signed zero")
    expect_error(issue199_portable_compare(1e-310 + 5e-324, 1e-310,
        "TCGABRCAobj", path = p), "subnormal support/values")
    expect_error(issue199_portable_compare(-Inf, Inf, "TCGABRCAobj", path = p),
        "nonfinite masks/signs")
})

test_that("compact historical files retain pinned payloads and original provenance", {
    manifest <- read.delim(test_path("fixtures", "historical-cb7a7c70",
        "manifest.tsv"), stringsAsFactors = FALSE)
    # Stage 3C changes storage only. Original file hashes remain provenance;
    # compact bytes and independent frozen reconstructed values are both checked.
    expected <- c(OSLO2EMIT0obj.rds = "e3a859333650c120905d1c1611401a84",
        TCGABRCAobj.rds = "f6ba54fa6777e45355a6879a22dd161e")
    original_md5 <- c(OSLO2EMIT0obj.rds = "32b39cceec5075f8eab246f9ac71aecb",
        TCGABRCAobj.rds = "5139e1f0a00cbd76d6869f2f8fbbbc4a")
    original_sha256 <- c(
        OSLO2EMIT0obj.rds = "43393c2bf06f0ea44b262d5525b75154ab99162fac5de898433556b6805450ec",
        TCGABRCAobj.rds = "b94309c6f915d5f5fcb18368f3876f281b2111acf50eeb25b167f201f499a45d")
    for (nm in names(expected)) {
        path <- test_path("fixtures", "historical-cb7a7c70", nm)
        row <- match(nm, manifest$file)
        expect_identical(unname(tools::md5sum(path)), unname(expected[nm]))
        expect_identical(manifest$md5[row], unname(expected[nm]))
        expect_identical(manifest$original_md5[row], unname(original_md5[nm]))
        expect_identical(manifest$original_sha256[row], unname(original_sha256[nm]))
    }
    expect_identical(dim(historical_example("OSLO2EMIT0obj")$se_obj), c(356L, 102L))
    old_tcga <- historical_example("TCGABRCAobj")
    expect_identical(dim(old_tcga$se_obj), c(5000L, 100L))
    expect_identical(dim(old_tcga$data_input$se_NC), c(37L, 100L))
    expect_identical(dim(old_tcga$data_input$se_SSP), c(178L, 100L))
})

test_that("expanded public examples retain original observations and metadata", {
    sizes <- list(OSLO2EMIT0obj = c(569L, 102L), TCGABRCAobj = c(5350L, 100L))
    for (nm in names(sizes)) {
        obj <- current_example(nm)
        old <- historical_example(nm)
        expect_identical(names(obj), c("se_obj", "data_input", "res"), info = nm)
        expect_identical(names(obj$data_input), c("se_NC", "se_SSP"), info = nm)
        expect_identical(names(obj$res), c("res_subtypes", "results"), info = nm)
        expect_identical(class(obj$se_obj), class(old$se_obj), info = nm)
        expect_identical(dim(obj$se_obj), sizes[[nm]], info = nm)
        first <- seq_len(nrow(old$se_obj))
        expect_identical(SummarizedExperiment::assay(obj$se_obj)[first, , drop = FALSE],
            SummarizedExperiment::assay(old$se_obj), info = nm)
        expect_identical(SummarizedExperiment::rowData(obj$se_obj)[first, ],
            SummarizedExperiment::rowData(old$se_obj), info = nm)
        old_metadata <- S4Vectors::metadata(old$se_obj)
        retained_metadata <- S4Vectors::metadata(obj$se_obj)
        expect_identical(setdiff(names(retained_metadata), names(old_metadata)),
            "issue199_reconstruction", info = nm)
        retained_metadata$issue199_reconstruction <- NULL
        # Removing the only named member may leave names=character(0) on list().
        if (!length(retained_metadata)) retained_metadata <- list()
        expect_identical(retained_metadata, old_metadata, info = nm)
        expect_identical(SummarizedExperiment::colData(obj$se_obj),
            SummarizedExperiment::colData(old$se_obj), info = nm)
        expect_identical(colnames(obj$se_obj), colnames(old$se_obj), info = nm)
        expect_false(anyDuplicated(rownames(obj$se_obj)) > 0L, info = nm)
        expect_false(anyDuplicated(colnames(obj$se_obj)) > 0L, info = nm)
        ids <- as.character(SummarizedExperiment::rowData(obj$se_obj)$ENTREZID)
        added <- ids[-first]
        expect_false(anyNA(added), info = nm)
        expect_false(any(added == ""), info = nm)
        expect_false(anyDuplicated(added) > 0L, info = nm)
        expect_identical(length(intersect(added, ids[first])), 0L, info = nm)
        expect_identical(length(added), if (nm == "OSLO2EMIT0obj") 213L else 350L)
        if (nm == "OSLO2EMIT0obj") expect_false(anyDuplicated(ids) > 0L)
        # TCGA's original missing/shared Entrez annotations remain unchanged;
        # only appended model-gene identifiers must be complete and unique.
        values <- as.matrix(SummarizedExperiment::assay(obj$se_obj))
        expect_true(all(is.finite(values)), info = nm)
        if (nm == "TCGABRCAobj") {
            expect_true(all(values >= 0 & values == floor(values)))
        }
    }
})

test_that("default Mapping reconstructs current examples without changing RNG", {
    withr::local_seed(199)
    for (nm in c("OSLO2EMIT0obj", "TCGABRCAobj")) {
        obj <- current_example(nm)
        before <- obj$se_obj
        rng <- .Random.seed
        mapped <- suppressMessages(Mapping(obj$se_obj,
            RawCounts = nm == "TCGABRCAobj", verbose = FALSE))
        expect_identical(.Random.seed, rng, info = nm)
        expect_identical(obj$se_obj, before, info = nm)
        expect_true(issue199_portable_compare(mapped, obj$data_input, nm,
            path = "/data_input"), info = nm)
        for (method in names(mapped)) {
            expect_identical(dim(mapped[[method]]),
                c(if (method == "se_NC") 50L else 339L, ncol(obj$se_obj)), info = nm)
            expect_identical(colnames(mapped[[method]]), colnames(obj$se_obj), info = nm)
            expect_identical(SummarizedExperiment::colData(mapped[[method]]),
                SummarizedExperiment::colData(obj$se_obj), info = nm)
        }
    }
})

test_that("expanded current AUTO results reproduce the accepted complete objects", {
    reference <- issue199_calls()
    withr::local_options(BreastSubtypeR.verbose = FALSE)
    withr::local_seed(199)
    for (nm in names(reference)) {
        obj <- current_example(nm)
        rng <- .Random.seed
        fresh <- issue199_fresh(nm)
        expect_identical(.Random.seed, rng, info = nm)
        expect_true(issue199_portable_compare(fresh$data_input, obj$data_input,
            nm, path = "/data_input"), info = nm)
        # Only measured component paths admit bounded elementwise roundoff;
        # all other values, structure and scientific decisions remain exact.
        # #233 expands only ssBC AUTO sample axes. Check new NA rows/axes,
        # then compare every retained payload against the unchanged fixture.
        compatible <- auto_ssbc_historical_view(fresh$res, obj$res)
        expect_true(issue199_portable_compare(compatible, obj$res, nm,
            path = "/res"), info = nm)
        issue199_expect_decisions(compatible, obj$res)
        methods <- setdiff(names(obj$res$res_subtypes), "entropy")
        expect_identical(obj$res$res_subtypes[, methods],
            reference[[nm]]$C[, methods], info = nm)
        expect_identical(fresh$res$res_subtypes, reference[[nm]]$C, info = nm)
        expect_identical(obj$res$res_subtypes$entropy, reference[[nm]]$C$entropy,
            info = nm)
        expect_identical(length(obj$res$results), 8L, info = nm)
        expect_identical(sum(is.na(obj$res$res_subtypes$ssBC.v2)),
            if (nm == "OSLO2EMIT0obj") 8L else 4L, info = nm)
        calls <- obj$res$res_subtypes[, methods, drop = FALSE]
        entropy <- apply(calls, 1L, function(x) {
            x <- x[!is.na(x)]
            if (!length(x)) return(NA_real_)
            p <- as.numeric(table(x)) / length(x)
            -sum(p * log2(p))
        })
        expect_equal(unname(entropy), obj$res$res_subtypes$entropy,
            tolerance = 1e-12, info = nm)
    }
})

test_that("OSLO expansion changes neither mapped values nor current method outputs", {
    withr::local_options(BreastSubtypeR.verbose = FALSE)
    old <- historical_example("OSLO2EMIT0obj")
    obj <- current_example("OSLO2EMIT0obj")
    b <- issue199_fresh("OSLO2EMIT0obj", historical = TRUE)
    c <- issue199_fresh("OSLO2EMIT0obj")
    expect_identical(b, c)
    methods <- setdiff(names(old$res$res_subtypes), "entropy")
    expect_identical(obj$res$res_subtypes[, methods],
        old$res$res_subtypes[, methods])
    expect_equal(obj$res$res_subtypes$entropy, old$res$res_subtypes$entropy,
        tolerance = 1e-12)
    expect_identical(b$res$res_subtypes[, methods],
        issue199_calls()$OSLO2EMIT0obj$B[, methods])
    expect_identical(b$res$res_subtypes, issue199_calls()$OSLO2EMIT0obj$B)
    # Historical stored ROR-P predates #198. It remains testable in the
    # historical suite and is deliberately not relabelled as current output.
})

test_that("TCGA expanded counts reproduce the 64 documented dataset-induced calls", {
    withr::local_options(BreastSubtypeR.verbose = FALSE)
    b <- issue199_fresh("TCGABRCAobj", historical = TRUE)
    c <- issue199_fresh("TCGABRCAobj")
    reference <- issue199_calls()$TCGABRCAobj
    methods <- setdiff(names(reference$B), "entropy")
    expect_identical(b$res$res_subtypes[, methods], reference$B[, methods])
    expect_identical(c$res$res_subtypes[, methods], reference$C[, methods])
    expect_identical(b$res$res_subtypes, reference$B)
    expect_identical(c$res$res_subtypes, reference$C)
    expect_true(issue199_portable_compare(b$data_input,
        historical_example("TCGABRCAobj")$data_input, "TCGABRCAobj",
        "B_mapping_vs_frozen_A", "/data_input"))
    rows <- lapply(methods, function(method) {
        a <- b$res$res_subtypes[[method]]
        z <- c$res$res_subtypes[[method]]
        a[is.na(a)] <- "<NA>"
        z[is.na(z)] <- "<NA>"
        changed <- which(a != z)
        data.frame(sample_id = rownames(reference$B)[changed],
            method = rep(method, length(changed)),
            old_call = a[changed], new_call = z[changed], stringsAsFactors = FALSE)
    })
    actual <- do.call(rbind, rows)
    expected <- read.delim(test_path("fixtures", "issue199-TCGA-B-C-changes.tsv"),
        stringsAsFactors = FALSE, check.names = FALSE)
    sort_rows <- function(x) {
        x <- x[order(x$method, x$sample_id), , drop = FALSE]
        rownames(x) <- NULL
        x
    }
    expect_identical(sort_rows(actual), sort_rows(expected))
    expect_identical(nrow(actual), 64L)
    expect_identical(length(unique(actual$sample_id)), 40L)
    counts <- c(genefu.robust = 8L, ssBC = 4L, ssBC.v2 = 6L, cIHC = 5L,
        cIHC.itr = 6L, PCAPAM50 = 8L, AIMS = 16L, sspbc = 11L)
    expect_identical(as.integer(table(actual$method)[names(counts)]), unname(counts))
})

test_that("expanded TCGA has complete proliferation panels with independent ROR arithmetic", {
    se <- current_example("TCGABRCAobj")$data_input$se_NC
    p18 <- c("ANLN", "CEP55", "ORC6L", "CCNE1", "EXO1", "PTTG1",
        "CDC20", "KIF2C", "RRM2", "CDC6", "KNTC2", "TYMS", "CDCA1",
        "MELK", "UBE2C", "CENPF", "MKI67", "UBE2T")
    expect_identical(length(intersect(ror198_panel, rownames(se))), 11L)
    expect_identical(length(intersect(p18, rownames(se))), 18L)
    res <- suppressMessages(BS_parker(se, calibration = "Internal",
        internal = "medianCtr", Subtype = TRUE, hasClinical = FALSE))
    ror198_expect_reference(res)
    expect_identical(as.character(res$score.ROR$patientID), colnames(se))
    expect_true(all(is.finite(res$score.ROR[["ROR-P (Subtype + Prolif)"]])))
})
