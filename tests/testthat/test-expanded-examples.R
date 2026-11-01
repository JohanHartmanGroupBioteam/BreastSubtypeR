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
        expect_identical(mapped, obj$data_input, info = nm)
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
        expect_identical(fresh$data_input, obj$data_input, info = nm)
        # A strict current-runtime replay; any future intentional correction must
        # receive a separate versioned baseline, never overwrite the old oracle.
        expect_identical(fresh$res, obj$res, info = nm)
        methods <- setdiff(names(obj$res$res_subtypes), "entropy")
        expect_identical(obj$res$res_subtypes[, methods],
            reference[[nm]]$C[, methods], info = nm)
        expect_equal(obj$res$res_subtypes$entropy, reference[[nm]]$C$entropy,
            tolerance = 1e-12, info = nm)
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
    expect_equal(b$res$res_subtypes$entropy, reference$B$entropy, tolerance = 1e-12)
    expect_identical(b$data_input, historical_example("TCGABRCAobj")$data_input)
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
