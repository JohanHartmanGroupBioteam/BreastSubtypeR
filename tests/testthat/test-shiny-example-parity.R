# Current Shiny examples are projections of the expanded public objects.
# The preserved historical Shiny ancillary OSLO fields are not model inputs.
shiny199_object <- function(name) {
    e <- new.env(parent = emptyenv())
    utils::data(list = name, package = "BreastSubtypeR", envir = e)
    e[[name]]
}
shiny199_server <- function() {
    skip_if_not_installed("shiny")
    skip_if_not_installed("bslib")
    e <- new.env(parent = asNamespace("shiny"))
    eval(parse(file = system.file("ShinyBreastSubtypeR", "server.R",
                                 package = "BreastSubtypeR"),
               encoding = "UTF-8", keep.source = FALSE), envir = e)
    e$server
}
shiny199_upload <- function(name) {
    path <- shiny199_checked_path(name)
    data.frame(name = name, size = file.info(path)$size,
               type = "text/plain", datapath = path)
}
shiny199_fields_equal <- function(actual, expected) {
    expect_true(all(names(expected) %in% names(actual)))
    for (field in names(expected)) {
        expect_identical(as.character(actual[[field]]),
                         as.character(expected[[field]]), info = field)
    }
}

# Pinned table bytes: data-raw/issue199/manifests/shiny-tables.tsv and
# line-ending-corrections-v1.tsv (canonical LF / documented historical CRLF).
# MD5 is an integrity check using base R, as in the historical fixture tests.
shiny199_checked_path <- function(name) {
    fingerprints <- list(
        "OSLO2EMIT0_GEX_log2.FPKM.txt" = c("a0a58a7ace7f8bc34cc91782fa58b16f",
            "3f963c4265baa3d481759d2b306ed1af"),
        "OSLO2EMIT0_anno.txt" = c("e4195ac2cebd9b582455f19083295aa6",
            "cc444d6b6ede27a081b5d20f816eb480"),
        "OSLO2EMIT0_clinical.txt" = c("540e76446c86f50c6c08bb0c60f33bc6",
            "102064e190a9674db8e53424e4bae677"),
        "TCGABRCAsub_GEX_RawCounts.txt" = c("1b8fb97f4ba04bdaeded2049ff4d7d58",
            "eb61b6eaf77994d335a533fd68dd9e78"),
        "TCGABRCAsub_anno.txt" = c("9353c305399ebfb6045510ec3a156799",
            "143c26e5a28df2d81641c2f33096d302"),
        "TCGABRCAsub_clinical.txt" = c("76dc3e9b7a8848c2481685476915c7ae",
            "6a279eb59eba50f3db7bfe48d81b446b"))
    path <- system.file("RshinyTest", name, package = "BreastSubtypeR")
    observed <- if (nzchar(path) && file.exists(path))
        unname(tools::md5sum(path)) else "<missing>"
    if (!name %in% names(fingerprints) || !observed %in% fingerprints[[name]])
        stop("Shiny table bytes differ from pinned provenance: ", name,
             "; observed MD5=", observed, "; ", shiny199_runtime(), call. = FALSE)
    path
}
shiny199_numeric_fingerprint <- function(x) {
    path <- tempfile("shiny-example-profile-")
    on.exit(unlink(path))
    writeBin(as.double(x), path, size = 8L, endian = "big")
    unname(tools::md5sum(path))
}
shiny199_runtime <- function() {
    sprintf("R=%s; OS=%s; arch=%s; sizeof.longdouble=%s; long.double=%s",
            as.character(getRversion()), unname(Sys.info()["sysname"]),
            R.version$arch, .Machine$sizeof.longdouble,
            capabilities("long.double")[[1L]])
}
shiny199_oslo_profile <- function(gex, expected) {
    # Entire 569 x 102 arrays, in column-major binary64 big-endian order.
    # Reference: shipped binary; alternative: original native artifact
    # e1cf864bfa6d23646a313c7798e0a4b41430a3aa672421760e21091e099a990e,
    # baseline/candidate collect/decimal-OSLO2EMIT0obj/complete-arrays.rds.
    # Captured natively on R 4.6.2; runtime eligibility is version-independent.
    # The historical f6dbdbb/3755d3d capture is NOT a run of this source.
    # No single-token probe, measured-error tolerance or OS-only exemption.
    fingerprint <- if (is.double(gex)) shiny199_numeric_fingerprint(gex) else
        paste0("<unavailable: type=", typeof(gex), ">")
    reference <- shiny199_numeric_fingerprint(expected)
    fail <- function(reason) stop(reason, "; observed fingerprint=", fingerprint,
        "; reference fingerprint=", reference, "; ", shiny199_runtime(), call. = FALSE)
    if (!identical(dim(gex), c(569L, 102L)) ||
        !identical(typeof(gex), "double") ||
        !identical(dimnames(gex), dimnames(expected)))
        fail(sprintf("OSLO input structure changed: dim=%s; type=%s; dimnames_match=%s",
             paste(dim(gex), collapse = "x"), typeof(gex),
             identical(dimnames(gex), dimnames(expected))))
    if (!identical(reference, "35f0427b506fa811a82ff98e6a7367da"))
        fail("OSLO binary reference changed")
    if (identical(fingerprint, "35f0427b506fa811a82ff98e6a7367da"))
        return("exact")
    native <- identical(unname(Sys.info()["sysname"]), "Darwin") &&
        R.version$arch %in% c("aarch64", "arm64") &&
        .Machine$sizeof.longdouble == 8L &&
        !isTRUE(capabilities("long.double")[[1L]])
    if (native && identical(fingerprint,
                            "30f29ac10f58103e923033d01facadc1"))
        return("native-full-table")
    fail("Unknown or inconsistent OSLO decimal-conversion profile")
}

test_that("the six Shiny tables retain provenance, schemas, IDs and associations", {
    for (name in c("TCGABRCAobj", "OSLO2EMIT0obj")) {
        obj <- shiny199_object(name)
        raw <- identical(name, "TCGABRCAobj")
        prefix <- if (raw) "TCGABRCAsub" else "OSLO2EMIT0"
        path <- function(suffix) shiny199_checked_path(paste0(prefix, suffix))
        gex <- as.matrix(utils::read.delim(path(if (raw)
            "_GEX_RawCounts.txt" else "_GEX_log2.FPKM.txt"),
            row.names = 1L, check.names = FALSE))
        storage.mode(gex) <- "double"
        expected <- as.matrix(SummarizedExperiment::assay(obj$se_obj))
        storage.mode(expected) <- "double"
        expect_identical(dimnames(gex), dimnames(expected))
        if (raw) {
            expect_identical(gex, expected)
        } else {
            expect_true(shiny199_oslo_profile(gex, expected) %in%
                        c("exact", "native-full-table"))
        }
        expect_identical(dim(gex), if (raw) c(5350L, 100L) else c(569L, 102L))
        expect_false(anyDuplicated(rownames(gex)) > 0L)
        expect_false(anyDuplicated(colnames(gex)) > 0L)
        expect_true(all(is.finite(gex)))
        if (raw) expect_true(all(gex >= 0 & gex == floor(gex)))
        ann <- utils::read.delim(path("_anno.txt"), check.names = FALSE)
        clinical <- utils::read.delim(path("_clinical.txt"), check.names = FALSE)
        parsed <- BreastSubtypeR:::.read_shiny_expression(
            path(if (raw) "_GEX_RawCounts.txt" else "_GEX_log2.FPKM.txt"),
            "\t", clinical$PatientID, ann$probe, raw)
        expect_identical(parsed, gex)
        expect_identical(ann$probe, rownames(gex))
        expect_identical(clinical$PatientID, colnames(gex))
        expect_identical(names(ann), if (raw)
            c("probe", "SYMBOL", "ENTREZID", "Length") else c("probe", "SYMBOL", "ENTREZID"))
        expect_identical(names(clinical), if (raw)
            c("PatientID", "ER", "PR", "HER2") else
            c("PatientID", "ER", "PR", "HER2", "TSIZE", "NODE", "IHC",
              "Prosigna.Subtype", "Prosigna.ROR", "MKI67.corrected",
              "Ki67.estimated", "IHC.estimated"))
        shiny199_fields_equal(ann, as.data.frame(SummarizedExperiment::rowData(obj$se_obj)))
        shiny199_fields_equal(clinical, as.data.frame(SummarizedExperiment::colData(obj$se_obj)))
    }
})

# Only this block may skip for decimal conversion. Every other block runs.
test_that("OSLO text is bitwise identical to its frozen binary on the exact profile", {
    obj <- shiny199_object("OSLO2EMIT0obj")
    path <- shiny199_checked_path("OSLO2EMIT0_GEX_log2.FPKM.txt")
    gex <- as.matrix(utils::read.delim(path, row.names = 1L, check.names = FALSE))
    expected <- as.matrix(SummarizedExperiment::assay(obj$se_obj))
    storage.mode(gex) <- storage.mode(expected) <- "double"
    if (shiny199_oslo_profile(gex, expected) == "native-full-table")
        skip(paste("OSLO frozen-binary exactness only: verified complete",
                   "58,038-cell archived native decimal-conversion profile;",
                   "all provenance, same-input, TCGA and output checks still run"))
    expect_identical(gex, expected)
})

test_that("installed Shiny examples map, calculate AUTO and download candidate calls", {
    withr::local_preserve_seed()
    original_mapping <- BreastSubtypeR::Mapping
    mapping_input <- mapping_args <- NULL
    testthat::local_mocked_bindings(Mapping = function(se_obj, ...) {
        mapping_input <<- se_obj
        mapping_args <<- list(...)
        original_mapping(se_obj = se_obj, ...)
    }, .package = "BreastSubtypeR")
    for (name in c("OSLO2EMIT0obj", "TCGABRCAobj")) {
        obj <- shiny199_object(name)
        raw <- identical(name, "TCGABRCAobj")
        shiny::testServer(shiny199_server(), {
            session$setInputs(BSmethod = "AUTO Mode", k_subtypes = "5",
                              is_raw_counts = if (raw) "raw" else "norm",
                              hasClinical = FALSE, export_kind = "calls")
            if (raw) {
                session$setInputs(
                    GEX = shiny199_upload("TCGABRCAsub_GEX_RawCounts.txt"),
                    clinic = shiny199_upload("TCGABRCAsub_clinical.txt"),
                    anno = shiny199_upload("TCGABRCAsub_anno.txt"))
            } else {
                session$setInputs(load_example_step1 = 1)
                expect_identical(dim(reactive_files$GEX), c(569L, 102L))
            }
            RNGkind("Mersenne-Twister", "Inversion", "Rejection")
            set.seed(1992026)
            session$setInputs(map = 1)
            expect_null(reactive_files$error)
            expect_true(mapping_current())
            expect_identical(reactive_files$RawCounts, raw)
            expected_input <- as.matrix(SummarizedExperiment::assay(obj$se_obj))
            storage.mode(expected_input) <- "double"
            prefix <- if (raw) "TCGABRCAsub" else "OSLO2EMIT0"
            path <- function(suffix) shiny199_checked_path(paste0(prefix, suffix))
            parsed_input <- as.matrix(utils::read.delim(path(if (raw)
                "_GEX_RawCounts.txt" else "_GEX_log2.FPKM.txt"),
                row.names = 1L, check.names = FALSE))
            storage.mode(parsed_input) <- "double"
            expect_identical(SummarizedExperiment::assay(mapping_input), parsed_input)
            expect_identical(reactive_files$GEX, parsed_input)
            if (raw || identical(shiny199_oslo_profile(parsed_input, expected_input), "exact"))
                expect_identical(SummarizedExperiment::assay(mapping_input), expected_input)
            ann <- utils::read.delim(path("_anno.txt"), check.names = FALSE)
            clinical <- utils::read.delim(path("_clinical.txt"), check.names = FALSE)
            rownames(ann) <- ann$probe
            rownames(clinical) <- clinical$PatientID
            expect_identical(as.data.frame(SummarizedExperiment::rowData(mapping_input)), ann)
            expect_identical(as.data.frame(SummarizedExperiment::colData(mapping_input)), clinical)
            expect_identical(mapping_args,
                             list(RawCounts = raw, impute = FALSE, verbose = TRUE))
            # Exact actual-server-input comparison, including its metadata.
            direct <- do.call(original_mapping, c(list(se_obj = mapping_input), mapping_args))
            # Keep the frozen Mapping rules, but supply their original binary
            # measurements, not a potentially different text-converted matrix.
            binary_input <- mapping_input
            SummarizedExperiment::assay(binary_input) <- expected_input
            binary_mapped <- do.call(original_mapping,
                c(list(se_obj = binary_input), mapping_args))
            mapped <- reactive_files$data_input
            expect_identical(mapped, direct)
            expect_identical(names(mapped), c("se_NC", "se_SSP"))
            for (component in names(mapped)) {
                expect_true(issue199_portable_compare(
                    SummarizedExperiment::assay(binary_mapped[[component]]),
                    SummarizedExperiment::assay(obj$data_input[[component]]), name,
                    path = paste0("/data_input/", component, "/assay")))
                expect_identical(rownames(mapped[[component]]), rownames(obj$data_input[[component]]))
                expect_identical(colnames(mapped[[component]]), colnames(obj$data_input[[component]]))
                shiny199_fields_equal(
                    as.data.frame(SummarizedExperiment::colData(mapped[[component]])),
                    as.data.frame(SummarizedExperiment::colData(obj$data_input[[component]])))
            }
            expect_true(auto_requirements()$ready)
            set.seed(20260831)
            session$setInputs(run = 1)
            expect_null(reactive_files$error)
            prediction <- current_prediction()
            expect_false(is.null(prediction))
            expected_calls <- data.frame(PatientID = rownames(obj$res$res_subtypes),
                                         obj$res$res_subtypes, check.names = FALSE)
            expect_identical(prediction$calls, expected_calls)
            expect_identical(prediction$plot, obj$res$res_subtypes)
            expect_false(prediction$want_4)
            expect_null(prediction$full)
            export_path <- output$download
            expect_true(file.exists(export_path))
            export <- utils::read.delim(export_path, comment.char = "#",
                                       check.names = FALSE)
            shiny199_fields_equal(export[, setdiff(names(export), "entropy"), drop = FALSE],
                                  expected_calls[, setdiff(names(expected_calls), "entropy"), drop = FALSE])
            if ("entropy" %in% names(export)) expect_equal(export$entropy, expected_calls$entropy,
                                                         tolerance = 1e-14)
            expect_identical(names(export), names(expected_calls))
            expect_equal(nrow(export), ncol(obj$se_obj))
            expect_match(readLines(export_path)[1L], "BreastSubtypeR", fixed = TRUE)
            expect_match(readLines(export_path)[2L], "Classes: 5", fixed = TRUE)
        })
    }
})

