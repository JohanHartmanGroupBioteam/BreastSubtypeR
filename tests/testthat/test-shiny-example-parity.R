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
    path <- system.file("RshinyTest", name, package = "BreastSubtypeR")
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

test_that("the six Shiny tables retain schemas, IDs and exact candidate values", {
    for (name in c("OSLO2EMIT0obj", "TCGABRCAobj")) {
        obj <- shiny199_object(name)
        raw <- identical(name, "TCGABRCAobj")
        prefix <- if (raw) "TCGABRCAsub" else "OSLO2EMIT0"
        exdir <- system.file("RshinyTest", package = "BreastSubtypeR")
        path <- function(suffix) file.path(exdir, paste0(prefix, suffix))
        gex <- as.matrix(utils::read.delim(path(if (raw)
            "_GEX_RawCounts.txt" else "_GEX_log2.FPKM.txt"),
            row.names = 1L, check.names = FALSE))
        storage.mode(gex) <- "double"
        expected <- as.matrix(SummarizedExperiment::assay(obj$se_obj))
        storage.mode(expected) <- "double"
        expect_identical(gex, expected)
        expect_identical(dim(gex), if (raw) c(5350L, 100L) else c(569L, 102L))
        expect_false(anyDuplicated(rownames(gex)) > 0L)
        expect_false(anyDuplicated(colnames(gex)) > 0L)
        expect_true(all(is.finite(gex)))
        if (raw) expect_true(all(gex >= 0 & gex == floor(gex)))
        ann <- utils::read.delim(path("_anno.txt"), check.names = FALSE)
        clinical <- utils::read.delim(path("_clinical.txt"), check.names = FALSE)
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

test_that("installed Shiny examples map, calculate AUTO and download candidate calls", {
    withr::local_preserve_seed()
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
            mapped <- reactive_files$data_input
            expect_identical(names(mapped), c("se_NC", "se_SSP"))
            for (component in names(mapped)) {
                expect_true(issue199_portable_compare(
                    SummarizedExperiment::assay(mapped[[component]]),
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

