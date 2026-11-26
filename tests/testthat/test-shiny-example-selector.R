test_that("packaged example configuration rejects missing or inconsistent tables", {
    read_example <- BreastSubtypeR:::.read_shiny_example
    expect_error(read_example("other"), "supported packaged example")
    folder <- withr::local_tempdir()
    expect_error(read_example("tcga", folder), "files are missing")
    files <- system.file("RshinyTest", package = "BreastSubtypeR")
    paths <- list.files(files, "^TCGABRCAsub_", full.names = TRUE)
    expect_true(all(file.copy(paths, folder)))
    path <- file.path(folder, "TCGABRCAsub_anno.txt")
    anno <- read.delim(path)
    anno$Length[1L] <- 0
    write.table(anno, path, sep = "\t", quote = FALSE, row.names = FALSE)
    expect_error(read_example("tcga", folder), "positive, finite Length")
    anno$Length[1L] <- 1
    write.table(anno[nrow(anno):1L, ], path, sep = "\t", quote = FALSE, row.names = FALSE)
    expect_error(read_example("tcga", folder), "dimensions or identifiers")
})

test_that("confirmed examples enforce Mapping mode and isolate each upload generation", {
    skip_if_not_installed("shiny")
    skip_if_not_installed("bslib")
    env <- new.env(parent = asNamespace("shiny"))
    path <- system.file("ShinyBreastSubtypeR", "server.R", package = "BreastSubtypeR")
    eval(parse(file = path, encoding = "UTF-8", keep.source = FALSE), env)
    dialog <- character()
    env$showModal <- function(ui, ...) dialog <<- as.character(ui)
    failure <- new.env()
    reader <- BreastSubtypeR:::.read_shiny_example
    original <- BreastSubtypeR::Mapping
    seen <- list()
    testthat::local_mocked_bindings(.read_shiny_example = function(id) {
        if (is.null(failure$directory)) reader(id) else reader(id, failure$directory)
    }, Mapping = function(se_obj, RawCounts, impute, verbose) {
        seen[[length(seen) + 1L]] <<- list(dim(se_obj), RawCounts, impute, verbose)
        original(se_obj, RawCounts = RawCounts, impute = impute, verbose = verbose)
    }, .package = "BreastSubtypeR")
    shiny::testServer(env$server, {
        session$setInputs(is_raw_counts = "raw", BSmethod = "PAM50", k_subtypes = "5",
                          calibration = "Internal", internal = "medianCtr")
        session$setInputs(load_example = 1)
        first <- dialog
        session$setInputs(load_example_step1 = 1)
        expect_identical(dialog, first)
        expect_match(dialog, 'value="oslo" checked', fixed = TRUE)
        expect_match(dialog, "TCGA-BRCA", fixed = TRUE)
        expect_null(reactive_files$source)
        for (id in c("oslo", "tcga", "tcga", "oslo")) {
            raw <- id == "tcga"
            session$setInputs(example_choice = id, confirm_example = length(seen) + 1L)
            expect_null(current_input_coverage())
            expect_null(output$plotSection)
            for (kind in c("calls", "full")) {
                session$setInputs(export_kind = kind)
                expect_error(output$download, "No current prediction")
            }
            session$setInputs(is_raw_counts = if (raw) "norm" else "raw")
            session$setInputs(map = length(seen) + 1L)
            expected_dim <- if (raw) c(5350L, 100L) else c(569L, 102L)
            expect_identical(seen[[length(seen)]], list(expected_dim, raw, FALSE, TRUE))
            expect_true(mapping_current())
            session$setInputs(is_raw_counts = if (raw) "raw" else "norm")
            expect_true(mapping_current())
            session$setInputs(run = length(seen))
            for (kind in c("calls", "full")) {
                session$setInputs(export_kind = kind)
                tbl <- read.delim(output$download, comment.char = "#")
                expect_identical(tbl$PatientID, colnames(reactive_files$GEX))
            }
            session$setInputs(load_example = 2, example_choice = if (raw) "oslo" else "tcga")
            expect_identical(reactive_files$example_id, id)
            expect_false(is.null(current_prediction()))
        }
        # A completion from a previous file binding cannot replace this example.
        file <- data.frame(name = "old.txt", datapath = tempfile())
        session$setInputs(GEX_3 = file)
        expect_identical(reactive_files$source, "example")
        session$setInputs(upload_started = "GEX_3")
        expect_identical(reactive_files$source, "example")
        session$setInputs(upload_started = "GEX_4")
        expect_identical(reactive_files$source, "upload")
        expect_error(output$download, "No current prediction")
        session$setInputs(GEX_4 = file)
        expect_identical(reactive_files$source, "upload")
        expect_true(raw_mode()) # restored the original manual preference
        expect_null(reactive_files$anno)
        expect_null(reactive_files$clinic)
        session$setInputs(map = 8)
        expect_match(reactive_files$error, "upload expression, clinical and annotation")
        expect_length(seen, 4L)
        session$setInputs(example_choice = "tcga", confirm_example = 6)
        expect_identical(reactive_files$source, "example")
        # A missing installation file during replacement follows the same
        # empty-state policy as a malformed table. The real loader raises it.
        failure$directory <- tempfile()
        session$setInputs(example_choice = "oslo", confirm_example = 7)
        expect_null(reactive_files$source)
        expect_null(reactive_files$GEX)
        expect_null(reactive_files$anno)
        expect_null(reactive_files$clinic)
        expect_null(reactive_files$RawCounts)
        expect_match(reactive_files$error, "files are missing")
        expect_error(output$download, "No current prediction")
        failure$directory <- NULL
        session$setInputs(example_choice = "tcga", confirm_example = 8)
        expect_null(reactive_files$error)
        expect_identical(reactive_files$example_id, "tcga")
    })
})
