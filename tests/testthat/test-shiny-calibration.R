test_that("PAM50 custom medians belong to the current Shiny dataset", {
    skip_if_not_installed("shiny")
    skip_if_not_installed("bslib")
    app <- new.env(parent = asNamespace("shiny"))
    path <- system.file("ShinyBreastSubtypeR/server.R", package = "BreastSubtypeR")
    eval(parse(path, encoding = "UTF-8"), envir = app)
    # Simulate an unavailable example installation without changing its files.
    missing_example <- FALSE
    app$system.file <- function(...) {
        if (missing_example) return(tempfile())
        base::system.file(...)
    }
    data("BreastSubtypeRobj", package = "BreastSubtypeR")
    meta <- function(path) data.frame(name = basename(path), datapath = path,
                                     size = file.size(path), type = "text/plain")
    write_reference <- function(column, ext) {
        md <- BreastSubtypeRobj$medians[, c("X", column)]
        names(md)[2L] <- "Given.mdns"
        path <- tempfile(fileext = ext)
        sep <- if (ext == ".csv") "," else "\t"
        write.table(md, path, sep = sep, quote = FALSE, row.names = FALSE)
        list(file = meta(path), values = read.delim(path, sep = sep, check.names = FALSE))
    }
    a <- write_reference("nCounter", ".csv")
    b <- write_reference("RNAseq.V2", ".txt")
    expect_false(identical(a$values, b$values))

    cohort_paths <- function(id) {
        prefix <- if (id == "tcga") "TCGABRCAsub" else "OSLO2EMIT0"
        expression <- if (id == "tcga") "_GEX_RawCounts.txt" else "_GEX_log2.FPKM.txt"
        filenames <- paste0(prefix, c(expression, "_anno.txt", "_clinical.txt"))
        vapply(setNames(filenames, c("GEX", "anno", "clinic")), function(f)
            system.file("RshinyTest", f, package = "BreastSubtypeR", mustWork = TRUE),
            character(1))
    }
    # Read the shipped tables independently of the installed server.
    direct <- lapply(c(oslo = "oslo", tcga = "tcga"), function(id) {
        paths <- cohort_paths(id)
        gex <- as.matrix(read.delim(paths["GEX"], row.names = 1L, check.names = FALSE))
        storage.mode(gex) <- "double"
        anno <- read.delim(paths["anno"], check.names = FALSE)
        clinic <- read.delim(paths["clinic"], check.names = FALSE)
        rownames(anno) <- anno$probe
        rownames(clinic) <- clinic$PatientID
        se <- SummarizedExperiment::SummarizedExperiment(list(counts = gex),
            rowData = S4Vectors::DataFrame(anno), colData = S4Vectors::DataFrame(clinic))
        BreastSubtypeR::Mapping(se, RawCounts = id == "tcga", impute = FALSE,
                                verbose = FALSE)$se_NC
    })
    parker <- BreastSubtypeR::BS_parker
    captured <- new.env()
    testthat::local_mocked_bindings(BS_parker = function(...) {
        captured$args <- list(...)
        captured$result <- do.call(parker, captured$args)
        captured$result
    }, .package = "BreastSubtypeR")

    shiny::testServer(app$server, {
        tick <- 0L
        next_click <- function() {
            tick <<- tick + 1L
            tick
        }
        example <- function() session$setInputs(load_example_step1 = next_click())
        map <- function() session$setInputs(map = next_click())
        set_controls <- function(values) {
            ns <- calibration_ns()
            names(values) <- vapply(names(values), ns, character(1))
            for (id in names(values))
                expect_match(output$pam50_calibration$html, paste0('id="', id, '"'), fixed = TRUE)
            do.call(session$setInputs, values)
        }
        custom <- function(ref) set_controls(list(calibration = "External",
                                                   external = "Given.mdns", medians = ref$file))
        upload <- function(id) {
            records <- lapply(cohort_paths(id), meta)
            for (id_name in names(records)) {
                expect_match(output[[paste0("upload_", id_name)]]$html,
                             paste0('id="', upload_id(id_name), '"'), fixed = TRUE)
            }
            names(records) <- vapply(names(records), upload_id, character(1))
            do.call(session$setInputs, records)
            session$setInputs(is_raw_counts = if (id == "tcga") "raw" else "norm")
            expect_identical(reactive_files$source, "upload")
            expect_identical(raw_mode(), id == "tcga")
        }
        blocked <- function() {
            expect_null(current_prediction())
            for (kind in c("calls", "full")) {
                session$setInputs(export_kind = kind)
                expect_error(output$download, "No current prediction")
            }
        }
        compare <- function(id, ref = NULL, calibration = "Internal", internal = "medianCtr",
                            external = NULL) {
            args <- list(se_obj = direct[[id]], hasClinical = FALSE, Subtype = FALSE)
            if (!is.null(ref)) {
                calibration <- "External"
                external <- "Given.mdns"
            }
            args$calibration <- calibration
            variant <- "parker.original"
            label <- calibration
            if (calibration == "Internal") {
                args$internal <- internal
                variant <- switch(internal, medianCtr = "parker.original", meanCtr = "genefu.scale",
                                  qCtr = "genefu.robust")
                label <- sprintf("Internal (%s)", internal)
            } else if (calibration == "External") {
                args$external <- external
                if (!is.null(ref)) args$medians <- ref$values
                label <- sprintf("External (%s)", if (is.null(ref)) external else "Custom (Given.mdns)")
            }
            header <- sprintf("# Method: PAM50 (variant: %s, calibration: %s) | Classes: 5 | ROR: OFF",
                              variant, label)
            captured$args <- captured$result <- NULL
            set.seed(274)
            session$setInputs(run = next_click())
            expect_identical(captured$args, args)
            set.seed(274)
            result <- do.call(parker, args)
            expect_identical(captured$result, result)
            pred <- current_prediction()
            expected <- app$.build_single_exports(result, FALSE)
            expect_identical(pred$calls, expected$calls)
            expect_identical(pred$full, expected$full)
            for (kind in c("calls", "full")) {
                session$setInputs(export_kind = kind)
                body <- readLines(output$download)
                expect_identical(body[2L], header)
                serialized <- capture.output(write.table(expected[[kind]], sep = "\t",
                                                         quote = FALSE, row.names = FALSE))
                expect_identical(body[!startsWith(body, "#")], serialized)
            }
        }
        reset <- function() {
            blocked()
            expect_null(reference_file())
            expect_identical(calibration_choice("calibration"), "Internal")
            expect_identical(calibration_choice("internal"), "medianCtr")
        }
        refused <- function() {
            captured$args <- NULL
            session$setInputs(run = next_click())
            expect_null(captured$args)
            expect_match(reactive_files$error, "Upload custom medians")
            blocked()
        }
        session$setInputs(BSmethod = "PAM50", k_subtypes = "5", hasClinical = FALSE)
        example()
        map()
        custom(a)
        compare("oslo", a)
        compare("oslo", a)
        set_controls(list(calibration = "None"))
        blocked()
        set_controls(list(calibration = "External"))
        expect_identical(reference_file(), a$file)
        map()
        blocked()
        compare("oslo", a)

        # The method can hide the upload. Its dataset still owns the file.
        session$setInputs(BSmethod = "AIMS")
        blocked()
        expect_identical(reference_file(), a$file)
        old_medians <- calibration_ns()("medians")
        upload("tcga")
        session$setInputs(BSmethod = "PAM50")
        reset()
        map()
        compare("tcga")
        # A late completion targets the old binding. Even a replay on the
        # current binding cannot claim the old temporary file.
        do.call(session$setInputs, setNames(list(a$file), old_medians))
        expect_false(is.null(current_prediction()))
        replay <- a$file
        replay$datapath <- file.path(dirname(replay$datapath), ".", basename(replay$datapath))
        set_controls(list(calibration = "External", external = "Given.mdns", medians = replay))
        refused()
        custom(b)
        compare("tcga", b)

        # Replacing a reference immediately blocks the previous export but
        # preserves mapping. A new file in this dataset is independently usable.
        session$setInputs(upload_started = calibration_ns()("medians"))
        expect_true(mapping_current())
        expect_null(reference_file())
        refused()
        custom(write_reference("RNAseq.V2", ".txt"))
        compare("tcga", b)
        old_medians <- calibration_ns()("medians")
        session$setInputs(upload_started = old_medians)
        example()
        reset()
        do.call(session$setInputs, setNames(list(write_reference("nCounter", ".csv")$file), old_medians))
        map()
        compare("oslo")

        # Reloading the same example also replaces the dataset.
        custom(write_reference("nCounter", ".csv"))
        example()
        reset()
        map()
        compare("oslo")
        for (choice in c("meanCtr", "qCtr")) {
            set_controls(list(internal = choice))
            blocked()
            compare("oslo", internal = choice)
        }
        set_controls(list(calibration = "None"))
        blocked()
        compare("oslo", calibration = "None")
        set_controls(list(calibration = "External", external = "nCounter"))
        compare("oslo", calibration = "External", external = "nCounter")
        upload("tcga")
        blocked()
        map()
        compare("tcga", calibration = "External", external = "nCounter")
        custom(write_reference("RNAseq.V2", ".txt"))
        session$setInputs(is_raw_counts = "norm")
        expect_false(raw_mode())
        reset()
        session$setInputs(is_raw_counts = "raw")
        expect_true(raw_mode())
        map()
        compare("tcga")

        # A failed/interrupted primary upload cannot leave the other cohort
        # runnable with the old reference. Other slots survive a single edit.
        custom(write_reference("RNAseq.V2", ".txt"))
        compare("tcga", b)
        old_clinic <- upload_id("clinic")
        session$setInputs(upload_started = old_clinic)
        reset()
        expect_null(uploads()$clinic)
        expect_false(is.null(uploads()$GEX))
        map()
        expect_false(mapping_current())
        expect_match(reactive_files$error, "upload expression")
        # Complete the single-slot edit, then replace all three uploaded tables.
        corrected <- tempfile(fileext = ".txt")
        expect_true(file.copy(cohort_paths("tcga")[["clinic"]], corrected))
        do.call(session$setInputs, setNames(list(meta(corrected)), old_clinic))
        map()
        custom(write_reference("RNAseq.V2", ".txt"))
        compare("tcga", b)
        upload("oslo")
        reset()
        map()
        compare("oslo")

        custom(write_reference("nCounter", ".csv"))
        compare("oslo", a)
        old_medians <- calibration_ns()("medians")
        missing_example <<- TRUE
        # Queue a reference completion in the same flush as the failed load.
        late <- write_reference("RNAseq.V2", ".txt")$file
        do.call(session$setInputs, c(list(load_example_step1 = next_click()), setNames(list(late), old_medians)))
        reset()
        expect_null(reactive_files$source)
        do.call(session$setInputs, setNames(list(meta(cohort_paths("tcga")[["clinic"]])), old_clinic))
        expect_null(uploads()$clinic)
        expect_null(reactive_files$source)
        missing_example <<- FALSE
        example()
        map()
        compare("oslo")
        set_controls(list(calibration = "External", external = "Given.mdns", medians = late))
        refused()
    })
})
