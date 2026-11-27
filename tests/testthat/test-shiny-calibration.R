test_that("PAM50 custom medians belong to the current Shiny dataset", {
    app <- new.env(parent = asNamespace("shiny"))
    path <- system.file("ShinyBreastSubtypeR/server.R", package = "BreastSubtypeR")
    eval(parse(path, encoding = "UTF-8"), envir = app)
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

    # A fresh map for each cohort, independent of the live reactive inputs.
    direct <- lapply(c(oslo = "oslo", tcga = "tcga"), function(id) {
        x <- BreastSubtypeR:::.read_shiny_example(id)
        se <- SummarizedExperiment::SummarizedExperiment(list(counts = x$GEX),
            rowData = S4Vectors::DataFrame(x$anno), colData = S4Vectors::DataFrame(x$clinic))
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
        example <- function(id) session$setInputs(example_choice = id, confirm_example = next_click())
        map <- function() session$setInputs(map = next_click())
        set_controls <- function(values, ns = calibration_ns()) {
            names(values) <- vapply(names(values), ns, character(1))
            do.call(session$setInputs, values)
        }
        custom <- function(ref) set_controls(list(calibration = "External",
                                                   external = "Given.mdns", medians = ref$file))
        upload <- function(id) {
            spec <- BreastSubtypeR:::.shiny_examples[[id]]
            paths <- c(GEX = spec$expression, clinic = paste0(spec$prefix, "_clinical.txt"),
                       anno = paste0(spec$prefix, "_anno.txt"))
            records <- lapply(paths, function(f) meta(system.file("RshinyTest", f, package = "BreastSubtypeR")))
            names(records) <- vapply(names(records), upload_id, character(1))
            do.call(session$setInputs, records)
            session$setInputs(is_raw_counts = if (id == "tcga") "raw" else "norm")
        }
        blocked <- function() {
            expect_null(current_prediction())
            for (kind in c("calls", "full")) {
                session$setInputs(export_kind = kind)
                expect_error(output$download, "No current prediction")
            }
        }
        compare <- function(id, ref = NULL) {
            args <- list(se_obj = direct[[id]], hasClinical = FALSE, Subtype = FALSE)
            if (is.null(ref)) {
                args$calibration <- "Internal"
                args$internal <- "medianCtr"
            } else {
                args$calibration <- "External"
                args$external <- "Given.mdns"
                args$medians <- ref$values
            }
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
                expect_match(body[2L], if (is.null(ref)) "Internal (medianCtr)" else "Custom (Given.mdns)",
                             fixed = TRUE)
                serialized <- capture.output(write.table(expected[[kind]], sep = "\t",
                                                         quote = FALSE, row.names = FALSE))
                expect_identical(body[!startsWith(body, "#")], serialized)
            }
        }
        check_replacement <- function(id) {
            blocked()
            map()
            compare(id)
        }
        session$setInputs(BSmethod = "PAM50", k_subtypes = "5", hasClinical = FALSE)
        example("oslo")
        map()
        custom(a)
        compare("oslo", a)
        ns_a <- calibration_ns()
        previous <- current_prediction()
        session$setInputs(load_example = 1, example_choice = "tcga")
        expect_identical(current_prediction(), previous)
        expect_identical(reference_file(), a$file)
        # Hide the custom file before switching, as in the baseline trace.
        set_controls(list(calibration = "Internal"))
        example("tcga")
        check_replacement("tcga")
        set_controls(list(calibration = "External", external = "Given.mdns", medians = a$file), ns_a)
        expect_false(is.null(current_prediction()))
        # Even replaying the retired path against the current control cannot use it.
        replay <- a$file
        replay$datapath <- file.path(dirname(replay$datapath), ".", basename(replay$datapath))
        set_controls(list(calibration = "External", external = "Given.mdns", medians = replay))
        captured$args <- NULL
        session$setInputs(run = next_click())
        expect_null(captured$args)
        expect_match(reactive_files$error, "Upload custom medians")
        blocked()
        custom(b)
        compare("tcga", b)
        example("oslo")
        check_replacement("oslo")

        # A genuine re-upload receives a new temporary path, even for the same contents.
        a <- write_reference("nCounter", ".csv")
        custom(a)
        compare("oslo", a)
        example("oslo")
        check_replacement("oslo")
        a <- write_reference("nCounter", ".csv")
        custom(a)
        compare("oslo", a)
        upload("tcga")
        check_replacement("tcga")
        b <- write_reference("RNAseq.V2", ".txt")
        custom(b)
        compare("tcga", b)
        upload("oslo")
        check_replacement("oslo")
        a <- write_reference("nCounter", ".csv")
        custom(a)
        compare("oslo", a)
        example("tcga")
        check_replacement("tcga")
        b <- write_reference("RNAseq.V2", ".txt")
        custom(b)
        compare("tcga", b)
        upload("oslo")
        check_replacement("oslo")

        custom(write_reference("nCounter", ".csv"))
        session$setInputs(run = next_click())
        expect_false(is.null(current_prediction()))
        session$setInputs(upload_started = calibration_ns()("medians"))
        blocked()
        expect_true(mapping_current())
        b <- write_reference("RNAseq.V2", ".txt")
        custom(b)
        compare("oslo", b)
        map()
        compare("oslo", b)
        example("invalid")
        blocked()
        expect_null(reactive_files$source)
        expect_type(reactive_files$error, "character")
        example("oslo")
        map()
        compare("oslo")
    })
})
