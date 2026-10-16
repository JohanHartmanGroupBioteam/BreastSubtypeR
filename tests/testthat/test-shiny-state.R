# State transitions and download callbacks of the installed Shiny app.
shiny_app_env <- function() {
    skip_if_not_installed("shiny")
    skip_if_not_installed("bslib")
    env <- new.env(parent = asNamespace("shiny"))
    env$Mapping <- BreastSubtypeR::Mapping
    sys.source(system.file("ShinyBreastSubtypeR/server.R",
                           package = "BreastSubtypeR"), env)
    env
}
shiny_defaults <- function(session, method = "PAM50") {
    session$setInputs(BSmethod = method, k_subtypes = "5",
                     is_raw_counts = "norm", hasClinical = FALSE,
                     export_kind = "calls", calibration = "Internal",
                     internal = "medianCtr", external = "Given.mdns",
                     s = "ER.v2", ratio = 54 / 64, iteration = 10)
    session$setInputs(load_example_step1 = 1)
    session$setInputs(map = 1)
}
shiny_upload <- function(name) {
    path <- system.file("RshinyTest", name, package = "BreastSubtypeR")
    data.frame(name = name, size = file.info(path)$size,
               type = "text/plain", datapath = path)
}

test_that("iteration control preserves its historical default, minimum and step", {
    env <- shiny_app_env()
    withr::local_dir(system.file("ShinyBreastSubtypeR", package = "BreastSubtypeR"))
    sys.source(system.file("ShinyBreastSubtypeR/ui.R",
                           package = "BreastSubtypeR"), env)
    html <- as.character(env$ui)
    control <- regmatches(html, regexpr('<input[^>]*id="iteration"[^>]*>', html))
    expect_length(control, 1L)
    expect_match(control, 'type="number"', fixed = TRUE)
    expect_match(control, 'value="100"', fixed = TRUE)
    expect_match(control, 'min="10"', fixed = TRUE)
    expect_match(control, 'step="10"', fixed = TRUE)
})

for (changed in c("GEX", "clinic", "anno")) {
    local({
        id <- changed
        test_that(paste(id, "replacement invalidates mapping, plots and both exports"), {
            env <- shiny_app_env()
            shiny::testServer(env$server, {
                shiny_defaults(session)
                session$setInputs(run = 1)
                expect_equal(nrow(current_prediction()$calls), 102L)
                expect_true(file.exists(output$download))
                expect_match(output$download_control$html, 'id="download"', fixed = TRUE)
                files <- list(GEX = "OSLO2EMIT0_GEX_log2.FPKM.txt",
                              clinic = "OSLO2EMIT0_clinical.txt",
                              anno = "OSLO2EMIT0_anno.txt")
                do.call(session$setInputs, setNames(list(shiny_upload(files[[id]])), id))
                expect_null(reactive_files$data_input)
                expect_null(reactive_files$prediction)
                expect_null(current_prediction())
                expect_null(output$plotSection)
                expect_null(output$download_control)
                expect_match(output$export_selector$html, "Run subtyping before downloading results.",
                             fixed = TRUE)
                expect_error(output$download, "No current prediction")
                session$setInputs(export_kind = "full")
                expect_error(output$download, "No current prediction")
                session$setInputs(map = 2)
                expect_null(reactive_files$data_input)
                expect_match(reactive_files$error, "upload expression, clinical and annotation")
                expect_false(session$isClosed())
            })
        })
    })
}

test_that("example reload invalidates AUTO exports and fresh prediction restores them", {
    env <- shiny_app_env()
    shiny::testServer(env$server, {
        shiny_defaults(session, "AUTO Mode")
        expect_null(output$download_control)
        expect_match(output$export_selector$html, "Run subtyping before downloading results.",
                     fixed = TRUE)
        set.seed(148); session$setInputs(run = 1)
        original <- current_prediction()$calls
        expect_equal(nrow(original), 102L)
        expect_true(file.exists(output$download))
        expect_match(output$download_control$html, 'id="download"', fixed = TRUE)
        expect_match(output$export_selector$html, 'value="calls"', fixed = TRUE)
        expect_false(grepl('value="full"', output$export_selector$html, fixed = TRUE))
        expect_match(output$export_selector$html, "Full metrics are unavailable for this method.",
                     fixed = TRUE)
        expect_match(output$plotSection$html, "multi_plot")
        session$setInputs(export_kind = "full")
        expect_null(output$download_control)
        expect_error(output$download, "does not provide this export")
        session$setInputs(export_kind = "calls", load_example_step1 = 2)
        expect_false(mapping_current())
        expect_null(output$plotSection)
        expect_null(output$download_control)
        expect_error(output$download, "No current prediction")
        session$setInputs(map = 2)
        set.seed(148); session$setInputs(run = 2)
        expect_identical(current_prediction()$calls, original)
        expect_match(output$plotSection$html, "multi_plot")
        expect_true(file.exists(output$download))
        expect_match(output$download_control$html, 'id="download"', fixed = TRUE)
    })
})

test_that("prediction parameters invalidate results while retaining mapped inputs", {
    env <- shiny_app_env()
    shiny::testServer(env$server, {
        shiny_defaults(session)
        mapped <- reactive_files$data_input
        changes <- list(list(k_subtypes = "4"), list(hasClinical = TRUE),
                        list(internal = "meanCtr"), list(calibration = "None"),
                        list(BSmethod = "cIHC.itr"), list(ratio = 1.5),
                        list(iteration = 11), list(BSmethod = "ssBC"), list(s = "ER"))
        for (i in seq_along(changes)) {
            session$setInputs(run = i)
            expect_false(is.null(current_prediction()))
            do.call(session$setInputs, changes[[i]])
            expect_null(current_prediction())
            expect_null(output$plotSection)
            expect_null(output$download_control)
            expect_error(output$download, "No current prediction")
            expect_identical(reactive_files$data_input, mapped)
            expect_true(mapping_current())
        }
        session$setInputs(BSmethod = "PAM50", calibration = "External", external = "nCounter")
        session$setInputs(run = 20)
        expect_false(is.null(current_prediction()))
        session$setInputs(external = "Given.mdns")
        expect_null(current_prediction())
        session$setInputs(medians = data.frame(name = "medians.txt", datapath = tempfile()))
        expect_null(current_prediction())
        expect_identical(reactive_files$data_input, mapped)
    })
})

test_that("export choice and inactive method parameters preserve a valid prediction", {
    env <- shiny_app_env()
    shiny::testServer(env$server, {
        shiny_defaults(session)
        session$setInputs(run = 1)
        p <- current_prediction()
        expect_match(output$export_selector$html, 'value="calls"', fixed = TRUE)
        expect_match(output$export_selector$html, 'value="full"', fixed = TRUE)
        session$setInputs(export_kind = "full", ratio = 2, iteration = 20, s = "TN")
        expect_identical(current_prediction(), p)
        expect_match(output$download_control$html, 'id="download"', fixed = TRUE)
        full <- read.delim(output$download, comment.char = "#", check.names = FALSE)
        expect_equal(full, p$full, ignore_attr = TRUE)
        expect_match(readLines(output$download)[2], "Method: PAM50")
        expect_match(readLines(output$download)[3], "Export: full")
    })
})

test_that("data type changes and explicit remapping revoke predictions", {
    env <- shiny_app_env()
    shiny::testServer(env$server, {
        shiny_defaults(session)
        session$setInputs(run = 1)
        session$setInputs(map = 2)
        expect_true(mapping_current())
        expect_null(current_prediction())
        expect_error(output$download, "No current prediction")
        session$setInputs(run = 2)
        session$setInputs(is_raw_counts = "raw")
        expect_false(mapping_current())
        expect_null(current_prediction())
        expect_error(output$download, "No current prediction")
    })
})

test_that("provenance rejects retained records even if invalidation is bypassed", {
    env <- shiny_app_env()
    shiny::testServer(env$server, {
        shiny_defaults(session)
        session$setInputs(run = 1)
        old <- reactive_files$prediction
        session$setInputs(k_subtypes = "4")
        reactive_files$prediction <- old
        expect_null(current_prediction())
        expect_null(output$download_control)
        expect_error(export_result("calls"), "No current prediction")
        expect_error(export_result("full"), "No current prediction")
        expect_error(output$download, "No current prediction")
        session$setInputs(k_subtypes = "5", load_example_step1 = 2, map = 2)
        reactive_files$prediction <- old
        expect_null(current_prediction())
        expect_error(output$download, "No current prediction")
    })
})

test_that("bad ratios are visible recoverable errors and valid values above one work", {
    env <- shiny_app_env()
    shiny::testServer(env$server, {
        shiny_defaults(session, "cIHC.itr")
        session$setInputs(run = 1)
        mapped <- reactive_files$data_input
        bad <- list(NULL, numeric(0), "bad", "1", NA_real_, NaN, Inf, -Inf, 0, -1, c(1, 2))
        for (i in seq_along(bad)) {
            session$setInputs(ratio = bad[[i]])
            session$setInputs(run = i + 1)
            expect_match(reactive_files$error, "finite positive number")
            expect_match(output$analysis_error$html, "finite positive number")
            expect_false(session$isClosed())
            expect_null(current_prediction())
            expect_null(output$download_control)
            expect_error(output$download, "No current prediction")
            expect_identical(reactive_files$data_input, mapped)
        }
        session$setInputs(ratio = 1.5)
        session$setInputs(run = 100)
        expect_null(reactive_files$error)
        expect_null(output$analysis_error)
        expect_false(session$isClosed())
        p <- current_prediction()
        expect_equal(p$parameters$ratio, 1.5)
        expect_identical(p$calls$PatientID, colnames(mapped$se_NC))
        expect_match(output$plotSection$html, "pie1")
        expect_match(output$plotSection$html, "heat2")
        expect_match(output$download_control$html, 'id="download"', fixed = TRUE)
        session$setInputs(export_kind = "full")
        expect_true(file.exists(output$download))
    })
})

test_that("invalid iterations and predictor errors discard previous results", {
    env <- shiny_app_env()
    shiny::testServer(env$server, {
        shiny_defaults(session, "cIHC.itr")
        session$setInputs(run = 1)
        mapped <- reactive_files$data_input
        bad <- list(NULL, numeric(0), "10", NA_real_, NaN, Inf, -Inf,
                    9, 0, -1, 10.5, c(10, 100), .Machine$integer.max + 1)
        for (i in seq_along(bad)) {
            session$setInputs(iteration = bad[[i]])
            session$setInputs(run = i + 1)
            expect_match(reactive_files$error, "whole number from 10 to 2147483647", fixed = TRUE)
            expect_match(output$analysis_error$html, "whole number from 10 to 2147483647", fixed = TRUE)
            expect_false(session$isClosed())
            expect_null(output$plotSection)
            expect_null(output$download_control)
            expect_error(output$download, "No current prediction")
            expect_identical(reactive_files$data_input, mapped)
            expect_true(mapping_current())
        }
        session$setInputs(BSmethod = "PAM50", calibration = "Internal", internal = "medianCtr")
        session$setInputs(run = 50)
        expect_false(is.null(current_prediction()))
        session$setInputs(internal = "invalid-calibration")
        session$setInputs(run = 51)
        expect_match(reactive_files$error, "Invalid internal calibration")
        expect_null(current_prediction())
        expect_false(session$isClosed())
        expect_error(output$download, "No current prediction")
        session$setInputs(internal = "medianCtr")
        session$setInputs(run = 52)
        expect_false(is.null(current_prediction()))
        expect_match(output$download_control$html, 'id="download"', fixed = TRUE)
    })
})

test_that("valid iteration counts retain direct predictor calls and numeric metrics", {
    env <- shiny_app_env()
    shiny::testServer(env$server, {
        shiny_defaults(session, "cIHC.itr")
        mapped <- reactive_files$data_input
        for (iterations in c(10, 11, 100)) {
            session$setInputs(iteration = iterations)
            set.seed(148)
            expected <- BreastSubtypeR::BS_cIHC.itr(
                se_obj = mapped$se_NC, iteration = iterations, ratio = 54 / 64,
                Subtype = FALSE, hasClinical = FALSE)
            set.seed(148); session$setInputs(run = iterations)
            p <- current_prediction()
            expect_null(reactive_files$error)
            expect_false(session$isClosed())
            expect_identical(reactive_files$data_input, mapped)
            expect_identical(p$parameters$iteration, iterations)
            expect_identical(p$calls$PatientID, expected$BS.all$PatientID)
            expect_identical(p$calls$Call_5class, expected$BS.all$BS)
            scores <- as.data.frame(expected$score.ROR)
            numeric_columns <- names(scores)[vapply(scores, is.numeric, logical(1))]
            expect_gt(length(numeric_columns), 0L)
            expect_true(all(numeric_columns %in% names(p$full)))
            ids <- if ("PatientID" %in% names(scores)) scores$PatientID else rownames(scores)
            for (column in numeric_columns) {
                expect_identical(p$full[[column]], scores[[column]][match(p$full$PatientID, ids)])
            }
            expect_match(output$download_control$html, 'id="download"', fixed = TRUE)
            session$setInputs(export_kind = "full")
            exported <- read.delim(output$download, comment.char = "#", check.names = FALSE)
            expect_equal(exported, p$full, ignore_attr = TRUE)
        }
    })
})

test_that("actual upload tables map and replace the example without sharing sessions", {
    env <- shiny_app_env()
    shiny::testServer(env$server, {
        shiny_defaults(session)
        session$setInputs(GEX = shiny_upload("OSLO2EMIT0_GEX_log2.FPKM.txt"),
                         clinic = shiny_upload("OSLO2EMIT0_clinical.txt"),
                         anno = shiny_upload("OSLO2EMIT0_anno.txt"))
        session$setInputs(map = 2)
        expect_true(mapping_current())
        expect_identical(reactive_files$source, "upload")
        session$setInputs(run = 1)
        expect_equal(nrow(current_prediction()$calls), 102L)
        session$setInputs(load_example_step1 = 2)
        expect_identical(reactive_files$source, "example")
        expect_null(current_prediction())
        session$setInputs(map = 3)
        session$setInputs(run = 2)
        expect_equal(nrow(current_prediction()$calls), 102L)
    })
    shiny::testServer(env$server, {
        session$setInputs(BSmethod = "PAM50", k_subtypes = "5")
        expect_null(reactive_files$data_input)
        expect_null(current_prediction())
        expect_error(output$download, "No current prediction")
    })
})


test_that("a failed rerun with unchanged parameters revokes both export callbacks", {
    env <- shiny_app_env()
    shiny::testServer(env$server, {
        shiny_defaults(session)
        session$setInputs(run = 1)
        expect_false(is.null(current_prediction()))
        fail_rerun <- function() {
            testthat::local_mocked_bindings(
                BS_parker = function(...) stop("Deterministic predictor failure"),
                .package = "BreastSubtypeR")
            session$setInputs(run = 2)
        }
        fail_rerun()
        expect_match(reactive_files$error, "Deterministic predictor failure")
        expect_null(current_prediction())
        expect_null(output$plotSection)
        expect_null(output$download_control)
        expect_true(mapping_current())
        expect_false(session$isClosed())
        expect_error(output$download, "No current prediction")
        session$setInputs(export_kind = "full")
        expect_error(output$download, "No current prediction")
    })
})

test_that("full metrics never fall back to calls when unavailable", {
    env <- shiny_app_env()
    shiny::testServer(env$server, {
        shiny_defaults(session)
        session$setInputs(run = 1)
        p <- reactive_files$prediction
        p$full <- NULL
        reactive_files$prediction <- p
        session$flushReact()
        expect_true(file.exists(output$download))
        expect_match(output$download_control$html, 'id="download"', fixed = TRUE)
        expect_false(grepl('value="full"', output$export_selector$html, fixed = TRUE))
        expect_match(output$export_selector$html, "Full metrics are unavailable for this method.",
                     fixed = TRUE)
        session$setInputs(export_kind = "full")
        expect_null(output$download_control)
        expect_error(output$download, "does not provide this export")
        expect_error(export_result("full"), "does not provide this export")
        session$setInputs(export_kind = "calls")
        expect_match(output$download_control$html, 'id="download"', fixed = TRUE)
        expect_true(file.exists(output$download))
    })
})


test_that("mapping errors stay connected and cannot fall back to previous data", {
    env <- shiny_app_env()
    shiny::testServer(env$server, {
        shiny_defaults(session)
        session$setInputs(run = 1)
        bad_file <- tempfile(fileext = ".txt")
        writeLines(c("wrong_column", "sample"), bad_file)
        bad_upload <- data.frame(name = "bad-clinical.txt", size = file.info(bad_file)$size,
                                 type = "text/plain", datapath = bad_file)
        session$setInputs(GEX = shiny_upload("OSLO2EMIT0_GEX_log2.FPKM.txt"),
                         clinic = bad_upload, anno = shiny_upload("OSLO2EMIT0_anno.txt"))
        session$setInputs(map = 2)
        expect_false(mapping_current())
        expect_match(output$analysis_error$html, "PatientID")
        expect_false(session$isClosed())
        expect_error(output$download, "No current prediction")
        session$setInputs(clinic = shiny_upload("OSLO2EMIT0_clinical.txt"))
        session$setInputs(map = 3)
        session$setInputs(run = 2)
        expect_true(mapping_current())
        expect_equal(nrow(current_prediction()$calls), 102L)
        expect_null(reactive_files$error)
    })
})

test_that("replacing custom medians invalidates predictions without remapping", {
    env <- shiny_app_env()
    shiny::testServer(env$server, {
        shiny_defaults(session)
        data("BreastSubtypeRobj", package = "BreastSubtypeR")
        medians <- BreastSubtypeRobj$medians[, c("X", "nCounter")]
        names(medians)[2] <- "Given.mdns"
        path <- tempfile(fileext = ".txt")
        write.table(medians, path, sep = "\t", quote = FALSE, row.names = FALSE)
        upload <- data.frame(name = "medians.txt", size = file.info(path)$size,
                             type = "text/plain", datapath = path)
        session$setInputs(calibration = "External", external = "Given.mdns", medians = upload)
        session$setInputs(run = 1)
        expect_false(is.null(current_prediction()))
        mapped <- reactive_files$data_input
        upload$name <- "replacement-medians.txt"
        session$setInputs(medians = upload)
        expect_null(current_prediction())
        expect_error(output$download, "No current prediction")
        expect_identical(reactive_files$data_input, mapped)
        session$setInputs(run = 2)
        expect_false(is.null(current_prediction()))
    })
})

test_that("a second live session cannot clear or inherit the first prediction", {
    env <- shiny_app_env()
    check_sessions <- function() {
        first <- shiny::MockShinySession$new()
        second <- shiny::MockShinySession$new()
        on.exit({first$close(); second$close()})
        shiny::withReactiveDomain(first, env$server(first$input, first$output, first))
        shiny_defaults(first)
        first$setInputs(run = 1)
        download_for <- function(s) {
            shiny::withReactiveDomain(s, shiny::isolate(s$getOutput("download")))
        }
        download_control_for <- function(s) {
            shiny::withReactiveDomain(s, shiny::isolate(s$getOutput("download_control")))
        }
        original <- readLines(download_for(first))
        expect_match(download_control_for(first)$html, 'id="download"', fixed = TRUE)
        shiny::withReactiveDomain(second, env$server(second$input, second$output, second))
        second$setInputs(BSmethod = "PAM50", k_subtypes = "5")
        expect_null(download_control_for(second))
        expect_error(download_for(second), "No current prediction")
        expect_null(download_control_for(second))
        expect_match(download_control_for(first)$html, 'id="download"', fixed = TRUE)
        second$setInputs(load_example_step1 = 1)
        expect_error(download_for(second), "No current prediction")
        current <- readLines(download_for(first))
        expect_identical(current[-3L], original[-3L])
        expect_match(current[3L], "^# Export: calls \\| Date: .+ \\| R: ")
        expect_match(current[3L], paste0(" | R: ", getRversion()), fixed = TRUE)
        expect_false(first$isClosed())
        expect_false(second$isClosed())
    }
    check_sessions()
})
