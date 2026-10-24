# State transitions and download callbacks of the installed Shiny app.
shiny_app_env <- function() {
    skip_if_not_installed("shiny")
    skip_if_not_installed("bslib")
    env <- new.env(parent = asNamespace("shiny"))
    server_path <- system.file("ShinyBreastSubtypeR/server.R",
                               package = "BreastSubtypeR")
    eval(parse(file = server_path, encoding = "UTF-8", keep.source = FALSE),
         envir = env)
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

test_that("mapping resolves from the installed namespace without injection", {
    skip_if_not_installed("shiny")
    skip_if_not_installed("bslib")
    ## test_package() attaches BreastSubtypeR. Use a separate vanilla process
    ## so its search path cannot mask an unqualified app function.
    script <- tempfile(fileext = ".R")
    log <- tempfile(fileext = ".log")
    on.exit(unlink(c(script, log)), add = TRUE)
    writeLines(c(
        paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
        paste0("defaults <- ", paste(deparse(shiny_defaults), collapse = "\n")),
        'stopifnot(!"package:BreastSubtypeR" %in% search(),',
        '    !exists("Mapping", .GlobalEnv, inherits = TRUE))',
        'env <- new.env(parent = asNamespace("shiny"))',
        'sys.source(system.file("ShinyBreastSubtypeR/server.R",',
        '    package = "BreastSubtypeR"), env)',
        'stopifnot(!exists("Mapping", env, inherits = FALSE))',
        'shiny::testServer(env$server, {',
        '    defaults(session)',
        '    stopifnot(mapping_current(), is.null(reactive_files$error),',
        '        ncol(reactive_files$data_input$se_NC) == 102L)',
        '})'
    ), script)
    rscript <- file.path(R.home("bin"),
        if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
    status <- system2(rscript, c("--vanilla", shQuote(script)),
        stdout = log, stderr = log)
    expect_identical(status, 0L, info = paste(readLines(log), collapse = "\n"))
})

test_that("silent aborts clear stale results but preserve meaningful errors", {
    env <- shiny_app_env()
    shiny::testServer(env$server, {
        shiny_defaults(session)
        session$setInputs(run = 1)
        expect_false(is.null(current_prediction()))
        ## Inject conditions into the mapping/prediction error handler.
        ## This isolates Shiny semantics; it is not a normal upload workflow.
        silent <- tryCatch(shiny::req(NULL), error = identity)
        expect_s3_class(silent, "shiny.silent.error")
        report_error(silent)
        expect_null(reactive_files$error)
        expect_null(output$analysis_error)
        expect_null(current_prediction())
        expect_error(export_result("calls"), "No current prediction")

        validation <- tryCatch(shiny::validate(shiny::need(FALSE,
            "Annotation IDs do not overlap.")), error = identity)
        expect_s3_class(validation, "shiny.silent.error")
        report_error(validation)
        session$flushReact()
        expect_identical(reactive_files$error, "Annotation IDs do not overlap.")
        expect_match(output$analysis_error$html, "Annotation IDs do not overlap.",
            fixed = TRUE)
        expect_null(current_prediction())
        expect_error(export_result("calls"), "No current prediction")

        report_error(simpleError("Ordinary prediction error."))
        expect_identical(reactive_files$error, "Ordinary prediction error.")
        session$setInputs(run = 2)
        expect_null(reactive_files$error)
        expect_false(is.null(current_prediction()))
        expect_equal(nrow(export_result("calls")$calls), 102L)
    })
})

test_that("a req abort inside prediction clears the result and can recover", {
    env <- shiny_app_env()
    shiny::testServer(env$server, {
        shiny_defaults(session, method = "AIMS")
        session$setInputs(run = 1)
        previous <- current_prediction()$calls
        mapped <- reactive_files$data_input
        ## Deliberately corrupt internal mapped state to reach req(se_SSP).
        ## Successful Mapping() does not normally produce this state.
        missing <- mapped
        missing$se_SSP <- NULL
        missing$x_SSP <- NULL
        reactive_files$data_input <- missing
        session$setInputs(run = 2)
        expect_null(reactive_files$error)
        expect_null(output$analysis_error)
        expect_null(current_prediction())
        expect_error(output$download, "No current prediction")
        reactive_files$data_input <- mapped
        session$setInputs(run = 3)
        expect_null(reactive_files$error)
        expect_identical(current_prediction()$calls, previous)
        expect_true(file.exists(output$download))
    })
})

test_that("input coverage uses the shared report for the current mapped inputs", {
    env <- shiny_app_env()
    shiny::testServer(env$server, {
        expect_null(current_input_coverage())
        expect_null(output$input_coverage)
        shiny_defaults(session)
        mapped <- reactive_files$data_input
        expected <- BreastSubtypeR::InputCoverage(mapped)
        expect_identical(current_input_coverage()$data, expected)
        expect_null(current_input_coverage()$error)
        html <- output$input_coverage$html
        expect_match(html, "Input coverage", fixed = TRUE)
        expect_match(html, '<th scope="col">Method/model</th>', fixed = TRUE)
        expect_match(html, '<th scope="row">NC-based (PAM50)</th>', fixed = TRUE)
        expect_match(html, "Genes represented", fixed = TRUE)
        expect_match(html, "Genes absent", fixed = TRUE)
        expect_match(html, "Rules evaluable", fixed = TRUE)
        expect_match(html, "N/A", fixed = TRUE)
        for (i in seq_len(nrow(expected))) {
            expect_match(html, expected$model[i], fixed = TRUE)
            expect_match(html, sprintf("%s/%s", expected$genes_present[i],
                                      expected$genes_total[i]), fixed = TRUE)
        }
        expect_match(html, "NC-based methods share the 50-gene PAM50 subtype reference; AIMS and both SSPBC models are SSP-based. These rows describe input coverage, not which methods were selected or run.",
                     fixed = TRUE)
        expect_match(html, "Descriptive coverage of processed expression inputs. Coverage does not establish prediction reliability or change method eligibility.",
                     fixed = TRUE)
        expect_identical(reactive_files$data_input, mapped)
        expect_null(current_prediction())
        session$setInputs(BSmethod = "AUTO Mode", k_subtypes = "4")
        expect_identical(current_input_coverage()$data, expected)
        expect_identical(output$input_coverage$html, html)
    })
})

test_that("input coverage displays rule ranges across processed samples", {
    env <- shiny_app_env()
    shiny::testServer(env$server, {
        shiny_defaults(session)
        mapped <- reactive_files$data_input
        SummarizedExperiment::assay(mapped$se_SSP)[, 1L] <- NA_real_
        reactive_files$data_input <- mapped
        session$flushReact()
        expected <- BreastSubtypeR::InputCoverage(mapped)
        expect_identical(current_input_coverage()$data, expected)
        for (i in which(!is.na(expected$rules_total))) {
            expect_equal(expected$rules_evaluable_min[i], 0L)
            expect_gt(expected$rules_evaluable_max[i], 0L)
            range <- sprintf("%s-%s/%s per sample", expected$rules_evaluable_min[i],
                             expected$rules_evaluable_max[i], expected$rules_total[i])
            expect_match(output$input_coverage$html, range, fixed = TRUE)
        }
        expect_identical(reactive_files$data_input, mapped)
    })
})

test_that("unavailable coverage is informative without changing prediction state", {
    env <- shiny_app_env()
    shiny::testServer(env$server, {
        shiny_defaults(session)
        mapped <- reactive_files$data_input
        key <- reactive_files$mapped_key
        fail_coverage <- function() {
            testthat::local_mocked_bindings(
                InputCoverage = function(data_input) stop("Deterministic coverage failure"),
                .package = "BreastSubtypeR")
            reactive_files$mapped_key <- NULL
            session$flushReact()
            expect_null(current_input_coverage())
            reactive_files$mapped_key <- key
            session$flushReact()
            expect_true(mapping_current())
            expect_null(current_input_coverage()$data)
            expect_match(output$input_coverage$html,
                         "Input coverage is unavailable: Deterministic coverage failure", fixed = TRUE)
            expect_false(grepl("Genes represented", output$input_coverage$html, fixed = TRUE))
            expect_null(reactive_files$error)
            session$setInputs(run = 1)
            expect_false(is.null(current_prediction()))
            expect_identical(reactive_files$data_input, mapped)
            expect_match(output$download_control$html, 'id="download"', fixed = TRUE)
            expect_false(session$isClosed())
        }
        fail_coverage()
        reactive_files$mapped_key <- NULL
        session$flushReact()
        expect_null(current_input_coverage())
        reactive_files$mapped_key <- key
        session$flushReact()
        expect_identical(current_input_coverage()$data, BreastSubtypeR::InputCoverage(mapped))
    })
})

test_that("coverage provenance rejects retained mapping from an older generation", {
    env <- shiny_app_env()
    shiny::testServer(env$server, {
        shiny_defaults(session)
        mapped <- reactive_files$data_input
        key <- reactive_files$mapped_key
        expect_false(is.null(current_input_coverage()))
        session$setInputs(load_example_step1 = 2)
        reactive_files$data_input <- mapped
        reactive_files$mapped_key <- key
        session$flushReact()
        expect_false(mapping_current())
        expect_null(current_input_coverage())
        expect_null(output$input_coverage)
        session$setInputs(map = 2)
        expect_identical(current_input_coverage()$data,
                         BreastSubtypeR::InputCoverage(reactive_files$data_input))
    })
})

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
                expect_null(current_input_coverage())
                expect_null(output$input_coverage)
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
                expect_null(output$input_coverage)
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
        expect_null(output$input_coverage)
        expect_null(output$plotSection)
        expect_null(output$download_control)
        expect_error(output$download, "No current prediction")
        session$setInputs(map = 2)
        expect_identical(current_input_coverage()$data,
                         BreastSubtypeR::InputCoverage(reactive_files$data_input))
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
        coverage <- current_input_coverage()
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
            expect_identical(current_input_coverage(), coverage)
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
        expect_identical(current_input_coverage()$data,
                         BreastSubtypeR::InputCoverage(reactive_files$data_input))
        expect_null(current_prediction())
        expect_error(output$download, "No current prediction")
        session$setInputs(run = 2)
        session$setInputs(is_raw_counts = "raw")
        expect_false(mapping_current())
        expect_null(current_input_coverage())
        expect_null(output$input_coverage)
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
        expect_match(output$plotSection$html, "height:720px", fixed = TRUE)
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
        expect_identical(current_input_coverage()$data,
                         BreastSubtypeR::InputCoverage(reactive_files$data_input))
        expect_identical(reactive_files$source, "upload")
        session$setInputs(run = 1)
        expect_equal(nrow(current_prediction()$calls), 102L)
        session$setInputs(load_example_step1 = 2)
        expect_identical(reactive_files$source, "example")
        expect_null(output$input_coverage)
        expect_null(current_prediction())
        session$setInputs(map = 3)
        session$setInputs(run = 2)
        expect_equal(nrow(current_prediction()$calls), 102L)
    })
    shiny::testServer(env$server, {
        session$setInputs(BSmethod = "PAM50", k_subtypes = "5")
        expect_null(reactive_files$data_input)
        expect_null(current_input_coverage())
        expect_null(output$input_coverage)
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
        expect_null(current_input_coverage())
        expect_null(output$input_coverage)
        expect_match(output$analysis_error$html, "PatientID")
        expect_false(session$isClosed())
        expect_error(output$download, "No current prediction")
        session$setInputs(clinic = shiny_upload("OSLO2EMIT0_clinical.txt"))
        session$setInputs(map = 3)
        session$setInputs(run = 2)
        expect_true(mapping_current())
        expect_identical(current_input_coverage()$data,
                         BreastSubtypeR::InputCoverage(reactive_files$data_input))
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
        coverage_for <- function(s) {
            shiny::withReactiveDomain(s, shiny::isolate(s$getOutput("input_coverage")))
        }
        original <- readLines(download_for(first))
        original_coverage <- coverage_for(first)
        expect_match(original_coverage$html, "Input coverage", fixed = TRUE)
        expect_match(download_control_for(first)$html, 'id="download"', fixed = TRUE)
        shiny::withReactiveDomain(second, env$server(second$input, second$output, second))
        second$setInputs(BSmethod = "PAM50", k_subtypes = "5")
        expect_null(coverage_for(second))
        expect_null(download_control_for(second))
        expect_error(download_for(second), "No current prediction")
        expect_null(download_control_for(second))
        expect_match(download_control_for(first)$html, 'id="download"', fixed = TRUE)
        second$setInputs(load_example_step1 = 1)
        expect_error(download_for(second), "No current prediction")
        expect_null(coverage_for(second))
        expect_identical(coverage_for(first), original_coverage)
        second$setInputs(is_raw_counts = "norm", map = 1)
        expect_match(coverage_for(second)$html, "Input coverage", fixed = TRUE)
        expect_identical(coverage_for(first), original_coverage)
        second$setInputs(is_raw_counts = "raw")
        expect_null(coverage_for(second))
        expect_identical(coverage_for(first), original_coverage)
        current <- readLines(download_for(first))
        expect_identical(current[-3L], original[-3L])
        expect_match(current[3L], "^# Export: calls \\| Date: .+ \\| R: ")
        expect_match(current[3L], paste0(" | R: ", getRversion()), fixed = TRUE)
        expect_false(first$isClosed())
        expect_false(second$isClosed())
    }
    check_sessions()
})

test_that("the presentation preserves controls, defaults and conditional settings", {
    env <- shiny_app_env()
    app_dir <- system.file("ShinyBreastSubtypeR", package = "BreastSubtypeR")
    withr::local_dir(app_dir)
    sys.source(file.path(app_dir, "ui.R"), env)
    html <- as.character(env$ui)
    for (id in c("load_example_step1", "about", "GEX", "is_raw_counts", "clinic",
                 "anno", "map", "BSmethod", "k_subtypes", "hasClinical", "calibration",
                 "external", "medians", "internal", "iteration", "ratio", "s", "run",
                 "auto_chip", "gex_help", "clin_help", "anno_help", "input_coverage",
                 "auto_preflight", "method_help", "calib_help", "analysis_error",
                 "plotSection", "export_selector", "download_control")) {
        expect_equal(length(gregexpr(paste0('id="', id, '"'), html, fixed = TRUE)[[1L]]),
                     1L, info = id)
        expect_match(html, paste0('id="', id, '"'), fixed = TRUE)
    }
    select_contract <- list(
        BSmethod = c("AUTO Mode", "PAM50", "cIHC", "cIHC.itr", "PCAPAM50", "ssBC", "AIMS", "sspbc"),
        calibration = c("None", "External", "Internal"),
        external = c("Given.mdns", "nCounter", "RNAseq.Freeze.20120907", "totalRNA.FFPE.20151111",
                     "RNAseq.V2", "RNAseq.V1", "GC.4x44Kcustom", "Agilent_244K",
                     "commercial_1x44k_postMeanCollapse_WashU",
                     "commercial_4x44k_postMeanCollapse_WashU_v2",
                     "htp1.5_WU_update", "arrayTrain_postMeanCollapse"),
        internal = c("medianCtr", "meanCtr", "qCtr"), s = c("ER", "ER.v2", "TN", "TN.v2"))
    defaults <- c(BSmethod = "AUTO Mode", calibration = "Internal", external = "Given.mdns",
                  internal = "medianCtr", s = "ER.v2")
    for (id in names(select_contract)) {
        block <- regmatches(html, regexpr(paste0('(?s)<select[^>]*id="', id,
                                                 '"[^>]*>.*?</select>'), html, perl = TRUE))
        options <- regmatches(block, gregexpr("<option[^>]*>", block))[[1L]]
        expect_identical(sub('.*value="([^"]*)".*', "\\1", options), select_contract[[id]])
        expect_match(options[grepl(" selected", options, fixed = TRUE)],
                     paste0('value="', defaults[[id]], '"'), fixed = TRUE)
    }
    for (id in c("is_raw_counts", "k_subtypes")) {
        controls <- regmatches(html, gregexpr(paste0('<input[^>]*name="', id, '"[^>]*>'), html))[[1L]]
        expect_identical(sub('.*value="([^"]*)".*', "\\1", controls),
                         if (id == "is_raw_counts") c("norm", "raw") else c("5", "4"))
        expect_match(controls[grepl("checked=", controls, fixed = TRUE)],
                     if (id == "is_raw_counts") 'value="norm"' else 'value="5"', fixed = TRUE)
    }
    clinical <- regmatches(html, regexpr('<input[^>]*id="hasClinical"[^>]*>', html))
    expect_false(grepl("checked", clinical, fixed = TRUE))
    ratio <- regmatches(html, regexpr('<input[^>]*id="ratio"[^>]*>', html))
    expect_match(ratio, 'value="0.84375"', fixed = TRUE)
    expect_match(ratio, 'step="any"', fixed = TRUE)
    for (id in c("GEX", "clinic", "anno", "medians")) {
        upload <- regmatches(html, regexpr(paste0('<input[^>]*id="', id, '"[^>]*>'), html))
        expect_match(upload, 'type="file"', fixed = TRUE)
        expect_match(upload, '.csv,.txt"', fixed = TRUE)
    }
    conditions <- regmatches(html, gregexpr('data-display-if="[^"]*"', html))[[1L]]
    conditions <- gsub("&#39;", "'", sub('data-display-if="([^"]*)"', "\\1", conditions), fixed = TRUE)
    expect_setequal(conditions, c(
        "input.BSmethod == 'PAM50' || input.BSmethod == 'cIHC' || input.BSmethod == 'cIHC.itr' || input.BSmethod == 'PCAPAM50' || input.BSmethod == 'ssBC'",
        "input.BSmethod == 'PAM50'", "input.calibration == 'External'",
        "input.external == 'Given.mdns'", "input.calibration == 'Internal'",
        "input.BSmethod == 'cIHC.itr'", "input.BSmethod == 'ssBC'"))
})

test_that("the interface identifies its workflow and provides accessible scientific help", {
    env <- shiny_app_env()
    app_dir <- system.file("ShinyBreastSubtypeR", package = "BreastSubtypeR")
    withr::local_dir(app_dir)
    sys.source(file.path(app_dir, "ui.R"), env)
    html <- as.character(env$ui)
    for (text in c("iBreastSubtypeR", "For research use", "About &amp; Citation",
                   "The application runs locally in your R session.", "1. Data and preprocessing",
                   "2. Analysis settings", "Results and downloads", "tab-separated text (TSV)",
                   "AUTO - select methods for this cohort", "Required: PatientID.",
                   "Required: probe and ENTREZID.", "gene Length in base pairs.")) {
        expect_match(html, text, fixed = TRUE)
    }
    expect_match(html, 'alt="BreastSubtypeR logo"', fixed = TRUE)
    expect_match(html, 'aria-label="Workflow"', fixed = TRUE)
    for (help in c("Expression format and processing", "Clinical fields and research ROR",
                   "Annotation format", "Method guide")) {
        expect_match(html, paste0("<summary>", help, "</summary>"), fixed = TRUE)
    }
    expect_false(grepl("reduces bias", html, fixed = TRUE))
    expect_false(grepl("Bioconductor-ready", html, fixed = TRUE))
    expect_match(html, paste("selects methods according to implemented receptor-annotation",
                            "and subgroup-size rules."), fixed = TRUE)
    expect_match(html, paste("It does not establish complete model-input adequacy",
                            "or clinical validity."), fixed = TRUE)
    expect_false(grepl("cohort composition and method requirements", html, fixed = TRUE))
    expect_false(grepl("Editor's Choice", html, fixed = TRUE))
})

test_that("the citation modal retains the publication and links its distinction", {
    env <- shiny_app_env()
    modal_html <- NULL
    env$showModal <- function(ui, ...) modal_html <<- as.character(ui)
    shiny::testServer(env$server, {
        session$setInputs(BSmethod = "AUTO Mode", k_subtypes = "5", about = 1)
    })
    expect_match(modal_html, "Please cite:", fixed = TRUE)
    expect_match(modal_html, "Yang Q., Hartman J., Sifakis E.G.", fixed = TRUE)
    expect_match(modal_html, ">Editor's Choice (2025)</a> - ", fixed = TRUE)
    expect_match(modal_html, paste0(
        "href='https://academic.oup.com/nargab/pages/editors-choice-articles",
        "?login=false#:~:text=BreastSubtypeR'"), fixed = TRUE)
    expect_match(modal_html, "https://doi.org/10.1093/nargab/lqaf131", fixed = TRUE)
    expect_match(modal_html, ">NAR Genomics and Bioinformatics</a>.</p>", fixed = TRUE)
})

test_that("AUTO awaits current preprocessing without a false preflight error", {
    env <- shiny_app_env()
    shiny::testServer(env$server, {
        session$setInputs(BSmethod = "AUTO Mode", k_subtypes = "5",
                          is_raw_counts = "norm")
        expect_false(mapping_current())
        expect_false(auto_requirements()$ready)
        expect_match(output$auto_chip$html, 'class="chip auto pending"', fixed = TRUE)
        expect_match(output$auto_chip$html, "AUTO: awaiting preprocessing", fixed = TRUE)
        expect_match(output$auto_chip$html, "Run Step 1", fixed = TRUE)
        expect_false(grepl("blocked|alert-danger|auto ready", output$auto_chip$html))
        expect_null(output$auto_preflight)
        expect_null(output$method_help)
        expect_null(output$analysis_error)

        shiny_defaults(session, "AUTO Mode")
        expect_true(mapping_current())
        expect_true(auto_requirements()$ready)
        expect_match(output$auto_chip$html, 'class="chip auto ready"', fixed = TRUE)
        expect_match(output$auto_preflight$html, "AUTO preflight: ready", fixed = TRUE)
        expect_match(output$auto_preflight$html, "ER+: 84", fixed = TRUE)
        expect_match(output$auto_preflight$html, "ER-: 18", fixed = TRUE)
        expect_match(output$auto_preflight$html, "HER2+: 8", fixed = TRUE)
        expect_match(output$auto_preflight$html, "HER2-: 94", fixed = TRUE)
        expect_null(output$method_help)

        session$setInputs(clinic = shiny_upload("OSLO2EMIT0_clinical.txt"))
        expect_false(mapping_current())
        expect_match(output$auto_chip$html, "AUTO: awaiting preprocessing", fixed = TRUE)
        expect_false(grepl("ER\\+|HER2|auto ready", output$auto_chip$html))
        expect_null(output$auto_preflight)
        expect_null(current_input_coverage())
        expect_null(current_prediction())
        session$setInputs(load_example_step1 = 2)
        session$setInputs(map = 2)
        expect_true(mapping_current())
        expect_match(output$auto_preflight$html, "AUTO preflight: ready", fixed = TRUE)
    })
})

test_that("mapped AUTO annotation failures retain blocking and recover after replacement", {
    env <- shiny_app_env()
    clinical <- shiny_upload("OSLO2EMIT0_clinical.txt")
    ph <- read.delim(clinical$datapath, check.names = FALSE)
    ph$HER2 <- NULL
    path <- tempfile(fileext = ".txt")
    withr::defer(unlink(path))
    write.table(ph, path, sep = "\t", quote = FALSE, row.names = FALSE)
    incomplete <- clinical
    incomplete$name <- "clinical-without-HER2.txt"
    incomplete$datapath <- path
    incomplete$size <- file.info(path)$size
    shiny::testServer(env$server, {
        shiny_defaults(session, "AUTO Mode")
        session$setInputs(GEX = shiny_upload("OSLO2EMIT0_GEX_log2.FPKM.txt"),
                          clinic = incomplete, anno = shiny_upload("OSLO2EMIT0_anno.txt"))
        expect_match(output$auto_chip$html, "AUTO: awaiting preprocessing", fixed = TRUE)
        expect_null(output$auto_preflight)
        session$setInputs(map = 2)
        expect_true(mapping_current())
        expect_null(reactive_files$error)
        expect_false(auto_requirements()$ready)
        expect_match(output$auto_chip$html, 'class="chip auto blocked"', fixed = TRUE)
        expect_match(output$auto_preflight$html, "AUTO preflight: action needed", fixed = TRUE)
        expect_match(output$auto_preflight$html, "missing: HER2", fixed = TRUE)
        session$setInputs(run = 1)
        expect_match(reactive_files$error, "AUTO needs ER and HER2", fixed = TRUE)
        expect_null(current_prediction())

        session$setInputs(clinic = clinical)
        expect_false(mapping_current())
        expect_match(output$auto_chip$html, "AUTO: awaiting preprocessing", fixed = TRUE)
        expect_null(output$auto_preflight)
        session$setInputs(map = 3)
        expect_true(mapping_current())
        expect_true(auto_requirements()$ready)
        expect_null(reactive_files$error)
        expect_match(output$auto_preflight$html, "AUTO preflight: ready", fixed = TRUE)
    })
})

test_that("pending AUTO presentation preserves mapping errors and recovery", {
    env <- shiny_app_env()
    shiny::testServer(env$server, {
        shiny_defaults(session, "AUTO Mode")
        session$setInputs(GEX = shiny_upload("OSLO2EMIT0_GEX_log2.FPKM.txt"),
                          clinic = shiny_upload("OSLO2EMIT0_clinical.txt"),
                          anno = shiny_upload("OSLO2EMIT0_anno.txt"))
        session$setInputs(is_raw_counts = "raw")
        expect_null(output$auto_preflight)
        session$setInputs(map = 2)
        expect_false(mapping_current())
        expect_match(reactive_files$error, "annotation lacks", fixed = TRUE)
        expect_match(reactive_files$error, "Length", fixed = TRUE)
        expect_match(output$analysis_error$html, "Analysis needs attention", fixed = TRUE)
        expect_match(output$analysis_error$html, "annotation lacks", fixed = TRUE)
        expect_match(output$auto_chip$html, "AUTO: awaiting preprocessing", fixed = TRUE)
        expect_null(output$auto_preflight)
        expect_null(current_prediction())
        session$setInputs(is_raw_counts = "norm")
        session$setInputs(map = 3)
        expect_true(mapping_current())
        expect_true(auto_requirements()$ready)
        expect_null(reactive_files$error)
        expect_null(output$analysis_error)
        expect_match(output$auto_preflight$html, "AUTO preflight: ready", fixed = TRUE)
    })
})

test_that("removing generic AUTO help preserves other classifier help", {
    env <- shiny_app_env()
    shiny::testServer(env$server, {
        for (method in c("PAM50", "AIMS", "sspbc")) {
            session$setInputs(BSmethod = method, k_subtypes = "5")
            expect_match(output$method_help$html,
                         if (method == "sspbc") "SSPBC" else method, fixed = TRUE)
            expect_null(output$auto_chip)
            expect_null(output$auto_preflight)
        }
        session$setInputs(BSmethod = "AUTO Mode")
        expect_null(output$method_help)
    })
})

test_that("scroll hints accompany wide figures only and are narrow-screen scoped", {
    env <- shiny_app_env()
    app_dir <- system.file("ShinyBreastSubtypeR", package = "BreastSubtypeR")
    withr::local_dir(app_dir)
    sys.source(file.path(app_dir, "ui.R"), env)
    html <- as.character(env$ui)
    css <- htmltools::renderTags(env$ui)$head
    expect_match(css, ".plot-scroll-hint { display: none;", fixed = TRUE)
    expect_match(css, "@media (max-width: 575px) {\n      .plot-scroll-hint { display: block; }",
                 fixed = TRUE)
    expect_match(css, "#multi_plot, #heat2 { min-width: 760px; }", fixed = TRUE)
    expect_false(grepl("Swipe horizontally", html, fixed = TRUE))
    shiny::testServer(env$server, {
        shiny_defaults(session)
        expect_null(output$plotSection)
        session$setInputs(run = 1)
        plot_html <- output$plotSection$html
        expect_equal(length(gregexpr("Swipe horizontally", plot_html, fixed = TRUE)[[1L]]), 1L)
        expect_match(plot_html, "Expression and subtype calls", fixed = TRUE)
        expect_match(plot_html, "Swipe horizontally to view the full plot.", fixed = TRUE)
    })
})
