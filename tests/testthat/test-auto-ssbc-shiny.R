test_that("installed Shiny AUTO and single ssBC calls, plots and downloads remain aligned", {
    skip_if_not_installed("shiny")
    skip_if_not_installed("bslib")
    env <- new.env(parent = asNamespace("shiny"))
    eval(parse(system.file("ShinyBreastSubtypeR/server.R", package = "BreastSubtypeR"),
        encoding = "UTF-8", keep.source = FALSE), env)
    e <- new.env(parent = emptyenv())
    utils::data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = e)
    input_data <- lapply(e$OSLO2EMIT0obj$data_input, function(x) {
        x <- x[, ncol(x):1, drop = FALSE]
        SummarizedExperiment::colData(x)$PatientID <- factor(colnames(x),
            levels = c("unused", rev(colnames(x))))
        SummarizedExperiment::colData(x)$ER[c(1, 3)] <- NA_character_
        SummarizedExperiment::colData(x)$HER2[c(4, 6)] <- NA_character_
        x
    })
    withr::local_seed(233)
    for (method in c("AUTO Mode", "ER", "ER.v2")) for (four in c(FALSE, TRUE)) {
        auto <- method == "AUTO Mode"
        expected <- suppressMessages(if (auto) BS_Multi(input_data, Subtype = four,
            hasClinical = FALSE) else expect_experimental_ror_pc(
                BS_ssBC(input_data$se_NC, s = method, Subtype = four, hasClinical = TRUE)))
        shiny::testServer(env$server, {
            session$setInputs(BSmethod = if (auto) method else "ssBC",
                s = if (auto) "ER" else method, k_subtypes = if (four) "4" else "5",
                is_raw_counts = "norm", hasClinical = TRUE, export_kind = "calls")
            # Prepared public input avoids repeating Mapping in this wrapper test.
            reactive_files$data_input <- input_data
            reactive_files$mapped_key <- input_key()
            session$flushReact()
            session$setInputs(run = 1)
            expect_null(reactive_files$error)
            prediction <- current_prediction()
            expect_false(is.null(prediction))
            expect_identical(prediction$calls$PatientID, colnames(input_data$se_NC))
            if (auto) {
                tab <- expected[[if (four) "res_subtypes.Subtype" else "res_subtypes"]]
                expect_identical(prediction$calls[, -1, drop = FALSE], tab)
                expect_null(prediction$full)
                expect_false(is.null(output$multi_plot))
            } else {
                expect_identical(prediction$calls[[2]], expected$BS.all[[if (four) "BS.Subtype" else "BS"]])
                scores <- expected$score.ROR[match(prediction$full$PatientID,
                    rownames(expected$score.ROR)), , drop = FALSE]
                for (n in intersect(names(scores), names(prediction$full)))
                    expect_identical(prediction$full[[n]], scores[[n]])
                expect_false(is.null(output$heat2))
                expect_false(is.null(output$pie1))
            }
            downloaded <- utils::read.delim(output$download, comment.char = "#", check.names = FALSE)
            exported_calls <- prediction$calls
            # TSVs deliberately omit R row names; PatientID is an explicit column.
            rownames(exported_calls) <- NULL
            expect_equal(downloaded, exported_calls, tolerance = 1e-12)
            session$setInputs(export_kind = "full")
            if (auto) expect_error(export_result("full")) else {
                downloaded <- utils::read.delim(output$download, comment.char = "#", check.names = FALSE)
                exported_full <- prediction$full
                rownames(exported_full) <- NULL
                # Allow only decimal text round-trip error, not scientific drift.
                expect_equal(downloaded, exported_full, tolerance = 1e-12)
            }
        })
    }
})
