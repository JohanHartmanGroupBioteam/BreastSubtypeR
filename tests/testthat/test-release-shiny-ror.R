r143_shiny_env <- function() {
    skip_if_not_installed("shiny")
    skip_if_not_installed("bslib")
    env <- new.env(parent = asNamespace("shiny"))
    eval(parse(system.file("ShinyBreastSubtypeR/server.R", package = "BreastSubtypeR"),
        encoding = "UTF-8", keep.source = FALSE), env)
    env
}
r143_shiny_defaults <- function(session) {
    session$setInputs(BSmethod = "PAM50", k_subtypes = "5",
        is_raw_counts = "norm", hasClinical = TRUE, export_kind = "calls",
        calibration = "Internal", internal = "medianCtr",
        s = "ER", iteration = 3, ratio = 54/64)
}

test_that("release Shiny notice follows selected complete results and run context", {
    env <- r143_shiny_env()
    notice <- BreastSubtypeR:::.experimental_ror_pc_notice()
    fields <- BreastSubtypeR:::.experimental_ror_pc_fields()
    shiny::testServer(env$server, {
        r143_shiny_defaults(session)
        reactive_files$data_input <- r143_input()
        r143_capture(session$setInputs(run = 1))
        expect_false(is.null(reactive_files$download_full))
        expect_null(output$experimental_ror_pc_notice)
        calls <- readLines(output$download)
        expect_length(grep("^# ", calls), 3L)
        expect_false(any(grepl(notice, calls, fixed = TRUE)))
        session$setInputs(export_kind = "full")
        html <- output$experimental_ror_pc_notice$html
        expect_match(html, notice, fixed = TRUE)
        for (label in unname(fields)) expect_match(html, label, fixed = TRUE)
        full <- readLines(output$download)
        expect_length(grep("^# ", full), 3L)
        expect_match(full[2L], notice, fixed = TRUE)
        exported <- read.delim(output$download, comment.char = "#", check.names = FALSE)
        expect_identical(names(exported), names(reactive_files$download_full))
        expect_identical(full[-seq_len(3L)], capture.output(utils::write.table(
            reactive_files$download_full, sep = "\t", quote = FALSE, row.names = FALSE)))
        # Switching display mode does not discard valid results.
        session$setInputs(export_kind = "calls")
        expect_null(output$experimental_ror_pc_notice)
        session$setInputs(export_kind = "full")
        expect_match(output$experimental_ror_pc_notice$html, notice, fixed = TRUE)
        # A changed analysis setting hides a stale notice without restructuring
        # the release's existing storage/download behavior.
        session$setInputs(internal = "qCtr")
        expect_null(output$experimental_ror_pc_notice)
        r143_capture(session$setInputs(run = 2))
        expect_match(output$experimental_ror_pc_notice$html, notice, fixed = TRUE)
        session$setInputs(hasClinical = FALSE)
        expect_null(output$experimental_ror_pc_notice)
        session$setInputs(run = 3)
        expect_null(output$experimental_ror_pc_notice)
        expect_false(any(grepl(notice, readLines(output$download), fixed = TRUE)))
        session$setInputs(BSmethod = "AUTO Mode")
        expect_null(output$experimental_ror_pc_notice)
        session$setInputs(run = 4)
        expect_length(grep("^# ", readLines(output$download)), 2L)
        expect_false(any(grepl(notice, readLines(output$download), fixed = TRUE)))
    })
})

test_that("release Shiny requires both fields and clears a notice before failed runs", {
    env <- r143_shiny_env()
    fields <- names(BreastSubtypeR:::.experimental_ror_pc_fields())
    shiny::testServer(env$server, {
        r143_shiny_defaults(session)
        reactive_files$data_input <- r143_input()
        r143_capture(session$setInputs(run = 1))
        session$setInputs(export_kind = "full")
        expect_false(is.null(output$experimental_ror_pc_notice))
        saved <- reactive_files$download_full
        reactive_files$download_full <- saved[, setdiff(names(saved), fields[2]), drop = FALSE]
        session$flushReact()
        expect_null(output$experimental_ror_pc_notice)
        reactive_files$download_full <- saved[FALSE, , drop = FALSE]
        session$flushReact()
        expect_null(output$experimental_ror_pc_notice)
        reactive_files$download_full <- saved
        session$flushReact()
        expect_false(is.null(output$experimental_ror_pc_notice))
        # req() aborts this new run before prediction; no old notice survives.
        reactive_files$data_input <- NULL
        session$setInputs(run = 2)
        expect_null(output$experimental_ror_pc_notice)
    })
})

test_that("release ssBC retains unrelated warning notifications", {
    env <- r143_shiny_env()
    shown <- character()
    env$showNotification <- function(ui, ...) shown <<- c(shown, as.character(ui))
    fixture_data <- r143_input()
    result <- r143_capture(suppressMessages(BS_ssBC(fixture_data$se_NC,
        s = "ER", hasClinical = TRUE)))$value
    local_mocked_bindings(BS_ssBC = function(...) {
        warning("release unrelated warning", call. = FALSE)
        BreastSubtypeR:::.warn_experimental_ror_pc(result)
        result
    }, .package = "BreastSubtypeR")
    shiny::testServer(env$server, {
        r143_shiny_defaults(session)
        session$setInputs(BSmethod = "ssBC", export_kind = "full")
        reactive_files$data_input <- fixture_data
        session$setInputs(run = 1)
        expect_match(output$experimental_ror_pc_notice$html,
            BreastSubtypeR:::.experimental_ror_pc_notice(), fixed = TRUE)
    })
    expect_true("release unrelated warning" %in% shown)
    expect_false(any(grepl("Experimental 18-gene ROR-PC", shown, fixed = TRUE)))
})
