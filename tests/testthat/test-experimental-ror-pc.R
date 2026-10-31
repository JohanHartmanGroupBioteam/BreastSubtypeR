## #221: signal at public return boundaries, never in the arithmetic helper.
ror221_capture <- function(expr) {
    warnings <- list()
    value <- withCallingHandlers(force(expr), warning = function(w) {
        warnings[[length(warnings) + 1L]] <<- w
        invokeRestart("muffleWarning")
    })
    list(value = value, warnings = warnings)
}

ror221_count <- function(x) sum(vapply(x$warnings, inherits, logical(1),
    what = "bsr_experimental_ror_pc_warning"))

ror221_fixture <- function() {
    env <- new.env()
    data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = env)
    env$OSLO2EMIT0obj$data_input
}

ror221_run <- function(method, se, clinical, subtype) {
    switch(method,
        parker = BS_parker(se, hasClinical = clinical, Subtype = subtype),
        cIHC = BS_cIHC(se, hasClinical = clinical, Subtype = subtype),
        cIHC.itr = BS_cIHC.itr(se, iteration = 3L,
            hasClinical = clinical, Subtype = subtype),
        PCAPAM50 = BS_PCAPAM50(se, hasClinical = clinical, Subtype = subtype),
        ssBC = BS_ssBC(se, s = "ER", hasClinical = clinical, Subtype = subtype))
}

test_that("each NC public return warns once without changing its result or RNG", {
    input <- ror221_fixture()$se_NC
    snapshot <- serialize(input, NULL)
    withr::local_preserve_seed()
    withr::local_options(BreastSubtypeR.verbose = FALSE)
    for (method in c("parker", "cIHC", "cIHC.itr", "PCAPAM50", "ssBC")) {
        for (clinical in c(FALSE, TRUE)) for (subtype in c(FALSE, TRUE)) {
            info <- paste(method, clinical, subtype)
            set.seed(221)
            before <- .Random.seed
            actual <- ror221_capture(suppressMessages(
                ror221_run(method, input, clinical, subtype)))
            expect_identical(.Random.seed, before, info = info)
            expect_identical(ror221_count(actual), as.integer(clinical), info = info)
            expect_length(actual$warnings, as.integer(clinical))
            if (clinical) {
                expect_identical(conditionMessage(actual$warnings[[1L]]), ror221_notice)
                expect_null(conditionCall(actual$warnings[[1L]]))
                expect_true(all(ror221_fields %in% names(actual$value$score.ROR)))
            } else expect_false(any(ror221_fields %in% names(actual$value$score.ROR)))

            # Remove only signaling, rerun the same unmodified calculation,
            # and require byte-identical complete objects (including all axes,
            # legacy keys, types, numerical scores, groups and attributes).
            baseline <- testthat::with_mocked_bindings(
                suppressMessages(ror221_run(method, input, clinical, subtype)),
                .warn_experimental_ror_pc = function(result) invisible(result),
                .package = "BreastSubtypeR")
            expect_identical(actual$value, baseline, info = info)
            expect_identical(.Random.seed, before, info = info)
        }
    }
    expect_identical(serialize(input, NULL), snapshot)
})

test_that("experimental warning is scoped to each call and absent from internal RORgroup", {
    se <- ror221_fixture()$se_NC
    first <- expect_experimental_ror_pc(suppressMessages(
        BS_parker(se, hasClinical = TRUE)))
    second <- expect_experimental_ror_pc(suppressMessages(
        BS_parker(se, hasClinical = TRUE)))
    expect_identical(first, second)
    internal <- ror221_capture(BreastSubtypeR:::RORgroup(first$outList,
        as.data.frame(SummarizedExperiment::colData(se)), hasClinical = TRUE))
    expect_length(internal$warnings, 0L)
    expect_identical(internal$value, first$score.ROR)
    withr::with_preserve_seed({
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv)
        expect_experimental_ror_pc(suppressMessages(
            BS_cIHC.itr(se, iteration = 2L, hasClinical = TRUE)))
        expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    })
})

test_that("public signaling follows actual nonempty fields instead of intent flags", {
    se <- ror221_fixture()$se_NC
    base <- data.frame(patientID = colnames(se), check.names = FALSE)
    base[[ror221_fields[1]]] <- NA_real_
    base[[ror221_fields[2]]] <- NA_character_
    cases <- list(no_result = NULL, no_scores = list(BS.all = data.frame()),
        empty = list(score.ROR = base[FALSE, , drop = FALSE]),
        one_field = list(score.ROR = base[, -3, drop = FALSE]),
        parker_only = list(score.ROR = data.frame(
            "ROR-PC (Subtype + Clinic + Prolif)" = 50,
            "ROR-PC Group (Subtype + Clinic + Prolif)" = "med", check.names = FALSE)),
        experimental = list(score.ROR = base))
    for (name in names(cases)) {
        mock_result <- cases[[name]]
        actual <- testthat::with_mocked_bindings(
            ror221_capture(BS_parker(se, hasClinical = TRUE)),
            makeCalls.parker = function(...) mock_result,
            .package = "BreastSubtypeR")
        expect_identical(actual$value, mock_result)
        expect_identical(ror221_count(actual), as.integer(name == "experimental"))
    }
})

test_that("BS_Multi aggregates real nested calculations without changing results", {
    input <- ror221_fixture()
    withr::local_options(BreastSubtypeR.verbose = FALSE)
    methods <- c("parker.original", "PCAPAM50", "cIHC.itr")
    for (clinical in c(FALSE, TRUE)) {
        actual <- ror221_capture(suppressMessages(
            BS_Multi(input, methods = methods, hasClinical = clinical)))
        expect_identical(ror221_count(actual), as.integer(clinical))
        expect_length(actual$warnings, as.integer(clinical))
        expect_named(actual$value$results, methods)
        for (method in methods) {
            expect_identical(all(ror221_fields %in%
                names(actual$value$results[[method]]$score.ROR)), clinical)
        }
        baseline <- testthat::with_mocked_bindings(
            suppressMessages(BS_Multi(input, methods = methods, hasClinical = clinical)),
            .warn_experimental_ror_pc = function(result) invisible(result),
            .package = "BreastSubtypeR")
        expect_identical(actual$value, baseline)
    }
})

ror221_ssp <- function(se_obj, ...) {
    ids <- colnames(se_obj)
    list(cl = matrix("LumA", nrow = length(ids), ncol = 1L,
        dimnames = list(ids, "Subtype")))
}

test_that("BS_Multi does not warn for SSP-only or failed experimental outputs", {
    input <- ror221_fixture()
    local_mocked_bindings(BS_AIMS = ror221_ssp, BS_sspbc = ror221_ssp,
        BS_PCAPAM50 = function(...) stop("synthetic PCAPAM50 failure"),
        .package = "BreastSubtypeR")
    ssp <- ror221_capture(suppressMessages(BS_Multi(input,
        methods = c("AIMS", "sspbc"), hasClinical = TRUE)))
    expect_length(ssp$warnings, 0L)
    failed <- ror221_capture(suppressMessages(BS_Multi(input,
        methods = c("PCAPAM50", "AIMS"), hasClinical = TRUE)))
    expect_identical(ror221_count(failed), 0L)
    expect_length(failed$warnings, 1L)
    expect_match(conditionMessage(failed$warnings[[1L]]),
        "PCAPAM50 failed in this iteration: synthetic PCAPAM50 failure", fixed = TRUE)
    expect_null(failed$value$results$PCAPAM50)
})

test_that("NODE and component warnings survive experimental aggregation", {
    input <- ror221_fixture()
    SummarizedExperiment::colData(input$se_NC)$NODE[1L] <- -1
    one <- ror221_capture(suppressMessages(BS_parker(input$se_NC, hasClinical = TRUE)))
    expect_identical(ror221_count(one), 1L)
    expect_length(one$warnings, 2L)
    expect_match(conditionMessage(one$warnings[[1L]]), "finite nonnegative integer")
    expect_true(is.na(one$value$score.ROR[[ror221_fields[2L]]][1L]))
    local_mocked_bindings(BS_PCAPAM50 = function(...) stop("synthetic failure"),
        .package = "BreastSubtypeR")
    multi <- ror221_capture(suppressMessages(BS_Multi(input,
        methods = c("parker.original", "PCAPAM50"), hasClinical = TRUE)))
    expect_identical(ror221_count(multi), 1L)
    expect_length(multi$warnings, 3L)
    messages <- vapply(multi$warnings, conditionMessage, character(1))
    expect_true(any(grepl("finite nonnegative integer", messages)))
    expect_true(any(grepl("PCAPAM50 failed in this iteration: synthetic failure", messages,
        fixed = TRUE)))
})

test_that("warn=2 remains explicit and selective handlers preserve unrelated warnings", {
    input <- ror221_fixture()
    withr::local_options(warn = 2)
    expect_error(suppressMessages(BS_parker(input$se_NC, hasClinical = TRUE)),
        "Experimental 18-gene ROR-PC", fixed = TRUE)
    # Nested experimental warnings must be muffled before warn=2 conversion;
    # only the final BS_Multi condition reaches this caller's handler.
    result <- expect_experimental_ror_pc(suppressMessages(BS_Multi(input,
        methods = c("parker.original", "PCAPAM50"), hasClinical = TRUE)))
    expect_false(is.null(result$results$PCAPAM50))
    expect_experimental_ror_pc(suppressMessages(BS_parker(input$se_NC)), n = 0L)
    SummarizedExperiment::colData(input$se_NC)$NODE[1L] <- -1
    expect_error(expect_experimental_ror_pc(suppressMessages(
        BS_parker(input$se_NC, hasClinical = TRUE))), "finite nonnegative integer")
})
