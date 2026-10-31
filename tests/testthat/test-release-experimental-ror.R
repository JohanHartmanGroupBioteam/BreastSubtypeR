test_that("release public boundaries warn once without altering results or RNG", {
    input <- r143_input()
    before_input <- serialize(input, NULL)
    withr::local_preserve_seed()
    withr::local_options(BreastSubtypeR.verbose = FALSE)
    for (method in c("BS_parker", "BS_cIHC", "BS_cIHC.itr",
                     "BS_PCAPAM50", "BS_ssBC")) {
        for (clinical in c(FALSE, TRUE)) for (subtype in c(FALSE, TRUE)) {
            set.seed(143)
            before_seed <- .Random.seed
            actual <- r143_capture(suppressMessages(
                r143_run(method, input, clinical, subtype)))
            expect_identical(.Random.seed, before_seed)
            expect_identical(r143_count(actual), as.integer(clinical))
            expect_length(actual$warnings, as.integer(clinical))
            expected <- testthat::with_mocked_bindings(
                suppressMessages(r143_run(method, input, clinical, subtype)),
                .warn_experimental_ror_pc = function(result) invisible(result),
                .package = "BreastSubtypeR")
            expect_identical(actual$value, expected)
            expect_identical(.Random.seed, before_seed)
            if (clinical) {
                expect_identical(conditionMessage(actual$warnings[[1]]),
                    BreastSubtypeR:::.experimental_ror_pc_notice())
                expect_null(conditionCall(actual$warnings[[1]]))
            }
        }
    }
    expect_identical(serialize(input, NULL), before_input)
    # A condition must not create a seed where none existed.
    if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    r143_capture(suppressMessages(BS_cIHC.itr(input$se_NC,
        iteration = 2L, hasClinical = TRUE)))
    expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
})

test_that("release conditions follow actual fields, including missing values", {
    input <- r143_input()
    table <- data.frame(patientID = "S1", check.names = FALSE)
    table[[r143_fields[1]]] <- NA_real_
    table[[r143_fields[2]]] <- NA_character_
    cases <- list(NULL, list(score.ROR = table[FALSE, ]),
        list(score.ROR = table[, -3]), list(score.ROR = table))
    for (i in seq_along(cases)) {
        result <- cases[[i]]
        actual <- testthat::with_mocked_bindings(
            r143_capture(BS_parker(input$se_NC, hasClinical = TRUE)),
            makeCalls.parker = function(...) result, .package = "BreastSubtypeR")
        expect_identical(actual$value, result)
        expect_identical(r143_count(actual), as.integer(i == 4L))
    }
    one <- r143_capture(suppressMessages(BS_parker(input$se_NC, hasClinical = TRUE)))
    two <- r143_capture(suppressMessages(BS_parker(input$se_NC, hasClinical = TRUE)))
    expect_identical(r143_count(two), 1L)
    expect_identical(one$value, two$value)
    internal <- r143_capture(BreastSubtypeR:::RORgroup(one$value$outList,
        as.data.frame(SummarizedExperiment::colData(input$se_NC)), hasClinical = TRUE))
    expect_length(internal$warnings, 0L)
    expect_identical(internal$value, one$value$score.ROR)
})

test_that("release Multi aggregates only experimental conditions", {
    input <- r143_input()
    for (clinical in c(FALSE, TRUE)) {
        actual <- r143_capture(suppressMessages(BS_Multi(input,
            methods = c("parker.original", "PCAPAM50", "cIHC.itr"),
            hasClinical = clinical)))
        expect_identical(r143_count(actual), as.integer(clinical))
        expect_length(actual$warnings, as.integer(clinical))
        expected <- testthat::with_mocked_bindings(
            suppressMessages(BS_Multi(input,
                methods = c("parker.original", "PCAPAM50", "cIHC.itr"),
                hasClinical = clinical)),
            .warn_experimental_ror_pc = function(result) invisible(result),
            .package = "BreastSubtypeR")
        expect_identical(actual$value, expected)
    }
    ssp <- function(se_obj, ...) list(cl = matrix("LumA", ncol = 1,
        nrow = ncol(se_obj), dimnames = list(colnames(se_obj), "Subtype")))
    local_mocked_bindings(BS_AIMS = ssp, BS_sspbc = ssp,
        BS_PCAPAM50 = function(...) stop("release sentinel error"),
        .package = "BreastSubtypeR")
    only_ssp <- r143_capture(suppressMessages(BS_Multi(input,
        methods = c("AIMS", "sspbc"), hasClinical = TRUE)))
    expect_length(only_ssp$warnings, 0L)
    failure <- r143_capture(suppressMessages(BS_Multi(input,
        methods = c("PCAPAM50", "AIMS"), hasClinical = TRUE)))
    expect_identical(r143_count(failure), 0L)
    expect_length(failure$warnings, 1L)
    expect_match(conditionMessage(failure$warnings[[1]]),
        "PCAPAM50 failed in this iteration: release sentinel error", fixed = TRUE)
    mixed <- r143_capture(suppressMessages(BS_Multi(input,
        methods = c("parker.original", "PCAPAM50"), hasClinical = TRUE)))
    expect_identical(r143_count(mixed), 1L)
    expect_length(mixed$warnings, 2L)
})

test_that("release NODE handling and warn=2 are preserved apart from new notice", {
    input <- r143_input()
    SummarizedExperiment::colData(input$se_NC)$NODE[1:5] <- c(-1, .5, Inf, NA, 4)
    actual <- r143_capture(suppressMessages(BS_parker(input$se_NC, hasClinical = TRUE)))
    # Release 1.4.2 has no development invalid-NODE warning/validation.
    expect_length(actual$warnings, 1L)
    expect_identical(r143_count(actual), 1L)
    expected <- testthat::with_mocked_bindings(
        suppressMessages(BS_parker(input$se_NC, hasClinical = TRUE)),
        .warn_experimental_ror_pc = function(result) invisible(result),
        .package = "BreastSubtypeR")
    expect_identical(actual$value, expected)
    withr::local_options(warn = 2)
    expect_error(suppressMessages(BS_parker(input$se_NC, hasClinical = TRUE)),
        "Experimental 18-gene ROR-PC", fixed = TRUE)
    actual <- withCallingHandlers(suppressMessages(BS_Multi(input,
        methods = c("parker.original", "PCAPAM50"), hasClinical = TRUE)),
        bsr_experimental_ror_pc_warning = function(w) invokeRestart("muffleWarning"))
    expect_false(is.null(actual$results$PCAPAM50))
    testthat::local_mocked_bindings(makeCalls.parker = function(...) {
        warning("unrelated sentinel")
    }, .package = "BreastSubtypeR")
    expect_error(withCallingHandlers(BS_parker(input$se_NC, hasClinical = TRUE),
        bsr_experimental_ror_pc_warning = function(w) invokeRestart("muffleWarning")),
        "unrelated sentinel", fixed = TRUE)
})
