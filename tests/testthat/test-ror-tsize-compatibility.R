test_that("Shiny infers logical TSIZE only for a whole logical column", {
    path <- withr::local_tempfile()
    for (sep in c(",", "\t")) for (quoted in c(FALSE, TRUE)) {
        upload <- function(size) {
            field <- if (quoted) paste0('"', size, '"') else size
            writeLines(c(paste("PatientID", "TSIZE", sep = sep),
                         paste(seq_along(field), field, sep = sep)), path)
            BreastSubtypeR:::.read_shiny_metadata(path, sep, "PatientID")$TSIZE
        }
        for (labels in list(c("TRUE", "FALSE"), c("T", "F"), c("TRUE", "F"))) {
            expect_identical(upload(labels), c(TRUE, FALSE))
            expect_identical(upload(c(labels, "", "NA", "NaN")),
                             c(TRUE, FALSE, NA, NA, NA))
        }
        mixed <- c("0", "1", "TRUE", "F", "2", "25", "1.00000000000000001",
                   "1e-999", "true", "unknown")
        expect_identical(upload(mixed), mixed)
        decoded <- BreastSubtypeR:::.decode_tsize(upload(mixed))
        expect_identical(decoded$value, c(0, 1, rep(NA_real_, 8)))
        expect_identical(decoded$invalid, c(FALSE, FALSE, rep(TRUE, 8)))
    }
})

test_that("exact binary labels retain numeric whitespace and genuine missingness", {
    decode <- BreastSubtypeR:::.decode_tsize
    expect_identical(decode(c("\v1", "0\f", "0x2p-1", "10e-1"))$value, c(1, 0, 1, 1))
    expect_identical(decode(addNA(factor(c("1", NA, "25")))),
                     list(value = c(1, NA_real_, NA_real_), invalid = c(FALSE, FALSE, TRUE)))
    expect_identical(decode(c("NA", "NaN", "", NA_character_))$invalid,
                     c(TRUE, TRUE, TRUE, FALSE))
})

test_that("labelled numeric storage survives an unavailable vctrs double cast", {
    skip_if_not_installed("vctrs")
    # Exercise the class reported on a saved numeric column. Its available
    # cast methods depend on which namespaces the caller has loaded.
    x <- vctrs::new_vctr(c(1, 0, 25, NA_real_), class = "haven_labelled")
    expect_identical(BreastSubtypeR:::.decode_tsize(x),
                     list(value = c(1, 0, NA_real_, NA_real_),
                          invalid = c(FALSE, FALSE, TRUE, FALSE)))
    opaque <- vctrs::new_vctr(c(1, 0, NA_real_), class = "opaque_clinical_measure")
    expect_identical(BreastSubtypeR:::.decode_tsize(opaque),
                     list(value = rep(NA_real_, 3), invalid = c(TRUE, TRUE, FALSE)))
})

test_that("haven value labels preserve their numeric codes", {
    skip_if_not_installed("haven")
    x <- haven::labelled(c(0, 1, 2, NA_real_), c(small = 0, large = 1))
    expect_identical(BreastSubtypeR:::.decode_tsize(x),
                     list(value = c(0, 1, NA_real_, NA_real_),
                          invalid = c(FALSE, FALSE, TRUE, FALSE)))
})

test_that("integer64 TSIZE preserves exact 0 and 1 even beside a large integer", {
    skip_if_not_installed("bit64")
    x <- bit64::as.integer64(c("0", "1", "9223372036854775807", NA_character_))
    expect_silent(actual <- BreastSubtypeR:::.decode_tsize(x))
    expect_identical(actual, list(value = c(0, 1, NA_real_, NA_real_),
                                  invalid = c(FALSE, FALSE, TRUE, FALSE)))
})

test_that("ssBC accepts binary raw values and excludes the singleton from TSIZE warnings", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = environment())
    se <- OSLO2EMIT0obj$data_input$se_NC[, c(1, 2, 5, 7)]
    se$ER <- c("ER+", "ER-", "ER+", "ER+")
    se$NODE <- rep(4, 4)
    withr::local_seed(382)
    for (four in c(FALSE, TRUE)) {
        score <- function() expect_experimental_ror_pc(
            BS_ssBC(se, "ER", Subtype = four, hasClinical = TRUE))
        se$TSIZE <- c(1, 0, 0, NA_real_)
        expect_warning(reference <- score(), "has n=1")
        se$TSIZE <- as.raw(c(1, 25, 0, 2))
        before <- .Random.seed
        conditions <- list()
        expect_warning(actual <- withCallingHandlers(score(),
            bsr_invalid_tsize_warning = function(w) {
                conditions[[length(conditions) + 1L]] <<- w
                invokeRestart("muffleWarning")
            }), "has n=1")
        expect_identical(.Random.seed, before)
        expect_length(conditions, 1L)
        expect_identical(conditions[[1L]]$patient_ids, colnames(se)[4])
        expect_identical(actual, reference)
        expect_identical(actual$score.ROR[[ror221_fields[2]]],
                         c("high", NA_character_, "high", NA_character_))
    }
})

# This function is also serialized into a fresh R session, so it uses only
# base R and explicit package names. The binary control carries genuine NA
# wherever the classed input must withhold clinical ROR.
check_public_tsize <- function(value, expected_codes, invalid, route) {
    e <- new.env(parent = emptyenv())
    utils::data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = e)
    d <- lapply(e$OSLO2EMIT0obj$data_input, function(z) z[, 12:1])
    z <- d$se_NC
    ids <- colnames(z)
    z$ER <- rep("ER+", 12)
    z$NODE <- rep(4, 12)
    rng <- get(".Random.seed", globalenv())
    run <- function(x) {
        z$TSIZE <- x
        original <- SummarizedExperiment::colData(z)
        d$se_NC <- z
        warnings <- list()
        out <- suppressMessages(withCallingHandlers(switch(route,
            parker = BreastSubtypeR::BS_parker(z, hasClinical = TRUE),
            ssBC = BreastSubtypeR::BS_ssBC(z, s = "ER", hasClinical = TRUE),
            multi = BreastSubtypeR::BS_Multi(d, methods = c("parker.original", "ssBC"),
                hasClinical = TRUE)), warning = function(w) {
                    warnings[[length(warnings) + 1L]] <<-
                        list(class = class(w)[1], ids = w$patient_ids, message = conditionMessage(w))
                    invokeRestart("muffleWarning")
                }))
        stopifnot(identical(original, SummarizedExperiment::colData(z)),
                  identical(rng, get(".Random.seed", globalenv())))
        list(out = out, warnings = warnings)
    }
    # The classed input goes first so the numeric control cannot prime dispatch.
    got <- run(value)
    control <- run(expected_codes)
    is_size <- vapply(got$warnings, function(w) w$class == "bsr_invalid_tsize_warning", TRUE)
    stopifnot(identical(got$out, control$out),
              identical(got$warnings[!is_size], control$warnings),
              sum(is_size) == as.integer(length(invalid) > 0L))
    if (length(invalid)) {
        warned <- got$warnings[is_size][[1L]]$ids
        stopifnot(setequal(warned, ids[invalid]), !anyDuplicated(warned))
    }
    scores <- c("ROR-C (Subtype + Clinic)", "ROR-PC (Subtype + Clinic + Prolif)",
                "ROR-PC (Subtype + Clinic + Prolif.Subtype)")
    groups <- sub(" (", " Group (", scores, fixed = TRUE)
    for (method in if (route == "multi") got$out$results else list(got$out)) {
        r <- method$score.ROR
        stopifnot(identical(rownames(r), ids))
        empty <- is.na(expected_codes)
        for (f in scores) stopifnot(is.double(r[[f]]), all(is.na(r[[f]][empty])))
        for (f in groups) stopifnot(is.character(r[[f]]), all(is.na(r[[f]][empty])))
        stopifnot(all(r[[groups[3]]][!empty] == "high"))
    }
    TRUE
}

test_that("stored character codes work through public clinical NC and multi-method calls", {
    withr::local_seed(314)
    exercise <- function(wrap) {
        valid <- rep(c("0", "1", " 1 ", "1e0"), 3)
        mixed <- valid
        mixed[5:8] <- c("2", "25", "T1", NA_character_)
        codes <- rep(c(0, 1, 1, 1), 3)
        for (route in c("parker", "ssBC", "multi")) {
            expect_true(check_public_tsize(wrap(valid), codes, integer(), route))
            expect_true(check_public_tsize(wrap(mixed), replace(codes, 5:8, NA_real_),
                                           5:7, route))
        }
    }
    skip_if_not_installed("vctrs")
    exercise(function(x) vctrs::new_vctr(x, class = "clinical_text"))
    skip_if_not_installed("haven")
    exercise(function(x) haven::labelled(x, c(small = "0", large = "1")))
    exercise(function(x) haven::labelled_spss(x, labels = c(small = "0", large = "1")))
})

test_that("serialized integer64 survives first-use row alignment without attaching bit64", {
    skip_if_not_installed("bit64")
    config <- withr::local_tempfile()
    script <- withr::local_tempfile(fileext = ".R")
    log <- withr::local_tempfile()
    x <- bit64::as.integer64(rep(c("0", "1", NA, "2", "-1", "9223372036854775807"), 2))
    check <- check_public_tsize
    environment(check) <- baseenv()
    saveRDS(list(x = x, libs = .libPaths(), check = check), config)
    writeLines(c(
        'a <- commandArgs(TRUE)',
        'cfg <- readRDS(a[1]); .libPaths(cfg$libs)',
        'suppressPackageStartupMessages(library(BreastSubtypeR))',
        'stopifnot(!isNamespaceLoaded("bit64"))',
        'set.seed(314)',
        'stopifnot(cfg$check(cfg$x, rep(c(0,1,NA_real_,NA_real_,NA_real_,NA_real_),2),',
        '                    c(4:6,10:12), a[2]))',
        'stopifnot(isNamespaceLoaded("bit64"), !"package:bit64" %in% search())'
    ), script)
    exe <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
    withr::local_envvar(c(R_TESTS = ""))
    for (route in c("parker", "ssBC", "multi")) {
        status <- system2(exe, c("--vanilla", shQuote(script), shQuote(config), route),
                          stdout = log, stderr = log)
        expect_identical(status, 0L, info = paste(route, paste(readLines(log), collapse = "\n")))
    }
})

test_that("unsupported numeric classes cannot contribute binary storage values", {
    for (x in list(as.Date(c(0, 1, NA), origin = "1970-01-01"),
                   as.difftime(c(0, 1, NA), units = "secs"),
                   c(0+0i, 1+0i, NA_complex_),
                   structure(c(0, 1, NA_real_), class = "unknown_measure"))) {
        expect_identical(BreastSubtypeR:::.decode_tsize(x),
                         list(value = rep(NA_real_, 3), invalid = c(TRUE, TRUE, FALSE)))
    }
})

test_that("logical and raw binary values retain their existing AsIs representation support", {
    for (x in list(I(c(FALSE, TRUE)), I(as.raw(0:1)))) {
        expect_identical(BreastSubtypeR:::.decode_tsize(x),
                         list(value = c(0, 1), invalid = c(FALSE, FALSE)))
    }
})

test_that("row selection cannot turn an unknown numeric class into binary TSIZE", {
    withr::local_seed(314)
    x <- structure(c(0, 1, 2, NA, 1, 0, 1, 25, 0, 1, 0, 1), class = "unknown_measure")
    for (route in c("parker", "ssBC", "multi"))
        expect_true(check_public_tsize(x, rep(NA_real_, 12), which(!is.na(x)), route))
})

test_that("SPSS missing codes and ranges score identically in fresh cold and warm sessions", {
    skip_if_not_installed("haven")
    char <- haven::labelled_spss(c("0", "1", "25", "T1", NA, "1e0"),
                                 na_values = c("1", "25"))
    num <- haven::labelled_spss(c(0, 1, 25, 2, NA_real_, NaN), na_values = c(1, 25))
    range <- haven::labelled_spss(c(0, 1, 2, NA_real_, NaN, 25), na_range = c(0, 1))
    cases <- list(
        list(value = rep(char, 2), codes = rep(c(0, NA, NA, NA, NA, 1), 2), bad = c(4, 10)),
        list(value = rep(num, 2), codes = rep(c(0, NA_real_, NA, NA, NA, NA), 2), bad = c(4, 10)),
        list(value = rep(range, 2), codes = rep(NA_real_, 12), bad = c(3, 6, 9, 12)))
    cases[[4L]] <- cases[[1L]]
    cases[[4L]]$value <- I(cases[[1L]]$value)
    check <- check_public_tsize
    environment(check) <- baseenv()
    config <- withr::local_tempfile()
    saveRDS(list(cases = cases, libs = .libPaths(), check = check), config)
    script <- withr::local_tempfile(fileext = ".R")
    log <- withr::local_tempfile()
    writeLines(c(
        'args <- commandArgs(TRUE)',
        'cfg <- readRDS(args[1L]); .libPaths(cfg$libs)',
        'suppressPackageStartupMessages(library(BreastSubtypeR))',
        'stopifnot(!isNamespaceLoaded("haven"))',
        'if (args[2L] == "absent") {',
        '    data("OSLO2EMIT0obj", package = "BreastSubtypeR")',
        '    z <- OSLO2EMIT0obj$data_input$se_NC[, 12:1]',
        '    genomic <- function() suppressMessages(BS_parker(z, hasClinical = FALSE))',
        '    set.seed(314); reference <- genomic(); seed <- .Random.seed',
        '    z$TSIZE <- cfg$cases[[4L]]$value',
        '    .libPaths(.Library, include.site = FALSE)',
        '    if (length(find.package("haven", quiet = TRUE))) quit(status = 77L)',
        '    msg <- tryCatch({BS_parker(z, hasClinical = TRUE); ""}, error = conditionMessage)',
        '    stopifnot(grepl("Install haven", msg), identical(genomic(), reference),',
        '              identical(.Random.seed, seed), !isNamespaceLoaded("haven"))',
        '    quit(status = 0L)',
        '}',
        'if (args[2L] == "warm") loadNamespace("haven")',
        'set.seed(314)',
        'case <- cfg$cases[[as.integer(args[3L])]]',
        'stopifnot(cfg$check(case$value, case$codes, case$bad, "parker"),',
        '          isNamespaceLoaded("haven"), !"package:haven" %in% search())'
    ), script)
    exe <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
    withr::local_envvar(c(R_TESTS = ""))
    for (i in seq_along(cases)) for (mode in c("cold", "warm")) {
        status <- system2(exe, c("--vanilla", shQuote(script), shQuote(config), mode, i),
                          stdout = log, stderr = log)
        expect_identical(status, 0L, info = paste(readLines(log), collapse = "\n"))
    }
    status <- system2(exe, c("--vanilla", shQuote(script), shQuote(config), "absent"),
                      stdout = log, stderr = log)
    skip_if(identical(status, 77L), "haven is installed in R's base library")
    expect_identical(status, 0L, info = paste(readLines(log), collapse = "\n"))
})
