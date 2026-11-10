# Each reference result uses binary TSIZE throughout. Only these six cells
# per unavailable ID are replaced when constructing expected public returns.
tsize_scores <- c("ROR-C (Subtype + Clinic)", "ROR-PC (Subtype + Clinic + Prolif)",
    "ROR-PC (Subtype + Clinic + Prolif.Subtype)")
tsize_groups <- sub(" (", " Group (", tsize_scores, fixed = TRUE)

tsize_expected <- function(result, ids) {
    s <- result$score.ROR
    if (is.null(s) || !all(tsize_scores %in% names(s))) return(result)
    at <- rownames(s) %in% ids
    for (field in tsize_scores) s[[field]][at] <- NA_real_
    for (field in tsize_groups) s[[field]][at] <- NA_character_
    result$score.ROR <- s
    result
}

expect_tsize_warning <- function(expr, ids) {
    conditions <- list()
    value <- withCallingHandlers(force(expr), bsr_invalid_tsize_warning = function(w) {
        conditions[[length(conditions) + 1L]] <<- w
        invokeRestart("muffleWarning")
    })
    expect_length(conditions, as.integer(length(ids) > 0L))
    if (length(ids)) {
        w <- conditions[[1L]]
        expect_setequal(w$patient_ids, ids)
        expect_length(w$patient_ids, length(ids))
        expect_match(conditionMessage(w), "TSIZE", fixed = TRUE)
        if (length(ids) < 4L) expect_match(conditionMessage(w), ids[1L], fixed = TRUE)
    }
    value
}

tsize_input <- function() {
    e <- new.env(parent = emptyenv())
    utils::data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = e)
    e$OSLO2EMIT0obj
}

test_that("the 12-tumor research ROR reference masks only invalid size codes", {
    raw <- tsize_input()$se_obj[, 1:12]
    z <- suppressMessages(Mapping(raw, RawCounts = FALSE, method = "max",
        impute = TRUE, verbose = FALSE))$se_NC
    id <- "OSLO2EMIT0.003"
    # Independent rank-correlation/proliferation oracle, retained ROR-01
    # audit values for this patient at binary T = 0 and 1, respectively.
    reference <- rbind(c(53.349203286208606, 64.716899984172997, 73.387721444930975),
        c(74.687450345032133, 80.199954101820055, 79.593049144930987))
    for (four in c(FALSE, TRUE)) {
        for (size in 0:1) {
            z$TSIZE <- rep(size, ncol(z))
            binary <- expect_experimental_ror_pc(ror198_run("parker.original", z, subtype = four))
            scores <- unlist(binary$score.ROR[id, tsize_scores], use.names = FALSE)
            expect_lt(max(abs(scores - reference[size + 1L, ])), 1e-10)
        }
        z$TSIZE[match(id, colnames(z))] <- 25
        result <- expect_tsize_warning(expect_experimental_ror_pc(
            ror198_run("parker.original", z, subtype = four)), id)
        expect_identical(result, tsize_expected(binary, id))
    }
})

test_that("resolved haven and serialized legacy TSIZE preserve public NC results", {
    skip_if_not_installed("haven")
    withr::local_seed(217)
    input <- lapply(tsize_input()$data_input, function(z) z[, 70:1, drop = FALSE])
    ids <- colnames(input$se_NC)
    input$se_NC$PatientID <- factor(ids, levels = rev(ids))
    input$se_NC$ER <- rep(c("ER+", "ER-"), c(40, 30))
    input$se_NC$HER2 <- rep("HER2-", length(ids))
    input$se_NC$TN <- rep("nonTN", length(ids))
    input$se_NC$NODE <- rep(4, length(ids))

    binary <- rep(0:1, 35)
    stored <- binary
    stored[3:6] <- c(9, NA_real_, NaN, 2)
    declared <- haven::labelled_spss(stored, labels = c(small = 0, large = 1), na_values = 9)
    oracle <- binary
    oracle[3:6] <- NA_real_
    legacy <- structure(replace(stored, 3, NA_real_), class = "haven_labelled",
                        labels = c(small = 0, large = 1))
    sizes <- list(vector = haven::zap_missing(declared),
                  frame = haven::zap_missing(data.frame(TSIZE = declared))$TSIZE,
                  legacy = unserialize(serialize(legacy, NULL, version = 2)))

    for (method in c("parker.original", "cIHC", "cIHC.itr", "PCAPAM50", "ssBC", "manual", "AUTO")) {
        d <- input
        # Route AUTO to ssBC.v2 and the SSP methods on the entire cohort.
        if (method == "AUTO") d$se_NC$ER <- rep("ER+", length(ids))
        multi <- method %in% c("manual", "AUTO")
        methods <- if (method == "manual") c("PCAPAM50", "ssBC") else "AUTO"
        run <- function(size, warnings, four, clinical) {
            d$se_NC$TSIZE <- size
            meta <- serialize(SummarizedExperiment::colData(d$se_NC), NULL)
            rng <- .Random.seed
            call <- function() {
                if (multi) return(BS_Multi(d, methods = methods, Subtype = four, hasClinical = clinical))
                ror198_run(method, d$se_NC, subtype = four, clinical = clinical)
            }
            out <- expect_tsize_warning(expect_experimental_ror_pc(suppressMessages(call()), clinical),
                                        warnings)
            expect_identical(serialize(SummarizedExperiment::colData(d$se_NC), NULL), meta)
            expect_identical(.Random.seed, rng)
            out
        }
        for (four in c(FALSE, TRUE)) for (clinical in c(FALSE, TRUE)) {
            control <- run(oracle, character(), four, clinical)
            if (method == "AUTO") expect_true("ssBC.v2" %in% names(control$results))
            if (clinical) {
                for (res in if (multi) control$results else list(control)) {
                    if (is.null(res$score.ROR)) next
                    expect_identical(res, tsize_expected(res, ids[3:6]))
                }
                if (method == "parker.original")
                    expect_true(all(control$score.ROR[[tsize_groups[3]]][-c(3:6)] == "high"))
            }
            for (size in sizes) {
                actual <- run(size, if (clinical) ids[6] else character(), four, clinical)
                expect_identical(actual, control)
            }
        }
    }
})

test_that("TSIZE is decoded exactly and independently of NODE", {
    z <- tsize_input()$data_input$se_NC[, 12:1]
    ids <- colnames(z)
    z$PatientID <- factor(ids, levels = rev(ids))
    # Genuine missing values and malformed but nonmissing tokens have
    # identical numerical availability but distinct warning behavior.
    numeric_sizes <- c(0, 1, -1, 0.5, 2, 25, Inf, -Inf, 1 + .Machine$double.eps,
        .Machine$double.xmin, NA_real_, NaN)
    text_sizes <- c(" 1 ", "0", "1.0", "1e0", "T1", "25mm", "NA", "NaN", "",
        "1.00000000000000001", "1e-999", NA_character_)
    text_control <- c(1, 0, 1, 1, rep(0, 8))
    for (node in c(0, 1, 3, 4, NA_real_, -1)) {
        z$NODE <- rep(node, length(ids))
        run <- function(values, bad) {
            z$TSIZE <- values
            expr <- function() expect_tsize_warning(expect_experimental_ror_pc(
                ror198_run("parker.original", z)), ids[bad])
            if (!is.na(node) && node < 0) {
                expect_warning(result <- expr(), "NODE must give")
                return(result)
            }
            expr()
        }
        binary <- run(c(0, 1, rep(0, 10)), integer())
        expect_identical(run(numeric_sizes, 3:10), tsize_expected(binary, ids[3:12]))
        binary <- run(text_control, integer())
        for (values in list(text_sizes, factor(text_sizes, levels = rev(unique(text_sizes[!is.na(text_sizes)]))))) {
            expect_identical(run(values, 5:11), tsize_expected(binary, ids[5:12]))
        }
        expect_identical(run(as.logical(text_control), integer()), binary)
        expect_identical(run(as.integer(text_control), integer()), binary)
        expect_identical(run(rep(25, 12), seq_along(ids)), tsize_expected(binary, ids))
    }
})

test_that("standalone clinical routes count only patients evaluated by each method", {
    z <- tsize_input()$data_input$se_NC[, 70:1]
    ids <- colnames(z)
    n <- length(ids)
    z$PatientID <- factor(ids, levels = rev(ids))
    z$ER <- rep(c("ER+", "ER-"), each = n / 2)
    z$HER2 <- rep("HER2-", n)
    z$HER2[c(2, 3, 40)] <- "HER2+" # id 40 is a singleton for ER.v2
    z$TN <- rep("nonTN", n)
    z$TN[1:4] <- "TN"              # only the first four are evaluated by TN.v2
    z$NODE <- rep(4, n)
    bad <- c(1L, 40L, 60L)
    withr::local_seed(901)
    for (method in c("parker.original", "cIHC", "cIHC.itr", "PCAPAM50", "ER", "TN", "ER.v2", "TN.v2")) {
        expected_ids <- switch(method, ER.v2 = ids[c(1, 60)], TN.v2 = ids[1], ids[bad])
        for (clinical in c(FALSE, TRUE)) {
            run <- function(input, warnings) {
                do_run <- function() {
                    if (method %in% c("ER", "TN", "ER.v2", "TN.v2"))
                        return(BS_ssBC(input, s = method, Subtype = TRUE, hasClinical = clinical))
                    ror198_run(method, input, clinical = clinical, iterations = 2L)
                }
                call <- function() expect_tsize_warning(expect_experimental_ror_pc(do_run(), clinical),
                    warnings)
                if (method == "ER.v2") {
                    expect_warning(result <- call(), "has n=1")
                    return(result)
                }
                call()
            }
            z$TSIZE <- rep(1, n)
            rng <- .Random.seed
            reference <- run(z, character())
            z$TSIZE[bad] <- 25
            before <- SummarizedExperiment::colData(z)
            result <- run(z, if (clinical) expected_ids else character())
            expect_identical(SummarizedExperiment::colData(z), before)
            expect_identical(.Random.seed, rng)
            expect_identical(result, tsize_expected(reference, ids[bad]))
        }
    }
})

test_that("manual and AUTO report invalid TSIZE once across eligible clinical methods", {
    input <- lapply(tsize_input()$data_input, function(z) z[, 70:1, drop = FALSE])
    z <- input$se_NC
    ids <- colnames(z)
    z$PatientID <- factor(ids, levels = rev(ids))
    z$ER <- c(rep("ER+", 40), rep("ER-", 30))
    z$HER2 <- c(rep("HER2-", 3), rep(NA_character_, 37), "HER2+", rep(NA_character_, 29))
    z$NODE <- rep(4, 70)
    bad <- c(1, 41)
    withr::local_seed(717)
    for (methods in list(c("PCAPAM50", "ssBC.v2"), "AUTO")) {
        run <- function(input, warnings) {
            expect_warning(result <- expect_tsize_warning(expect_experimental_ror_pc(
                suppressMessages(BS_Multi(input, methods = methods, hasClinical = TRUE))),
                warnings), "has n=1")
            result
        }
        z$TSIZE <- rep(0, 70)
        input$se_NC <- z
        rng <- .Random.seed
        expected <- run(input, character())
        input$se_NC$TSIZE[bad] <- -1
        actual <- run(input, ids[bad])
        expected$results <- lapply(expected$results, tsize_expected, ids = ids[bad])
        expect_identical(actual, expected)
        expect_identical(.Random.seed, rng)
    }
})

test_that("invalid TSIZE honors warning-as-error even in the PCAPAM50 component", {
    input <- tsize_input()$data_input
    input$se_NC$TSIZE <- rep(25, ncol(input$se_NC))
    withr::local_options(list(warn = 2))
    expect_error(suppressMessages(BS_Multi(input, methods = c("PCAPAM50", "sspbc"), hasClinical = TRUE)),
        "TSIZE")
})
