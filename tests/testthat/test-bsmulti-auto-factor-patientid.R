# Native AUTO restoration tests using public packaged expressions.
m202_predictor_values <- function(predictors) {
    lapply(predictors, function(x) {
        x$clinical$PatientID <- as.character(x$clinical$PatientID)
        x
    })
}
m202_valid_ids <- function(x, expected = x) {
    stopifnot(is.character(x), length(x) == length(expected), !anyNA(x), !anyDuplicated(x),
        all(nzchar(x)), !anyNA(expected), !anyDuplicated(expected), setequal(x, expected))
    invisible(TRUE)
}
m202_capture_run <- function(input, subtype = TRUE, clinical = TRUE, methods = "AUTO") {
    traces <- list()
    predictors <- list()
    selection <- NULL
    native_multi <- BreastSubtypeR::BS_Multi
    native_predict <- BreastSubtypeR::BS_ssBC
    native_select <- getFromNamespace("get_methods", "BreastSubtypeR")
    record <- function(tag, e) {
        method <- get("method", e)
        name <- if (method == "ssBC")
            "res_ssBC"
        else "res_ssBC.v2"
        traces[[paste(method, tag, sep = "/")]] <<- list(result = get(name, e)$BS.all, na_rows = get("na_df",
            e), pheno_ids = as.character(get("pheno", e)$PatientID))
    }
    hits <- 0L
    walk <- function(x) {
        if (missing(x))
            return(quote(expr = ))
        if (!is.call(x))
            return(x)
        if (identical(x[[1]], as.name("<-")) && any(vapply(c("res_ssBC", "res_ssBC.v2"), function(n) identical(x[[2]],
            substitute(A$BS.all, list(A = as.name(n)))), logical(1)))) {
            tag <- if (identical(x[[3]][[1]], as.name("rbind")))
                "before"
            else "after"
            hits <<- hits + 1L
            return(bquote({
                .(x)
                .record(.(tag), environment())
            }))
        }
        as.call(lapply(as.list(x), walk))
    }
    body(native_multi) <- walk(body(native_multi))
    stopifnot(hits == 4L)
    environment(native_multi) <- list2env(list(.record = record), parent = environment(native_multi))
    predictor <- function(se_obj, s, Subtype, hasClinical) {
        actual <- native_predict(se_obj, s = s, Subtype = Subtype, hasClinical = hasClinical)
        predictors[[s]] <<- list(ids = colnames(se_obj), expression = SummarizedExperiment::assay(se_obj),
            clinical = as.data.frame(SummarizedExperiment::colData(se_obj)), result = actual)
        actual
    }
    selector <- function(pheno) {
        selection <<- native_select(pheno)
        selection
    }
    messages <- warnings <- character()
    result <- withCallingHandlers(testthat::with_mocked_bindings(native_multi(input, methods = methods,
        Subtype = subtype, hasClinical = clinical), BS_ssBC = predictor, get_methods = selector,
        .package = "BreastSubtypeR"), message = function(m) {
        messages <<- c(messages, conditionMessage(m))
        invokeRestart("muffleMessage")
    }, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    list(result = result, traces = traces, predictors = predictors, selection = selection,
        messages = messages, warnings = warnings)
}
m202_encode_input <- function(input, encoding, universe = colnames(input$se_NC)) {
    ids <- colnames(input$se_NC)
    m202_valid_ids(ids)
    levels <- switch(encoding, character = NULL, default = sort(universe), reversed = rev(universe),
        unused = c(paste0("unused-", 1:5), universe), shuffled = withr::with_seed(202, sample(universe)))
    for (n in names(input)) if (inherits(input[[n]], "SummarizedExperiment")) {
        stopifnot(identical(colnames(input[[n]]), ids))
        SummarizedExperiment::colData(input[[n]])$PatientID <- if (encoding == "character")
            ids
        else factor(ids, levels = levels)
    }
    input
}
m202_make_cases <- function(original) {
    subset_input <- function(x, idx) lapply(x, function(z) if (inherits(z, "SummarizedExperiment"))
        z[, idx, drop = FALSE]
    else z)
    annotate <- function(x, er, her2, tn = "nonTN") {
        for (n in names(x)) if (inherits(x[[n]], "SummarizedExperiment")) {
            SummarizedExperiment::colData(x[[n]])$ER <- er
            SummarizedExperiment::colData(x[[n]])$HER2 <- her2
            SummarizedExperiment::colData(x[[n]])$TN <- rep(tn, length.out = ncol(x[[n]]))
        }
        x
    }
    both <- annotate(subset_input(original, 1:30), rep(c("ER+", "ER-"), c(20, 10)), rep(c("HER2+",
        "HER2-", "HER2+", "HER2-"), c(10, 10, 5, 5)))
    complete <- annotate(subset_input(original, 1:17), rep("ER+", 17), rep(c("HER2+", "HER2-"),
        c(9, 8)))
    tn <- annotate(subset_input(original, 1:18), rep("ER-", 18), rep("HER2-", 18), "TN")
    make <- function(name, x, subtype = TRUE, clinical = TRUE, methods = "AUTO", universe = colnames(x$se_NC)) list(name = name,
        input = x, subtype = subtype, clinical = clinical, methods = methods, universe = universe)
    cases <- list(make("original", original), make("reversed", subset_input(original, ncol(original$se_NC):1)))
    for (s in c(FALSE, TRUE)) for (c in c(FALSE, TRUE)) cases[[length(cases) + 1L]] <- make(paste0("both-",
        s, "-", c), both, s, c)
    cases <- c(cases, list(make("subset-retained", subset_input(both, c(1:17, 21:30)), universe = colnames(both$se_NC)),
        make("complete-auto", complete), make("tn-auto", tn), make("manual", both, methods = c("ssBC",
            "ssBC.v2")), make("manual-nonclinical", both, FALSE, FALSE, c("ssBC", "ssBC.v2")),
        make("tn-manual", tn, methods = c("ssBC", "ssBC.v2"))))
    names(cases) <- vapply(cases, function(x) x$name, character(1))
    cases
}

test_that("native AUTO restoration keeps eligible calls with unused factor levels", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR")
    cs <- m202_make_cases(OSLO2EMIT0obj$data_input)[["both-TRUE-TRUE"]]
    reference <- m202_capture_run(m202_encode_input(cs$input, "character"))
    actual <- m202_capture_run(m202_encode_input(cs$input, "unused"))
    ids <- colnames(cs$input$se_NC)
    expect_length(actual$warnings, 0L)
    expect_identical(actual$selection, reference$selection)
    expect_identical(m202_predictor_values(actual$predictors),
        m202_predictor_values(reference$predictors))
    for (method in c("ssBC", "ssBC.v2")) {
        before <- actual$traces[[paste0(method, "/before")]]
        m202_valid_ids(as.character(before$result$PatientID), ids)
        expect_identical(before, reference$traces[[paste0(method, "/before")]])
        expect_identical(actual$result$results[[method]]$BS.all$PatientID, ids)
        expect_identical(actual$result$results[[method]]$BS.all,
            reference$result$results[[method]]$BS.all)
        expect_identical(actual$result$res_subtypes[[method]],
            reference$result$res_subtypes[[method]])
        expect_identical(actual$result$res_subtypes.Subtype[[method]],
            reference$result$res_subtypes.Subtype[[method]])
    }
})

m202_assert_restoration <- function(run, input, subtype) {
    ids <- colnames(input$se_NC)
    expect_identical(rownames(run$result$res_subtypes), ids)
    if (subtype) expect_identical(rownames(run$result$res_subtypes.Subtype), ids)
    expect_length(run$warnings, 0L)
    for (route in names(run$predictors)) {
        method <- if (route %in% c("ER", "TN")) "ssBC" else "ssBC.v2"
        predictor <- run$predictors[[route]]
        eligible <- predictor$ids
        m202_valid_ids(eligible)
        expect_true(all(eligible %in% ids))
        m202_valid_ids(as.character(predictor$result$BS.all$PatientID), eligible)
        expect_identical(predictor$expression,
            SummarizedExperiment::assay(input$se_NC[, eligible, drop = FALSE]))
        expect_identical(predictor$clinical,
            as.data.frame(SummarizedExperiment::colData(input$se_NC[, eligible, drop = FALSE])))
        calls <- run$result$results[[method]]$BS.all
        m202_valid_ids(as.character(calls$PatientID), ids)
        expect_identical(calls$PatientID, ids)
        expect_identical(rownames(calls), ids)
        for (column in c("BS", if (subtype) "BS.Subtype")) {
            expected <- rep(NA_character_, length(ids))
            expected[match(predictor$result$BS.all$PatientID, ids)] <-
                as.character(predictor$result$BS.all[[column]])
            observed <- if (column == "BS") run$result$res_subtypes[[method]] else
                run$result$res_subtypes.Subtype[[method]]
            expect_identical(observed, expected)
            expect_identical(as.character(calls[[column]]), expected)
            expect_identical(is.na(observed), is.na(expected))
        }
        if (length(eligible) < length(ids)) {
            before <- run$traces[[paste0(method, "/before")]]
            m202_valid_ids(as.character(before$na_rows$PatientID), setdiff(ids, eligible))
            expect_true(all(is.na(before$na_rows[, -1, drop = FALSE])))
            expect_identical(run$traces[[paste0(method, "/after")]]$result, calls)
        }
    }
}

test_that("AUTO factor restoration preserves modes, clinical values and predictor results", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR")
    cases <- m202_make_cases(OSLO2EMIT0obj$data_input)
    cases <- cases[grepl("^both-", names(cases))]
    for (cs in cases) {
        reference <- m202_capture_run(cs$input, cs$subtype, cs$clinical)
        m202_assert_restoration(reference, cs$input, cs$subtype)
        for (encoding in c("default", "reversed", "unused", "shuffled")) {
            input <- m202_encode_input(cs$input, encoding)
            saved <- input
            actual <- m202_capture_run(input, cs$subtype, cs$clinical)
            expect_identical(input, saved)
            expect_identical(actual$selection, reference$selection)
            expect_identical(actual$result, reference$result)
            expect_identical(actual$messages, reference$messages)
            expect_identical(m202_predictor_values(actual$predictors),
                m202_predictor_values(reference$predictors))
            m202_assert_restoration(actual, input, cs$subtype)
        }
    }
})

test_that("manual, TN and complete AUTO paths preserve routing and caller state", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR")
    cases <- m202_make_cases(OSLO2EMIT0obj$data_input)
    cases <- cases[c("subset-retained", "complete-auto", "tn-auto", "manual",
        "manual-nonclinical", "tn-manual")]
    old_kind <- RNGkind()
    old_seed <- get0(".Random.seed", .GlobalEnv, inherits = FALSE)
    locale <- Sys.getlocale()
    withr::defer({
        do.call(RNGkind, as.list(old_kind))
        if (is.null(old_seed)) {
            if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
                rm(".Random.seed", envir = .GlobalEnv)
        } else assign(".Random.seed", old_seed, .GlobalEnv)
    })
    for (j in seq_along(cases)) {
        cs <- cases[[j]]
        for (name in names(cs$input)) if (inherits(cs$input[[name]], "SummarizedExperiment"))
            cs$input[[name]] <- cs$input[[name]][, ncol(cs$input[[name]]):1, drop = FALSE]
        reference <- m202_capture_run(cs$input, cs$subtype, cs$clinical, cs$methods)
        input <- m202_encode_input(cs$input, "unused", cs$universe)
        saved <- input
        RNGkind(if (j %% 2L) "Mersenne-Twister" else "L'Ecuyer-CMRG")
        if (j %% 2L) set.seed(202) else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv)
        before <- get0(".Random.seed", .GlobalEnv, inherits = FALSE)
        before_kind <- RNGkind()
        actual <- m202_capture_run(input, cs$subtype, cs$clinical, cs$methods)
        expect_identical(get0(".Random.seed", .GlobalEnv, inherits = FALSE), before)
        expect_identical(RNGkind(), before_kind)
        expect_identical(Sys.getlocale(), locale)
        expect_identical(input, saved)
        expect_identical(actual$result, reference$result)
        expect_identical(actual$selection, reference$selection)
        expect_identical(names(actual$predictors), names(reference$predictors))
        if (grepl("^tn-", cs$name)) expect_identical(names(actual$predictors), c("TN", "TN.v2"))
        m202_assert_restoration(actual, input, cs$subtype)
        plain <- suppressMessages(BS_Multi(input, methods = cs$methods,
            Subtype = cs$subtype, hasClinical = cs$clinical))
        expect_identical(plain, actual$result)
    }
})
