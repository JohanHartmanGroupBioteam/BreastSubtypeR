cihc_itr_ids <- function(actual, expected) {
    complete <- is.character(actual) && length(actual) == length(expected) &&
        !anyNA(actual) && !anyDuplicated(actual) && all(nzchar(actual)) &&
        !anyNA(expected) && !anyDuplicated(expected) && setequal(actual, expected)
    expect_true(complete)
    if (!complete) stop("Incomplete sample identities in iterative cIHC comparison")
}

cihc_itr_outputs <- function(result, ids, subtype, iterations) {
    axes <- list(result$BS.all$PatientID, as.character(result$score.ROR$patientID),
        names(result$outList$predictions), rownames(result$outList$distances),
        rownames(result$outList$dist.RORSubtype), colnames(result$outList$testData),
        rownames(result$BS.itr.keep))
    if (subtype) axes <- c(axes, list(names(result$outList$predictions.FourSubtype),
        rownames(result$BS.itr.keep.Subtype)))
    for (axis in axes) {
        cihc_itr_ids(axis, ids)
        expect_identical(axis, ids)
    }
    consensus <- function(calls) names(which.max(table(calls)))
    expect_false(anyNA(result$BS.itr.keep))
    expect_identical(ncol(result$BS.itr.keep), iterations)
    expect_identical(apply(result$BS.itr.keep, 1, consensus), result$outList$predictions)
    expect_identical(unname(result$outList$predictions), result$BS.all$BS)
    expect_false(anyNA(result$BS.all$BS))
    if (subtype) {
        expect_false(anyNA(result$BS.itr.keep.Subtype))
        expect_identical(ncol(result$BS.itr.keep.Subtype), iterations)
        expect_identical(apply(result$BS.itr.keep.Subtype, 1, consensus),
            result$outList$predictions.FourSubtype)
        expect_identical(unname(result$outList$predictions.FourSubtype),
            result$BS.all$BS.Subtype)
        expect_false(anyNA(result$BS.all$BS.Subtype))
    }
}

cihc_itr_capture <- function(input, subtype = TRUE, clinical = TRUE,
                             policy = "target", ratio = 54 / 64, iterations = 2L) {
    captured <- list()
    native <- getFromNamespace("makeCalls_ihc.iterative", "BreastSubtypeR")
    record <- function(e) {
        captured[[e$itr]] <<- list(input = e$.input_expression,
            clinical = get("df.cln", e), plan = get("plan", e),
            intended = as.character(e$keep), selected = e$mbal.ihc,
            mdns = e$mdns, medians = e$medians.all, output = e$out)
    }
    found <- 0L
    walk <- function(x) {
        if (missing(x)) return(quote(expr = ))
        if (!is.call(x)) return(x)
        if (identical(x[[1]], as.name("function")) &&
            identical(names(x[[2]]), "itr")) {
            found <<- found + 1L
            x[[3]] <- bquote({
                .input_expression <- mat
                on.exit(.record(environment()), add = TRUE)
                .(x[[3]])
            })
            return(x)
        }
        as.call(lapply(as.list(x), walk))
    }
    body(native) <- walk(body(native))
    stopifnot(found == 1L)
    environment(native) <- list2env(list(.record = record), parent = environment(native))
    result <- testthat::with_mocked_bindings(
        expect_experimental_ror_pc(BS_cIHC.itr(input, Subtype = subtype, hasClinical = clinical,
            er_subsample = policy, ratio = ratio, iteration = iterations, seed = 118),
            n = as.integer(clinical)),
        makeCalls_ihc.iterative = native, .package = "BreastSubtypeR")
    list(result = result, captured = captured)
}

cihc_itr_centering <- function(run, input, iterations) {
    ids <- colnames(input)
    expression <- SummarizedExperiment::assay(input)
    clinical <- SummarizedExperiment::colData(input)
    expect_identical(names(run$captured), paste0("itr.", seq_len(iterations)))
    for (z in run$captured) {
        cihc_itr_ids(colnames(z$input), ids)
        cihc_itr_ids(as.character(z$clinical$PatientID), ids)
        expect_identical(z$input, expression)
        expect_false(anyNA(z$intended))
        expect_identical(anyDuplicated(z$intended), 0L)
        expect_true(all(z$intended %in% ids))
        expect_identical(colnames(z$selected), z$intended)
        expect_identical(z$selected, expression[, z$intended, drop = FALSE])
        expect_identical(z$mdns, apply(z$selected, 1, median, na.rm = TRUE))
        matched <- match(ids, as.character(z$clinical$PatientID))
        expect_false(anyNA(matched))
        for (column in c("ER", "TSIZE", "NODE"))
            expect_identical(z$clinical[[column]][matched], clinical[[column]])
        selected_er <- clinical$ER[match(z$intended, ids)]
        if (z$plan$group == "none") {
            cihc_itr_ids(z$intended, ids)
        } else {
            expect_identical(sum(selected_er == z$plan$group), z$plan$size)
            other <- if (z$plan$group == "ER+") "ER-" else "ER+"
            expect_identical(z$intended[selected_er == other], ids[clinical$ER == other])
        }
    }
}

test_that("native iterative cIHC selects each centering matrix by biological ID", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR")
    input <- OSLO2EMIT0obj$data_input$se_NC
    ids <- colnames(input)
    reference <- expect_silent(cihc_itr_capture(input))
    plain <- expect_silent(expect_experimental_ror_pc(BS_cIHC.itr(input, iteration = 2L,
        Subtype = TRUE, hasClinical = TRUE)))
    expect_identical(reference$result, plain)
    cihc_itr_outputs(reference$result, ids, TRUE, 2L)
    cihc_itr_centering(reference, input, 2L)
    SummarizedExperiment::colData(input)$PatientID <- factor(ids, levels = rev(ids))
    actual <- expect_silent(cihc_itr_capture(input))
    cihc_itr_outputs(actual$result, ids, TRUE, 2L)
    cihc_itr_centering(actual, input, 2L)
    expect_identical(actual$result, reference$result)
    for (i in seq_len(2L)) {
        expect_identical(actual$captured[[i]]$intended, reference$captured[[i]]$intended)
        expect_identical(actual$captured[[i]]$medians, reference$captured[[i]]$medians)
    }
})

test_that("iterative factor IDs preserve both policies and all output modes", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR")
    original <- OSLO2EMIT0obj$data_input$se_NC
    all_ids <- colnames(original)
    er <- SummarizedExperiment::colData(original)$ER
    subsets <- list(seq_along(all_ids), rev(seq_along(all_ids)),
        c(which(er == "ER+")[1:12], which(er == "ER-")),
        c(which(er == "ER+")[1:9], which(er == "ER-")))
    shuffled <- withr::with_seed(201, sample(all_ids))
    for (j in seq_along(subsets)) {
        character_input <- original[, subsets[[j]], drop = FALSE]
        ids <- colnames(character_input)
        ratio <- if (j <= 2L) 54 / 64 else 0.5
        for (policy in c("target", "larger")) {
            at_target <- j == 4L && policy == "target"
            for (subtype in c(FALSE, TRUE)) for (clinical in c(FALSE, TRUE)) {
                invoke <- function(input) cihc_itr_capture(input, subtype, clinical, policy, ratio)
                if (at_target) {
                    expect_message(reference <- invoke(character_input), "already matches")
                } else reference <- expect_silent(invoke(character_input))
                cihc_itr_outputs(reference$result, ids, subtype, 2L)
                cihc_itr_centering(reference, character_input, 2L)
                encodings <- list(factor(ids), factor(ids, levels = rev(all_ids)),
                    factor(ids, levels = c("unused-before", all_ids, "unused-after")),
                    factor(ids, levels = shuffled))
                for (encoding in encodings) {
                    input <- character_input
                    SummarizedExperiment::colData(input)$PatientID <- encoding
                    saved <- input
                    if (at_target) {
                        expect_message(actual <- invoke(input), "already matches")
                    } else actual <- expect_silent(invoke(input))
                    expect_identical(input, saved)
                    cihc_itr_outputs(actual$result, ids, subtype, 2L)
                    cihc_itr_centering(actual, input, 2L)
                    expect_identical(actual$result, reference$result)
                    expect_identical(is.na(actual$result$score.ROR),
                        is.na(reference$result$score.ROR))
                    expect_identical(actual$captured, reference$captured)
                }
            }
        }
    }
})

test_that("iterative cIHC rejects malformed IDs before sampling", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR")
    original <- OSLO2EMIT0obj$data_input$se_NC
    ids <- colnames(original)
    cases <- list(missing = replace(ids, 1, NA_character_),
        empty = replace(ids, 1, ""), duplicated = replace(ids, 1, ids[2]),
        unmatched = replace(ids, 1, "not-in-assay"))
    patterns <- c(missing = "must not contain missing or empty",
        empty = "must not contain missing or empty",
        duplicated = "must be unique; duplicated", unmatched = "not found among assay column names")
    for (name in names(cases)) for (as_factor in c(FALSE, TRUE)) {
        input <- original
        SummarizedExperiment::colData(input)$PatientID <-
            if (as_factor) factor(cases[[name]]) else cases[[name]]
        expect_error(BS_cIHC.itr(input, iteration = 2L), patterns[[name]])
    }
})

test_that("iterative ID conversion preserves RNG kinds and absent caller state", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR")
    input <- OSLO2EMIT0obj$data_input$se_NC
    ids <- colnames(input)
    SummarizedExperiment::colData(input)$PatientID <-
        factor(ids, levels = c("unused", rev(ids)))
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
    for (kind in c("Mersenne-Twister", "L'Ecuyer-CMRG")) {
        RNGkind(kind, "Inversion", "Rejection")
        for (present in c(FALSE, TRUE)) for (policy in c("target", "larger")) {
            if (present) set.seed(201) else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
                rm(".Random.seed", envir = .GlobalEnv)
            before <- get0(".Random.seed", .GlobalEnv, inherits = FALSE)
            before_kind <- RNGkind()
            actual <- expect_silent(expect_experimental_ror_pc(BS_cIHC.itr(input, iteration = 2L,
                Subtype = TRUE, hasClinical = TRUE, er_subsample = policy)))
            expect_identical(exists(".Random.seed", .GlobalEnv, inherits = FALSE), present)
            expect_identical(get0(".Random.seed", .GlobalEnv, inherits = FALSE), before)
            expect_identical(RNGkind(), before_kind)
            expect_identical(Sys.getlocale(), locale)
            cihc_itr_outputs(actual, ids, TRUE, 2L)
        }
    }
})
