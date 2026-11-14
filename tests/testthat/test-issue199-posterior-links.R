# Synthetic posteriors use an independent which.max selection oracle.
posterior199_object <- function(n = 3L) {
    ids <- if (n) paste0("patient /", seq_len(n)) else character()
    make <- function(k) {
        p <- rbind(c(.5, .25, .125, .0625, .0625),
                   c(.75, .25, 0, 0, 0),
                   c(.125, .5, .25, .0625, .0625))
        p <- p[rep(seq_len(3L), length.out = n), , drop = FALSE]
        dimnames(p) <- list(NULL, c("Basal", "Her2", "LumA", "LumB", "Normal"))
        j <- vapply(seq_len(n), function(i) which.max(p[i, ]), integer(1))
        list(cl = matrix(colnames(p)[j], n, 1L, dimnames = list(ids, k)),
             prob = matrix(p[cbind(seq_len(n), j)], n, 1L, dimnames = list(ids, k)),
             all.probs = setNames(list(p), k),
             rules.matrix = matrix(rep(c(0, 1), length.out = 2L * n), 2L, n,
                 dimnames = list(c("a<b", "c<d"), NULL)),
             data.used = matrix(1, 2L, n, dimnames = list(c("a", "b"), ids)))
    }
    list(res_subtypes = data.frame(AIMS = rep("synthetic", n),
             sspbc = rep("synthetic", n), row.names = ids),
         results = list(AIMS = make("20"), sspbc = make("24")), unlisted = 1)
}

posterior199_relink <- function(x, method) {
    o <- x$results[[method]]
    for (k in names(o$all.probs)) {
        p <- o$all.probs[[k]]
        for (i in seq_len(nrow(p))) {
            if (all(is.na(p[i, ]))) {
                o$prob[i, k] <- NA_real_
                o$cl[i, k] <- NA_character_
            } else {
                j <- which.max(p[i, ])
                o$prob[i, k] <- p[i, j]
                o$cl[i, k] <- colnames(p)[j]
            }
        }
    }
    x$results[[method]] <- o
    x
}

# Capture all decision expectations, including deliberate failures, separately.
posterior199_decisions <- function(actual, expected) {
    reporter <- testthat::ListReporter$new()
    reporter$start_file("posterior-link-control")
    testthat::with_reporter(reporter, testthat::test_that("decision layer", {
        issue199_expect_decisions(actual, expected)
    }))
    result <- as.data.frame(reporter$get_results())
    nrow(result) == 1L && result$passed > 0L &&
        !any(result$failed | result$error | result$warning | result$skipped)
}
posterior199_compare <- function(actual, expected) {
    issue199_portability_unit(actual, expected, "OSLO2EMIT0obj", path = "/res")
}
posterior199_accept <- function(actual, expected = actual) {
    expect_true(posterior199_compare(actual, expected))
    expect_true(posterior199_decisions(actual, expected))
}

test_that("complete posterior links accept approved and bounded coherent results", {
    reference <- posterior199_object()
    posterior199_accept(reference)
    for (method in names(reference$results)) {
        actual <- reference
        actual$results[[method]]$all.probs[[1L]][1L, 1L] <- .5 + 2^-53
        actual <- posterior199_relink(actual, method)
        expect_false(identical(actual$results[[method]]$prob,
            reference$results[[method]]$prob, num.eq = FALSE))
        posterior199_accept(actual, reference)
    }
    for (dataset in c("OSLO2EMIT0obj", "TCGABRCAobj")) {
        e <- new.env(parent = emptyenv())
        utils::data(list = dataset, package = "BreastSubtypeR", envir = e)
        x <- e[[dataset]]$res
        expect_true(issue199_portable_compare(x, x, dataset, path = "/res"))
        expect_true(posterior199_decisions(x, x))
        for (method in c("AIMS", "sspbc")) {
            expect_true(issue199_posterior_link(x$results[[method]], method))
            actual <- x
            p <- actual$results[[method]]$all.probs[[1L]]
            j <- which.max(p[1L, ])
            p[1L, j] <- p[1L, j] - abs(p[1L, j]) * .Machine$double.eps
            actual$results[[method]]$all.probs[[1L]] <- p
            actual <- posterior199_relink(actual, method)
            expect_false(identical(actual$results[[method]]$prob,
                x$results[[method]]$prob, num.eq = FALSE))
            expect_true(issue199_portable_compare(actual, x, dataset, path = "/res"))
            expect_true(posterior199_decisions(actual, x))
        }
    }
})

test_that("posterior links support first ties, whole masks, singletons and empty comparisons", {
    for (n in c(0L, 1L, 3L)) {
        x <- posterior199_object(n)
        posterior199_accept(x)
        for (method in names(x$results)) {
            if (n) {
                x$results[[method]]$all.probs[[1L]][1L, ] <- c(.5, .5, 0, 0, 0)
                x <- posterior199_relink(x, method)
                expect_identical(x$results[[method]]$cl[1L, 1L], "Basal")
            }
        }
        posterior199_accept(x)
        for (method in names(x$results)) {
            if (n) {
                x$results[[method]]$all.probs[[1L]][1L, ] <- NA_real_
                x$results[[method]]$rules.matrix[, 1L] <- NA_real_
                x <- posterior199_relink(x, method)
            }
        }
        posterior199_accept(x)
    }
    # Zero rows are a helper-only comparison contract, not a public classifier call.
})

test_that("selected posterior values and calls must match each parent bit exactly", {
    reference <- posterior199_object()
    for (method in names(reference$results)) {
        for (kind in c("ulp", "nonmaximum", "class", "prob_na", "cl_na",
                       "prob_nan", "prob_inf")) {
            actual <- reference
            o <- actual$results[[method]]
            if (kind == "ulp") o$prob[1L, 1L] <- .5 + 2^-53
            if (kind == "nonmaximum") o$prob[1L, 1L] <- .25
            if (kind == "class") o$cl[1L, 1L] <- "LumB"
            if (kind == "prob_na") o$prob[1L, 1L] <- NA_real_
            if (kind == "cl_na") o$cl[1L, 1L] <- NA_character_
            if (kind == "prob_nan") o$prob[1L, 1L] <- NaN
            if (kind == "prob_inf") o$prob[1L, 1L] <- Inf
            actual$results[[method]] <- o
            expect_error(posterior199_compare(actual, reference), "selected probability/class link",
                info = paste(method, kind))
            expect_false(posterior199_decisions(actual, reference), info = paste(method, kind))
            # The reference is validated as strictly as the actual object.
            expect_error(posterior199_compare(reference, actual), "selected probability/class link")
            expect_false(posterior199_decisions(reference, actual))
        }
        actual <- reference
        actual$results[[method]]$all.probs[[1L]][1L, ] <- 0
        actual <- posterior199_relink(actual, method)
        actual$results[[method]]$prob[1L, 1L] <- -0
        expect_error(issue199_posterior_link(actual$results[[method]], method),
            "selected probability/class link")
    }
})

test_that("posterior ties, ranks and first argmax remain exact across objects", {
    for (method in c("AIMS", "sspbc")) {
        for (kind in c("make_tie", "break_tie", "rank", "argmax")) {
            reference <- posterior199_object()
            p <- reference$results[[method]]$all.probs[[1L]]
            p[1L, ] <- if (kind == "break_tie") c(.5, .5, 0, 0, 0) else
                if (kind == "rank") c(.5, .25, .25 - 2^-55, 0, 0) else
                    c(.5, .5 - 2^-54, 0, 0, 0)
            reference$results[[method]]$all.probs[[1L]] <- p
            reference <- posterior199_relink(reference, method)
            actual <- reference
            if (kind == "make_tie") p[1L, 2L] <- .5
            if (kind == "break_tie") p[1L, 2L] <- .5 - 2^-54
            if (kind == "rank") p[1L, 2:3] <- rev(p[1L, 2:3])
            if (kind == "argmax") p[1L, 1:2] <- rev(p[1L, 1:2])
            actual$results[[method]]$all.probs[[1L]] <- p
            actual <- posterior199_relink(actual, method)
            expect_error(posterior199_compare(actual, reference),
                "posterior ranks/ties/argmax|unlisted value must be exact", info = kind)
            expect_false(posterior199_decisions(actual, reference), info = kind)
        }
        for (value in c(2^-1074, 1e-300, -0)) {
            reference <- posterior199_object()
            actual <- reference
            actual$results[[method]]$all.probs[[1L]][2L, 3L] <- value
            expect_error(posterior199_compare(actual, reference),
                "zero support/signs|signed zero")
            if (value != 0) expect_false(posterior199_decisions(actual, reference))
        }
    }
})

test_that("partial missing rows, NaN and infinite posteriors are refused", {
    for (method in c("AIMS", "sspbc")) {
        for (kind in c("partial", "mask_only", "nan", "inf", "masked_nan", "masked_call")) {
            reference <- posterior199_object()
            actual <- reference
            o <- actual$results[[method]]
            if (kind == "partial") o$all.probs[[1L]][1L, 2L] <- NA_real_
            if (kind == "mask_only") o$all.probs[[1L]][1L, ] <- NA_real_
            if (kind == "nan") o$all.probs[[1L]][1L, 2L] <- NaN
            if (kind == "inf") o$all.probs[[1L]][1L, 2L] <- Inf
            if (kind %in% c("masked_nan", "masked_call")) {
                o$all.probs[[1L]][1L, ] <- NA_real_
                o$prob[1L, 1L] <- if (kind == "masked_nan") NaN else NA_real_
                o$cl[1L, 1L] <- if (kind == "masked_call") "Basal" else NA_character_
            }
            actual$results[[method]] <- o
            expect_error(posterior199_compare(actual, reference),
                "posterior rows|selected probability/class link", info = kind)
            expect_false(posterior199_decisions(actual, reference), info = kind)
        }
        reference <- posterior199_object()
        actual <- reference
        actual$results[[method]]$all.probs[[1L]][1L, ] <- NA_real_
        actual <- posterior199_relink(actual, method)
        expect_error(posterior199_compare(actual, reference), "derived posterior|unlisted value")
        expect_false(posterior199_decisions(actual, reference))
    }
})

test_that("posterior models, types, sample and class axes are structural contracts", {
    reference <- posterior199_object()
    for (method in names(reference$results)) {
        for (kind in c("prob_vector", "prob_integer", "cl_numeric", "transpose",
                       "prob_rows", "posterior_rows", "rules_count", "rules_names",
                       "posterior_names", "class_empty", "class_duplicate", "class_na",
                       "class_absent", "model_empty", "model_duplicate", "model_order")) {
            actual <- reference
            o <- actual$results[[method]]
            if (kind == "prob_vector") o$prob <- as.vector(o$prob)
            if (kind == "prob_integer") storage.mode(o$prob) <- "integer"
            if (kind == "cl_numeric") o$cl <- matrix(0, nrow(o$cl), ncol(o$cl),
                dimnames = dimnames(o$cl))
            if (kind == "transpose") o$prob <- t(o$prob)
            if (kind == "prob_rows") o$prob <- o$prob[-1L, , drop = FALSE]
            if (kind == "posterior_rows") o$all.probs[[1L]] <- o$all.probs[[1L]][-1L, , drop = FALSE]
            if (kind == "rules_count") o$rules.matrix <- o$rules.matrix[, -1L, drop = FALSE]
            if (kind == "rules_names") colnames(o$rules.matrix) <- rev(rownames(o$prob))
            if (kind == "posterior_names") rownames(o$all.probs[[1L]]) <- rev(rownames(o$prob))
            if (kind == "class_empty") colnames(o$all.probs[[1L]])[1L] <- ""
            if (kind == "class_duplicate") colnames(o$all.probs[[1L]])[1L] <- "Her2"
            if (kind == "class_na") colnames(o$all.probs[[1L]])[1L] <- NA_character_
            if (kind == "class_absent") colnames(o$all.probs[[1L]]) <- NULL
            if (kind == "model_empty") names(o$all.probs) <- ""
            if (kind == "model_duplicate") o$all.probs <- c(o$all.probs, o$all.probs)
            if (kind == "model_order") names(o$all.probs) <- "other"
            actual$results[[method]] <- o
            # Self-comparison must not make a malformed complete result valid.
            expect_error(posterior199_compare(actual, actual),
                "structure|axes|axis|names", info = paste(method, kind))
            expect_false(posterior199_decisions(actual, actual), info = paste(method, kind))
        }
        for (field in c("prob", "cl", "all.probs", "rules.matrix")) {
            for (rename in c(FALSE, TRUE)) {
                actual <- reference
                if (rename) names(actual$results[[method]])[
                    names(actual$results[[method]]) == field] <- paste0(field, "_renamed")
                else actual$results[[method]][[field]] <- NULL
                expect_error(posterior199_compare(actual, actual), "posterior fields")
                expect_false(posterior199_decisions(actual, actual))
            }
        }
        actual <- reference
        actual$results[[method]][c("prob", "cl")] <- NULL
        expect_error(posterior199_compare(actual, actual), "posterior fields")
        expect_false(posterior199_decisions(actual, actual))
        actual <- reference
        actual$results[[method]]$all.probs[[1L]] <-
            actual$results[[method]]$all.probs[[1L]][, -5L, drop = FALSE]
        expect_error(posterior199_compare(actual, reference), "type/class/length")
        expect_false(posterior199_decisions(actual, reference))
        actual <- reference
        attr(actual$results[[method]]$prob, "extra") <- "changed"
        expect_error(posterior199_compare(actual, reference), "attributes/names/dimensions")
        expect_false(posterior199_decisions(actual, reference))
        actual <- reference
        o <- actual$results[[method]]
        o$prob <- o$prob[3:1, , drop = FALSE]
        o$cl <- o$cl[3:1, , drop = FALSE]
        o$all.probs <- lapply(o$all.probs, function(p) p[3:1, , drop = FALSE])
        o$rules.matrix <- o$rules.matrix[, 3:1, drop = FALSE]
        actual$results[[method]] <- o
        expect_error(posterior199_compare(actual, reference), "attributes/names/dimensions")
        expect_false(posterior199_decisions(actual, reference))
        actual <- reference
        actual$results[[method]]$all.probs[[1L]] <-
            actual$results[[method]]$all.probs[[1L]][3:1, , drop = FALSE]
        expect_error(posterior199_compare(actual, reference), "selected probability/class link")
        expect_false(posterior199_decisions(actual, reference))
    }
})

test_that("every posterior model is linked and original numeric safeguards remain active", {
    reference <- posterior199_object()
    for (method in names(reference$results)) {
        o <- reference$results[[method]]
        o$prob <- cbind(o$prob, extra = o$prob[, 1L])
        o$cl <- cbind(o$cl, extra = o$cl[, 1L])
        o$all.probs$extra <- o$all.probs[[1L]]
        reference$results[[method]] <- o
    }
    posterior199_accept(reference)
    for (method in names(reference$results)) {
        actual <- reference
        actual$results[[method]]$prob <- actual$results[[method]]$prob[, 2:1, drop = FALSE]
        actual$results[[method]]$cl <- actual$results[[method]]$cl[, 2:1, drop = FALSE]
        actual$results[[method]]$all.probs <- actual$results[[method]]$all.probs[2:1]
        expect_error(posterior199_compare(actual, reference), "attributes/names/dimensions")
        expect_false(posterior199_decisions(actual, reference))
        actual <- reference
        actual$results[[method]]$prob[1L, "extra"] <- .5 + 2^-53
        expect_error(posterior199_compare(actual, reference), "selected probability/class link")
        expect_false(posterior199_decisions(actual, reference))
        actual <- reference
        actual$results[[method]]$all.probs$extra[1L, 1L] <- .5 + 2^-53
        actual <- posterior199_relink(actual, method)
        expect_error(posterior199_compare(actual, reference), "unlisted value must be exact")
        actual <- reference
        actual$results[[method]]$all.probs[[1L]][1L, 1L] <- .5 + 2^-40
        actual <- posterior199_relink(actual, method)
        expect_true(posterior199_decisions(actual, reference))
        expect_error(posterior199_compare(actual, reference), "numeric bound exceeded")
        actual <- reference
        actual$results[[method]]$rules.matrix[1L, 1L] <- 1
        expect_error(posterior199_compare(actual, reference), "unlisted value must be exact")
        expect_false(posterior199_decisions(actual, reference))
        actual$results[[method]]$rules.matrix[1L, 1L] <- NA_real_
        expect_error(posterior199_compare(actual, reference), "unlisted value must be exact")
        expect_false(posterior199_decisions(actual, reference))
    }
    actual <- reference; actual$unlisted <- 1 + .Machine$double.eps
    expect_error(posterior199_compare(actual, reference), "unlisted value must be exact")
    actual <- reference; actual$results$AIMS$data.used <- NULL
    expect_error(posterior199_compare(actual, actual), "unused.*policy path")
    actual <- reference; actual$results <- c(actual$results, actual$results["AIMS"])
    expect_error(posterior199_compare(actual, actual), "multiply visited.*policy path")
})

test_that("intentionally partial generic objects retain their existing comparison contract", {
    x <- posterior199_object()
    for (method in names(x$results))
        x$results[[method]] <- x$results[[method]][c("all.probs", "data.used")]
    expect_true(posterior199_compare(x, x))
    actual <- x
    actual$results$AIMS$all.probs[[1L]][1L, 1L] <- .5 + 2^-53
    expect_true(posterior199_compare(actual, x))
    actual$results$AIMS$all.probs[[1L]][1L, 1L] <- .5 + 2^-40
    expect_error(posterior199_compare(actual, x), "numeric bound exceeded")
})

test_that("posterior helper comparisons leave seed and RNG kind unchanged", {
    withr::local_preserve_seed()
    set.seed(199)
    for (n in c(0L, 1L, 3L)) {
        x <- posterior199_object(n)
        for (method in names(x$results)) {
            if (n) {
                x$results[[method]]$all.probs[[1L]][1L, ] <- c(.5, .5, 0, 0, 0)
                x <- posterior199_relink(x, method)
            }
        }
        before <- .Random.seed; kind <- RNGkind()
        posterior199_accept(x)
        expect_identical(.Random.seed, before)
        expect_identical(RNGkind(), kind)
    }
})
