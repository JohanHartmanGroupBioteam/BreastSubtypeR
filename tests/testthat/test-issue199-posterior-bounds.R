b32_compare <- function(actual, reference, dataset = "OSLO2EMIT0obj") {
    issue199_portable_compare(actual, reference, dataset, path = "/res")
}

b32_relink <- function(parent) {
    p <- parent$all.probs[[1L]]
    for (i in seq_len(nrow(p))) {
        j <- which.max(p[i, ])
        parent$prob[i, 1L] <- if (length(j)) p[i, j] else NA_real_
        parent$cl[i, 1L] <- if (length(j)) colnames(p)[j] else NA_character_
    }
    parent
}

test_that("scientific comparison requires the complete pinned fixture at either root", {
    dataset <- "OSLO2EMIT0obj"
    fixture <- issue199_ror_approved(dataset)
    ref <- fixture$res
    expect_true(issue199_portable_compare(fixture, fixture, dataset))
    expect_true(b32_compare(ref, ref))
    x <- fixture
    x$se_obj <- NULL
    expect_false(identical(x, fixture))
    expect_error(issue199_portable_compare(x, x, dataset), "reference provenance")

    x <- ref
    x$results <- x$results[c("AIMS", "sspbc")]
    x$res_subtypes <- x$res_subtypes[1:3, c("AIMS", "sspbc"), drop = FALSE]
    rownames(x$res_subtypes) <- paste0("changed-", 1:3)
    for (method in names(x$results)) {
        o <- x$results[[method]]
        o$all.probs[[1L]] <- o$all.probs[[1L]][1:3, , drop = FALSE]
        for (field in c("prob", "cl")) {
            o[[field]] <- o[[field]][1:3, , drop = FALSE]
            rownames(o[[field]]) <- rownames(x$res_subtypes)
        }
        o$rules.matrix <- o$rules.matrix[, 1:3, drop = FALSE]
        o$data.used <- o$data.used[, 1:3, drop = FALSE]
        colnames(o$data.used) <- rownames(x$res_subtypes)
        o$BS.all <- o$BS.all[1:3, , drop = FALSE]
        o$BS.all$PatientID <- rownames(x$res_subtypes)
        x$results[[method]] <- o
    }
    expect_true(issue199_portability_unit(x, x, dataset, path = "/res"))
    expect_error(b32_compare(x, x), "reference provenance")
    x <- ref
    x$results$AIMS <- NULL
    expect_false("AIMS" %in% names(x$results))
    expect_error(b32_compare(x, ref), "Missing scientific posterior")
    expect_error(b32_compare(x, x), "reference provenance")
    x <- ref
    names(x$results$sspbc$all.probs) <- "23"
    expect_identical(names(x$results$sspbc$all.probs), "23")
    expect_error(b32_compare(x, ref), "posterior structure/axes")
    expect_error(b32_compare(x, x), "reference provenance")
    expect_error(issue199_portable_compare(1, 1, dataset,
        path = "/res/results/AIMS/all.probs/20"), "full parent context")
    x <- ref
    p <- x$results$AIMS$all.probs[[1L]]
    j <- which.min(p[1L, ])
    x$results$AIMS$all.probs[[1L]][1L, j] <- p[1L, j] * (1 + 2^-50)
    expect_false(identical(x, ref, num.eq = FALSE))
    expect_true(b32_compare(x, ref))
})

test_that("changed installed fixture and model bytes fail their provenance checks", {
    file <- tempfile()
    writeBin(as.raw(1:3), file)
    on.exit(unlink(file))
    check <- function(reader, basename, ...) {
        # Bind only the copied reader to a private file resolver.
        env <- new.env(parent = environment(reader))
        env$system.file <- function(..., package) {
            args <- list(...)
            if (identical(args[[2L]], basename)) return(file)
            base::system.file(..., package = package)
        }
        environment(reader) <- env
        reader(...)
    }
    expect_error(check(issue199_ror_approved, "OSLO2EMIT0obj.RData", "OSLO2EMIT0obj"),
        "file provenance")
    for (file_name in c("AIMSmodel.RData", "sspbc.models.RData"))
        expect_error(check(issue199_posterior_models, file_name), "model provenance")
})


test_that("both approved datasets route only normal posterior values to the model ceilings", {
    withr::local_seed(321)
    seed <- .Random.seed
    kind <- RNGkind()
    for (dataset in c("OSLO2EMIT0obj", "TCGABRCAobj")) {
        ref <- issue199_ror_approved(dataset)$res
        for (method in c("AIMS", "sspbc")) {
            q <- ref$results[[method]]$all.probs[[1L]]
            ceiling <- if (method == "AIMS") 2^-35 else 2^-33
            # Exercise both an unselected coordinate and the exact B2 link.
            for (j in c(which.min(q[1L, ]), which.max(q[1L, ]))) {
                expect_gt(q[1L, j], .Machine$double.xmin)
                for (fraction in c(0.5, 1.01)) {
                    x <- ref
                    p <- q
                    p[1L, j] <- q[1L, j]*(1 - fraction*ceiling)
                    gap <- abs(p[1L, j] - q[1L, j])/q[1L, j]
                    expect_gt(gap, 0)
                    expect_identical(unname(gap > ceiling), fraction > 1)
                    x$results[[method]]$all.probs[[1L]] <- p
                    x$results[[method]] <- b32_relink(x$results[[method]])
                    expect_true(issue199_posterior_link(x$results[[method]], method))
                    if (fraction < 1) expect_true(b32_compare(x, ref, dataset))
                    else expect_error(b32_compare(x, ref, dataset), "posterior bound exceeded")
                }
            }
        }
    }
    expect_identical(.Random.seed, seed)
    expect_identical(RNGkind(), kind)
})

test_that("rules, axes, selected probabilities and rank changes reach distinct safeguards", {
    ref <- issue199_ror_approved("OSLO2EMIT0obj")$res
    original <- ref$results$AIMS
    compare <- function(o) {
        x <- ref
        x$results$AIMS <- o
        b32_compare(x, ref)
    }
    a <- original
    a$rules.matrix[1L, 1L] <- 1 - a$rules.matrix[1L, 1L]
    expect_false(identical(a$rules.matrix, original$rules.matrix))
    expect_error(compare(a), "posterior rules")
    a <- original
    a$prob[1L, 1L] <- a$prob[1L, 1L] + 2^-53
    expect_false(identical(a$prob, original$prob, num.eq = FALSE))
    expect_error(compare(a), "selected probability/class link")
    a <- original
    colnames(a$all.probs[[1L]])[1L] <- "renamed"
    expect_error(compare(a), "posterior structure/axes")
    a <- original
    row_order <- c(2L, 1L, 3:nrow(a$prob))
    a$cl <- a$cl[row_order, , drop = FALSE]
    a$prob <- a$prob[row_order, , drop = FALSE]
    a$all.probs[[1L]] <- a$all.probs[[1L]][row_order, , drop = FALSE]
    expect_false(identical(rownames(a$prob), rownames(original$prob)))
    expect_error(compare(a), "posterior structure/axes")
    a <- original
    a$all.probs$extra <- a$all.probs[[1L]]
    expect_error(compare(a), "posterior structure/axes")
    a <- original
    a$extra <- TRUE
    expect_error(compare(a), "type/class/length|attributes/names/dimensions")
    q <- original$all.probs[[1L]]
    order <- order(q[1L, ], decreasing = TRUE)
    for (change in c("swap", "first tie", "nonmaximum tie")) {
        p <- q
        if (change == "swap") p[1L, order[1:2]] <- q[1L, order[2:1]]
        if (change == "first tie") p[1L, order[2L]] <- q[1L, order[1L]]
        if (change == "nonmaximum tie") p[1L, order[4L]] <- q[1L, order[3L]]
        expect_false(identical(p, q, num.eq = FALSE))
        a <- original
        a$all.probs[[1L]] <- p
        a <- b32_relink(a)
        expect_true(issue199_posterior_link(a, change))
        expect_error(compare(a), "posterior ranks/ties/argmax")
    }
})

test_that("scientific posterior masks, signed zeros and subnormal values remain exact", {
    dataset <- "TCGABRCAobj"
    ref <- issue199_ror_approved(dataset)$res
    original <- ref$results$sspbc
    q <- original$all.probs[[1L]]
    tiny <- which(q > 0 & q < .Machine$double.xmin)
    zero <- which(q == 0)[1L]
    expect_length(tiny, 1L)
    negative_zero <- readBin(as.raw(c(0,0,0,0,0,0,0,128)), "double", 1L, endian = "little")
    expect_identical(1/negative_zero, -Inf)
    for (change in c("zero to tiny", "tiny to zero", "signed", "tiny ulp", "tiny to normal",
                     "partial NA", "whole NA", "NaN", "Inf")) {
        p <- q
        if (change == "zero to tiny") p[zero] <- 2^-1074
        if (change == "tiny to zero") p[tiny] <- 0
        if (change == "signed") p[zero] <- negative_zero
        if (change == "tiny ulp") p[tiny] <- q[tiny] + 2^-1074
        if (change == "tiny to normal") p[tiny] <- .Machine$double.xmin
        if (change == "partial NA") p[1L, 1L] <- NA_real_
        if (change == "whole NA") p[1L, ] <- NA_real_
        if (change == "NaN") p[1L, 1L] <- NaN
        if (change == "Inf") p[1L, 1L] <- Inf
        expect_false(identical(q, p, num.eq = FALSE, single.NA = FALSE))
        x <- ref
        x$results$sspbc$all.probs[[1L]] <- p
        reason <- switch(change, "zero to tiny" =, "tiny to zero" = "zero support",
            signed = "signed zero", "tiny ulp" =, "tiny to normal" = "subnormal support/values",
            "whole NA" = "nonfinite masks/signs", "posterior rows")
        if (change == "whole NA") x$results$sspbc <- b32_relink(x$results$sspbc)
        if (reason != "posterior rows") expect_true(issue199_posterior_link(x$results$sspbc, change))
        expect_error(b32_compare(x, ref, dataset), reason)
    }
})

test_that("each numerical guard refuses an isolated violation", {
    expect_error(issue199_posterior_x(38), "G-X")
    terms <- list(L = matrix(c(-4, -3, -2, -1, 0), 1L), magnitude = matrix(4, 1L, 5L))
    q <- matrix(c(.02, .04, .1, .24, .6), 1L)
    expect_type(issue199_posterior_domain(terms, q, "AIMS"), "list")
    for (method in c("AIMS", "sspbc")) {
        bad <- terms
        bad$magnitude[1L] <- if (method == "AIMS") 521 else 1521
        expect_error(issue199_posterior_domain(bad, q, method), "G-S1")
        bad <- terms
        bad$L[1L] <- if (method == "AIMS") -501 else -1521
        expect_error(issue199_posterior_domain(bad, q, method), "G-S2")
    }
    bad <- terms
    bad$L[1L, ] <- c(-log(.Machine$double.xmax), -100, -90, -80, 0)
    expect_error(issue199_posterior_domain(bad, q, "sspbc"), "G-O")
    near <- q
    near[1L] <- q[2L]*(1 - 2^-21)
    expect_error(issue199_posterior_domain(terms, near, "AIMS"), "G-R")
})

b32_expect_oracle <- function(observed, calculated) {
    if (!identical(observed, calculated, num.eq = FALSE, single.NA = FALSE))
        stop("same-runtime posterior oracle mismatch", call. = FALSE)
    invisible(TRUE)
}

test_that("fresh public classifiers match the independent reciprocal oracle exactly", {
    models <- issue199_posterior_models()
    withr::local_seed(322)
    for (dataset in c("OSLO2EMIT0obj", "TCGABRCAobj")) {
        fixture <- issue199_ror_approved(dataset)
        for (method in names(models)) {
            model <- models[[method]]
            nb <- model$one.vs.all.tsp[[model$k]]
            se <- fixture$data_input$se_SSP
            fresh <- if (method == "AIMS") BS_AIMS(se) else BS_sspbc(se, ssp.name = "ssp.pam50")
            rules <- fresh$rules.matrix
            expect_true(identical(rules, fixture$res$results[[method]]$rules.matrix, num.eq = FALSE))
            seed <- .Random.seed
            kind <- RNGkind()
            want <- issue199_posterior_oracle(nb, rules)
            expect_true(b32_expect_oracle(fresh$all.probs[[as.character(model$k)]], want))
            expect_identical(.Random.seed, seed)
            expect_identical(RNGkind(), kind)
            expect_identical(issue199_posterior_oracle(nb, rules[rev(seq_len(nrow(rules))), , drop = FALSE]), want)
            L <- issue199_posterior_terms(nb, rules)$L
            alternative <- t(apply(L, 1L, function(v) {
                e <- exp(v - max(v))
                e/sum(e)
            }))
            dimnames(alternative) <- dimnames(want)
            expect_false(identical(alternative, want, num.eq = FALSE))
            expect_error(b32_expect_oracle(alternative, want), "same-runtime posterior oracle mismatch")
        }
    }
})

test_that("model perturbations separate provenance, arithmetic and cross-runtime sensitivity", {
    models <- issue199_posterior_models()
    dataset <- "OSLO2EMIT0obj"
    ref <- issue199_ror_approved(dataset)$res
    for (method in names(models)) {
        model <- models[[method]]
        nb <- model$one.vs.all.tsp[[model$k]]
        nb$isnumeric <- setNames(rep(TRUE, length(nb$tables)), names(nb$tables))
        o <- ref$results[[method]]
        rules <- o$rules.matrix
        q <- o$all.probs[[1L]]
        nd <- data.frame(t(rules))
        pristine <- issue199_posterior_oracle(nb, rules)
        ceiling <- if (method == "AIMS") 2^-35 else 2^-33
        for (change in c("prior", "mean", "SD", "floor")) {
            sizes <- if (change == "floor") 1e-7 else c(1e-12, 1e-9)
            for (size in sizes) {
                copy <- nb
                if (change == "prior") copy$apriori[1L] <- if (size == 1e-12)
                    nb$apriori[1L]*(1+size) else nb$apriori[1L]+1L
                if (change %in% c("mean", "SD")) {
                    j <- if (change == "mean") 1L else 2L
                    copy$tables[[1L]][1L, j] <- nb$tables[[1L]][1L, j]*(1+size)
                }
                threshold <- if (change == "floor") .001*(1+size) else .001
                mutated <- predict(copy, nd, type = "raw", threshold = threshold)
                if (change != "floor") {
                    expect_false(identical(copy, nb))
                    expect_true(b32_expect_oracle(mutated, issue199_posterior_oracle(copy, rules)))
                }
                if (method == "AIMS" && change == "floor") {
                    expect_true(b32_expect_oracle(mutated, pristine))
                    next
                }
                expect_false(identical(mutated, pristine, num.eq = FALSE))
                expect_error(b32_expect_oracle(mutated, pristine), "same-runtime posterior oracle mismatch")
                # Transplant one normal, nonmaximal coordinate to isolate E10.
                eligible <- which(pristine > 1e-200 & pristine < .25 & mutated > 0,
                    arr.ind = TRUE)
                ratios <- mutated[eligible]/pristine[eligible]
                j <- which.max(abs(ratios - 1))
                cell <- eligible[j, , drop = FALSE]
                x <- o
                x$all.probs[[1L]][cell] <- q[cell]*ratios[j]
                x <- b32_relink(x)
                actual <- x$all.probs[[1L]][cell]
                gap <- abs(actual-q[cell])/max(actual, q[cell])
                expect_gt(gap, 0)
                expect_true(issue199_posterior_link(x, method))
                check <- function() issue199_posterior_method(x, o, model, method, dataset)
                if (size == 1e-12) {
                    expect_lt(gap, ceiling)
                    expect_type(check(), "list")
                } else {
                    expect_gt(gap, ceiling)
                    expect_error(check(), "posterior bound exceeded")
                }
            }
        }
    }
})

test_that("missing evidence and the density floor follow the original primitive order", {
    model <- issue199_posterior_models()$sspbc
    nb <- model$one.vs.all.tsp[[24L]]
    nb$tables <- list(x = cbind(rep(0, 5), c(0, -1, .5, .5, .5)),
        y = cbind(rep(.2, 5), rep(.5, 5)))
    nb$isnumeric <- c(x = TRUE, y = TRUE)
    rules <- rbind(x = c(0, 1, NA_real_), y = c(NA_real_, 0, NA_real_))
    nd <- data.frame(t(rules))
    raw <- predict(nb, nd, type = "raw")
    want <- issue199_posterior_oracle(nb, rules)
    expect_true(identical(raw[1:2, , drop = FALSE], want[1:2, , drop = FALSE], num.eq = FALSE))
    expect_true(all(is.na(want[3L, ])) && !any(is.nan(want[3L, ])))
    expect_true(all(is.finite(raw[3L, ])))
    mutant <- predict(nb, nd, type = "raw", threshold = .0010000001)
    expect_false(identical(mutant[2L, ], raw[2L, ], num.eq = FALSE))
    expect_error(b32_expect_oracle(mutant[2L, ], want[2L, ]), "same-runtime posterior oracle mismatch")
})
