rule_test_models <- function() {
    env <- new.env(parent = emptyenv())
    data(list = c("AIMSmodel", "sspbc.models"),
        package = "BreastSubtypeR", envir = env)
    list(AIMS = env$AIMSmodel, ssp.pam50 = env$sspbc.models$ssp.pam50,
        ssp.subtype = env$sspbc.models$ssp.subtype)
}

rule_test_pairs <- function(model) {
    strsplit(model$all.pairs, "<", fixed = TRUE)
}

rule_test_input <- function(model) {
    pairs <- rule_test_pairs(model)
    genes <- unique(unlist(pairs, use.names = FALSE))
    x <- matrix(NA_real_, length(genes), 6L, dimnames = list(genes,
        c("missing/z", "complete B", "all FALSE", "one FALSE", "partial", "one TRUE")))
    x[, 2L] <- 2 + sin(seq_along(genes))
    x[, 3L] <- 1
    # Select two genes that make exactly one of the actual rules evaluable.
    unique_pair <- which(vapply(pairs, function(pair) {
        sum(vapply(pairs, function(p) all(p %in% pair), logical(1))) == 1L
    }, logical(1)))[1L]
    pair <- pairs[[unique_pair]]
    x[pair, 4L] <- c(2, 1)
    x[seq_len(ceiling(length(genes) / 2)), 5L] <- 1
    x[pair, 5L] <- c(1, 2)
    x[pair, 6L] <- c(1, 2)
    x
}

rule_test_se <- function(x) {
    SummarizedExperiment::SummarizedExperiment(
        assays = list(expression = x),
        colData = S4Vectors::DataFrame(PatientID = colnames(x),
            row.names = colnames(x)))
}

rule_test_wrapper <- function(x, key) {
    rownames(x) <- sub("^e", "", rownames(x))
    se <- rule_test_se(x)
    suppressMessages(if (key == "AIMS") BS_AIMS(se) else BS_sspbc(se, key))
}

# Independent native posterior reference: no zero-evidence masking. This is
# also the original prediction for every supported sample, including FALSEs.
rule_test_native <- function(rules, model, key) {
    lapply(model$k, function(k) {
        fit <- model$one.vs.all.tsp[[k]]
        if (key == "AIMS") {
            fit$isnumeric <- as.list(setNames(rep(TRUE, length(fit$tables)),
                names(fit$tables)))
        }
        stats::predict(fit, data.frame(t(rules)), type = "raw")
    })
}

expect_rule_predictions <- function(result, model, key, samples) {
    rules <- result$rules.matrix
    native <- rule_test_native(rules, model, key)
    expect_identical(names(result), c("cl", "prob", "all.probs",
        "rules.matrix", "data.used", "EntrezID.used"))
    expect_identical(dim(result$cl), c(length(samples), length(model$k)))
    expect_identical(dimnames(result$prob), dimnames(result$cl))
    expect_identical(rownames(result$cl), samples)
    expect_identical(colnames(result$cl), as.character(model$k))
    expect_identical(names(result$all.probs), as.character(model$k))
    expect_identical(dim(rules), c(length(model$all.pairs), length(samples)))
    expect_identical(rownames(rules), model$all.pairs)
    expect_identical(colnames(result$data.used), samples)
    for (i in seq_along(model$k)) {
        k <- model$k[i]
        feature_rows <- match(names(model$one.vs.all.tsp[[k]]$tables),
            make.names(rownames(rules)))
        unsupported <- colSums(!is.na(rules[feature_rows, , drop = FALSE])) == 0
        posterior <- native[[i]]
        calls <- colnames(posterior)[max.col(posterior, ties.method = "first")]
        probability <- apply(posterior, 1, max)
        calls[unsupported] <- NA_character_
        probability[unsupported] <- NA_real_
        posterior[unsupported, ] <- NA_real_
        expect_identical(unname(result$cl[, i]), calls)
        expect_identical(unname(result$prob[, i]), unname(probability))
        expect_identical(result$all.probs[[i]], posterior)
    }
}

for (key in c("AIMS", "ssp.pam50", "ssp.subtype")) {
    test_that(paste(key, "masks only zero-evidence samples, including singly"), {
        model <- rule_test_models()[[key]]
        snapshot <- serialize(model, NULL)
        x <- rule_test_input(model)
        result <- rule_test_wrapper(x, key)
        direct <- suppressMessages(if (key == "AIMS") {
            BreastSubtypeR:::.apply.nbc_AIMS(x, rownames(x), model)
        } else {
            BreastSubtypeR:::.apply.nbc(x, rownames(x), model, zero.rules.na = TRUE)
        })
        expect_identical(result, direct)
        expected_rules <- t(vapply(rule_test_pairs(model), function(pair) {
            as.numeric(x[pair[1L], ] < x[pair[2L], ])
        }, numeric(ncol(x))))
        rownames(expected_rules) <- model$all.pairs
        expect_identical(result$rules.matrix, expected_rules)
        count <- colSums(!is.na(result$rules.matrix))
        expect_equal(count[c(1, 2, 3, 4, 6)],
            c(0, length(model$all.pairs), length(model$all.pairs), 1, 1))
        expect_gt(count[5], 1)
        expect_lt(count[5], length(model$all.pairs))
        expect_equal(sum(result$rules.matrix[, 3:4], na.rm = TRUE), 0)
        expect_rule_predictions(result, model, key, colnames(x))
        for (j in seq_len(ncol(x))) {
            single <- rule_test_wrapper(x[, j, drop = FALSE], key)
            expect_rule_predictions(single, model, key, colnames(x)[j])
            expect_identical(single$cl, result$cl[j, , drop = FALSE])
            expect_identical(single$prob, result$prob[j, , drop = FALSE])
            expect_identical(single$rules.matrix,
                result$rules.matrix[, j, drop = FALSE])
            expect_identical(single$all.probs[[1]],
                result$all.probs[[1]][j, , drop = FALSE])
        }
        expect_identical(serialize(model, NULL), snapshot)
        expect_identical(serialize(rule_test_models()[[key]], NULL), snapshot)
    })

    test_that(paste(key, "handles absent genes without losing matrix dimensions"), {
        model <- rule_test_models()[[key]]
        x <- rule_test_input(model)[, c(2L, 3L), drop = FALSE]
        irrelevant <- x[1:2, , drop = FALSE]
        rownames(irrelevant) <- c("not-a-model-gene", "also-absent")
        inputs <- list(no_model_genes = irrelevant, empty = x[FALSE, , drop = FALSE],
            one_gene = x[1L, , drop = FALSE],
            absent_pair_partners = x[c(1L, 3L), , drop = FALSE])
        for (input in inputs) {
            for (cols in list(1L, 2:1)) {
                selected <- input[, cols, drop = FALSE]
                result <- rule_test_wrapper(selected, key)
                expect_true(all(is.na(result$rules.matrix)))
                expect_rule_predictions(result, model, key, colnames(selected))
            }
        }
        # Absent rows with one remaining evaluable rule must retain native calls.
        pair <- rownames(x)[!is.na(rule_test_input(model)[, 4L])]
        result <- rule_test_wrapper(x[pair, , drop = FALSE], key)
        expect_equal(colSums(!is.na(result$rules.matrix)), c(1, 1))
        expect_rule_predictions(result, model, key, colnames(x))
    })
}

test_that("other SSPBC models retain their native zero-evidence behavior", {
    env <- new.env(parent = emptyenv())
    data("sspbc.models", package = "BreastSubtypeR", envir = env)
    for (key in setdiff(names(env$sspbc.models), c("ssp.pam50", "ssp.subtype"))) {
        model <- env$sspbc.models[[key]]
        x <- rule_test_input(model)[, c(1L, 2L, 3L), drop = FALSE]
        result <- rule_test_wrapper(x, key)
        native <- rule_test_native(result$rules.matrix, model, key)
        expect_identical(unname(result$all.probs), native)
        for (i in seq_along(model$k)) {
            expect_identical(unname(result$cl[, i]),
                colnames(native[[i]])[max.col(native[[i]], ties.method = "first")])
            expect_identical(unname(result$prob[, i]),
                unname(apply(native[[i]], 1, max)))
        }
    }
})

test_that("BS_Multi propagates per-model missing calls and native details", {
    models <- rule_test_models()
    pairs <- lapply(models, function(model) {
        lapply(rule_test_pairs(model), function(p) sub("^e", "", p))
    })
    genes <- unique(unlist(pairs, use.names = FALSE))
    # A real PAM50 rule whose two genes support neither AIMS nor ssp.subtype.
    unique_pair <- which(vapply(pairs$ssp.pam50, function(pair) {
        !any(vapply(c(pairs$AIMS, pairs$ssp.subtype),
            function(p) all(p %in% pair), logical(1)))
    }, logical(1)))[1L]
    expect_false(is.na(unique_pair))
    x <- matrix(1, length(genes), 4L, dimnames = list(genes,
        c("unsupported z", "FALSE a", "pam50 only", "complete b")))
    x[, 1L] <- NA_real_
    x[, 3L] <- NA_real_
    x[pairs$ssp.pam50[[unique_pair]], 3L] <- c(1, 2)
    x[, 4L] <- 2 + sin(seq_along(genes))
    for (cols in list(1:4, 1L, 3L, 4L)) {
        se <- rule_test_se(x[, cols, drop = FALSE])
        input <- list(se_NC = se[, rev(seq_len(ncol(se))), drop = FALSE], se_SSP = se)
        for (subtype in c(FALSE, TRUE)) {
            result <- suppressMessages(BS_Multi(input,
                methods = c("sspbc", "AIMS"), Subtype = subtype))
            samples <- colnames(input$se_NC)
            expect_identical(rownames(result$res_subtypes), samples)
            for (key in names(models)) {
                direct <- rule_test_wrapper(x[, cols, drop = FALSE], key)
                method <- if (key == "AIMS") "AIMS" else "sspbc"
                if (key != "ssp.subtype") {
                    expect_identical(result$results[[method]][names(direct)], direct)
                    expect_identical(result$res_subtypes[[method]],
                        unname(direct$cl[samples, 1L]))
                }
                if (subtype && key != "ssp.pam50") {
                    calls <- unname(direct$cl[samples, 1L])
                    if (key == "AIMS") calls[which(calls == "Normal")] <- NA_character_
                    expect_identical(result$res_subtypes.Subtype[[method]], calls)
                }
            }
            if (1L %in% cols) {
                expect_true(all(is.na(result$res_subtypes["unsupported z", ])))
                if (subtype) expect_true(all(is.na(
                    result$res_subtypes.Subtype["unsupported z", ])))
            }
            if (3L %in% cols) {
                expect_true(is.na(result$res_subtypes["pam50 only", "AIMS"]))
                expect_false(is.na(result$res_subtypes["pam50 only", "sspbc"]))
                expect_equal(result$res_subtypes["pam50 only", "entropy"], 0)
                if (subtype) expect_true(all(is.na(
                    result$res_subtypes.Subtype["pam50 only", ])))
            }
        }
    }
})
