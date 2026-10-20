sspbc_ror_data <- function() {
    result <- new.env(parent = emptyenv())
    data("sspbc.models", package = "BreastSubtypeR", envir = result)
    data("sspbc.models.fullname", package = "BreastSubtypeR", envir = result)
    result
}

sspbc_ror_public_input <- function() {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = environment())
    OSLO2EMIT0obj$data_input$se_SSP
}

sspbc_ror_complete_input <- function() {
    model <- sspbc_ror_data()$sspbc.models[["ssp.ror"]]
    genes <- unique(unlist(strsplit(model$all.pairs, "<", fixed = TRUE)))
    # Deterministic positive expression for every model gene; no clinical
    # interpretation is attached to these artificial samples.
    expr <- outer(seq_along(genes), seq_len(3), function(i, j) {
        ((i * (2 * j + 1)) %% 97) + 1
    })
    dimnames(expr) <- list(
        sub("^e", "", genes),
        c("patient 03/A-B", "patient.01", "patient 02")
    )
    SummarizedExperiment::SummarizedExperiment(assays = list(expression = expr))
}

sspbc_native_result <- function(input, model) {
    suppressMessages(BreastSubtypeR:::applySSP(
        gex = SummarizedExperiment::assay(input),
        id = rownames(input), id.type = "EntrezGene", ssp.name = model
    ))
}

test_that("both packaged model collections contain the fitted SSP-ROR model", {
    models <- sspbc_ror_data()
    expect_true("ssp.ror" %in% names(models$sspbc.models))
    fullname <- paste0(
        "Training_Run19081Genes_noNorm_SSP.scaled.ROR.tot.asT0.c005.",
        "Fcc15_5x5foldCV.num.rules.50_21.selRules.AIMS.GS.RData"
    )
    expect_true(fullname %in% names(models$sspbc.models.fullname))
    short <- models$sspbc.models[["ssp.ror"]]
    full <- models$sspbc.models.fullname[[fullname]]
    expect_identical(short$k, 21)
    expect_length(short$all.pairs, 361)
    expect_identical(short$all.pairs, full$all.pairs)
    expect_identical(short$k, full$k)
    expect_identical(short$selected.pairs.list, full$selected.pairs.list)
    # The short-name model additionally carries isnumeric metadata.
    for (field in c("apriori", "tables", "levels", "call")) {
        expect_identical(
            short$one.vs.all.tsp[[21]][[field]],
            full$one.vs.all.tsp[[21]][[field]], info = field
        )
    }
    expect_identical(
        names(short$one.vs.all.tsp[[21]]$apriori),
        sprintf("c%03d", seq(5, 90, by = 5))
    )
})

test_that("the public SSP-ROR example exactly preserves native results", {
    input <- sspbc_ror_public_input()
    result <- suppressMessages(BS_sspbc(input, ssp.name = "ssp.ror"))
    expect_identical(result, sspbc_native_result(input, "ssp.ror"))
    expect_identical(rownames(result$cl), colnames(input))
    expect_identical(rownames(result$prob), colnames(input))
    expect_identical(colnames(result$data.used), colnames(input))
    # This reduced public example does not cover all SSP-ROR rules. Preserve
    # the native missingness instead of imputing or changing the fitted model.
    expect_true(anyNA(result$rules.matrix))
    expect_true(all(colSums(!is.na(result$rules.matrix)) == 75L))
})

test_that("SSP-ROR retains native types and shapes for one or more samples", {
    input <- sspbc_ror_complete_input()
    model <- sspbc_ror_data()$sspbc.models[["ssp.ror"]]
    labels <- sprintf("c%03d", seq(5, 90, by = 5))
    for (indices in list(2L, c(3L, 1L, 2L))) {
        subset <- input[, indices, drop = FALSE]
        result <- suppressMessages(BS_sspbc(subset, ssp.name = "ssp.ror"))
        n <- ncol(subset)
        expect_identical(result, sspbc_native_result(subset, "ssp.ror"))
        expect_type(result, "list")
        expect_identical(names(result), c(
            "cl", "prob", "all.probs", "rules.matrix",
            "data.used", "EntrezID.used"
        ))
        expect_true(is.matrix(result$cl))
        expect_type(result$cl, "character")
        expect_identical(dim(result$cl), c(n, 1L))
        expect_identical(dimnames(result$cl), list(colnames(subset), "21"))
        expect_true(all(result$cl %in% labels))
        expect_true(is.matrix(result$prob))
        expect_type(result$prob, "double")
        expect_identical(dimnames(result$prob), dimnames(result$cl))
        expect_identical(dim(result$prob), c(n, 1L))
        expect_true(all(is.finite(result$prob)))
        expect_true(all(result$prob >= 0 & result$prob <= 1))

        expect_type(result$all.probs, "list")
        expect_identical(names(result$all.probs), "21")
        probabilities <- result$all.probs[["21"]]
        expect_true(is.matrix(probabilities))
        expect_type(probabilities, "double")
        expect_identical(dim(probabilities), c(n, 18L))
        expect_identical(dimnames(probabilities), list(NULL, labels))
        expect_true(all(is.finite(probabilities)))
        expect_equal(unname(rowSums(probabilities)), rep(1, n))
        expect_equal(unname(result$prob[, 1]), apply(probabilities, 1, max))
        expect_identical(
            unname(result$cl[, 1]),
            labels[max.col(probabilities, ties.method = "first")]
        )

        expect_true(is.matrix(result$rules.matrix))
        expect_type(result$rules.matrix, "double")
        expect_identical(dim(result$rules.matrix), c(361L, n))
        expect_identical(
            dimnames(result$rules.matrix), list(model$all.pairs, NULL)
        )
        expect_true(all(result$rules.matrix %in% c(0, 1)))
        expect_true(is.matrix(result$data.used))
        expect_type(result$data.used, "double")
        expect_identical(dim(result$data.used), c(296L, n))
        expect_identical(colnames(result$data.used), colnames(subset))
        expect_type(result$EntrezID.used, "character")
        expect_identical(result$EntrezID.used, paste0("e", rownames(subset)))
        expect_identical(
            unname(result$data.used),
            unname(SummarizedExperiment::assay(subset))
        )
    }
})

test_that("SSP-ROR sample order agrees between individual and batch calls", {
    input <- sspbc_ror_complete_input()[, c(3L, 1L, 2L), drop = FALSE]
    batch <- suppressMessages(BS_sspbc(input, ssp.name = "ssp.ror"))
    for (i in seq_len(ncol(input))) {
        single <- suppressMessages(BS_sspbc(
            input[, i, drop = FALSE], ssp.name = "ssp.ror"
        ))
        expect_identical(single$cl, batch$cl[i, , drop = FALSE])
        expect_identical(single$prob, batch$prob[i, , drop = FALSE])
        expect_identical(
            single$all.probs[["21"]],
            batch$all.probs[["21"]][i, , drop = FALSE]
        )
        expect_identical(
            single$rules.matrix, batch$rules.matrix[, i, drop = FALSE]
        )
        expect_identical(single$data.used, batch$data.used[, i, drop = FALSE])
        expect_identical(single$EntrezID.used, batch$EntrezID.used)
    }
})

test_that("SSPBC rejects unknown and partial model names", {
    input <- sspbc_ror_complete_input()[, 1, drop = FALSE]
    for (name in c("unknown", "ssp.ROR", "ssp.r", "")) {
        expect_error(
            BS_sspbc(input, ssp.name = name),
            "Specified ssp.name not among names in sspbc.models",
            fixed = TRUE
        )
    }
})

test_that("existing SSPBC subtype models retain baseline output", {
    input <- sspbc_ror_public_input()[, c(7L, 2L, 5L), drop = FALSE]
    # Labels verified at unmodified main:
    # 67d982cb324ec376ea599929eb41837a50fd09b3.
    expected <- list(
        ssp.pam50 = c("LumA", "Basal", "Normal"),
        ssp.subtype = c("LumA", "Basal", "LumA")
    )
    for (name in names(expected)) {
        result <- suppressMessages(BS_sspbc(input, ssp.name = name))
        expect_identical(result, sspbc_native_result(input, name))
        expect_identical(unname(result$cl[, 1]), expected[[name]])
        expect_identical(rownames(result$cl), colnames(input))
        single <- suppressMessages(BS_sspbc(
            input[, 2L, drop = FALSE], ssp.name = name
        ))
        expect_identical(single, sspbc_native_result(
            input[, 2L, drop = FALSE], name
        ))
        expect_identical(single$cl, result$cl[2L, , drop = FALSE])
    }
    expect_identical(
        suppressMessages(BS_sspbc(input)),
        suppressMessages(BS_sspbc(input, ssp.name = "ssp.pam50"))
    )
})
