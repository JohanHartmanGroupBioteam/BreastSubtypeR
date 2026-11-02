sspbc_ror_data <- function() {
    result <- new.env(parent = emptyenv())
    data("sspbc.models", package = "BreastSubtypeR", envir = result)
    data("sspbc.models.fullname", package = "BreastSubtypeR", envir = result)
    result
}

sspbc_ror_public_object <- function(name = "OSLO2EMIT0obj") {
    env <- new.env(parent = emptyenv())
    data(list = name, package = "BreastSubtypeR", envir = env)
    env[[name]]
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

test_that("subtype-only inputs preserve native incomplete-rule results", {
    # Default mapped inputs deliberately retain only the subtype signature.
    # Full SSP-ROR coverage requires the opt-in input tested separately below.
    input <- sspbc_ror_public_object()$data_input$se_SSP
    result <- suppressMessages(BS_sspbc(input, ssp.name = "ssp.ror"))
    expect_identical(result, sspbc_native_result(input, "ssp.ror"))
    expect_identical(rownames(result$cl), colnames(input))
    expect_identical(rownames(result$prob), colnames(input))
    expect_identical(colnames(result$data.used), colnames(input))
    # Derive the available rules from the supplied genes, without a historical
    # rule-count oracle tied to an older packaged dataset.
    pairs <- strsplit(sspbc_ror_data()$sspbc.models$ssp.ror$all.pairs, "<", fixed = TRUE)
    available <- vapply(pairs, function(pair) {
        all(pair %in% paste0("e", rownames(input)))
    }, logical(1))
    expect_true(any(available))
    expect_true(any(!available))
    expect_true(anyNA(result$rules.matrix))
    expect_identical(unname(!is.na(result$rules.matrix)),
        matrix(rep(available, ncol(input)), ncol = ncol(input)))
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

test_that("existing SSPBC subtype models retain current native output", {
    input <- sspbc_ror_public_object()$data_input$se_SSP[, c(7L, 2L, 5L), drop = FALSE]
    for (name in c("ssp.pam50", "ssp.subtype")) {
        result <- suppressMessages(BS_sspbc(input, ssp.name = name))
        expect_identical(result, sspbc_native_result(input, name))
        model <- sspbc_ror_data()$sspbc.models[[name]]
        labels <- names(model$one.vs.all.tsp[[model$k]]$apriori)
        expect_true(all(result$cl %in% labels))
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

test_that("expanded public datasets supply complete native SSP-ROR outputs", {
    model <- sspbc_ror_data()$sspbc.models$ssp.ror
    genes <- unique(unlist(strsplit(model$all.pairs, "<", fixed = TRUE)))
    ids <- sub("^e", "", genes)
    labels <- sprintf("c%03d", seq(5, 90, by = 5))
    components <- c("cl", "prob", "all.probs", "rules.matrix",
        "data.used", "EntrezID.used")
    for (dataset in c("OSLO2EMIT0obj", "TCGABRCAobj")) {
        object <- sspbc_ror_public_object(dataset)
        raw <- identical(dataset, "TCGABRCAobj")
        source <- object$se_obj
        expect_identical(dim(source), if (raw) c(5350L, 100L) else c(569L, 102L))
        mapped <- Mapping(source, RawCounts = raw, include_ssp_ror = TRUE,
            verbose = FALSE)
        input <- mapped$se_SSP_ROR
        expect_identical(dim(input), c(296L, ncol(source)))
        expect_identical(rownames(input), ids)
        expect_identical(colnames(input), colnames(source))
        expect_identical(SummarizedExperiment::colData(input),
            SummarizedExperiment::colData(source))
        batch <- suppressMessages(BS_sspbc(input, ssp.name = "ssp.ror"))
        # Full cohort, full reordered cohort and singleton: all six native
        # components must agree, including the positional rule/posterior axes.
        for (indices in list(seq_len(ncol(input)), rev(seq_len(ncol(input))), 7L)) {
            subset <- input[, indices, drop = FALSE]
            result <- suppressMessages(BS_sspbc(subset, ssp.name = "ssp.ror"))
            native <- sspbc_native_result(subset, "ssp.ror")
            expect_identical(names(result), components)
            for (component in components) {
                expect_identical(result[[component]], native[[component]],
                    info = paste(dataset, component))
            }
            expect_identical(result$cl, batch$cl[indices, , drop = FALSE])
            expect_identical(result$prob, batch$prob[indices, , drop = FALSE])
            expect_identical(result$all.probs[["21"]],
                batch$all.probs[["21"]][indices, , drop = FALSE])
            expect_identical(result$rules.matrix,
                batch$rules.matrix[, indices, drop = FALSE])
            expect_identical(result$data.used, batch$data.used[, indices, drop = FALSE])
            expect_identical(result$EntrezID.used, batch$EntrezID.used)
            expect_identical(dimnames(result$cl), list(colnames(subset), "21"))
            expect_identical(dimnames(result$prob), dimnames(result$cl))
            expect_identical(colnames(result$data.used), colnames(subset))
            expect_identical(result$EntrezID.used, genes)
            expect_identical(dim(result$rules.matrix), c(361L, length(indices)))
            expect_true(all(colSums(!is.na(result$rules.matrix)) == 361L))
            expect_true(all(result$rules.matrix %in% c(0, 1)))
            expect_true(all(result$cl %in% labels))
            probabilities <- result$all.probs[["21"]]
            expect_identical(dimnames(probabilities), list(NULL, labels))
            expect_identical(dim(probabilities), c(length(indices), 18L))
            expect_true(all(is.finite(probabilities)))
            expect_true(all(probabilities >= 0 & probabilities <= 1))
            expect_equal(unname(rowSums(probabilities)), rep(1, length(indices)))
            expect_equal(unname(result$prob[, 1]), apply(probabilities, 1, max))
            expect_identical(unname(result$cl[, 1]),
                labels[max.col(probabilities, ties.method = "first")])
            # Executable counterpart of the documented table/CSV recipe.
            table <- data.frame(SampleID = rownames(result$cl),
                ROR_bin = result$cl[, "21"],
                Bin_probability = result$prob[, "21"], probabilities,
                row.names = NULL, check.names = FALSE)
            expect_identical(table$SampleID, colnames(subset))
            expect_identical(table$ROR_bin, unname(result$cl[, "21"]))
            expect_identical(table$Bin_probability, unname(result$prob[, "21"]))
            expect_identical(names(table), c("SampleID", "ROR_bin", "Bin_probability", labels))
            expect_identical(unname(as.matrix(table[labels])), unname(probabilities))
        }
    }
})
