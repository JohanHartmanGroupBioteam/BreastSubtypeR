# These synthetic cohorts reproduce the NCP-03 investigation at ffab0e4.
ncp03_auto_input <- function(pos, neg, measured) {
    model <- new.env(parent = emptyenv())
    utils::data("BreastSubtypeRobj", envir = model, package = "BreastSubtypeR")
    registry <- model$BreastSubtypeRobj$genes.signature
    centroids <- as.matrix(model$BreastSubtypeRobj$centroid)
    patients <- seq_len(pos + neg)
    ids <- sprintf("N%02d", rev(patients))
    clinical <- S4Vectors::DataFrame(PatientID = factor(ids, levels = c("unused", rev(ids))),
        ER = rep(c("ER+", "ER-"), c(pos, neg)), HER2 = rep("HER2-", length(ids)),
        TSIZE = patients %% 2, NODE = patients %% 3, row.names = ids)
    expr <- outer(seq_len(nrow(registry)), patients,
        function(g, s) 7 + sin(g * s * .19) + cos(g + s * .23))
    pam <- match(rownames(centroids), registry$Symbol)
    subtype <- 1L + (patients - 1L) %% ncol(centroids)
    expr[pam, ] <- expr[pam, ] + 2 * centroids[, subtype]
    expr[pam[seq_along(pam) > measured], 3L] <- NA_real_
    dimnames(expr) <- list(as.character(registry$EntrezGene.ID), ids)
    nc <- expr[pam, , drop = FALSE]
    rownames(nc) <- rownames(centroids)
    se <- function(x) SummarizedExperiment::SummarizedExperiment(
        assays = list(expression = x), colData = clinical)
    list(input = list(se_NC = se(nc), se_SSP = se(2 ^ expr[registry$SSP_based == "Yes", ])),
        centroids = centroids)
}

# All 33 patients calibrate each iteration in the 15+/18- cohort.
# In the 15+/25- cohort the default ratio draws 18 ER- and keeps every ER+.
# Independent Spearman votes use only the complete peers for prediction, and
# retain the original sparse patient in every calibration draw that selects it.
ncp03_auto_votes <- function(se, centroids, four) {
    x <- SummarizedExperiment::assay(se)
    pos <- which(se$ER == "ER+")
    neg <- which(se$ER == "ER-")
    draws <- withr::with_seed(118L, replicate(100L,
        if (length(neg) == 18L) c(pos, neg) else c(sample(neg, 18L), pos),
        simplify = FALSE))
    classes <- colnames(centroids)
    if (four) classes <- classes[1:4]
    ids <- colnames(x)
    votes <- vapply(draws, function(keep) {
        medians <- apply(x[, keep, drop = FALSE], 1L, stats::median, na.rm = TRUE)
        centered <- sweep(x[, -3L, drop = FALSE], 1L, medians)
        rho <- stats::cor(centroids[, classes, drop = FALSE], centered, method = "spearman")
        result <- stats::setNames(rep(NA_character_, length(ids)), ids)
        result[-3L] <- classes[apply(rho, 2L, which.max)]
        result
    }, stats::setNames(character(length(ids)), ids))
    colnames(votes) <- paste0("itr.", seq_len(ncol(votes)))
    votes
}

ncp03_multi <- function(input, four, clinical, methods = "AUTO") {
    # PCAPAM50 cannot fit these sparse mixed cohorts. Its warning is independent
    # of the experimental research ROR notice aggregated by BS_Multi.
    expect_warning(result <- suppressMessages(expect_experimental_ror_pc(
        BS_Multi(input, methods = methods, Subtype = four, hasClinical = clinical),
        n = as.integer(clinical))), "^PCAPAM50 failed in this iteration:")
    result
}

ncp03_expect_summary <- function(summary, votes, ids, other) {
    expect_identical(rownames(summary), ids)
    expect_identical(summary[names(other)], other)
    expect_identical(summary$cIHC.itr, unname(apply(votes, 1L, function(row) {
        counts <- table(row)
        if (length(counts)) names(counts)[which.max(counts)] else NA_character_
    })))
    calls <- summary[setdiff(names(summary), "entropy")]
    entropy <- apply(calls, 1L, function(row) {
        frequencies <- prop.table(table(row))
        -sum(frequencies * log2(frequencies))
    })
    expect_equal(summary$entropy, unname(entropy), tolerance = 1e-14)
    expect_identical(which(is.na(summary$cIHC.itr)), 3L)
}

for (design in list(c(18L, 0L), c(18L, 1L), c(25L, 1L))) {
    for (four in c(FALSE, TRUE)) {
        test_that(sprintf("AUTO preserves peers with %d ER-, %d measured genes and Subtype=%s",
            design[1L], design[2L], four), {
            withr::local_options(list(BreastSubtypeR.verbose = FALSE))
            data <- ncp03_auto_input(15L, design[1L], design[2L])
            input <- data$input
            ids <- colnames(input$se_NC)
            clinical <- design[1L] == 18L
            methods <- c(if (clinical) c("parker.original", "genefu.scale"), "genefu.robust",
                "ssBC", "ssBC.v2", "cIHC", "cIHC.itr", "PCAPAM50", "AIMS", "sspbc")
            other_methods <- setdiff(methods, "cIHC.itr")
            actual <- ncp03_multi(input, four, clinical)
            expect_identical(names(actual$results), methods)
            # All patients satisfy both ssBC clinical-subgroup gates here.
            native <- ncp03_multi(input, four, clinical, other_methods)
            for (method in other_methods) {
                if (method %in% c("ssBC", "ssBC.v2")) {
                    expect_auto_ssbc_axes(actual$results[[method]], native$results[[method]], ids)
                } else {
                    expect_identical(actual$results[[method]], native$results[[method]])
                }
            }

            itr <- actual$results$cIHC.itr
            expect_identical(itr$BS.all$PatientID, ids)
            five_votes <- ncp03_auto_votes(input$se_NC, data$centroids, FALSE)
            expect_identical(itr$BS.itr.keep, five_votes)
            ncp03_expect_summary(actual$res_subtypes, five_votes, ids, native$res_subtypes[other_methods])
            expect_identical(itr$BS.all$BS, actual$res_subtypes$cIHC.itr)
            if (four) {
                four_votes <- ncp03_auto_votes(input$se_NC, data$centroids, TRUE)
                expect_identical(itr$BS.itr.keep.Subtype, four_votes)
                ncp03_expect_summary(actual$res_subtypes.Subtype, four_votes, ids,
                    native$res_subtypes.Subtype[other_methods])
                expect_identical(itr$BS.all$BS.Subtype, actual$res_subtypes.Subtype$cIHC.itr)
            } else {
                expect_null(actual$res_subtypes.Subtype)
            }
            expect_identical(as.character(actual$res_subtypes$AIMS[3L]), "Her2")
            expect_identical(as.character(actual$res_subtypes$sspbc[3L]), "Her2")
        })
    }
}

test_that("SSP-only AUTO retains both SSP calls for a patient with no PAM50 data", {
    input <- ncp03_auto_input(5L, 3L, 0L)$input
    for (four in c(FALSE, TRUE)) {
        auto <- suppressMessages(BS_Multi(input, Subtype = four))
        manual <- suppressMessages(BS_Multi(input, methods = c("AIMS", "sspbc"), Subtype = four))
        expect_identical(auto, manual)
        expect_identical(names(auto$results), c("AIMS", "sspbc"))
        expect_identical(rownames(auto$res_subtypes), colnames(input$se_NC))
        expect_identical(as.character(auto$res_subtypes["N06", c("AIMS", "sspbc")]), c("Her2", "Her2"))
        expect_identical(auto$res_subtypes["N06", "entropy"], 0)
    }
})
