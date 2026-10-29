## Proliferation panel coverage and nonfinite ROR outputs in RORgroup() (#206).
## The synthetic 'out' object fixes the centroid distances, so removing or
## perturbing proliferation genes changes only the proliferation-dependent
## columns and every expectation can be computed by hand from the published
## coefficients: ROR-P = 100 * (-0.0951505484 * LumA + 0.3385116381 * p + 0.35) / 0.85
## with Basal/Her2/LumB distances at zero, and the additional score
## 54.7690 * (-0.3172 * LumA + 0.1981 * p18 + 0.1133 * TSIZE + 0.8826).

prolif11 <- c("CCNB1", "UBE2C", "BIRC5", "KNTC2", "CDC20", "PTTG1", "RRM2",
    "MKI67", "TYMS", "CEP55", "CDCA1")
prolif18 <- c("ANLN", "CEP55", "ORC6L", "CCNE1", "EXO1", "PTTG1", "CDC20",
    "KIF2C", "RRM2", "CDC6", "KNTC2", "TYMS", "CDCA1", "MELK", "UBE2C",
    "CENPF", "MKI67", "UBE2T")
score_cols <- c("ROR-S (Subtype Only)", "ROR-P (Subtype + Prolif)",
    "ROR-C (Subtype + Clinic)", "ROR-PC (Subtype + Clinic + Prolif)",
    "ROR-PC (Subtype + Clinic + Prolif.Subtype)")
group_cols <- c("ROR-S Group (Subtype Only)", "ROR-P Group (Subtype + Prolif)",
    "ROR-C Group (Subtype + Clinic)", "ROR-PC Group (Subtype + Clinic + Prolif)",
    "ROR-PC Group (Subtype + Clinic + Prolif.Subtype)")
prolif_scores <- c("ROR-P (Subtype + Prolif)", "ROR-PC (Subtype + Clinic + Prolif)")
prolif_groups <- c("ROR-P Group (Subtype + Prolif)",
    "ROR-PC Group (Subtype + Clinic + Prolif)")
add_score <- "ROR-PC (Subtype + Clinic + Prolif.Subtype)"
add_group <- "ROR-PC Group (Subtype + Clinic + Prolif.Subtype)"

## Four samples: CCNB1 and BIRC5 belong only to the 11-gene panel, ANLN and
## CENPF only to the 18-gene panel, UBE2C and MKI67 to both.
make_panel_input <- function() {
    ids <- paste0("P", 1:4)
    genes <- c("ESR1", "ERBB2", "CCNB1", "BIRC5", "UBE2C", "MKI67", "ANLN", "CENPF")
    testData <- rbind(
        ESR1 = c(0.5, 0.4, 0.3, 0.2),
        ERBB2 = c(-0.5, -0.4, -0.3, -0.2),
        CCNB1 = c(-1.0, 1.0, 0.0, 0.5),
        BIRC5 = c(-0.5, 0.5, 0.5, 0.5),
        UBE2C = c(-1.5, 1.5, -0.5, 0.5),
        MKI67 = c(-1.0, 1.0, 0.0, 0.5),
        ANLN = c(0.2, 0.4, -0.6, 0.8),
        CENPF = c(-0.4, 0.6, 0.1, -0.3)
    )
    colnames(testData) <- ids
    centroids <- matrix(0, nrow = length(genes), ncol = 5,
        dimnames = list(genes, c("Basal", "Her2", "LumA", "LumB", "Normal")))
    distances <- matrix(0, nrow = 4, ncol = 5,
        dimnames = list(ids, colnames(centroids)))
    distances[, "LumA"] <- c(-0.9, 0.9, 0, 0.3)
    distances[, "Normal"] <- c(0.1, 0.2, 0.3, 0.4)
    out <- list(
        predictions = setNames(c("Basal", "LumA", "Normal", "LumB"), ids),
        predictions.FourSubtype = setNames(c("Basal", "LumA", "LumB", "LumB"), ids),
        testData = testData,
        distances = distances,
        dist.RORSubtype = distances[, 1:4],
        centroids = centroids
    )
    df.cln <- data.frame(PatientID = ids, TSIZE = c(0, 1, 0, 1),
        NODE = c(0, 2, 4, NA), stringsAsFactors = FALSE)
    list(out = out, df.cln = df.cln)
}

run_ror <- function(input, hasClinical = TRUE, Subtype = TRUE) {
    conditions <- character()
    res <- withCallingHandlers(
        suppressMessages(BreastSubtypeR:::RORgroup(
            input$out, input$df.cln, Subtype = Subtype, hasClinical = hasClinical
        )),
        warning = function(w) {
            conditions <<- c(conditions, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    attr(res, "conditions") <- conditions
    res
}

expected_rorp <- function(lumA, p) 100 * (-0.0951505484 * lumA + 0.3385116381 * p + 0.35) / 0.85
expected_rorpc <- function(lumA, p, xT) {
    100 * (-0.090436516 * lumA + 0.131605960 * xT + 0.327259375 * p + 0.35) / 0.85
}
expected_add <- function(lumA, p18, xT) 54.7690 * (-0.3172 * lumA + 0.1981 * p18 + 0.1133 * xT + 0.8826)
parker_group <- function(unscaled, low, high) {
    ifelse(unscaled > high, "high", ifelse(unscaled < low, "low", "med"))
}
keep_genes <- function(input, genes) {
    input$out$testData <- input$out$testData[
        rownames(input$out$testData) %in% genes, , drop = FALSE]
    input
}
expect_valid_or_missing <- function(x, info = NULL) {
    expect_type(x, "character")
    expect_true(all(is.na(x) | x %in% c("low", "med", "high")), info = info)
    expect_false(any(x %in% c("NA", "NaN", "Inf", "-Inf")), info = info)
    expect_false(any(grepl("[0-9]", x)), info = info)
}
expect_structure <- function(res, input, info = NULL) {
    ids <- names(input$out$predictions)
    expect_identical(as.character(res$patientID), ids, info = info)
    expect_identical(rownames(res), ids, info = info)
    expect_identical(as.character(res$Call), unname(input$out$predictions), info = info)
    if ("Call.Subtype" %in% names(res)) {
        expect_identical(as.character(res$Call.Subtype),
            unname(input$out$predictions.FourSubtype), info = info)
    }
}

test_that("complete panels reproduce the published arithmetic for both panels", {
    input <- make_panel_input()
    res <- run_ror(input)
    expect_length(attr(res, "conditions"), 0)
    td <- input$out$testData
    p11 <- colMeans(td[c("CCNB1", "BIRC5", "UBE2C", "MKI67"), ])
    p18 <- colMeans(td[c("UBE2C", "MKI67", "ANLN", "CENPF"), ])
    lumA <- input$out$distances[, "LumA"]
    xT <- input$df.cln$TSIZE
    expect_equal(unname(res[["ROR-P (Subtype + Prolif)"]]), unname(expected_rorp(lumA, p11)))
    expect_equal(unname(res[["ROR-PC (Subtype + Clinic + Prolif)"]]),
        unname(expected_rorpc(lumA, p11, xT)))
    expect_equal(unname(res[[add_score]]), unname(expected_add(lumA, p18, xT)))
    expect_identical(unname(res[["ROR-P Group (Subtype + Prolif)"]]),
        unname(parker_group(-0.0951505484 * lumA + 0.3385116381 * p11, -0.25, 0.1)))
    for (g in group_cols) expect_valid_or_missing(res[[g]], info = g)
    expect_identical(unname(res[[add_group]]), c("med", "high", "high", NA))
    expect_structure(res, input)
})

test_that("one represented gene in either panel no longer drops dimensions", {
    input <- make_panel_input()
    one11 <- keep_genes(input, c("ESR1", "ERBB2", "CCNB1", "ANLN", "CENPF"))
    ## CCNB1 is the only 11-gene panel member left; ANLN and CENPF remain in
    ## the 18-gene panel. The 11-gene score is CCNB1 itself (mean of one gene).
    for (hasClinical in c(TRUE, FALSE)) {
        for (Subtype in c(TRUE, FALSE)) {
            info <- paste("hasClinical", hasClinical, "Subtype", Subtype)
            res <- run_ror(one11, hasClinical = hasClinical, Subtype = Subtype)
            expect_length(attr(res, "conditions"), 0)
            expect_identical(nrow(res), 4L, info = info)
            expect_structure(res, one11, info = info)
            expect_equal(unname(res[["ROR-P (Subtype + Prolif)"]]),
                unname(expected_rorp(one11$out$distances[, "LumA"],
                    one11$out$testData["CCNB1", ])), info = info)
            expect_identical("Call.Subtype" %in% names(res), Subtype, info = info)
            expect_identical(add_score %in% names(res), hasClinical, info = info)
        }
    }
    one18 <- keep_genes(input, c("ESR1", "ERBB2", "CCNB1", "BIRC5", "ANLN"))
    res <- run_ror(one18)
    expect_length(attr(res, "conditions"), 0)
    expect_identical(nrow(res), 4L)
    expect_equal(unname(res[[add_score]]), unname(expected_add(
        one18$out$distances[, "LumA"], one18$out$testData["ANLN", ],
        one18$df.cln$TSIZE)))
    for (g in group_cols) expect_valid_or_missing(res[[g]], info = g)
    ## one sample with complete genes
    single <- input
    single$out$testData <- input$out$testData[, 1, drop = FALSE]
    single$out$distances <- input$out$distances[1, , drop = FALSE]
    single$out$dist.RORSubtype <- input$out$dist.RORSubtype[1, , drop = FALSE]
    single$out$predictions <- input$out$predictions[1]
    single$out$predictions.FourSubtype <- input$out$predictions.FourSubtype[1]
    single$df.cln <- input$df.cln[1, ]
    full <- run_ror(input)
    res1 <- run_ror(single)
    expect_identical(nrow(res1), 1L)
    expect_identical(res1[, names(full)], full[1, names(full)])
})

test_that("partial panels average the represented genes without a threshold", {
    input <- make_panel_input()
    partial <- keep_genes(input, c("ESR1", "ERBB2", "CCNB1", "UBE2C", "ANLN"))
    res <- run_ror(partial)
    expect_length(attr(res, "conditions"), 0)
    td <- partial$out$testData
    p11 <- colMeans(td[c("CCNB1", "UBE2C"), ])
    p18 <- colMeans(td[c("UBE2C", "ANLN"), ])
    lumA <- partial$out$distances[, "LumA"]
    expect_equal(unname(res[["ROR-P (Subtype + Prolif)"]]), unname(expected_rorp(lumA, p11)))
    expect_equal(unname(res[[add_score]]), unname(expected_add(lumA, p18, partial$df.cln$TSIZE)))
    ## per-sample missing genes are averaged over the observed values
    na_one <- input
    na_one$out$testData["CCNB1", "P2"] <- NA
    na_one$out$testData["ANLN", "P3"] <- NaN
    res2 <- run_ror(na_one)
    full <- run_ror(input)
    p11 <- colMeans(input$out$testData[c("CCNB1", "BIRC5", "UBE2C", "MKI67"), ])
    p11["P2"] <- mean(input$out$testData[c("BIRC5", "UBE2C", "MKI67"), "P2"])
    expect_equal(unname(res2[["ROR-P (Subtype + Prolif)"]]),
        unname(expected_rorp(input$out$distances[, "LumA"], p11)))
    p18 <- colMeans(input$out$testData[c("UBE2C", "MKI67", "ANLN", "CENPF"), ])
    p18["P3"] <- mean(input$out$testData[c("UBE2C", "MKI67", "CENPF"), "P3"])
    expect_equal(unname(res2[[add_score]]), unname(expected_add(
        input$out$distances[, "LumA"], p18, input$df.cln$TSIZE)))
    other <- setdiff(names(full), c(prolif_scores, prolif_groups, add_score, add_group))
    expect_identical(res2[, other], full[, other])
    expect_identical(res2[c(1, 4), ], full[c(1, 4), ])
})

test_that("zero represented 11-gene panel genes give missing scores and groups", {
    input <- make_panel_input()
    zero11 <- keep_genes(input, c("ESR1", "ERBB2", "ANLN", "CENPF"))
    full <- run_ror(input)
    for (hasClinical in c(TRUE, FALSE)) {
        res <- run_ror(zero11, hasClinical = hasClinical)
        info <- paste("hasClinical", hasClinical)
        expect_length(attr(res, "conditions"), 0)
        expect_structure(res, zero11, info = info)
        for (s in intersect(prolif_scores, names(res))) {
            expect_type(res[[s]], "double")
            expect_true(all(is.na(res[[s]])), info = s)
            expect_false(any(is.nan(res[[s]])), info = s)
        }
        for (g in intersect(prolif_groups, names(res))) {
            expect_type(res[[g]], "character")
            expect_true(all(is.na(res[[g]])), info = g)
        }
        other <- setdiff(names(res), c(prolif_scores, prolif_groups, add_score, add_group))
        expect_identical(res[, other], full[, other], info = info)
    }
    res <- run_ror(zero11)
    expect_equal(unname(res[[add_score]]), unname(expected_add(
        zero11$out$distances[, "LumA"],
        colMeans(zero11$out$testData[c("ANLN", "CENPF"), ]), zero11$df.cln$TSIZE)))
    path <- tempfile(fileext = ".txt")
    utils::write.table(res, path, sep = "\t", quote = FALSE, row.names = FALSE)
    restored <- utils::read.delim(path, check.names = FALSE)
    unlink(path)
    expect_true(all(is.na(restored[["ROR-P Group (Subtype + Prolif)"]])))
    expect_false(any(restored[["ROR-P Group (Subtype + Prolif)"]] %in% c("NA", "NaN")))
})

test_that("zero represented 18-gene panel genes keep the NODE >= 4 override", {
    input <- make_panel_input()
    zero18 <- keep_genes(input, c("ESR1", "ERBB2", "CCNB1", "BIRC5"))
    full <- run_ror(input)
    res <- run_ror(zero18)
    expect_length(attr(res, "conditions"), 0)
    expect_true(all(is.na(res[[add_score]])))
    expect_false(any(is.nan(res[[add_score]])))
    ## NODE 0 and 1--3 are missing; NODE >= 4 with nonmissing TSIZE stays high;
    ## NA NODE stays missing.
    expect_identical(unname(res[[add_group]]), c(NA, NA, "high", NA))
    other <- setdiff(names(full), c(prolif_scores, prolif_groups, add_score, add_group))
    expect_identical(res[, other], full[, other])
    expect_equal(unname(res[["ROR-P (Subtype + Prolif)"]]), unname(expected_rorp(
        zero18$out$distances[, "LumA"], colMeans(zero18$out$testData[c("CCNB1", "BIRC5"), ]))))
    both <- keep_genes(input, c("ESR1", "ERBB2"))
    res0 <- run_ror(both)
    expect_true(all(is.na(res0[, c(prolif_scores, add_score)])))
    expect_identical(unname(res0[[add_group]]), c(NA, NA, "high", NA))
    for (g in prolif_groups) expect_true(all(is.na(res0[[g]])) && is.character(res0[[g]]))
})

test_that("per-sample NA, NaN, Inf and -Inf panel values give missing outputs for that sample only", {
    input <- make_panel_input()
    full <- run_ror(input)
    ## all 11-panel values missing, all 18-panel values missing, and single
    ## infinite values in panel-specific or shared genes
    cases <- list(
        list(gene = c("CCNB1", "BIRC5", "UBE2C", "MKI67"), value = NA_real_),
        list(gene = c("CCNB1", "BIRC5", "UBE2C", "MKI67"), value = NaN),
        list(gene = "CCNB1", value = Inf),
        list(gene = "BIRC5", value = -Inf),
        list(gene = c("UBE2C", "MKI67", "ANLN", "CENPF"), value = NA_real_),
        list(gene = "ANLN", value = Inf),
        list(gene = "CENPF", value = -Inf),
        list(gene = "UBE2C", value = Inf),
        list(gene = c("UBE2C", "MKI67"), value = c(Inf, -Inf))
    )
    for (sample in c("P1", "P3")) {
        for (case in cases) {
            perturbed <- input
            perturbed$out$testData[case$gene, sample] <- case$value
            res <- run_ror(perturbed)
            info <- paste(sample, paste(case$gene, collapse = "+"),
                paste(case$value, collapse = "/"))
            ## an infinite value makes every panel containing the gene missing;
            ## a missing value makes a panel missing only when no observed value
            ## remains, otherwise that panel is a changed but finite partial mean
            touched11 <- any(case$gene %in% prolif11)
            touched18 <- any(case$gene %in% prolif18)
            hits11 <- touched11
            hits18 <- touched18
            if (all(is.na(case$value))) {
                hits11 <- all(intersect(rownames(input$out$testData), prolif11) %in% case$gene)
                hits18 <- all(intersect(rownames(input$out$testData), prolif18) %in% case$gene)
            }
            changed_scores <- c(if (touched11) prolif_scores, if (touched18) add_score)
            changed_groups <- c(if (touched11) prolif_groups, if (touched18) add_group)
            expect_length(attr(res, "conditions"), 0)
            expect_structure(res, perturbed, info = info)
            affected_scores <- c(if (hits11) prolif_scores, if (hits18) add_score)
            affected_groups <- c(if (hits11) prolif_groups, if (hits18) add_group)
            i <- match(sample, res$patientID)
            for (s in affected_scores) {
                expect_true(is.na(res[[s]][i]) && !is.nan(res[[s]][i]), info = paste(info, s))
            }
            for (g in affected_groups) {
                expect_type(res[[g]], "character")
                expect_valid_or_missing(res[[g]], info = paste(info, g))
                if (g == add_group && sample == "P3") {
                    ## NODE >= 4 keeps the independent high-risk override
                    expect_identical(unname(res[[g]][i]), "high", info = info)
                } else {
                    expect_true(is.na(res[[g]][i]), info = paste(info, g))
                }
            }
            partial <- setdiff(changed_scores, affected_scores)
            for (s in partial) {
                expect_true(is.finite(res[[s]][i]), info = paste(info, s))
            }
            unaffected <- setdiff(names(full), c(changed_scores, changed_groups))
            expect_identical(res[, unaffected], full[, unaffected], info = info)
            expect_identical(res[-i, ], full[-i, ], info = info)
        }
    }
})

test_that("nonfinite numerical ROR scores never receive a risk category", {
    labels <- BreastSubtypeR:::.ror_pc_subtype_labels
    expect_identical(labels(c(A = NA_real_, B = NaN, C = -Inf, D = Inf, E = 50), 40, 60),
        c(A = NA, B = NA, C = NA, D = NA, E = "med"))
    expect_identical(BreastSubtypeR:::.finite_group(c("high", "low", "med", "x"),
        c(1, NA, NaN, Inf)), c("high", NA, NA, NA))
    expect_identical(BreastSubtypeR:::.finite_group(c(0.5, 0.2), c(0.5, 0.2)), c(0.5, 0.2))
    ## an infinite tumor-size covariate makes every clinical score missing
    input <- make_panel_input()
    input$df.cln$TSIZE <- c(Inf, 0, -Inf, 1)
    res <- run_ror(input)
    for (s in c("ROR-C (Subtype + Clinic)", "ROR-PC (Subtype + Clinic + Prolif)", add_score)) {
        expect_true(all(is.na(res[[s]][c(1, 3)])) && !any(is.nan(res[[s]][c(1, 3)])), info = s)
        expect_true(all(is.finite(res[[s]][c(2, 4)])), info = s)
    }
    for (g in c("ROR-C Group (Subtype + Clinic)", "ROR-PC Group (Subtype + Clinic + Prolif)")) {
        expect_identical(is.na(res[[g]]), c(TRUE, FALSE, TRUE, FALSE), info = g)
        expect_valid_or_missing(res[[g]], info = g)
    }
    expect_identical(unname(res[[add_group]]), c(NA, "high", "high", NA))
})

test_that("missing proliferation scores leave sample order, calls and RNG untouched", {
    input <- make_panel_input()
    input$out$testData[c("CCNB1", "BIRC5", "UBE2C", "MKI67"), "P2"] <- NA
    set.seed(206)
    seed_before <- .Random.seed
    res <- run_ror(input)
    expect_identical(.Random.seed, seed_before)
    expect_structure(res, input)
    shuffled <- input
    shuffled$df.cln <- shuffled$df.cln[c(3, 1, 4, 2), ]
    expect_identical(run_ror(shuffled), res)
})

test_that("complete OSLO2EMIT0 inputs reproduce the packaged results exactly", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = environment())
    se <- OSLO2EMIT0obj$data_input$se_NC
    expect_true(all(prolif11 %in% rownames(se)) && all(prolif18 %in% rownames(se)))
    set.seed(1)
    seed_before <- .Random.seed
    fresh <- suppressMessages(BS_cIHC(se, Subtype = FALSE, hasClinical = FALSE))
    expect_identical(.Random.seed, seed_before)
    packaged <- OSLO2EMIT0obj$res$results$cIHC$score.ROR
    expect_identical(fresh$score.ROR[, names(packaged)], packaged)
    expect_identical(fresh$BS.all, OSLO2EMIT0obj$res$results$cIHC$BS.all)
    fresh_ssbc <- suppressMessages(BS_ssBC(se, s = "ER", Subtype = FALSE, hasClinical = FALSE))
    packaged_ssbc <- OSLO2EMIT0obj$res$results$ssBC$score.ROR
    expect_identical(fresh_ssbc$score.ROR[, names(packaged_ssbc)], packaged_ssbc)
    clinical <- suppressMessages(BS_parker(se, calibration = "Internal",
        internal = "medianCtr", Subtype = TRUE, hasClinical = TRUE))
    expect_identical(as.character(clinical$score.ROR$patientID), colnames(se))
    expect_true(all(is.finite(as.matrix(clinical$score.ROR[, score_cols]))))
    for (g in group_cols) expect_true(all(clinical$score.ROR[[g]] %in% c("low", "med", "high")))
    ## explicit recomputation of the proliferation arithmetic from the output
    ol <- clinical$outList
    p11 <- colMeans(ol$testData[prolif11, ])
    rorp <- 100 * (-0.0009299747 * ol$distances[, 1] + 0.0692289192 * ol$distances[, 2] +
        -0.0951505484 * ol$distances[, 3] + 0.0493487685 * ol$distances[, 4] +
        0.3385116381 * p11 + 0.35) / 0.85
    expect_equal(unname(clinical$score.ROR[["ROR-P (Subtype + Prolif)"]]), unname(rorp))
})

test_that("incomplete TCGABRCA panels run as a documented compatibility limitation", {
    data("TCGABRCAobj", package = "BreastSubtypeR", envir = environment())
    se <- TCGABRCAobj$data_input$se_NC
    present11 <- intersect(prolif11, rownames(se))
    present18 <- intersect(prolif18, rownames(se))
    expect_identical(length(present11), 7L)
    expect_identical(length(present18), 8L)
    res <- suppressMessages(BS_parker(se, calibration = "Internal",
        internal = "medianCtr", Subtype = TRUE, hasClinical = FALSE))
    expect_identical(as.character(res$score.ROR$patientID), colnames(se))
    expect_identical(nrow(res$score.ROR), ncol(se))
    expect_true(all(is.finite(res$score.ROR[["ROR-P (Subtype + Prolif)"]])))
    expect_true(all(res$score.ROR[["ROR-P Group (Subtype + Prolif)"]] %in% c("low", "med", "high")))
    expect_identical(as.integer(table(res$score.ROR[["ROR-P Group (Subtype + Prolif)"]])[c("high", "low", "med")]),
        c(50L, 32L, 18L))
    ol <- res$outList
    p11 <- colMeans(ol$testData[present11, , drop = FALSE])
    rorp <- 100 * (-0.0009299747 * ol$distances[, 1] + 0.0692289192 * ol$distances[, 2] +
        -0.0951505484 * ol$distances[, 3] + 0.0493487685 * ol$distances[, 4] +
        0.3385116381 * p11 + 0.35) / 0.85
    expect_equal(unname(res$score.ROR[["ROR-P (Subtype + Prolif)"]]), unname(rorp))
    ## a sample with every represented panel gene missing is the only one affected
    assay_na <- SummarizedExperiment::assay(se)
    assay_na[present11, 1] <- NA
    se_na <- se
    SummarizedExperiment::assay(se_na) <- assay_na
    res_na <- suppressMessages(BS_parker(se_na, calibration = "Internal",
        internal = "medianCtr", Subtype = TRUE, hasClinical = FALSE))
    expect_true(is.na(res_na$score.ROR[["ROR-P (Subtype + Prolif)"]][1]))
    expect_true(is.na(res_na$score.ROR[["ROR-P Group (Subtype + Prolif)"]][1]))
    expect_type(res_na$score.ROR[["ROR-P Group (Subtype + Prolif)"]], "character")
    expect_identical(as.character(res_na$score.ROR$patientID), colnames(se))
})

expect_finite_or_na <- function(x, info = NULL) {
    expect_type(x, "double")
    expect_true(all(is.finite(x) | (is.na(x) & !is.nan(x))), info = info)
}

test_that("every exported ROR score is finite or NA and nonfinite results have NA groups", {
    helper <- BreastSubtypeR:::.finite_score
    expect_identical(helper(c(a = 1, b = NA, c = NaN, d = Inf, e = -Inf)),
        c(a = 1, b = NA_real_, c = NA_real_, d = NA_real_, e = NA_real_))
    expect_false(any(is.nan(helper(c(NaN, 1)))))
    expect_identical(helper(c(x = 0.5, y = -2)), c(x = 0.5, y = -2))
    input <- make_panel_input()
    full <- run_ror(input)
    parker_scores <- c("ROR-S (Subtype Only)", "ROR-P (Subtype + Prolif)",
        "ROR-C (Subtype + Clinic)", "ROR-PC (Subtype + Clinic + Prolif)")
    parker_groups <- c("ROR-S Group (Subtype Only)", "ROR-P Group (Subtype + Prolif)",
        "ROR-C Group (Subtype + Clinic)", "ROR-PC Group (Subtype + Clinic + Prolif)")
    ## 1. infinite tumor size: clinical scores NA, nonclinical untouched
    for (value in c(Inf, -Inf)) {
        tsize <- input
        tsize$df.cln$TSIZE[2] <- value
        res <- run_ror(tsize)
        info <- paste("TSIZE", value)
        for (s in c("ROR-C (Subtype + Clinic)", "ROR-PC (Subtype + Clinic + Prolif)", add_score)) {
            expect_finite_or_na(res[[s]], info = info)
            expect_true(is.na(res[[s]][2]), info = paste(info, s))
        }
        for (g in c("ROR-C Group (Subtype + Clinic)", "ROR-PC Group (Subtype + Clinic + Prolif)")) {
            expect_valid_or_missing(res[[g]], info = info)
            expect_true(is.na(res[[g]][2]), info = paste(info, g))
        }
        ## P2 has NODE = 2: the finite-score cutoffs give no category
        expect_true(is.na(res[[add_group]][2]), info = info)
        unaffected <- c("ROR-S (Subtype Only)", "ROR-S Group (Subtype Only)",
            "ROR-P (Subtype + Prolif)", "ROR-P Group (Subtype + Prolif)", "Call", "Call.Subtype")
        expect_identical(res[, unaffected], full[, unaffected], info = info)
        expect_identical(res[-2, ], full[-2, ], info = info)
    }
    ## 2. nonfinite proliferation scores: already NA before scaling
    prolif <- input
    prolif$out$testData["CCNB1", "P1"] <- Inf
    res <- run_ror(prolif)
    for (s in prolif_scores) expect_true(is.na(res[[s]][1]) && !is.nan(res[[s]][1]), info = s)
    for (g in prolif_groups) expect_true(is.na(res[[g]][1]), info = g)
    ## 3. extreme finite inputs that overflow only in the final scaling
    overflow <- input
    overflow$out$distances[2, "LumA"] <- -1e308 # unscaled scores finite, scaled Inf
    overflow$out$dist.RORSubtype[2, "LumA"] <- -1e308
    unscaled <- 0.04210193 * 0 + -0.35235561 * -1e308
    expect_true(is.finite(unscaled) && !is.finite(100 * (unscaled + 0.35) / 0.85))
    res <- run_ror(overflow)
    for (s in c(parker_scores, add_score)) {
        expect_finite_or_na(res[[s]], info = s)
        expect_true(is.na(res[[s]][2]), info = s)
    }
    for (g in c(parker_groups, add_group)) {
        expect_valid_or_missing(res[[g]], info = g)
        expect_true(is.na(res[[g]][2]), info = g)
    }
    expect_identical(res[-2, ], full[-2, ])
    overflow_p <- input
    overflow_p$out$testData[c("CCNB1", "BIRC5", "UBE2C", "MKI67"), "P4"] <- 1e308
    res <- run_ror(overflow_p)
    expect_true(is.finite(mean(rep(1e308, 4))))
    for (s in prolif_scores) expect_true(is.na(res[[s]][4]) && !is.nan(res[[s]][4]), info = s)
    for (g in prolif_groups) expect_true(is.na(res[[g]][4]), info = g)
    expect_identical(res[-4, ], full[-4, ])
    expect_identical(res[, c("ROR-S (Subtype Only)", "ROR-S Group (Subtype Only)")],
        full[, c("ROR-S (Subtype Only)", "ROR-S Group (Subtype Only)")])
    overflow_t <- input
    overflow_t$df.cln$TSIZE[3] <- 1e308
    res <- run_ror(overflow_t)
    for (s in c("ROR-C (Subtype + Clinic)", "ROR-PC (Subtype + Clinic + Prolif)", add_score)) {
        expect_true(is.na(res[[s]][3]) && !is.nan(res[[s]][3]), info = s)
    }
    for (g in c("ROR-C Group (Subtype + Clinic)", "ROR-PC Group (Subtype + Clinic + Prolif)")) {
        expect_true(is.na(res[[g]][3]), info = g)
    }
    ## 6. NODE >= 4 (P3) keeps the independent override with a missing score
    expect_identical(unname(res[[add_group]][3]), "high")
    ## 4. ordinary finite results: every exported score finite, groups valid
    for (s in c(parker_scores, add_score)) expect_true(all(is.finite(full[[s]])), info = s)
})

test_that("finite exact-equality scores keep the deferred #172 text after the finiteness check", {
    input <- make_panel_input()
    ## zero proliferation values and LumA = -1.0509660919621142 give an unscaled
    ## ROR-P of exactly 0.1, the upper ROR-P cutoff
    input$out$testData[c("CCNB1", "BIRC5", "UBE2C", "MKI67"), "P2"] <- 0
    input$out$distances[2, "LumA"] <- -1.0509660919621142
    input$out$dist.RORSubtype[2, "LumA"] <- -1.0509660919621142
    unscaled <- -0.0951505484 * -1.0509660919621142 + 0.3385116381 * 0
    skip_if_not(identical(unscaled, 0.1), "platform arithmetic does not reach exact equality")
    res <- run_ror(input)
    expect_identical(res[["ROR-P Group (Subtype + Prolif)"]][2], as.character(unscaled))
    expect_identical(res[["ROR-P (Subtype + Prolif)"]][2], 100 * (unscaled + 0.35) / 0.85)
    expect_identical(BreastSubtypeR:::.finite_group(c("0.1", "high", NA), c(52.9, 70, 1)),
        c("0.1", "high", NA))
})

test_that("missing additional scores are NA rather than NaN with the NODE >= 4 override kept", {
    input <- make_panel_input()
    input$out$dist.RORSubtype[1, 1] <- NaN
    input$out$dist.RORSubtype[3, 1] <- NaN
    res <- run_ror(input)
    expect_true(is.na(res[[add_score]][1]) && !is.nan(res[[add_score]][1]))
    expect_true(is.na(res[[add_score]][3]) && !is.nan(res[[add_score]][3]))
    expect_identical(unname(res[[add_group]]), c(NA, "high", "high", NA))
    other <- setdiff(names(res), c(add_score, add_group))
    expect_identical(res[, other], run_ror(make_panel_input())[, other])
})
