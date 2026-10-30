test_that("complete PAM50 ROR-P and original ROR-PC match the Parker reference", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = environment())
    se <- OSLO2EMIT0obj$data_input$se_NC
    expect_identical(dim(se), c(50L, 102L))
    res <- BS_parker(se, calibration = "Internal", internal = "medianCtr",
        Subtype = TRUE, hasClinical = TRUE)
    clinical <- as.data.frame(SummarizedExperiment::colData(se))
    ror198_expect_reference(res, clinical)
    ## Independently executed original distributed script on the complete
    ## packaged input (same official archive/hash as the test oracle).
    at <- match(c("OSLO2EMIT0.001", "OSLO2EMIT0.068"), res$score.ROR$patientID)
    expect_false(anyNA(at))
    expect_equal(res$score.ROR[[ror198_changed[1]]][at],
        c(38.5170379202764, 86.9930315059914), tolerance = 1e-12)
    expect_equal(res$score.ROR[[ror198_changed[3]]][at],
        c(39.1454520393945, 100.380820648158), tolerance = 1e-12)
    expect_identical(res$score.ROR[[ror198_changed[2]]][at], c("med", "high"))
    expect_identical(res$score.ROR[[ror198_changed[4]]][at], c("med", "high"))

    ## Every input is resolved to represented PAM50 rows before ROR scoring.
    expanded <- SummarizedExperiment::SummarizedExperiment(
        assays = list(rbind(as.matrix(SummarizedExperiment::assay(se)),
            NOT_A_PAM50_GENE = seq_len(ncol(se)) * 1e6)),
        colData = SummarizedExperiment::colData(se))
    expect_identical(BS_parker(expanded, calibration = "Internal",
        internal = "medianCtr", Subtype = TRUE, hasClinical = TRUE), res)
})

test_that("the ROR standardization is within each sample and never mutates its input", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = environment())
    se <- OSLO2EMIT0obj$data_input$se_NC
    clinical <- as.data.frame(SummarizedExperiment::colData(se))
    out <- BS_parker(se, calibration = "Internal", Subtype = TRUE)$outList
    snapshot <- serialize(out, NULL)
    set.seed(198)
    rng <- .Random.seed
    score <- BreastSubtypeR:::RORgroup(out, clinical, Subtype = TRUE, hasClinical = TRUE)
    expect_identical(serialize(out, NULL), snapshot)
    expect_identical(.Random.seed, rng)
    shifted <- out
    shifted$testData <- sweep(sweep(out$testData, 2,
        seq(0.5, 3, length.out = ncol(se)), "*"), 2, seq_len(ncol(se)), "+")
    shifted_score <- BreastSubtypeR:::RORgroup(shifted, clinical,
        Subtype = TRUE, hasClinical = TRUE)
    for (nm in ror198_changed[c(1, 3)]) {
        expect_equal(shifted_score[[nm]], score[[nm]], tolerance = 1e-12)
    }
    for (nm in ror198_changed[c(2, 4)]) expect_identical(shifted_score[[nm]], score[[nm]])
    ## Wrong-axis and 11-gene-only standardization cannot satisfy this oracle.
    p <- ror198_proliferation(out$testData)
    wrong_axis <- colMeans(t(scale(t(out$testData)))[ror198_panel, ])
    panel_only <- colMeans(scale(out$testData[ror198_panel, ]))
    expect_gt(max(abs(p - wrong_axis)), 0.1)
    expect_gt(max(abs(p - panel_only)), 0.1)
})

test_that("all NC methods use the final PAM50 matrix for original ROR only", {
    methods <- c("parker.original", "genefu.scale", "genefu.robust", "cIHC",
        "cIHC.itr", "PCAPAM50", "ssBC", "ssBC.v2")
    for (dataset in c("OSLO2EMIT0obj", "TCGABRCAobj")) {
        data(list = dataset, package = "BreastSubtypeR", envir = environment())
        se <- get(dataset)$data_input$se_NC
        # The packaged TCGA subset has no TSIZE/NODE; exercise clinical code
        # with explicit synthetic covariates in this temporary test input only.
        if (dataset == "TCGABRCAobj") {
            SummarizedExperiment::colData(se)$TSIZE <- rep(c(0, 1), length.out = ncol(se))
            SummarizedExperiment::colData(se)$NODE <- rep(c(0, 1, 2, 3, 4), length.out = ncol(se))
        }
        clinical <- as.data.frame(SummarizedExperiment::colData(se))
        for (method in methods) for (hasClinical in c(FALSE, TRUE)) for (Subtype in c(FALSE, TRUE)) {
            info <- paste(dataset, method, hasClinical, Subtype)
            set.seed(198)
            rng <- .Random.seed
            res <- suppressMessages(ror198_run(method, se, hasClinical, Subtype))
            expect_identical(.Random.seed, rng, info = info)
            expect_identical(res$BS.all$PatientID, colnames(se), info = info)
            expect_identical(as.character(res$score.ROR$patientID), colnames(se), info = info)
            expect_identical(colnames(res$outList$testData), colnames(se), info = info)
            ror198_expect_reference(res, if (hasClinical) clinical else NULL)
            ## Re-enter scoring with the very same classifier object; every
            ## field and internal axis must remain byte-for-byte intact.
            before <- serialize(res$outList, NULL)
            rescored <- BreastSubtypeR:::RORgroup(res$outList, clinical,
                hasClinical = hasClinical, Subtype = Subtype)
            expect_identical(serialize(res$outList, NULL), before, info = info)
            if (hasClinical) {
                p18 <- c("ANLN", "CEP55", "ORC6L", "CCNE1", "EXO1", "PTTG1",
                    "CDC20", "KIF2C", "RRM2", "CDC6", "KNTC2", "TYMS", "CDCA1",
                    "MELK", "UBE2C", "CENPF", "MKI67", "UBE2T")
                td <- res$outList$testData
                raw18 <- colMeans(td[rownames(td) %in% p18, , drop = FALSE], na.rm = TRUE)
                d <- res$outList$dist.RORSubtype
                tsize <- clinical$TSIZE[match(rownames(d), clinical$PatientID)]
                expected <- 54.7690 * (-0.0067 * d[, 1] + 0.4317 * d[, 2] -
                    0.3172 * d[, 3] + 0.4894 * d[, 4] + 0.1981 * raw18[rownames(d)] +
                    0.1133 * tsize + 0.8826)
                expect_equal(unname(res$score.ROR[["ROR-PC (Subtype + Clinic + Prolif.Subtype)"]]),
                    unname(expected), tolerance = 1e-12, info = info)
            }
        }
    }
})

test_that("all Parker calibration variants standardize after calibration", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = environment())
    se <- OSLO2EMIT0obj$data_input$se_NC
    clinical <- as.data.frame(SummarizedExperiment::colData(se))
    variants <- list(list(calibration = "None"),
        list(calibration = "Internal", internal = "medianCtr"),
        list(calibration = "Internal", internal = "meanCtr"),
        list(calibration = "Internal", internal = "qCtr"),
        list(calibration = "External", external = "RNAseq.V2"),
        list(calibration = "External", external = "Given.mdns", medians =
            data.frame(X = rownames(se), medians = apply(SummarizedExperiment::assay(se), 1, median))))
    for (args in variants) for (hasClinical in c(FALSE, TRUE)) for (Subtype in c(FALSE, TRUE)) {
        res <- do.call(BS_parker, c(list(se_obj = se, hasClinical = hasClinical,
            Subtype = Subtype), args))
        ror198_expect_reference(res, if (hasClinical) clinical else NULL)
    }
})

ror198_synthetic <- function() {
    data("BreastSubtypeRobj", package = "BreastSubtypeR", envir = environment())
    genes <- rownames(BreastSubtypeRobj$centroid)
    ids <- c("sample-z", "sample-a", "sample-m")
    td <- outer(seq_along(genes), 1:3, function(g, s) sin(g + s) + g / 20)
    dimnames(td) <- list(genes, ids)
    d <- matrix(rep(c(0.1, 0.2, -0.3, 0.4, -0.1), each = 3), nrow = 3,
        dimnames = list(ids, c("Basal", "Her2", "LumA", "LumB", "Normal")))
    out <- list(predictions = setNames(rep("LumB", 3), ids),
        predictions.FourSubtype = setNames(rep("LumB", 3), ids),
        testData = td, distances = d, dist.RORSubtype = d[, 1:4],
        centroids = BreastSubtypeRobj$centroid)
    list(out = out, clinical = data.frame(PatientID = ids, TSIZE = c(0, 1, 0),
        NODE = c(0, 2, 4)))
}

test_that("zero, one, partial and complete coverage retain dimensions and finite-or-NA outputs", {
    input <- ror198_synthetic()
    sets <- list(zero = character(), one_PAM50 = "CCNB1", one_nonprolif = "ESR1",
        one_prolif = c("CCNB1", "ESR1", "ERBB2"), partial = c("CCNB1", "MKI67", "ESR1", "ANLN"),
        all_prolif = ror198_panel, complete = rownames(input$out$testData),
        zero_prolif = setdiff(rownames(input$out$testData), ror198_panel))
    for (label in names(sets)) for (hasClinical in c(FALSE, TRUE)) for (Subtype in c(FALSE, TRUE)) {
        out <- input$out
        out$testData <- out$testData[sets[[label]], , drop = FALSE]
        snapshot <- serialize(out, NULL)
        res <- BreastSubtypeR:::RORgroup(out, input$clinical,
            hasClinical = hasClinical, Subtype = Subtype)
        expect_identical(serialize(out, NULL), snapshot, info = label)
        expect_identical(as.character(res$patientID), colnames(out$testData), info = label)
        ror198_expect_reference(list(outList = out, score.ROR = res),
            if (hasClinical) input$clinical else NULL)
        for (nm in intersect(ror198_changed[c(1, 3)], names(res))) {
            expect_true(all(is.finite(res[[nm]]) | (is.na(res[[nm]]) & !is.nan(res[[nm]]))), info = label)
        }
        if (label %in% c("zero", "one_PAM50", "one_nonprolif", "zero_prolif")) {
            expect_true(all(is.na(res[[ror198_changed[1]]])), info = label)
            expect_true(all(is.na(res[[ror198_changed[2]]])), info = label)
        } else expect_true(all(is.finite(res[[ror198_changed[1]]])), info = label)
    }
})

test_that("undefined sample standardization never supplies a fallback proliferation score", {
    input <- ror198_synthetic()
    td <- input$out$testData
    cases <- list(all_NA = rep(NA_real_, nrow(td)), all_NaN = rep(NaN, nrow(td)),
        constant = rep(3, nrow(td)), one_finite = c(3, rep(NA_real_, nrow(td) - 1)),
        overflow_SD = seq(-1, 1, length.out = nrow(td)) * 1e308)
    for (value in c(Inf, -Inf)) {
        sample <- td[, 2]
        sample[match("ESR1", rownames(td))] <- value # outside both proliferation panels
        cases[[as.character(value)]] <- sample
    }
    base <- BreastSubtypeR:::RORgroup(input$out, input$clinical, hasClinical = TRUE, Subtype = TRUE)
    for (label in names(cases)) {
        out <- input$out
        out$testData[, 2] <- cases[[label]]
        res <- BreastSubtypeR:::RORgroup(out, input$clinical, hasClinical = TRUE, Subtype = TRUE)
        for (nm in ror198_changed) {
            expect_true(is.na(res[[nm]][2]), info = paste(label, nm))
            if (is.numeric(res[[nm]])) expect_false(is.nan(res[[nm]][2]), info = label)
        }
        expect_identical(res[-2, ], base[-2, ], info = label)
        protected <- c("ROR-S (Subtype Only)", "ROR-S Group (Subtype Only)",
            "ROR-C (Subtype + Clinic)", "ROR-C Group (Subtype + Clinic)", "Call", "Call.Subtype")
        expect_identical(res[, protected], base[, protected], info = label)
        if (label %in% c("Inf", "-Inf")) {
            extra <- c("ROR-PC (Subtype + Clinic + Prolif.Subtype)",
                "ROR-PC Group (Subtype + Clinic + Prolif.Subtype)")
            expect_identical(res[, extra], base[, extra], info = label)
        }
    }
    ## No seed is created by scoring, including the undefined cases.
    withr::with_preserve_seed({
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv)
        invisible(BreastSubtypeR:::RORgroup(input$out, input$clinical, hasClinical = TRUE))
        expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    })
})
