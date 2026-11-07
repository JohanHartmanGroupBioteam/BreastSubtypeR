# All model targets in these two packaged inputs have one annotated probe.
# Join model ID -> probe label -> assay row independently of Mapping().
map01_example_oracle <- function(se, raw) {
    models <- new.env(parent = emptyenv())
    utils::data(list = c("BreastSubtypeRobj", "sspbc.models"),
        package = "BreastSubtypeR", envir = models)
    signature <- models$BreastSubtypeRobj$genes.signature
    pam50 <- models$BreastSubtypeRobj$genes.sig50
    pairs <- strsplit(models$sspbc.models$ssp.ror$all.pairs, "<", fixed = TRUE)
    ror_ids <- unique(sub("^e", "", unlist(pairs, use.names = FALSE)))
    annotation <- SummarizedExperiment::rowData(se)
    project <- function(matrix, ids, labels = ids) {
        ids <- as.character(ids)
        annotated <- as.character(annotation$ENTREZID)
        matches <- annotated %in% ids
        stopifnot(!anyDuplicated(annotated[matches]), all(ids %in% annotated))
        probes <- as.character(annotation$probe)[match(ids, annotated)]
        rows <- match(probes, rownames(matrix))
        stopifnot(!anyNA(rows))
        out <- matrix[rows, , drop = FALSE]
        rownames(out) <- as.character(labels)
        out
    }
    nc <- SummarizedExperiment::assay(se)
    if (raw) {
        counts <- edgeR::DGEList(counts = round(nc))
        keep <- rowMeans(edgeR::cpm(counts) > 0.1) >= 0.30
        subset <- counts[keep, , keep.lib.sizes = TRUE]
        counts$samples$norm.factors <- edgeR::calcNormFactors(subset,
            method = "upperquartile")$samples$norm.factors
        nc <- edgeR::cpm(counts, log = TRUE, prior.count = 1)
        linear <- edgeR::rpkm(counts, gene.length = as.numeric(annotation$Length),
            normalized.lib.sizes = FALSE)
    } else {
        linear <- 2^nc
    }
    list(se_NC = project(nc, pam50$EntrezGene.ID, pam50$Symbol),
        se_SSP = project(linear, signature$EntrezGene.ID[signature$SSP_based == "Yes"]),
        se_SSP_ROR = project(linear, ror_ids))
}

test_that("packaged log expression and raw counts map factor probes by label", {
    for (name in c("OSLO2EMIT0obj", "TCGABRCAobj")) {
        raw <- name == "TCGABRCAobj"
        se <- current_example(name)$se_obj
        factor_se <- se
        probes <- SummarizedExperiment::rowData(se)$probe
        SummarizedExperiment::rowData(factor_se)$probe <- factor(probes)
        before <- serialize(factor_se, NULL)
        expected <- map01_example_oracle(se, raw)
        expect_identical(unname(vapply(expected, nrow, integer(1))), c(50L, 339L, 296L))
        for (include in c(FALSE, TRUE)) {
            panels <- if (include) names(expected) else c("se_NC", "se_SSP")
            for (impute in c(FALSE, TRUE)) {
                reference <- Mapping(se, RawCounts = raw, impute = impute,
                    include_ssp_ror = include, verbose = FALSE)
                actual <- Mapping(factor_se, RawCounts = raw, impute = impute,
                    include_ssp_ror = include, verbose = FALSE)
                expect_identical(actual, reference)
                expect_identical(names(actual), panels)
                for (panel in panels) {
                    output <- actual[[panel]]
                    expect_identical(as.matrix(SummarizedExperiment::assay(output)),
                        expected[[panel]], info = paste(name, panel))
                    expect_identical(SummarizedExperiment::colData(output),
                        SummarizedExperiment::colData(se))
                }
            }
        }
        expect_identical(serialize(factor_se, NULL), before)
    }
})

test_that("OSLO factor mapping restores Parker research ROR and AUTO calls", {
    withr::local_options(BreastSubtypeR.verbose = FALSE)
    oslo <- current_example("OSLO2EMIT0obj")
    se <- oslo$se_obj
    character_input <- Mapping(se, verbose = FALSE)
    SummarizedExperiment::rowData(se)$probe <- factor(SummarizedExperiment::rowData(se)$probe)
    factor_input <- Mapping(se, verbose = FALSE)
    parker <- function(mapped) suppressMessages(suppressWarnings(
        BS_parker(mapped$se_NC, calibration = "None", hasClinical = TRUE),
        classes = "bsr_experimental_ror_pc_warning"))
    expected <- parker(character_input)
    actual <- parker(factor_input)
    expect_identical(actual$BS.all, expected$BS.all)
    call <- actual$BS.all$BS[match("OSLO2EMIT0.012", actual$BS.all$PatientID)]
    expect_identical(as.character(call), "Her2")
    expect_identical(actual$score.ROR, expected$score.ROR)

    withr::local_seed(20260831, .rng_kind = "Mersenne-Twister",
        .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
    auto <- suppressMessages(BS_Multi(factor_input, methods = "AUTO", Subtype = FALSE,
        hasClinical = FALSE, er_subsample = "target"))
    methods <- c("genefu.robust", "AIMS", "sspbc")
    expect_identical(auto$res_subtypes[, methods], oslo$res$res_subtypes[, methods])
})
