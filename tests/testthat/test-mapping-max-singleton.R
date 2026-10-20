mapping178_annotation <- function(x, ids = seq_len(nrow(x))) {
    data.frame(probe = rownames(x), ENTREZID = ids)
}

mapping178_input <- function(n = 3L) {
    env <- new.env()
    data("BreastSubtypeRobj", package = "BreastSubtypeR", envir = env)
    ref <- env$BreastSubtypeRobj
    ids <- rev(as.character(ref$genes.signature$EntrezGene.ID))
    samples <- c("sample 3", "s/1", "S.2")[seq_len(n)]
    counts <- outer(seq_along(ids), seq_len(n),
        function(i, j) 20 + (i * 37 + j * 83 + i * j * 17) %% 4000)
    dimnames(counts) <- list(paste0("probe_", ids), samples)
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts),
        rowData = S4Vectors::DataFrame(
            probe = rownames(counts), ENTREZID = ids,
            Length = 800 + (seq_along(ids) %% 7) * 137),
        colData = S4Vectors::DataFrame(PatientID = samples))
    list(se = se, ref = ref)
}

mapping178_map <- function(se, raw = FALSE, impute = FALSE) {
    suppressMessages(Mapping(se, RawCounts = raw, method = "max",
        impute = impute, verbose = FALSE))
}

mapping178_expected <- function(nc, ssp, input) {
    ref <- input$ref
    ids <- as.character(SummarizedExperiment::rowData(input$se)$ENTREZID)
    rownames(nc) <- rownames(ssp) <- ids
    nc <- nc[as.character(ref$genes.sig50$EntrezGene.ID), , drop = FALSE]
    rownames(nc) <- ref$genes.sig50$Symbol
    nc <- as.data.frame(nc)
    colnames(nc) <- colnames(input$se)
    ssp_ids <- as.character(ref$genes.signature$EntrezGene.ID[
        ref$genes.signature$SSP_based == "Yes"])
    list(nc = nc, ssp = ssp[ssp_ids, , drop = FALSE])
}

mapping178_expect_output <- function(actual, expected, se) {
    expect_identical(names(actual), c("se_NC", "se_SSP"))
    expect_identical(SummarizedExperiment::assay(actual$se_NC), expected$nc)
    expect_identical(SummarizedExperiment::assay(actual$se_SSP), expected$ssp)
    for (out in actual) {
        expect_s4_class(out, "SummarizedExperiment")
        expect_identical(colnames(out), colnames(se))
        expect_identical(SummarizedExperiment::colData(out),
            SummarizedExperiment::colData(se))
    }
}

test_that("max preserves the exact singleton-NA reproducer from issue 178", {
    x <- matrix(c(10, NA, 30, 4, 5, 6), nrow = 2, byrow = TRUE,
        dimnames = list(c("pA", "pB"), c("s1", "s2", "s3")))
    y <- mapping178_annotation(x, c(55165L, 100L))
    expected <- x
    rownames(expected) <- c("55165", "100")
    expect_identical(BreastSubtypeR:::duplicate_genes(x, y, "max"), expected)
})

test_that("max singleton values retain type, dimensions, identity and order", {
    x <- rbind(partial = c(2, NA, 4), missing = rep(NA_real_, 3),
        nan = c(NaN, 3, NaN), complete = c(7, 1, 8),
        nonfinite = c(Inf, -Inf, 5), negative = rep(-Inf, 3))
    colnames(x) <- c("sample 3", "s/1", "S.2")
    ids <- c(55165L, 100L, 20L, 9L, 300L, 4L)
    y <- mapping178_annotation(x, ids)[rev(seq_len(nrow(x))), ]
    for (rows in c(list(seq_len(nrow(x))), as.list(seq_len(nrow(x))))) {
        for (cols in list(seq_len(ncol(x)), 1L, 2L)) {
            selected <- x[rows, cols, drop = FALSE]
            expected <- selected
            rownames(expected) <- as.character(ids[rows])
            expect_identical(
                BreastSubtypeR:::duplicate_genes(selected, y, "max"), expected)
        }
    }
    integer_input <- matrix(1:2, nrow = 2,
        dimnames = list(c("a", "b"), "sample"))
    expected <- matrix(c(1, 2), nrow = 2,
        dimnames = list(c("9", "3"), "sample"))
    expect_identical(BreastSubtypeR:::duplicate_genes(integer_input,
        mapping178_annotation(integer_input, c(9L, 3L)), "max"), expected)
})

test_that("max combines a missing singleton with complete duplicate groups", {
    x <- rbind(singleton = c(1, NA, 3), low = c(10, 1, 1),
        high = c(2, 7, 7))
    colnames(x) <- c("sample 3", "s/1", "S.2")
    y <- mapping178_annotation(x, c(200L, 100L, 100L))
    expected <- x[c("singleton", "high"), , drop = FALSE]
    rownames(expected) <- c("200", "100")
    expect_identical(BreastSubtypeR:::duplicate_genes(x, y, "max"), expected)
})

test_that("max still selects one cohort-wide probe and keeps finite tie order", {
    x <- rbind(singleton = c(1, 2, 3), low = c(10, 1, 1),
        high = c(2, 7, 7), tied = c(8, 4, 4))
    colnames(x) <- c("sample 3", "s/1", "S.2")
    y <- mapping178_annotation(x, c(200L, 100L, 100L, 100L))
    for (probes in list(rownames(x), c("tied", "singleton", "high", "low"))) {
        selected <- x[probes, , drop = FALSE]
        expected <- if (probes[1] == "singleton") {
            x[c("singleton", "high"), , drop = FALSE]
        } else {
            x[c("tied", "singleton"), , drop = FALSE]
        }
        rownames(expected) <- if (probes[1] == "singleton") {
            c("200", "100")
        } else {
            c("100", "200")
        }
        expect_identical(BreastSubtypeR:::duplicate_genes(selected, y, "max"),
            expected)
    }
    # In one sample the same original maximum rule selects the largest value.
    expected <- matrix(c(1, 10), nrow = 2,
        dimnames = list(c("200", "100"), "sample 3"))
    expect_identical(BreastSubtypeR:::duplicate_genes(x[, 1, drop = FALSE],
        y, "max"), expected)
})

test_that("max duplicate-group missingness keeps baseline results and errors", {
    x <- rbind(a = c(1, NA, 3), b = c(2, 5, 4))
    colnames(x) <- c("sample 3", "s/1", "S.2")
    y <- mapping178_annotation(x, c(100L, 100L))
    expected <- x["b", , drop = FALSE]
    rownames(expected) <- "100"
    expect_identical(BreastSubtypeR:::duplicate_genes(x, y, "max"), expected)
    x["a", ] <- c(Inf, -Inf, 3) # NaN row sum still loses to the finite row.
    expect_identical(BreastSubtypeR:::duplicate_genes(x, y, "max"), expected)
    for (missing in list(NA_real_, NaN)) {
        x[,] <- missing
        expect_error(BreastSubtypeR:::duplicate_genes(x, y, "max"),
            "length of 'dimnames' \\[2\\] not equal to array extent")
        with_singleton <- rbind(x, complete = c(1, 2, 3))
        annotation <- mapping178_annotation(with_singleton, c(100L, 100L, 200L))
        expect_error(BreastSubtypeR:::duplicate_genes(
            with_singleton, annotation, "max"),
            "length of 'dimnames' \\[1\\] not equal to array extent")
    }
})

test_that("mean median iqr and stdev retain their existing missing-value behavior", {
    x <- rbind(a = c(1, 3, 5), b = c(9, 7, 5),
        partial = c(2, NA, 6), missing = rep(NA_real_, 3),
        nan = c(NaN, 4, 8))
    colnames(x) <- c("sample 3", "s/1", "S.2")
    y <- mapping178_annotation(x, c(200L, 200L, 100L, 300L, 400L))
    for (method in c("mean", "median", "iqr", "stdev")) {
        expected <- rbind(c(1, 3, 5), c(2, NA, 6),
            rep(NA_real_, 3), c(NaN, 4, 8))
        if (method %in% c("mean", "median")) expected[1, ] <- c(5, 5, 5)
        if (method == "mean") {
            expected[2, 2] <- NaN
            expected[3, ] <- NaN
        }
        if (method == "median") expected[4, 1] <- NA_real_
        dimnames(expected) <- list(c("200", "100", "300", "400"), colnames(x))
        expect_identical(BreastSubtypeR:::duplicate_genes(x, y, method), expected)
    }
})

test_that("nonraw Mapping retains missing singleton measurements in both outputs", {
    input <- mapping178_input()
    SummarizedExperiment::assay(input$se) <-
        log2(SummarizedExperiment::assay(input$se) + 1)
    target <- match("55165", SummarizedExperiment::rowData(input$se)$ENTREZID)
    for (values in list(c(10, NA, 30), rep(NA_real_, 3), c(NaN, 5, 6))) {
        SummarizedExperiment::assay(input$se)[target, ] <- values
        x <- SummarizedExperiment::assay(input$se)
        expected <- mapping178_expected(x, 2^x, input)
        mapping178_expect_output(mapping178_map(input$se), expected, input$se)
        for (cols in list(1:3, which(is.na(values))[1])) {
            one <- input$se[target, cols, drop = FALSE]
            result <- mapping178_map(one)
            nc <- as.data.frame(SummarizedExperiment::assay(one))
            rownames(nc) <- "CEP55"
            colnames(nc) <- colnames(one)
            ssp <- 2^SummarizedExperiment::assay(one)
            rownames(ssp) <- "55165"
            mapping178_expect_output(result, list(nc = nc, ssp = ssp), one)
        }
    }
})

test_that("nonraw Mapping preserves a missing singleton alongside duplicate genes", {
    input <- mapping178_input()
    se <- input$se
    SummarizedExperiment::assay(se) <- log2(SummarizedExperiment::assay(se) + 1)
    ids <- as.character(SummarizedExperiment::rowData(se)$ENTREZID)
    target <- match("55165", ids)
    alternative <- match(as.character(input$ref$genes.sig50$EntrezGene.ID[1]), ids)
    SummarizedExperiment::assay(se)[target, 2] <- NA_real_
    extra <- se[alternative, , drop = FALSE]
    rownames(extra) <- "alternative_probe"
    SummarizedExperiment::rowData(extra)$probe <- "alternative_probe"
    SummarizedExperiment::assay(extra) <- SummarizedExperiment::assay(extra) + 1
    mixed <- rbind(extra, se)
    expected_source <- SummarizedExperiment::assay(se)
    expected_source[alternative, ] <- SummarizedExperiment::assay(extra)[1, ]
    expected <- mapping178_expected(expected_source, 2^expected_source, input)
    mapping178_expect_output(mapping178_map(mixed), expected, mixed)
})

test_that("complete Mapping controls preserve both scales and impute settings", {
    for (n in c(1L, 3L)) {
        input <- mapping178_input(n)
        counts <- SummarizedExperiment::assay(input$se)
        for (raw in c(FALSE, TRUE)) {
            se <- input$se
            if (raw) {
                # Existing normalization on the full count matrix; no new policy.
                y <- edgeR::DGEList(counts = round(counts))
                keep <- rowMeans(edgeR::cpm(y) > 0.1) >= 0.30
                subset <- edgeR::calcNormFactors(y[keep, , keep.lib.sizes = TRUE],
                    method = "upperquartile", p = 0.75)
                y$samples$norm.factors <- subset$samples$norm.factors
                nc <- edgeR::cpm(y, normalized.lib.sizes = TRUE,
                    log = TRUE, prior.count = 1)
                ssp <- edgeR::rpkm(y,
                    gene.length = SummarizedExperiment::rowData(se)$Length,
                    normalized.lib.sizes = FALSE, log = FALSE)
            } else {
                nc <- log2(counts + 1)
                SummarizedExperiment::assay(se) <- nc
                ssp <- 2^nc
            }
            seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
            actual <- mapping178_map(se, raw = raw)
            expected <- mapping178_expected(nc, ssp, input)
            mapping178_expect_output(actual, expected, se)
            expect_identical(mapping178_map(se, raw = raw, impute = TRUE), actual)
            expect_identical(get0(".Random.seed", envir = .GlobalEnv,
                inherits = FALSE), seed)
        }
    }
})

test_that("nonraw imputation and early rejection of raw NA counts stay unchanged", {
    withr::local_seed(178)
    input <- mapping178_input()
    raw <- input$se
    target <- match("55165", SummarizedExperiment::rowData(raw)$ENTREZID)
    SummarizedExperiment::assay(raw)[target, 2] <- NA_real_
    for (impute in c(FALSE, TRUE)) {
        expect_error(mapping178_map(raw, raw = TRUE, impute = impute),
            "NA counts not allowed", fixed = TRUE)
    }
    se <- input$se
    SummarizedExperiment::assay(se) <- log2(SummarizedExperiment::assay(se) + 1)
    SummarizedExperiment::assay(se)[target, 2] <- NA_real_
    filled <- se
    SummarizedExperiment::assay(filled) <- BreastSubtypeR:::impute_missing(
        SummarizedExperiment::assay(se), verbose = FALSE)
    actual <- mapping178_map(se, impute = TRUE)
    expect_identical(actual, mapping178_map(filled))
    expect_false(anyNA(SummarizedExperiment::assay(actual$se_NC)))
    expect_false(anyNA(SummarizedExperiment::assay(actual$se_SSP)))
})
