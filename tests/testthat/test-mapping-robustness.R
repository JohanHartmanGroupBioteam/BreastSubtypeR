# Historical example input is pinned independently of public data; see fixtures/README.md.
## Regression tests for probe-to-gene collapsing and Mapping() with a single
## sample. duplicate_genes() used to drop the matrix dimensions of one-sample
## input and then fail in apply() ("dim(X) must have a positive length"), so
## Mapping() could not process a SummarizedExperiment with one sample.

make_probe_input <- function() {
    x <- matrix(
        c(
            1, 2, 3, # pA (gene 100)
            10, 20, 30, # pB (gene 100, larger)
            5, 5, 5, # pC (gene 200)
            7, 8, 9 # pD (no Entrez ID)
        ),
        nrow = 4,
        byrow = TRUE,
        dimnames = list(c("pA", "pB", "pC", "pD"), c("s1", "s2", "s3"))
    )
    y <- data.frame(
        probe = c("pA", "pB", "pC"),
        ENTREZID = c(100L, 100L, 200L),
        stringsAsFactors = FALSE
    )
    list(x = x, y = y)
}

test_that("duplicate_genes() returns a gene-by-sample matrix with dimnames", {
    input <- make_probe_input()
    res <- BreastSubtypeR:::duplicate_genes(input$x, input$y, "max")
    expect_true(is.matrix(res))
    expect_identical(dim(res), c(2L, 3L))
    expect_identical(rownames(res), c("100", "200"))
    expect_identical(colnames(res), c("s1", "s2", "s3"))
    ## "max" keeps the probe with the largest row sum
    expect_identical(unname(res["100", ]), c(10, 20, 30))
    expect_identical(unname(res["200", ]), c(5, 5, 5))
})

test_that("duplicate_genes() handles a single-sample matrix", {
    input <- make_probe_input()
    x1 <- input$x[, "s1", drop = FALSE]
    for (method in c("max", "mean", "median", "iqr")) {
        res <- BreastSubtypeR:::duplicate_genes(x1, input$y, method)
        expect_true(is.matrix(res), info = method)
        expect_identical(dim(res), c(2L, 1L), info = method)
        expect_identical(dimnames(res), list(c("100", "200"), "s1"), info = method)
        expect_identical(unname(res["200", "s1"]), 5, info = method)
    }
    res_max <- BreastSubtypeR:::duplicate_genes(x1, input$y, "max")
    expect_identical(unname(res_max["100", "s1"]), 10)
    ## SD is undefined for each duplicate probe with only one sample (#218).
    expect_error(BreastSubtypeR:::duplicate_genes(x1, input$y, "stdev"),
        class = "BreastSubtypeR_duplicate_undefined")
})

test_that("Mapping() accepts a SummarizedExperiment with one sample", {
    load_historical_example("OSLO2EMIT0obj")
    se <- OSLO2EMIT0obj$se_obj
    full <- suppressMessages(Mapping(
        se_obj = se, RawCounts = FALSE, method = "max",
        impute = TRUE, verbose = FALSE
    ))
    one <- suppressMessages(Mapping(
        se_obj = se[, 1], RawCounts = FALSE, method = "max",
        impute = TRUE, verbose = FALSE
    ))
    expect_identical(ncol(one$se_NC), 1L)
    expect_identical(ncol(one$se_SSP), 1L)
    expect_identical(nrow(one$se_NC), nrow(full$se_NC))
    expect_identical(nrow(one$se_SSP), nrow(full$se_SSP))
    expect_identical(colnames(one$se_NC), colnames(se)[1])
    expect_equal(
        SummarizedExperiment::assay(one$se_NC)[, 1],
        SummarizedExperiment::assay(full$se_NC)[, 1]
    )
    expect_equal(
        SummarizedExperiment::assay(one$se_SSP)[, 1],
        SummarizedExperiment::assay(full$se_SSP)[, 1]
    )
})
