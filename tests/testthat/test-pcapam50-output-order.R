## Check all sample-indexed outputs, with complete identities and structure
## asserted before matching or subsetting. Historical numeric fixtures retain
## expect_equal() tolerances for differences between supported environments.
pcapam50_order_assert <- function(ok, what) {
    expect_true(ok, info = what)
    if (!isTRUE(ok)) stop(what, call. = FALSE)
}

pcapam50_order_ids <- function(x, ids) {
    pcapam50_order_assert(
        is.character(x) && is.character(ids) && length(x) == length(ids) &&
            !anyNA(x) && !anyNA(ids) && !anyDuplicated(x) &&
            !anyDuplicated(ids) && all(nzchar(x)) && all(nzchar(ids)) &&
            setequal(x, ids),
        "Sample or gene identities must be complete, unique and one-to-one"
    )
}

pcapam50_order_structure <- function(res, ids, genes, subtype, clinical) {
    pcapam50_order_assert(
        identical(names(res), c("BS.all", "score.ROR", "mdns.fl", "outList")),
        "Required result components must be present"
    )
    ol <- res$outList
    pcapam50_order_assert(
        !anyDuplicated(names(ol)) &&
            all(c("predictions", "testData", "distances", "dist.RORSubtype",
                  "centroids") %in% names(ol)) &&
            identical(!is.null(ol$predictions.FourSubtype), subtype) &&
            identical(names(res$BS.all),
                      c("PatientID", "BS", if (subtype) "BS.Subtype")) &&
            !anyDuplicated(names(res$score.ROR)) &&
            "patientID" %in% names(res$score.ROR),
        "Required output names must be present exactly once"
    )
    axes <- list(
        res$BS.all$PatientID, as.character(res$score.ROR$patientID),
        names(ol$predictions), rownames(ol$distances),
        rownames(ol$dist.RORSubtype), colnames(ol$testData)
    )
    if (subtype) axes <- c(axes, list(names(ol$predictions.FourSubtype)))
    for (axis in axes) pcapam50_order_ids(axis, ids)
    n <- length(ids)
    g <- length(genes)
    pcapam50_order_assert(
        identical(dim(res$BS.all), c(n, if (subtype) 3L else 2L)) &&
            identical(dim(res$score.ROR),
                      c(n, 14L + as.integer(subtype) + 6L * clinical)) &&
            identical(dim(ol$testData), c(g, n)) &&
            identical(dim(ol$distances), c(n, 5L)) &&
            identical(dim(ol$dist.RORSubtype), c(n, 4L)) &&
            identical(dim(ol$centroids), c(g, 5L)) &&
            nrow(res$mdns.fl) == g,
        "Sample and gene dimensions must be preserved"
    )
    pcapam50_order_ids(rownames(ol$testData), genes)
    pcapam50_order_ids(rownames(ol$centroids), genes)
    pcapam50_order_ids(rownames(res$mdns.fl), genes)
    pcapam50_order_ids(colnames(ol$distances), colnames(ol$centroids))
    pcapam50_order_assert(
        !anyNA(colnames(ol$dist.RORSubtype)) &&
            !anyDuplicated(colnames(ol$dist.RORSubtype)) &&
            all(colnames(ol$dist.RORSubtype) %in% colnames(ol$centroids)) &&
            !anyDuplicated(colnames(res$mdns.fl)) &&
            "v1PAM.mdns" %in% colnames(res$mdns.fl) &&
            !anyNA(res$BS.all$BS) && !anyNA(ol$predictions),
        "Gene summaries, subtype names and calls must be complete"
    )
    axes
}

test_that("BS_PCAPAM50 returns all sample-indexed outputs in input order", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR")
    se <- OSLO2EMIT0obj$data_input$se_NC
    ids <- colnames(se)
    pcapam50_order_ids(ids, as.character(SummarizedExperiment::colData(se)$PatientID))
    genes <- rownames(se)
    pcapam50_order_ids(genes, genes)
    cd <- SummarizedExperiment::colData(se)
    pcapam50_order_assert(
        is.numeric(cd$TSIZE) && !anyNA(cd$TSIZE) && all(cd$TSIZE %in% 0:1) &&
            is.numeric(cd$NODE) && !anyNA(cd$NODE) && all(is.finite(cd$NODE)) &&
            all(cd$NODE >= 0) && all(cd$NODE == floor(cd$NODE)),
        "Packaged clinical covariates must be valid"
    )
    orders <- list(original = ids, reversed = rev(ids),
                   shuffled = withr::with_seed(20261022, sample(ids)))
    for (order in orders) {
        pcapam50_order_ids(order, ids)
        input <- se[, order, drop = FALSE]
        expect_identical(as.character(SummarizedExperiment::colData(input)$PatientID), order)
        for (subtype in c(FALSE, TRUE)) {
            ## Clinical and nonclinical calls are exercised for each order.
            res <- BS_PCAPAM50(input, Subtype = subtype,
                               hasClinical = subtype, seed = 118)
            axes <- pcapam50_order_structure(res, order, genes, subtype, subtype)
            for (axis in axes) expect_identical(axis, order)
        }
    }
})

test_that("output ordering retains samples with unknown ER", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR")
    se <- OSLO2EMIT0obj$data_input$se_NC
    ids <- colnames(se)
    pcapam50_order_ids(ids, as.character(SummarizedExperiment::colData(se)$PatientID))
    genes <- rownames(se)
    pcapam50_order_ids(genes, genes)
    SummarizedExperiment::colData(se)$ER[c(2L, 50L, 77L)] <- NA_character_
    expect_message(
        res <- BS_PCAPAM50(se, Subtype = TRUE, hasClinical = FALSE, seed = 118),
        "3 sample(s) without an evaluable ER/IHC status are excluded",
        fixed = TRUE
    )
    axes <- pcapam50_order_structure(res, ids, genes, TRUE, FALSE)
    for (axis in axes) expect_identical(axis, ids)
})

test_that("BS_PCAPAM50 preserves packaged OSLO outputs except original ROR correction", {
    data("OSLO2EMIT0obj", package = "BreastSubtypeR")
    frozen <- OSLO2EMIT0obj$res$results$PCAPAM50
    pcapam50_order_assert(!is.null(frozen), "The required packaged PCAPAM50 fixture is missing")
    se <- OSLO2EMIT0obj$data_input$se_NC
    ids <- colnames(se)
    genes <- rownames(se)
    pcapam50_order_ids(ids, as.character(SummarizedExperiment::colData(se)$PatientID))
    pcapam50_order_ids(genes, genes)
    fresh <- BS_PCAPAM50(
        se_obj = se,
        Subtype = FALSE,
        hasClinical = FALSE,
        seed = 118
    )
    fresh_axes <- pcapam50_order_structure(fresh, ids, genes, FALSE, FALSE)
    pcapam50_order_structure(frozen, ids, genes, FALSE, FALSE)
    for (axis in fresh_axes) expect_identical(axis, ids)
    pcapam50_order_assert(
        identical(colnames(fresh$score.ROR), colnames(frozen$score.ROR)) &&
            identical(colnames(fresh$outList$distances), colnames(frozen$outList$distances)) &&
            identical(colnames(fresh$outList$dist.RORSubtype), colnames(frozen$outList$dist.RORSubtype)) &&
            identical(dimnames(fresh$outList$centroids), dimnames(frozen$outList$centroids)),
        "Historical and fresh result dimensions and names must agree"
    )

    ## calls, aligned by PatientID
    bi <- match(ids, frozen$BS.all$PatientID)
    expect_false(anyNA(bi))
    expect_identical(fresh$BS.all$BS, frozen$BS.all$BS[bi])

    ## ROR table, aligned by PatientID
    si <- match(ids, as.character(frozen$score.ROR$patientID))
    expect_false(anyNA(si))
    expect_identical(colnames(fresh$score.ROR), colnames(frozen$score.ROR))
    ## Keep the retained fixture; #198 intentionally changes only these ROR
    ## columns. Assert their replacement independently instead of regenerating it.
    ror198_expect_reference(fresh)
    for (col in setdiff(colnames(fresh$score.ROR), ror198_changed)) {
        fresh_col <- fresh$score.ROR[[col]]
        frozen_col <- frozen$score.ROR[[col]][si]
        if (is.numeric(fresh_col)) {
            expect_equal(unname(fresh_col), unname(frozen_col), label = col)
        } else {
            expect_identical(
                as.character(fresh_col), as.character(frozen_col),
                label = col
            )
        }
    }

    ## nearest-centroid internals, aligned by sample name
    expect_identical(
        fresh$outList$predictions, frozen$outList$predictions[ids]
    )
    expect_equal(
        fresh$outList$distances, frozen$outList$distances[ids, , drop = FALSE]
    )
    expect_equal(
        fresh$outList$dist.RORSubtype,
        frozen$outList$dist.RORSubtype[ids, , drop = FALSE]
    )
    expect_equal(
        fresh$outList$testData, frozen$outList$testData[, ids, drop = FALSE]
    )
    expect_equal(fresh$outList$centroids, frozen$outList$centroids)

    ## medians: the frozen table predates two reference-platform columns, so
    ## compare the shared columns, which include the PCAPAM50 cohort medians
    shared <- intersect(colnames(fresh$mdns.fl), colnames(frozen$mdns.fl))
    expect_true("v1PAM.mdns" %in% shared)
    expect_equal(
        as.matrix(fresh$mdns.fl[genes, shared, drop = FALSE]),
        as.matrix(frozen$mdns.fl[genes, shared, drop = FALSE])
    )
})
