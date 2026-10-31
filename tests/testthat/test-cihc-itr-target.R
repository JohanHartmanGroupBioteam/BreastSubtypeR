## cIHC.itr ER subsampling: target-aware default (er_subsample = "target")
## and the behaviour of BreastSubtypeR <= 1.5.5 (er_subsample = "larger").
##
## Reference outputs in fixtures/cihc_itr_1.5.5.rds were produced by
## BreastSubtypeR 1.5.5 at commit 0dcc0c8a846290fcca3f266bd79248bff84fb58b,
## before this change, with fixtures/make-cihc-itr-1.5.5.R.

plan <- BreastSubtypeR:::.cihc_itr_plan
used <- function(p, n_pos, n_neg) {
    switch(p$group,
        "ER+" = c(p$size, n_neg),
        "ER-" = c(n_pos, p$size),
        none = c(n_pos, n_neg)
    )
}

test_that("the target plan subsamples the over-represented ER group", {
    expect_identical(plan(54, 64), list(group = "none", size = 0L))
    expect_identical(plan(40, 100), list(group = "ER-", size = 48L))
    expect_identical(plan(47, 53), list(group = "ER+", size = 45L))
    expect_identical(plan(100, 40), list(group = "ER+", size = 34L))
    expect_identical(plan(50, 50), list(group = "ER+", size = 43L))
    expect_identical(plan(1, 2), list(group = "none", size = 0L))
    expect_identical(plan(1, 5), list(group = "ER-", size = 2L))
    expect_identical(plan(2, 1), list(group = "ER+", size = 1L))
    expect_identical(plan(47, 53, ratio = 1), list(group = "ER-", size = 47L))
    expect_identical(plan(40, 100, ratio = 2), list(group = "ER-", size = 20L))
})

test_that("the larger plan is the rule of BreastSubtypeR <= 1.5.5", {
    expect_identical(plan(54, 64, er_subsample = "larger"),
        list(group = "ER-", size = 46L))
    expect_identical(plan(47, 53, er_subsample = "larger"),
        list(group = "ER-", size = 40L))
    expect_identical(plan(40, 100, er_subsample = "larger"),
        list(group = "ER-", size = 34L))
    expect_identical(plan(50, 50, er_subsample = "larger"),
        list(group = "ER+", size = 43L))
    expect_error(plan(47, 53, ratio = 2, er_subsample = "larger"),
        "please specify a ratio less than")
})

test_that("both plans select the same subset whenever n_ER+ >= n_ER-", {
    for (n_pos in seq_len(40)) {
        for (n_neg in seq_len(n_pos)) {
            expect_identical(
                used(plan(n_pos, n_neg), n_pos, n_neg),
                used(plan(n_pos, n_neg, er_subsample = "larger"), n_pos, n_neg)
            )
        }
    }
})

test_that("target subsets keep the limiting group and match the target", {
    for (r in c(54 / 64, 0.4, 1, 2.5)) {
        for (n_pos in c(1, 2, 5, 17, 54, 100)) {
            for (n_neg in c(1, 3, 18, 64, 100)) {
                p <- plan(n_pos, n_neg, ratio = r)
                u <- used(p, n_pos, n_neg)
                expect_true(all(u >= 1L) && u[1] <= n_pos && u[2] <= n_neg)
                if (p$group == "ER+") {
                    expect_identical(u[2], n_neg)
                    expect_lt(abs(u[1] - r * n_neg), 1)
                }
                if (p$group == "ER-") {
                    expect_identical(u[1], n_pos)
                    expect_lt(abs(u[2] - n_pos / r), 1)
                }
            }
        }
    }
})

test_that("ratio and group sizes are validated", {
    for (bad in list(54:64, 0, -1, NA_real_, Inf, NaN, "54/64", TRUE, numeric(0))) {
        expect_error(plan(84, 18, ratio = bad), "single positive number")
        expect_error(plan(84, 18, ratio = bad, er_subsample = "larger"),
            "single positive number")
    }
    expect_error(plan(0, 10), "requires both ER\\+ and ER- samples")
    expect_error(plan(10, 0, er_subsample = "larger"),
        "requires both ER\\+ and ER- samples")
})

test_that("extreme finite ratios give valid integer sample sizes", {
    tiny <- c(.Machine$double.xmin, 5e-324, 1e-300)
    huge <- c(.Machine$double.xmax, 1e300)
    for (r in c(tiny, huge)) {
        for (n in list(c(1, 1), c(2, 1), c(1, 2), c(1e6, 3), c(3, 1e6))) {
            p <- expect_silent(plan(n[1], n[2], ratio = r))
            expect_true(is.integer(p$size) && !is.na(p$size))
            u <- used(p, n[1], n[2])
            expect_true(all(u >= 1L) && u[1] <= n[1] && u[2] <= n[2])
        }
    }
    expect_identical(plan(100, 40, ratio = 1e-300), list(group = "ER+", size = 1L))
    expect_identical(plan(40, 100, ratio = 1e300), list(group = "ER-", size = 1L))
})

ref <- readRDS(test_path("fixtures", "cihc_itr_1.5.5.rds"))
data("OSLO2EMIT0obj", package = "BreastSubtypeR")
se <- OSLO2EMIT0obj$data_input$se_NC
parts <- c("BS.all", "score.ROR", "BS.itr.keep", "BS.itr.keep.Subtype")

run_case <- function(nm, er_subsample) {
    r <- ref$cihc_itr[[nm]]
    args <- list(se[, r$ids], Subtype = TRUE, hasClinical = TRUE,
        er_subsample = er_subsample)
    if (grepl("seed7", nm)) args$seed <- 7
    if (grepl("ratio1", nm)) args$ratio <- 1
    if (grepl("ratio2", nm)) args$ratio <- 2
    res <- withr::with_options(list(BreastSubtypeR.verbose = FALSE),
        expect_experimental_ror_pc(do.call(BS_cIHC.itr, args)))
    out <- res$outList
    confidence_reference <- vapply(seq_along(out$predictions), function(j) {
        suppressWarnings(round(1 - cor.test(out$testData[, j],
            out$centroids[rownames(out$testData), out$predictions[j]],
            method = "spearman")$p.value, 2))
    }, numeric(1))
    c(res[parts], list(testData = res$outList$testData,
        distances = res$outList$distances,
        confidence_reference = confidence_reference))
}
expect_as_1.5.5 <- function(new, old) {
    # calls, per-iteration calls and sample order exactly; numbers to tolerance
    expect_identical(new$BS.all, old$BS.all)
    expect_identical(new$BS.itr.keep, old$BS.itr.keep)
    expect_identical(new$BS.itr.keep.Subtype, old$BS.itr.keep.Subtype)
    # Retain the historical fixture, whose Confidence contains false zeros.
    # Keep original ROR-P/PC in the historical fixture as well; #198 changes
    # only these four columns. Test their new values on the final averaged
    # matrix while retaining every iteration/consensus/averaging assertion.
    other_scores <- setdiff(colnames(old$score.ROR), c("Confidence", ror198_changed))
    expect_identical(colnames(new$score.ROR), colnames(old$score.ROR))
    expect_equal(new$score.ROR[, other_scores, drop = FALSE],
        old$score.ROR[, other_scores, drop = FALSE], tolerance = 1e-12)
    expect_identical(new$score.ROR$Confidence, new$confidence_reference)
    clinical <- as.data.frame(SummarizedExperiment::colData(se))
    ror198_expect_reference(list(score.ROR = new$score.ROR,
        outList = list(testData = new$testData, distances = new$distances)), clinical)
    expect_equal(new$testData, old$testData, tolerance = 1e-12)
    expect_equal(new$distances, old$distances, tolerance = 1e-12)
}

test_that("er_subsample = 'larger' reproduces the retained 1.5.5 outputs", {
    for (nm in names(ref$cihc_itr)) {
        expect_as_1.5.5(run_case(nm, "larger"), ref$cihc_itr[[nm]])
    }
})

test_that("the default reproduces 1.5.5 where its selection is unchanged", {
    same <- c("full_84_18", "equal_18_18", "erpos_45_18", "erpos_45_18_seed7",
        "erpos_45_18_ratio1", "erpos_45_18_ratio2", "erneg_7_18_ratio1")
    for (nm in same) {
        expect_as_1.5.5(run_case(nm, "target"), ref$cihc_itr[[nm]])
    }
})

test_that("the default changes cIHC.itr only where ER- outnumbers ER+", {
    for (nm in c("erneg_7_18", "target_15_18", "near_16_18", "erneg_7_18_seed7")) {
        new <- run_case(nm, "target")
        old <- ref$cihc_itr[[nm]]
        expect_identical(new$BS.all$PatientID, old$BS.all$PatientID)
        expect_identical(rownames(new$score.ROR), rownames(old$score.ROR))
        expect_false(identical(new$BS.itr.keep, old$BS.itr.keep))
    }
    # at the 54:64 target (15 ER+ : 18 ER-), every iteration uses all samples
    new <- run_case("target_15_18", "target")
    expect_true(all(apply(new$BS.itr.keep, 1, function(x) length(unique(x)) == 1L)))
})

test_that("the default is reproducible for a given seed", {
    ids <- ref$cihc_itr$erneg_7_18$ids
    a <- withr::with_options(list(BreastSubtypeR.verbose = FALSE),
        BS_cIHC.itr(se[, ids], seed = 7))
    b <- withr::with_options(list(BreastSubtypeR.verbose = FALSE),
        BS_cIHC.itr(se[, ids], seed = 7))
    expect_identical(a[c("BS.all", "score.ROR", "BS.itr.keep")],
        b[c("BS.all", "score.ROR", "BS.itr.keep")])
})

test_that("a cohort at the target reports that all samples are used", {
    ids <- ref$cihc_itr$target_15_18$ids
    expect_message(
        withr::with_options(list(BreastSubtypeR.verbose = TRUE),
            BS_cIHC.itr(se[, ids])),
        "already matches the ER\\+ : ER- target"
    )
})

test_that("an absent ER group and an invalid ratio are errors", {
    er <- SummarizedExperiment::colData(se)$ER
    expect_error(BS_cIHC.itr(se[, er == "ER+"]), "requires both ER\\+ and ER-")
    expect_error(BS_cIHC.itr(se[, er == "ER-"], er_subsample = "larger"),
        "requires both ER\\+ and ER-")
    expect_error(BS_cIHC.itr(se, ratio = 54:64), "single positive number")
    expect_error(BS_cIHC.itr(se, er_subsample = "smaller"), "should be one of")
})

test_that("BS_Multi() changes only cIHC.itr and reproduces 1.5.5 with 'larger'", {
    methods_all <- c("parker.original", "genefu.scale", "genefu.robust", "ssBC",
        "ssBC.v2", "cIHC", "cIHC.itr", "PCAPAM50", "AIMS", "sspbc")
    din <- OSLO2EMIT0obj$data_input
    quiet <- function(expr) {
        withr::with_options(list(BreastSubtypeR.verbose = FALSE),
            suppressWarnings(suppressMessages(expr)))
    }
    keep <- c("res_subtypes", "res_subtypes.Subtype")

    full <- quiet(BS_Multi(din, methods = methods_all, Subtype = TRUE))
    expect_identical(full[keep], ref$multi$full_manual)
    auto <- quiet(BS_Multi(din, methods = "AUTO", Subtype = TRUE))
    expect_identical(auto[keep], ref$multi$full_auto)

    ids <- ref$multi$erneg_7_18_manual$ids
    sub <- list(se_NC = din$se_NC[, ids], se_SSP = din$se_SSP[, ids])
    old <- ref$multi$erneg_7_18_manual
    legacy <- quiet(BS_Multi(sub, methods = methods_all, Subtype = TRUE,
        er_subsample = "larger"))
    expect_identical(legacy[keep], old[keep])
    new <- quiet(BS_Multi(sub, methods = methods_all, Subtype = TRUE))
    others <- setdiff(methods_all, "cIHC.itr")
    for (k in keep) {
        expect_identical(new[[k]][, others], old[[k]][, others])
    }
})
