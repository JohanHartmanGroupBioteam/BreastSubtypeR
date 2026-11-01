# Historical example input is pinned independently of public data; see fixtures/README.md.
## The packaged OSLO2-EMIT0 example ships a frozen AUTO result
## (OSLO2EMIT0obj$res). A fresh AUTO run must reproduce it. The cohort has 84
## ER+ and 18 ER- tumors, and the ER- count equals the ER- minimum of 18, so
## this also guards the inclusive (>=) sample subsetting for ssBC/ssBC.v2: with
## a strict comparison the 18 ER- tumors are excluded from ssBC and come back
## NA, which does not match the stored result.

test_that("a fresh AUTO run reproduces the packaged frozen result", {
    load_historical_example("OSLO2EMIT0obj")

    old_options <- options(BreastSubtypeR.verbose = FALSE)
    on.exit(options(old_options), add = TRUE)
    old_rng_kind <- RNGkind()
    on.exit(do.call(RNGkind, as.list(old_rng_kind)), add = TRUE)
    RNGkind(
        kind = "Mersenne-Twister",
        normal.kind = "Inversion",
        sample.kind = "Rejection"
    )
    set.seed(20260831)

    fresh <- suppressMessages(BS_Multi(
        data_input = OSLO2EMIT0obj$data_input,
        methods = "AUTO",
        Subtype = FALSE,
        hasClinical = FALSE
    ))
    frozen <- OSLO2EMIT0obj$res
    methods <- setdiff(colnames(frozen$res_subtypes), "entropy")

    expect_identical(colnames(fresh$res_subtypes), colnames(frozen$res_subtypes))
    expect_identical(rownames(fresh$res_subtypes), rownames(frozen$res_subtypes))
    ## Discrete calls must match exactly; entropy is a floating-point quantity
    ## and is compared with the default testthat tolerance, because
    ## transcendental functions are not bit-identical across build platforms.
    expect_identical(
        fresh$res_subtypes[, methods, drop = FALSE],
        frozen$res_subtypes[, methods, drop = FALSE]
    )
    expect_equal(fresh$res_subtypes$entropy, frozen$res_subtypes$entropy)
    ## every ER- tumor receives an ssBC call (the ER- count equals the minimum)
    pheno <- as.data.frame(
        SummarizedExperiment::colData(OSLO2EMIT0obj$data_input$se_NC)
    )
    expect_identical(sum(pheno$ER == "ER-"), 18L)
    expect_false(anyNA(fresh$res_subtypes$ssBC))
})
