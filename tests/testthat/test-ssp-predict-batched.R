ssp_expect_all_cases <- function() {
    cases <- ssp_predict_cases()
    for (name in names(cases)) {
        for (subtype in c(FALSE, TRUE)) {
            args <- cases[[name]]$args
            args$Subtype <- subtype
            args$info <- paste(name, "Subtype =", subtype)
            do.call(expect_ssp_scalar_equivalence, args)
        }
    }
}

test_that("sspPredict matches the independent scalar predictor on complete and irregular panels", {
    ssp_expect_all_cases()
})

test_that("forcing all three scalar paths preserves complete values and conditions", {
    local_mocked_bindings(.can_batch_spearman = function(centres, samples) FALSE,
        .package = "BreastSubtypeR")
    ssp_expect_all_cases()
})

test_that("escalated correlation warnings stop at the same call as in the scalar predictor", {
    cases <- ssp_predict_cases()
    withr::local_options(list(warn = 2L))
    for (name in c("all_constant_sample", "ror_constant_sample",
                   "standardized_overflow", "late_alignment_error")) {
        for (subtype in c(FALSE, TRUE)) {
            args <- cases[[name]]$args
            args$Subtype <- subtype
            args$info <- paste(name, "warn = 2, Subtype =", subtype)
            outcome <- do.call(expect_ssp_scalar_equivalence, args)
            expect_true(outcome$failed, info = args$info)
            expect_true(any(vapply(outcome$trace, function(item)
                identical(item$type, "warning"), logical(1))), info = args$info)
        }
    }
})
