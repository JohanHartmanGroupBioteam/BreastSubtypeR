## Invalid method vectors must fail before the input promise is evaluated.
## This also prevents partial classifier execution for a late AUTO entry.
expect_bs_multi_method_error <- function(methods, expected, clean_call = FALSE) {
    input_used <- FALSE
    err <- tryCatch(
        BS_Multi(
            data_input = {
                input_used <- TRUE
                stop("data_input should not be evaluated")
            },
            methods = methods
        ),
        error = identity
    )
    expect_false(input_used)
    expect_s3_class(err, "error")
    expect_identical(conditionMessage(err), expected)
    if (clean_call) expect_null(conditionCall(err))
}

test_that("BS_Multi rejects mixed and repeated AUTO before reading input", {
    selections <- list(
        c("AUTO", "AIMS"),
        c("AIMS", "AUTO"),
        c("AIMS", "AUTO", "sspbc"),
        c("AUTO", "AUTO"),
        c("AUTO", "AUTO", "AUTO"),
        c("AUTO", "AUTO", "AIMS")
    )
    for (selected in selections) {
        expect_bs_multi_method_error(
            selected,
            "'AUTO' must be the only entry in 'methods'.",
            clean_call = TRUE
        )
    }
})

test_that("BS_Multi rejects duplicate manual methods with clear names", {
    selections <- list(
        c("AIMS", "AIMS"),
        c("AIMS", "AIMS", "sspbc"),
        c("AIMS", "AIMS", "AIMS"),
        c("sspbc", "AIMS", "sspbc", "AIMS", "AIMS")
    )
    repeated <- c("AIMS", "AIMS", "AIMS", "sspbc, AIMS")
    for (i in seq_along(selections)) {
        expect_bs_multi_method_error(
            selections[[i]],
            paste0("Duplicate method(s) specified: ", repeated[i]),
            clean_call = TRUE
        )
    }
})

test_that("BS_Multi preserves invalid-name validation precedence", {
    selections <- list(
        c("AUTO", "not-a-method"),
        c("AIMS", "AIMS", "not-a-method"),
        c("AIMS", NA_character_)
    )
    invalid <- c("not-a-method", "not-a-method", "NA")
    for (i in seq_along(selections)) {
        expect_bs_multi_method_error(
            selections[[i]],
            paste0("Invalid method(s) specified: ", invalid[i])
        )
    }
})

test_that("BS_Multi preserves the minimum manual-method requirement", {
    for (selected in list(character(), NULL, "AIMS", "sspbc")) {
        expect_bs_multi_method_error(
            selected,
            "Select at least two methods or set method to 'AUTO'."
        )
    }
})

test_that("single AUTO and distinct manual vectors pass method validation", {
    selections <- list(
        "AUTO",
        c("AIMS", "sspbc"),
        c("sspbc", "AIMS"),
        c("ssBC", "ssBC.v2", "AIMS")
    )
    for (selected in selections) {
        expect_error(
            suppressMessages(BS_Multi(
                data_input = stop("method validation passed"),
                methods = selected
            )),
            "method validation passed",
            fixed = TRUE
        )
    }
})

test_that("valid manual selections preserve order and native SSP calls", {
    data_env <- new.env(parent = emptyenv())
    data("OSLO2EMIT0obj", envir = data_env, package = "BreastSubtypeR")
    original <- data_env[["OSLO2EMIT0obj"]]$data_input
    input <- list(
        se_NC = original$se_NC[, c(3, 1, 2)],
        se_SSP = original$se_SSP[, c(3, 1, 2)]
    )
    samples <- colnames(input$se_NC)
    selected <- c("sspbc", "AIMS")
    out <- suppressMessages(BS_Multi(
        input, methods = selected, Subtype = TRUE
    ))
    reversed <- suppressMessages(BS_Multi(
        input, methods = rev(selected), Subtype = TRUE
    ))

    expect_identical(names(out$results), selected)
    expect_identical(names(reversed$results), rev(selected))
    expect_identical(names(out$res_subtypes), c(selected, "entropy"))
    expect_identical(names(out$res_subtypes.Subtype), c(selected, "entropy"))
    expect_identical(rownames(out$res_subtypes), samples)
    expect_identical(rownames(out$res_subtypes.Subtype), samples)
    expect_identical(
        reversed$res_subtypes[, names(out$res_subtypes), drop = FALSE],
        out$res_subtypes
    )
    expect_identical(
        reversed$res_subtypes.Subtype[
            , names(out$res_subtypes.Subtype), drop = FALSE
        ],
        out$res_subtypes.Subtype
    )
    expect_identical(reversed$results[selected], out$results)

    aims <- BS_AIMS(input$se_SSP)$cl
    ssp5 <- BS_sspbc(input$se_SSP, ssp.name = "ssp.pam50")$cl
    ssp4 <- BS_sspbc(input$se_SSP, ssp.name = "ssp.subtype")$cl
    expect_identical(out$results$AIMS$cl, aims)
    expect_identical(out$results$sspbc$cl, ssp5)
    expect_identical(out$res_subtypes$AIMS, unname(aims[samples, 1]))
    expect_identical(out$res_subtypes$sspbc, unname(ssp5[samples, 1]))
    aims4 <- unname(aims[samples, 1])
    aims4[which(aims4 == "Normal")] <- NA_character_
    expect_identical(out$res_subtypes.Subtype$AIMS, aims4)
    expect_identical(out$res_subtypes.Subtype$sspbc, unname(ssp4[samples, 1]))
})
