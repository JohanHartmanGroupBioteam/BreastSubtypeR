# Emulate the outer Hmisc display label without requiring Hmisc in the test library.
layer_display_label <- function(x) {
    attr(x, "label") <- "Synthetic size"
    class(x) <- c("labelled", class(x))
    x
}

layer_input <- function() {
    utils::data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = environment())
    d <- lapply(OSLO2EMIT0obj$data_input, function(z) z[, 24:1])
    z <- d$se_NC
    z$PatientID <- factor(colnames(z), levels = rev(colnames(z)))
    z$NODE <- rep(c(0, 2, 7), 8)
    z$ER <- rep(c("ER+", "ER-"), 12)
    z$HER2 <- rep(c("HER2-", "HER2-", "HER2+", "HER2+"), 6)
    d$se_NC <- z
    d
}

test_that("clinical routes stop before unverified declarations lose their class or attributes", {
    skip_if_not_installed("haven")
    d <- layer_input()
    v <- rep(c(0, 1), 12)
    labelled <- layer_display_label(haven::labelled(v))
    spss <- layer_display_label(haven::labelled_spss(v, na_values = 1))
    cases <- list()
    add <- function(x, key, declaration) {
        attr(x, key) <- declaration
        cases[[length(cases) + 1L]] <<- list(x = x, key = key)
    }
    for (code in 0:1) add(labelled, "na_values", code)
    add(labelled, "na_range", c(0, 1))
    add(labelled, "special.miss", list(obs = which(v == 1)))
    add(spss, "special.miss", 0)
    add(I(spss), "missing.codes", 0)
    add(v, "na_values", 1)
    add(factor(v), "missing.codes", 1L)
    add(S4Vectors::Rle(v), "special.miss", 1)
    add(structure(v, class = "unknown_measure"), "na_range", c(0, 1))
    for (case in cases) {
        d$se_NC$TSIZE <- case$x
        before <- SummarizedExperiment::colData(d$se_NC)
        expect_error(suppressMessages(BS_parker(d$se_NC, hasClinical = TRUE)), case$key)
        expect_error(suppressMessages(BS_Multi(d, methods = c("parker.original", "ssBC"),
                                               hasClinical = TRUE)), case$key)
        expect_identical(SummarizedExperiment::colData(d$se_NC), before)
    }
    d$se_NC$TSIZE <- v
    genomic <- function(d) suppressWarnings(suppressMessages(BS_parker(d$se_NC, hasClinical = FALSE)))
    ssp <- function(d) suppressWarnings(suppressMessages(BS_Multi(d, methods = c("AIMS", "sspbc"),
                                                                 hasClinical = TRUE)))
    ref <- genomic(d)
    ref_ssp <- ssp(d)
    d$se_NC$TSIZE <- cases[[5L]]$x
    expect_identical(genomic(d), ref)
    expect_identical(ssp(d), ref_ssp)
})

test_that("SPSS declarations require compatible storage, range structure and class order", {
    skip_if_not_installed("haven")
    decode <- BreastSubtypeR:::.decode_tsize
    spss <- haven::labelled_spss(c(0, 1, 25))
    for (v in list("1", NA_real_, factor(1), matrix(1))) {
        x <- spss
        attr(x, "na_values") <- v
        expect_error(decode(x), "na_values")
    }
    for (v in list(c("0", "1"), c(0, NaN), 1, 0:2, c(1, 0), c(1, 1), matrix(0:1))) {
        x <- spss
        attr(x, "na_range") <- v
        expect_error(decode(x), "na_range")
    }
    for (classes in list("haven_labelled_spss", c("opaque", class(spss)),
                         c("labelled", "labelled", class(spss)), c(class(spss)[1:3], "integer"))) {
        x <- spss
        class(x) <- classes
        expect_error(decode(x), "class stack")
    }
    x <- spss
    attr(x, "na_range") <- c(0, 1)
    attr(x, "na_values") <- c(2, 3)
    expect_identical(decode(x)$value, c(NA_real_, NA_real_, NA_real_))
    expect_identical(decode(x)$invalid, c(FALSE, FALSE, TRUE))
    for (values in list(numeric(), 1:4, c(0, 1, 2, 3))) {
        valid <- haven::labelled_spss(c(0, 1, 25), na_values = values)
        expect_identical(decode(valid), decode(as.numeric(haven::zap_missing(valid))))
    }
    for (classes in list(rev(class(spss)), c(class(spss), "AsIs"))) {
        x <- haven::labelled_spss(c(0, 1, 25), na_values = 1)
        class(x) <- classes
        expect_identical(decode(x), decode(c(0, NA_real_, 25)))
    }
})

test_that("declaration wrappers retain nested metadata and storage checks", {
    skip_if_not_installed("haven")
    decode <- BreastSubtypeR:::.decode_tsize
    x <- haven::labelled_spss(c(0, 1, NA, 25))
    for (key in c("na_values", "na_range")) {
        codes <- if (key == "na_values") 1 else c(0.5, 1.5)
        for (wrap in list(I, layer_display_label, haven::labelled)) {
            declaration <- wrap(codes)
            attr(x, key) <- declaration
            expect_identical(decode(x), decode(c(0, NA_real_, NA_real_, 25)))
            attr(declaration, "missing.codes") <- 1
            attr(x, key) <- declaration
            expect_error(decode(x), "missing.codes")
        }
        attr(x, key) <- structure(codes, class = c("labelled", "integer"))
        expect_error(decode(x), "class or storage")
        attr(x, key) <- NULL
    }
    r <- S4Vectors::Rle(c(0, 1, 0))
    S4Vectors::metadata(r) <- list(na.values = 1)
    expect_error(decode(r), "na.values")
    S4Vectors::metadata(r) <- list()
    S4Vectors::mcols(r)$missing.codes <- 1L
    expect_error(decode(r), "missing.codes")
    for (key in c("na_values2", "na_range2", "na.values", "na.range")) {
        y <- x
        attr(y, key) <- 1
        expect_error(decode(y), key, fixed = TRUE)
    }
})

test_that("layered scores and invalid patient IDs agree with a separate missingness oracle", {
    skip_if_not_installed("haven")
    withr::local_seed(590)
    d <- layer_input()
    v <- rep(c(0, 1, 25, NA_real_, NaN, 2), 4)
    run <- function(x) {
        z <- d$se_NC
        z$TSIZE <- x
        before <- serialize(SummarizedExperiment::colData(z), NULL)
        seed <- .Random.seed
        warnings <- list()
        result <- suppressMessages(withCallingHandlers(BS_parker(z, hasClinical = TRUE),
            warning = function(w) {
                warnings[[length(warnings) + 1L]] <<- w
                invokeRestart("muffleWarning")
            }))
        expect_identical(serialize(SummarizedExperiment::colData(z), NULL), before)
        expect_identical(.Random.seed, seed)
        list(output = result, warnings = warnings)
    }
    compare <- function(x, oracle) {
        invalid <- which(!is.na(oracle) & !oracle %in% 0:1)
        ref <- run(replace(oracle, invalid, NA_real_))
        actual <- run(x)
        size <- vapply(actual$warnings, inherits, TRUE, what = "bsr_invalid_tsize_warning")
        expect_identical(actual$output, ref$output)
        expect_identical(actual$warnings[!size], ref$warnings)
        expect_equal(sum(size), as.integer(length(invalid) > 0L))
        if (length(invalid))
            expect_identical(sort(actual$warnings[size][[1L]]$patient_ids), sort(colnames(d$se_NC)[invalid]))
    }
    compare(layer_display_label(haven::labelled(v, c(small = 0, large = 1))), v)
    for (s in list(haven::labelled_spss(v, na_values = c(1, 25)),
                   haven::labelled_spss(v, na_values = 0),
                   haven::labelled_spss(v, na_range = c(0, 1), na_values = 25),
                   haven::labelled_spss(as.character(v), na_range = c("0", "1")))) {
        oracle <- unclass(haven::zap_missing(s))
        if (is.character(oracle)) oracle[oracle == "NaN" & !is.na(oracle)] <- "2"
        for (x in list(s, layer_display_label(s), I(layer_display_label(s)), layer_display_label(I(s))))
            compare(x, oracle)
    }
    compare(I(structure(v, class = "unknown_measure")), ifelse(is.na(v), NA_real_, 2))
    compare(S4Vectors::Rle(v), v)
    if (requireNamespace("bit64", quietly = TRUE))
        compare(I(bit64::as.integer64(v)), v)
})
