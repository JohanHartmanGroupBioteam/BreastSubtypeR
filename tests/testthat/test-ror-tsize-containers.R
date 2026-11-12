container_input <- function() {
    e <- new.env(parent = emptyenv())
    utils::data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = e)
    lapply(e$OSLO2EMIT0obj$data_input, function(z) z[, 40:1])
}

container_run <- function(d, route, four = FALSE) {
    warnings <- list()
    result <- suppressMessages(withCallingHandlers(switch(route,
        parker = BS_parker(d$se_NC, hasClinical = TRUE, Subtype = four),
        cIHC = BS_cIHC(d$se_NC, hasClinical = TRUE, Subtype = four),
        iterative = BS_cIHC.itr(d$se_NC, iteration = 3L, hasClinical = TRUE, Subtype = four),
        PCA = BS_PCAPAM50(d$se_NC, hasClinical = TRUE, Subtype = four),
        ER = BS_ssBC(d$se_NC, "ER", hasClinical = TRUE, Subtype = four),
        ER.v2 = BS_ssBC(d$se_NC, "ER.v2", hasClinical = TRUE, Subtype = four),
        multi = BS_Multi(d, methods = c("ssBC", "ssBC.v2"), hasClinical = TRUE, Subtype = four)),
        warning = function(w) {
            warnings[[length(warnings) + 1L]] <<- list(class = class(w), ids = w$patient_ids,
                                                       message = conditionMessage(w))
            invokeRestart("muffleWarning")
        }))
    list(result = result, warnings = warnings)
}

test_that("Rle and layered annotations retain stored codes and missingness", {
    decode <- BreastSubtypeR:::.decode_tsize
    vals <- c(0, NA, 1, 2, 25)
    for (x in list(vals, as.integer(vals), as.character(vals), factor(vals, levels = c(25, 2, 1, 0)),
                   c(FALSE, NA, TRUE))) {
        expect_identical(decode(S4Vectors::Rle(x)), decode(x))
    }
    for (x in list(c(0, 1), 0:1, c(0, NA, 1), c(TRUE, FALSE), as.raw(0:1), c("0", "1"))) {
        expect_identical(decode(I(x)), decode(x))
    }
    x <- structure(vals, class = "unknown_measure")
    expect_identical(decode(I(x)), list(value = rep(NA_real_, 5), invalid = !is.na(vals)))
    expect_error(decode(S4Vectors::DataFrame(TSIZE = 0:1)), "Unsupported S4 TSIZE")
    expect_error(decode(S4Vectors::Rle(as.raw(0:1))), "Rle run values")

    for (v in list(as.integer(vals), c(vals, NaN))) {
        x <- structure(v, label = "Size", labels = c("1" = 0, "0" = 1),
                       class = c("labelled", class(v)))
        expect_identical(decode(x), decode(v))
        for (key in c("na_values", "na_range", "special.miss", "missing.codes")) {
            unsafe <- x
            attr(unsafe, key) <- 1
            expect_error(decode(unsafe), paste0(key, ".*missing codes to NA"))
        }
        class(x) <- c("labelled", "unknown_measure", class(v))
        expect_identical(decode(x)$invalid, !is.na(v))
    }
})

test_that("public NC and multi routes score supported containers by PatientID", {
    withr::local_seed(851)
    d <- container_input()
    z <- d$se_NC
    n <- ncol(z)
    z$ER <- rep(c("ER+", "ER-"), length.out = n)
    z$HER2 <- rep(c("HER2+", "HER2+", "HER2-", "HER2-"), length.out = n)
    z$NODE <- rep(c(0, 1, 4, 7), length.out = n)
    z$PatientID <- factor(colnames(z), levels = rev(colnames(z)))
    for (mixed in c(FALSE, TRUE)) {
        codes <- rep(0:1, length.out = n)
        if (mixed) codes[c(3, 7, 18)] <- c(NA_integer_, 2L, 25L)
        invalid <- which(!is.na(codes) & !codes %in% 0:1)
        wrapped <- list(S4Vectors::Rle(codes), I(as.double(codes)), I(codes),
                        structure(codes, class = c("labelled", "integer"),
                                  label = "Size", labels = c("small" = 1L, "large" = 0L)))
        z$TSIZE <- replace(codes, invalid, NA_integer_)
        d$se_NC <- z
        for (route in c("parker", "cIHC", "iterative", "PCA", "ER", "ER.v2", "multi")) {
            reference <- container_run(d, route)
            for (v in wrapped) {
                subject <- d
                subject$se_NC$TSIZE <- v
                before <- serialize(SummarizedExperiment::colData(subject$se_NC), NULL)
                seed <- .Random.seed
                actual <- container_run(subject, route)
                is_bad <- vapply(actual$warnings, function(w) "bsr_invalid_tsize_warning" %in% w$class, TRUE)
                expect_identical(actual$result, reference$result)
                expect_identical(actual$warnings[!is_bad], reference$warnings)
                expect_equal(sum(is_bad), as.integer(mixed))
                if (mixed) expect_setequal(actual$warnings[is_bad][[1L]]$ids, colnames(z)[invalid])
                expect_identical(serialize(SummarizedExperiment::colData(subject$se_NC), NULL), before)
                expect_identical(.Random.seed, seed)
            }
        }
    }
})

test_that("unknown storage is masked during cIHC, iterative cIHC and both PCAPAM50 stages", {
    withr::local_seed(822)
    d <- container_input()
    n <- ncol(d$se_NC)
    d$se_NC$NODE <- rep(c(0, 2, 4, 7), length.out = n)
    codes <- rep(0:1, length.out = n)
    codes[3L] <- NA_integer_
    stages <- list()
    original <- BreastSubtypeR:::RORgroup
    testthat::local_mocked_bindings(RORgroup = function(...) {
        out <- original(...)
        stages[[length(stages) + 1L]] <<- out
        out
    }, .package = "BreastSubtypeR")
    for (route in c("cIHC", "iterative", "PCA")) {
        d$se_NC$TSIZE <- rep(NA_real_, n)
        stages <- list()
        reference <- container_run(d, route, four = TRUE)
        expected_stages <- stages
        if (route == "PCA") expect_length(expected_stages, 2L)
        for (v in list(structure(codes, class = "unknown_measure"),
                       I(structure(codes, class = "unknown_measure")))) {
            stages <- list()
            d$se_NC$TSIZE <- v
            actual <- container_run(d, route, four = TRUE)
            expect_identical(stages, expected_stages)
            expect_identical(actual$result, reference$result)
            is_bad <- vapply(actual$warnings, function(w) "bsr_invalid_tsize_warning" %in% w$class, TRUE)
            expect_equal(sum(is_bad), 1)
            if (sum(is_bad) == 1L) {
                expect_setequal(actual$warnings[is_bad][[1L]]$ids, colnames(d$se_NC)[!is.na(codes)])
                expect_match(actual$warnings[is_bad][[1L]]$message, "unsupported column class")
            }
            expect_identical(actual$warnings[!is_bad], reference$warnings)
        }
    }
})

test_that("ssBC refuses excluded-row empty IDs while RORgroup keeps its alignment policy", {
    withr::local_seed(834)
    d <- container_input()
    n <- ncol(d$se_NC)
    d$se_NC$PatientID <- as.character(d$se_NC$PatientID)
    d$se_NC$TSIZE <- rep(0:1, length.out = n)
    d$se_NC$ER <- c(rep("ER+", n - 1L), "ER-")
    d$se_NC$HER2 <- rep("HER2-", n)
    d$se_NC$ER[5L] <- NA_character_
    d$se_NC$TSIZE[c(5L, n)] <- 25
    altered <- d
    altered$se_NC$PatientID[5L] <- ""
    for (route in c("ER", "ER.v2", "multi")) for (four in c(FALSE, TRUE)) {
        before <- serialize(d, NULL); rng <- .Random.seed
        expected <- container_baseline_run(d, route, four)
        actual <- container_run(d, route, four)
        expect_identical(actual, expected)
        expect_identical(serialize(d, NULL), before)
        expect_identical(.Random.seed, rng)
        expect_error(container_run(altered, route, four), class = "bsr_patientid_mismatch_error")
    }

    out <- suppressWarnings(BS_parker(d$se_NC, hasClinical = FALSE))$outList
    ids <- names(out$predictions)
    pheno <- data.frame(PatientID = ids, NODE = 4, TSIZE = 0)
    expected <- BreastSubtypeR:::RORgroup(out, pheno, hasClinical = TRUE)
    extra <- rbind(pheno[n:1, ], data.frame(PatientID = c("", NA, "unused", "unused"), NODE = 4, TSIZE = 25))
    expect_identical(BreastSubtypeR:::RORgroup(out, extra, hasClinical = TRUE), expected)
    for (malformed in list(rbind(pheno, pheno[1L, ]), pheno[-1L, ],
                           transform(pheno, PatientID = replace(PatientID, 1L, "")),
                           transform(pheno, PatientID = replace(PatientID, 1L, NA_character_)))) {
        expect_error(BreastSubtypeR:::RORgroup(out, malformed, hasClinical = TRUE), "exactly one clinical record")
    }
})

test_that("labelled missing declarations stop clinical scoring and leave genomic calls intact", {
    withr::local_seed(827)
    d <- container_input()
    x <- rep(0:1, length.out = ncol(d$se_NC))
    d$se_NC$TSIZE <- x
    genomic <- function(d) suppressMessages(BS_parker(d$se_NC, hasClinical = FALSE))
    reference <- genomic(d)
    for (key in c("na_values", "na_range", "special.miss", "missing.codes")) {
        v <- structure(x, class = c("labelled", "integer"), label = "Size")
        attr(v, key) <- 1L
        d$se_NC$TSIZE <- v
        original <- serialize(SummarizedExperiment::colData(d$se_NC), NULL)
        seed <- .Random.seed
        for (route in c("parker", "multi"))
            expect_error(container_run(d, route), paste0(key, ".*missing codes to NA"))
        expect_identical(genomic(d), reference)
        expect_identical(serialize(SummarizedExperiment::colData(d$se_NC), NULL), original)
        expect_identical(.Random.seed, seed)
    }
})
