# Fixture-only contract tests; synthetic columns below exercise arithmetic
# guards, not the admission/provenance of a full scientific comparison.
ror31_fresh <- local({
    cache <- new.env(parent = emptyenv())
    function(dataset) {
        if (!exists(dataset, cache, inherits = FALSE)) {
            obj <- issue199_ror_approved(dataset)
            value <- withr::with_preserve_seed({
                kind <- RNGkind()
                tryCatch({
                    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
                    set.seed(1992026)
                    mapped <- suppressMessages(Mapping(obj$se_obj,
                        RawCounts = dataset == "TCGABRCAobj", method = "max",
                        impute = TRUE, verbose = FALSE))
                    set.seed(20260831)
                    suppressMessages(BS_Multi(mapped, methods = "AUTO", Subtype = FALSE,
                        hasClinical = FALSE, er_subsample = "target"))
                }, finally = do.call(RNGkind, as.list(kind)))
            })
            assign(dataset, value, cache)
        }
        get(dataset, cache, inherits = FALSE)
    }
})

ror31_next <- function(x) {
    stopifnot(length(x) == 1L, is.finite(x), x != 0)
    b <- as.integer(writeBin(x, raw(), size = 8L, endian = "little"))
    carry <- if (x > 0) 1L else -1L
    for (j in seq_len(8L)) {
        t <- b[j] + carry
        b[j] <- t %% 256L
        carry <- if (t < 0) -1L else if (t > 255) 1L else 0L
        if (!carry) break
    }
    readBin(as.raw(b), "double", n = 1L, size = 8L, endian = "little")
}

ror31_formula_check <- function(result) {
    want <- issue199_ror_formula(result$outList$testData, result$outList$distances)$score
    if (!identical(result$score.ROR[["ROR-P (Subtype + Prolif)"]], want,
        num.eq = FALSE, single.NA = FALSE)) stop("same-runtime formula identity", call. = FALSE)
    invisible(TRUE)
}

ror31_method <- function(a, b = a, dataset = "TCGABRCAobj", method = "ssBC") {
    ref <- issue199_ror_approved(dataset)$res$results[[method]]
    issue199_ror_method(a, b, ref, dataset, method)
}

ror31_arithmetic <- function() {
    ref <- issue199_ror_approved("TCGABRCAobj")$res$results$ssBC$outList
    v <- matrix(rep(c(-1, 1), 25L), 50L, 1L,
        dimnames = list(rownames(ref$testData), "synthetic"))
    d <- matrix(0, 1L, 5L, dimnames = list("synthetic", colnames(ref$distances)))
    list(v = v, d = d)
}

test_that("ROR fixture provenance and all ten routes retain exactly 998 comparisons", {
    routes <- issue199_ror_routes()
    expect_identical(nrow(routes), 10L)
    expect_identical(anyDuplicated(routes), 0L)
    total <- 0L
    for (dataset in c("OSLO2EMIT0obj", "TCGABRCAobj")) {
        obj <- issue199_ror_approved(dataset)
        paths <- issue199_ror_context(obj$res, obj$res, dataset)
        expect_length(paths, 5L)
        expect_true(issue199_portable_compare(obj$res, obj$res, dataset, path = "/res"))
        for (method in routes$method[routes$dataset == dataset]) {
            o <- obj$res$results[[method]]
            n <- issue199_ror_method(o, o, o, dataset, method)
            total <- total + n
            expect_identical(colnames(o$outList$testData), as.character(o$score.ROR$patientID))
        }
    }
    expect_identical(total, 998L)
    expect_error(issue199_ror_approved("unlisted"), "Unknown ROR fixture")
})

test_that("ROR formula identity applies to fresh current-runtime results without changing RNG", {
    withr::local_seed(3112026)
    rng <- .Random.seed
    total <- 0L
    for (dataset in c("OSLO2EMIT0obj", "TCGABRCAobj")) {
        approved <- issue199_ror_approved(dataset)$res
        full <- ror31_fresh(dataset)
        expect_identical(.Random.seed, rng)
        fresh <- auto_ssbc_historical_view(full, approved)
        expect_true(issue199_portable_compare(fresh, approved, dataset, path = "/res"))
        expect_identical(fresh$res_subtypes, approved$res_subtypes)
        for (method in issue199_ror_routes()$method[seq_len(5L)]) {
            expect_true(ror31_formula_check(fresh$results[[method]]), info = paste(dataset, method))
            total <- total + ncol(fresh$results[[method]]$outList$testData)
        }
    }
    expect_identical(total, 998L)
    expect_identical(.Random.seed, rng)
})

test_that("ROR formula regression detects published-coefficient and oracle mutants", {
    approved <- issue199_ror_approved("TCGABRCAobj")$res
    fresh <- auto_ssbc_historical_view(ror31_fresh("TCGABRCAobj"), approved)
    original <- BreastSubtypeR:::RORgroup
    change <- function(expr, value) {
        if (missing(expr)) return(quote(expr = ))
        if (is.numeric(expr) && length(expr) == 1L && identical(expr, value)) return(expr + 1e-10)
        if (is.call(expr)) for (j in seq_along(expr)) expr[[j]] <- change(expr[[j]], value)
        expr
    }
    # Copy the function; never modify the package namespace or files.
    for (coefficient in c(.0009299747, .0692289192, .0951505484, .0493487685, .3385116381)) {
        mutant <- original
        body(mutant) <- change(body(mutant), coefficient)
        expect_false(identical(body(mutant), body(original)))
        o <- fresh$results$ssBC
        o$score.ROR <- mutant(o$outList, NULL, Subtype = FALSE, hasClinical = FALSE)
        expect_error(ror31_formula_check(o), "same-runtime formula identity")
    }
    oracle <- issue199_ror_formula
    body(oracle) <- change(body(oracle), .3385116381)
    o <- fresh$results$ssBC
    expect_false(identical(oracle(o$outList$testData, o$outList$distances)$score,
        o$score.ROR[["ROR-P (Subtype + Prolif)"]], num.eq = FALSE))
    expect_identical(BreastSubtypeR:::RORgroup, original)
})

test_that("ROR ceiling and independent parent rules have effective controls", {
    dataset <- "TCGABRCAobj"
    ref <- issue199_ror_approved(dataset)$res
    b <- ref$results$ssBC
    f <- "ROR-P (Subtype + Prolif)"
    a <- b
    a$score.ROR[[f]][1L] <- ror31_next(b$score.ROR[[f]][1L] + 2^-37)
    expect_gt(abs(a$score.ROR[[f]][1L] - b$score.ROR[[f]][1L]), 2^-37)
    expect_error(ror31_method(a, b), "fixture ROR bound exceeded")
    a <- b
    a$outList$testData["CCNB1", 1L] <- a$outList$testData["CCNB1", 1L] + 1e-10
    expect_error(ror31_method(a, b), "elementwise numeric bound exceeded")
    a <- b
    a$outList$testData["CCNB1", 1L] <- ror31_next(a$outList$testData["CCNB1", 1L])
    a$score.ROR[[f]] <- issue199_ror_formula(a$outList$testData, a$outList$distances)$score
    expect_true(is.integer(ror31_method(a, b)))
    actual <- ref; actual$results$ssBC <- a
    expect_true(issue199_portable_compare(actual, ref, dataset, path = "/res"))
    oslo <- issue199_ror_approved("OSLO2EMIT0obj")$res$results$ssBC
    a <- oslo; a$outList$testData[1L, 1L] <- ror31_next(a$outList$testData[1L, 1L])
    expect_error(ror31_method(a, oslo, "OSLO2EMIT0obj"), "exact OSLO parent")
})

test_that("ROR structural gates precede numerical and missing-score decisions", {
    b <- issue199_ror_approved("TCGABRCAobj")$res$results$ssBC
    f <- "ROR-P (Subtype + Prolif)"
    for (column in c(1L, 5L)) {
        a <- b
        a$outList$distances[1L, column] <- ror31_next(a$outList$distances[1L, column])
        expect_error(ror31_method(a, b), "five-distance identity")
    }
    for (missing in c(FALSE, TRUE)) {
        a <- b; expected <- b
        if (missing) { a$score.ROR[[f]][] <- NA_real_; expected$score.ROR[[f]][] <- NA_real_ }
        a$outList$testData <- a$outList$testData[c(2L, 1L, 3:50), , drop = FALSE]
        expect_error(ror31_method(a, expected), "gene/sample identities")
        a <- expected
        colnames(a$outList$testData)[1:2] <- rev(colnames(a$outList$testData)[1:2])
        expect_error(ror31_method(a, expected), "gene/sample identities")
        a <- expected; a$score.ROR$patientID[1L] <- a$score.ROR$patientID[2L]
        expect_error(ror31_method(a, expected), "score/prediction identities")
        a <- expected; a$outList$distances[1L, 5L] <- ror31_next(a$outList$distances[1L, 5L])
        expect_error(ror31_method(a, expected), "five-distance identity")
    }
    a <- b; a$outList$testData[1L, 1L] <- NA_real_
    expect_error(ror31_method(a, b), "nonfinite masks/signs")
    a <- b; a$outList$testData[1:2, 1L] <- rev(a$outList$testData[1:2, 1L])
    expect_error(ror31_method(a, b), "parent ranks/ties")
    a <- b; a$outList$testData[1L, 1L] <- a$outList$testData[2L, 1L]
    expect_error(ror31_method(a, b), "parent ranks/ties")
    a <- b; a$score.ROR[["ROR-P Group (Subtype + Prolif)"]][1L] <- "corrupt"
    expect_error(ror31_method(a, b), "ROR risk groups")
})

test_that("ROR numerical guards accept their interior and refuse the exterior", {
    q <- ror31_arithmetic()
    domain <- function(v = q$v, d = q$d) issue199_ror_domain(v, d, "synthetic guard")
    expect_type(domain(), "list")
    centered <- q$v - colMeans(q$v)
    sd0 <- attr(scale(centered), "scaled:scale")
    for (inside in c(TRUE, FALSE)) {
        direction <- if (inside) -1 else 1
        v <- centered * ((5/16)*(1-direction*2^-30)/sd0)
        if (inside) expect_type(domain(v), "list") else expect_error(domain(v), "SD >= 5/16")
        v <- centered + 1.5*sd0*(1+direction*2^-30)
        if (inside) expect_type(domain(v), "list") else expect_error(domain(v), "M1 <= 1.5 SD")
        v <- q$v * if (inside) 16 else 16*(1+2^-40)
        if (inside) expect_type(domain(v), "list") else expect_error(domain(v), "magnitude <= 16")
    }
    expect_error(domain(q$v[-1L, , drop = FALSE]), "50-gene domain")
    expect_error(domain(rbind(q$v, extra = 0)), "50-gene domain")
    v <- q$v; rownames(v)[match("CCNB1", rownames(v))] <- "not-panel"
    expect_error(domain(v), "11-panel domain")
    for (value in c(NA_real_, NaN, Inf, -Inf)) {
        v <- q$v; v[1L, 1L] <- value
        expect_error(domain(v), "finite parent domain")
    }
    p <- issue199_ror_formula(q$v, q$d)$raw
    coefficients <- c(-.0009299747, .0692289192, -.0951505484, .0493487685)
    cut <- c(-.25, .1)[which.min(abs(p-c(-.25, .1)))]
    for (gap in c(2^-39, 2^-41)) {
        d <- q$d; d[1L, 1:4] <- sign(coefficients)*(cut + gap - p)/sum(abs(coefficients))
        if (gap > 2^-40) expect_type(domain(d = d), "list") else
            expect_error(domain(d = d), "raw cutoff margin > 2\\^-40")
    }
})

test_that("ROR payload patterns and bare-leaf refusal cannot bypass the contract", {
    b <- issue199_ror_approved("TCGABRCAobj")$res$results$ssBC
    f <- "ROR-P (Subtype + Prolif)"
    for (value in c(NA_real_, NaN, Inf, -Inf)) {
        a <- b; a$score.ROR[[f]][1L] <- value
        expect_error(ror31_method(a, b), "nonfinite masks/signs")
        expected <- a
        expect_error(ror31_method(a, expected), "finite ROR domain")
    }
    a <- b; expected <- b
    a$score.ROR[[f]][1L] <- -0; expected$score.ROR[[f]][1L] <- 0
    expect_error(ror31_method(a, expected), "signed zero")
    for (method in c("genefu.robust", "ssBC", "ssBC.v2")) {
        path <- paste0("/res/results/", method, "/score.ROR/ROR-P (Subtype + Prolif)")
        expect_true(issue199_ror_policy_leaf_unit(1+2^-47, 1, "TCGABRCAobj", path = path))
        expect_error(issue199_ror_policy_leaf_unit(1+2^-44, 1, "TCGABRCAobj", path = path), "numeric bound exceeded")
        expect_error(issue199_portable_compare(1, 1, "TCGABRCAobj", path = path), "full parent context")
    }
    cut <- 100*(-.25+.35)/.85
    path <- "/res/results/ssBC/score.ROR/ROR-P (Subtype + Prolif)"
    expect_error(issue199_ror_policy_leaf_unit(cut+2^-49, cut, "TCGABRCAobj", path = path), "ROR cutoff sides")
    expect_error(issue199_ror_policy_leaf_unit(1, 1, "TCGABRCAobj",
        path = "/res/results/cIHC.itr/score.ROR/ROR-P (Subtype + Prolif)"), "Not a legacy")
})

test_that("ROR provenance, exclusions, entropy and route visits remain mandatory", {
    dataset <- "OSLO2EMIT0obj"
    ref <- issue199_ror_approved(dataset)$res
    a <- ref; a$results$ssBC <- NULL
    expect_error(issue199_portable_compare(a, ref, dataset, path = "/res"), "missing or duplicated ROR routes")
    a <- ref; a$res_subtypes$entropy[1L] <- a$res_subtypes$entropy[1L]+1
    expect_error(issue199_portable_compare(a, ref, dataset, path = "/res"), "subtype calls or entropy")
    a <- ref; k <- which(is.na(a$res_subtypes$ssBC.v2))[1L]; a$res_subtypes$ssBC.v2[k] <- "LumA"
    expect_error(issue199_portable_compare(a, ref, dataset, path = "/res"), "evaluated/excluded sample sets")
    a <- ref; a$results$ssBC$outList$testData[1L, 1L] <- 1
    expect_error(issue199_portable_compare(a, a, dataset, path = "/res"), "reference provenance")
    routes <- issue199_ror_context(ref, ref, dataset)
    expect_true(issue199_ror_visits(routes, routes))
    expect_error(issue199_ror_visits(routes[-1L], routes), "unused or multiply visited")
    expect_error(issue199_ror_visits(c(routes, routes[1L]), routes), "unused or multiply visited")
    for (case in list(c("cIHC.itr", "ROR-P (Subtype + Prolif)"),
                      c("ssBC", "ROR-S (Subtype Only)"))) {
        method <- case[1L]; field <- case[2L]
        frozen <- ref$results[[method]]$score.ROR
        a <- ref
        a$results[[method]]$score.ROR[[field]][1L] <- ror31_next(frozen[[field]][1L])
        # Data-frame replacement can reorder attributes. Keep this mutation
        # numerical only; the comparator's structural contract stays exact.
        attributes(a$results[[method]]$score.ROR) <- attributes(frozen)
        changed <- a$results[[method]]$score.ROR
        expect_identical(attributes(changed), attributes(frozen))
        expect_identical(names(attributes(changed)), names(attributes(frozen)))
        expect_identical(names(changed), names(frozen))
        expect_identical(class(changed), class(frozen))
        expect_identical(lapply(changed, class), lapply(frozen, class))
        expect_identical(rownames(changed), rownames(frozen))
        expect_identical(changed$patientID, frozen$patientID)
        expect_false(identical(changed[[field]][1L], frozen[[field]][1L],
            num.eq = FALSE, single.NA = FALSE))
        restored <- a
        restored$results[[method]]$score.ROR[[field]][1L] <- frozen[[field]][1L]
        attributes(restored$results[[method]]$score.ROR) <- attributes(frozen)
        expect_true(identical(restored, ref, num.eq = FALSE, single.NA = FALSE))
        expect_error(issue199_portable_compare(a, ref, dataset, path = "/res"),
            "unlisted value must be exact")
    }
})
