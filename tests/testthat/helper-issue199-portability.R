# Measured PR #228 Linux allowances; every other value and all structure stay exact.
# Historical TCGA B uses the independently verified Windows-B == frozen-A mapping
# reference bridge. OSLO historical A/B is deliberately outside this policy.
issue199_portability_policy <- function() {
    data.frame(
        scope = c(rep("C_vs_approved_C", 18L), "B_mapping_vs_frozen_A"),
        dataset = c(rep("OSLO2EMIT0obj", 5L), rep("TCGABRCAobj", 14L)),
        path = c(
            "/data_input/se_SSP/assay",
            "/res/results/AIMS/all.probs/20",
            "/res/results/AIMS/data.used",
            "/res/results/sspbc/all.probs/24",
            "/res/results/sspbc/data.used",
            "/data_input/se_NC/assay",
            "/res/results/genefu.robust/score.ROR/ROR-P (Subtype + Prolif)",
            "/res/results/genefu.robust/outList/testData",
            "/res/results/ssBC/score.ROR/ROR-P (Subtype + Prolif)",
            "/res/results/ssBC/outList/testData",
            "/res/results/ssBC.v2/score.ROR/ROR-P (Subtype + Prolif)",
            "/res/results/ssBC.v2/outList/testData",
            "/res/results/cIHC/outList/testData",
            "/res/results/cIHC.itr/score.ROR/ROR-P (Subtype + Prolif)",
            "/res/results/cIHC.itr/outList/testData",
            "/res/results/PCAPAM50/outList/testData",
            "/res/results/AIMS/all.probs/20",
            "/res/results/sspbc/all.probs/24",
            "/data_input/se_NC/assay"),
        atol = c(0, 0, 0, 0, 0, 2^-48, 2^-46, 2^-48, 2^-46, 2^-48, 2^-46, 2^-48, 2^-48, 2^-46, 2^-48, 2^-48, 0, 0, 2^-48),
        rtol = c(2^-51, 2^-44, 2^-51, 2^-44, 2^-51, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2^-44, 2^-45, 0),
        stringsAsFactors = FALSE)
}

# Exact within-object selection; masked rows are handled before first argmax.
issue199_posterior_link <- function(o, where) {
    fail <- function(reason) stop(paste(where, reason), call. = FALSE)
    fields <- c("prob", "cl", "all.probs", "rules.matrix")
    if (!is.list(o) || anyDuplicated(names(o)) || !all(fields %in% names(o)))
        fail("posterior fields missing or duplicated")
    ks <- names(o$all.probs)
    if (!is.list(o$all.probs) || !length(ks) ||
        length(ks) != length(o$all.probs) || anyNA(ks) ||
        any(ks == "") || anyDuplicated(ks) ||
        !is.matrix(o$prob) || !is.double(o$prob) ||
        !is.matrix(o$cl) || !is.character(o$cl) ||
        !identical(dim(o$prob), dim(o$cl)) ||
        !identical(dimnames(o$prob), dimnames(o$cl)) ||
        !identical(colnames(o$prob), ks))
        fail("prob/cl/all.probs names or structure")
    samples <- rownames(o$prob)
    if (!is.matrix(o$rules.matrix) || ncol(o$rules.matrix) != nrow(o$prob) ||
        (!is.null(colnames(o$rules.matrix)) &&
            !identical(colnames(o$rules.matrix), samples)))
        fail("rules.matrix sample axis")
    for (k in ks) {
        p <- o$all.probs[[k]]
        classes <- colnames(p)
        if (!is.matrix(p) || !is.double(p) || nrow(p) != nrow(o$prob) ||
            !length(classes) || anyNA(classes) || any(classes == "") ||
            anyDuplicated(classes) ||
            (!is.null(rownames(p)) && !identical(rownames(p), samples)))
            fail("all.probs class/sample axes or structure")
        n_na <- rowSums(is.na(p))
        masked <- n_na == ncol(p)
        if (any(is.nan(p)) || any(n_na != 0L & !masked) ||
            !all(is.finite(p[!masked, , drop = FALSE])))
            fail("posterior rows must be finite or entirely NA_real_")
        if (any(p < 0 | p > 1, na.rm = TRUE)) fail("posterior outside [0,1]")
        probability <- rep(NA_real_, nrow(p))
        call <- rep(NA_character_, nrow(p))
        if (any(!masked)) {
            j <- max.col(p[!masked, , drop = FALSE], ties.method = "first")
            probability[!masked] <- p[cbind(which(!masked), j)]
            call[!masked] <- classes[j]
        }
        if (!identical(unname(o$prob[, k]), probability, num.eq = FALSE) ||
            !identical(unname(o$cl[, k]), call))
            fail("selected probability/class link")
    }
    invisible(TRUE)
}

issue199_portable_compare <- function(actual, expected, dataset,
    scope = "C_vs_approved_C", path = "") {
    policy <- issue199_portability_policy()
    policy <- policy[policy$dataset == dataset & policy$scope == scope, , drop = FALSE]
    if (!nrow(policy)) stop("Unknown portability comparison scope", call. = FALSE)
    roots <- c("", "/data_input", "/res", "/data_input/se_NC/assay",
        "/data_input/se_SSP/assay", policy$path)
    if (!path %in% roots) stop("Unknown portability root path", call. = FALSE)
    policy <- policy[policy$path == path | startsWith(policy$path, paste0(path, "/")),
        , drop = FALSE]
    visited <- character()
    derived <- c("/res/results/AIMS/prob", "/res/results/sspbc/prob")
    linked <- character()
    fail <- function(p, reason) stop(paste(p, reason), call. = FALSE)
    ranks <- function(m) t(apply(m, 1L, function(v)
        rank(-v, ties.method = "min", na.last = "keep")))
    numeric_payload <- function(a, b, p, rule) {
        if (!identical(typeof(a), typeof(b)) ||
            !identical(attributes(a), attributes(b), num.eq = FALSE)) fail(p, "numeric type/attributes")
        if (!identical(is.na(a), is.na(b)) || !identical(is.nan(a), is.nan(b)) ||
            !identical(is.infinite(a), is.infinite(b)) ||
            !identical(a[is.infinite(a)], b[is.infinite(b)])) fail(p, "nonfinite masks/signs")
        finite <- is.finite(a)
        x <- a[finite]; y <- b[finite]
        if (!identical(x == 0, y == 0) || !identical(sign(x), sign(y)))
            fail(p, "zero support/signs")
        z <- which(x == 0)
        if (!identical(1 / x[z], 1 / y[z])) fail(p, "signed zero")
        tiny_x <- x != 0 & abs(x) < .Machine$double.xmin
        tiny_y <- y != 0 & abs(y) < .Machine$double.xmin
        if (!identical(tiny_x, tiny_y) || !identical(x[tiny_x], y[tiny_y]))
            fail(p, "subnormal support/values")
        if (!is.double(a)) fail(p, "allowances require double payloads")
        if (!all(abs(x - y) <= rule$atol + rule$rtol * pmax(abs(x), abs(y))))
            fail(p, "elementwise numeric bound exceeded")
    }
    walk <- function(a, b, p) {
        if (!identical(typeof(a), typeof(b)) || !identical(class(a), class(b)) ||
            length(a) != length(b)) fail(p, "type/class/length")
        if (isS4(a)) {
            slots <- methods::slotNames(a)
            if (!identical(slots, methods::slotNames(b))) fail(p, "S4 slots")
            ax <- attributes(a); bx <- attributes(b)
            ax[slots] <- NULL; bx[slots] <- NULL
            if (!identical(ax, bx, num.eq = FALSE)) fail(p, "S4 attributes")
            for (n in slots) walk(methods::slot(a, n), methods::slot(b, n),
                paste0(p, "/", n))
            return(invisible(NULL))
        }
        if (!identical(attributes(a), attributes(b), num.eq = FALSE)) fail(p, "attributes/names/dimensions")
        # Resolve the complete single-assay storage path, preserving every S4
        # container, slot, assay name/count, storage class and column attribute.
        assay_paths <- c(
            "/data_input/se_NC/assays/data/listData/counts",
            "/data_input/se_SSP/assays/data/listData/counts")
        i <- match(p, assay_paths)
        if (!is.na(i)) p <- c("/data_input/se_NC/assay", "/data_input/se_SSP/assay")[i]
        rule <- policy[policy$path == p, , drop = FALSE]
        if (nrow(rule)) {
            visited <<- c(visited, p)
            if (is.data.frame(a)) {
                if (!all(vapply(a, is.double, logical(1))) ||
                    !all(vapply(b, is.double, logical(1)))) fail(p, "assay double columns")
                for (j in seq_along(a)) numeric_payload(a[[j]], b[[j]], p, rule)
            } else {
                numeric_payload(a, b, p, rule)
            }
            if (p %in% c("/res/results/AIMS/all.probs/20", "/res/results/sspbc/all.probs/24")) {
                if (!is.matrix(a) || !all(a >= 0 & a <= 1, na.rm = TRUE) ||
                    !identical(ranks(a), ranks(b)) ||
                    !identical(max.col(a, ties.method = "first"),
                        max.col(b, ties.method = "first"))) fail(p, "posterior ranks/ties/argmax/bounds")
            }
            if (endsWith(p, "/score.ROR/ROR-P (Subtype + Prolif)")) {
                cuts <- 100 * (c(-.25, .1) + .35) / .85
                if (!identical(outer(a, cuts, "<"), outer(b, cuts, "<")) ||
                    !identical(outer(a, cuts, ">"), outer(b, cuts, ">")))
                    fail(p, "ROR cutoff sides")
            }
        } else if (is.list(a)) {
            # Existing generic fixtures omit all three selection/rule fields.
            # Any present field requires a complete parent in both objects.
            if (p %in% dirname(derived) && any(c("prob", "cl", "rules.matrix") %in%
                c(names(a), names(b)))) {
                issue199_posterior_link(a, paste(p, "actual"))
                issue199_posterior_link(b, paste(p, "reference"))
                linked <<- c(linked, p)
            }
            for (j in seq_along(a)) walk(a[[j]], b[[j]], paste0(p, "/",
                if (length(names(a))) names(a)[j] else j))
        } else if (p %in% derived) {
            if (!dirname(p) %in% linked || !identical(is.na(a), is.na(b)))
                fail(p, "derived posterior without verified link or with changed mask")
        } else if (!identical(a, b, num.eq = FALSE)) fail(p, "unlisted value must be exact")
        invisible(NULL)
    }
    walk(actual, expected, path)
    if (!identical(sort(unique(visited)), sort(policy$path)) || anyDuplicated(visited))
        fail(path, "unused or multiply visited portability policy path")
    invisible(TRUE)
}

issue199_expect_decisions <- function(actual, expected) {
    testthat::expect_identical(actual$res_subtypes, expected$res_subtypes)
    rank_rows <- function(x) t(apply(as.matrix(x), 1L, function(v)
        rank(-v, ties.method = "min", na.last = "keep")))
    compare_ranks <- function(a, b) {
        testthat::expect_identical(rank_rows(a), rank_rows(b))
        testthat::expect_identical(max.col(as.matrix(a), ties.method = "first"),
            max.col(as.matrix(b), ties.method = "first"))
    }
    for (method in names(actual$results)) {
        a <- actual$results[[method]]; b <- expected$results[[method]]
        if (method %in% c("AIMS", "sspbc")) {
            testthat::expect_true(issue199_posterior_link(a, paste(method, "actual")))
            testthat::expect_true(issue199_posterior_link(b, paste(method, "reference")))
            testthat::expect_identical(a$cl, b$cl)
            testthat::expect_true(identical(a$rules.matrix, b$rules.matrix, num.eq = FALSE))
            testthat::expect_true(identical(attributes(a$prob), attributes(b$prob), num.eq = FALSE))
            testthat::expect_identical(is.na(a$prob), is.na(b$prob))
            testthat::expect_true(identical(attributes(a$all.probs), attributes(b$all.probs), num.eq = FALSE))
            for (k in names(a$all.probs)) {
                testthat::expect_true(identical(attributes(a$all.probs[[k]]),
                    attributes(b$all.probs[[k]]), num.eq = FALSE))
                compare_ranks(a$all.probs[[k]], b$all.probs[[k]])
            }
        } else {
            for (n in c("distances", "dist.RORSubtype"))
                compare_ranks(a$outList[[n]], b$outList[[n]])
            testthat::expect_identical(a$score.ROR$Confidence, b$score.ROR$Confidence)
            for (n in grep("^ROR-[SP] ", names(a$score.ROR), value = TRUE)) {
                if (!is.numeric(a$score.ROR[[n]])) next
                cuts <- if (startsWith(n, "ROR-S "))
                    100 * (c(-.15, .1) + .35) / .85 else 100 * (c(-.25, .1) + .35) / .85
                testthat::expect_identical(outer(a$score.ROR[[n]], cuts, "<"),
                    outer(b$score.ROR[[n]], cuts, "<"))
                testthat::expect_identical(outer(a$score.ROR[[n]], cuts, ">"),
                    outer(b$score.ROR[[n]], cuts, ">"))
            }
            groups <- grep("^ROR-[SP] Group ", names(a$score.ROR), value = TRUE)
            testthat::expect_identical(a$score.ROR[, groups, drop = FALSE],
                b$score.ROR[, groups, drop = FALSE])
        }
    }
}
