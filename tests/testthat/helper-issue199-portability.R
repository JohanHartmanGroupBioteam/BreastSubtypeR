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

# These ten fixture routes supersede the old leaf tolerance only in complete
# expanded-example comparisons. All 19 policy rows remain unchanged.
issue199_ror_routes <- function() {
    data.frame(dataset = rep(c("OSLO2EMIT0obj", "TCGABRCAobj"), each = 5L),
        method = rep(c("genefu.robust", "ssBC", "ssBC.v2", "cIHC", "PCAPAM50"), 2L),
        stringsAsFactors = FALSE)
}

issue199_numeric_patterns <- function(a, b, p) {
    fail <- function(reason) stop(paste(p, reason), call. = FALSE)
    if (!identical(typeof(a), typeof(b)) ||
        !identical(attributes(a), attributes(b), num.eq = FALSE)) fail("numeric type/attributes")
    if (!identical(is.na(a), is.na(b)) || !identical(is.nan(a), is.nan(b)) ||
        !identical(is.infinite(a), is.infinite(b)) ||
        !identical(a[is.infinite(a)], b[is.infinite(b)])) fail("nonfinite masks/signs")
    x <- a[is.finite(a)]; y <- b[is.finite(b)]
    if (!identical(x == 0, y == 0) || !identical(sign(x), sign(y))) fail("zero support/signs")
    z <- which(x == 0)
    if (!identical(1 / x[z], 1 / y[z])) fail("signed zero")
    tx <- x != 0 & abs(x) < .Machine$double.xmin
    ty <- y != 0 & abs(y) < .Machine$double.xmin
    if (!identical(tx, ty) || !identical(x[tx], y[ty])) fail("subnormal support/values")
    if (!is.double(a)) fail("allowances require double payloads")
    invisible(TRUE)
}

issue199_ror_export_sides <- function(a, b, where) {
    cuts <- 100 * (c(-.25, .1) + .35) / .85
    if (!identical(outer(a, cuts, "<"), outer(b, cuts, "<")) ||
        !identical(outer(a, cuts, ">"), outer(b, cuts, ">")))
        stop(paste(where, "ROR cutoff sides"), call. = FALSE)
    invisible(TRUE)
}

# Unit tests of the three retained policy rows, never a scientific entry point.
issue199_ror_policy_leaf_unit <- function(actual, expected, dataset,
    scope = "C_vs_approved_C", path) {
    paths <- paste0("/res/results/", c("genefu.robust", "ssBC", "ssBC.v2"),
        "/score.ROR/ROR-P (Subtype + Prolif)")
    if (!identical(dataset, "TCGABRCAobj") || !identical(scope, "C_vs_approved_C") ||
        !path %in% paths) stop("Not a legacy ROR policy-leaf unit case", call. = FALSE)
    policy <- issue199_portability_policy()
    rule <- policy[policy$dataset == dataset & policy$scope == scope & policy$path == path, ]
    stopifnot(nrow(rule) == 1L)
    issue199_numeric_patterns(actual, expected, path)
    finite <- is.finite(actual)
    if (!all(abs(actual[finite] - expected[finite]) <= rule$atol +
        rule$rtol * pmax(abs(actual[finite]), abs(expected[finite]))))
        stop(paste(path, "elementwise numeric bound exceeded"), call. = FALSE)
    issue199_ror_export_sides(actual, expected, path)
    invisible(TRUE)
}

issue199_ror_approved <- function(dataset) {
    md5 <- c(OSLO2EMIT0obj = "04612ddb44d51c8b9faef21c7de892a0",
             TCGABRCAobj = "67359af00cb83c2bc0eb2b990b82012a")
    if (!dataset %in% names(md5)) stop("Unknown ROR fixture", call. = FALSE)
    file <- system.file("data", paste0(dataset, ".RData"), package = "BreastSubtypeR")
    if (!nzchar(file) || !identical(unname(tools::md5sum(file)), unname(md5[dataset])))
        stop(paste(dataset, "approved-example file provenance"), call. = FALSE)
    env <- new.env(parent = emptyenv())
    utils::data(list = dataset, package = "BreastSubtypeR", envir = env)
    env[[dataset]]
}

# Independent published formula: no production ROR/proliferation helper.
# scale() retains colMeans centering, sample SD and division order; mean() has
# R's correction pass. Never apply bitwise reconstruction across runtimes.
issue199_ror_formula <- function(v, d) {
    z <- scale(v)
    panel <- rownames(v) %in% c("CCNB1", "UBE2C", "BIRC5", "KNTC2", "CDC20",
        "PTTG1", "RRM2", "MKI67", "TYMS", "CEP55", "CDCA1")
    p <- apply(z[panel, , drop = FALSE], 2L, mean, na.rm = TRUE)
    raw <- -.0009299747 * d[, 1L] + .0692289192 * d[, 2L] +
        -.0951505484 * d[, 3L] + .0493487685 * d[, 4L] + .3385116381 * p
    list(raw = unname(raw), score = unname(100 * (raw + .35) / .85),
        sd = unname(attr(z, "scaled:scale")), m1 = unname(colMeans(abs(v))),
        panel = sum(panel))
}

# Arithmetic unit interface; full scientific admission additionally requires
# issue199_ror_context. No caller-supplied thresholds or fallback modes.
issue199_ror_domain <- function(v, d, where) {
    fail <- function(reason, j = 1L, observed = "")
        stop(paste(where, "sample", j, reason, observed), call. = FALSE)
    # Necessary runtime properties; not proof of compiler/rounding conformance.
    if (.Machine$double.base != 2 || .Machine$double.digits != 53 ||
        .Machine$double.rounding != 5 || .Machine$sizeof.pointer != 8 ||
        (!is.null(.Machine$longdouble.digits) && .Machine$longdouble.digits < 53))
        fail("unsupported arithmetic model")
    if (!is.matrix(v) || !is.double(v) || nrow(v) != 50L || !ncol(v))
        fail("50-gene domain")
    panel <- rownames(v) %in% c("CCNB1", "UBE2C", "BIRC5", "KNTC2", "CDC20",
        "PTTG1", "RRM2", "MKI67", "TYMS", "CEP55", "CDCA1")
    if (sum(panel) != 11L) fail("11-panel domain", observed = sum(panel))
    if (!all(is.finite(v))) fail("finite parent domain")
    if (any(abs(v) > 16)) fail("magnitude <= 16", observed = max(abs(v)))
    if (!is.matrix(d) || !is.double(d) || !identical(dim(d), c(ncol(v), 5L)) ||
        !all(is.finite(d)) || any(abs(d) > 1)) fail("five-distance domain")
    o <- issue199_ror_formula(v, d)
    bad <- which(!is.finite(o$sd) | o$sd < 5/16)
    if (length(bad)) fail("SD >= 5/16", bad[1L], o$sd[bad[1L]])
    bad <- which(!is.finite(o$m1) | o$m1 > 1.5 * o$sd)
    if (length(bad)) fail("M1 <= 1.5 SD", bad[1L], o$m1[bad[1L]])
    margin <- pmin(abs(o$raw + .25), abs(o$raw - .1))
    bad <- which(!is.finite(margin) | margin <= 2^-40)
    if (length(bad)) fail("raw cutoff margin > 2^-40", bad[1L], margin[bad[1L]])
    o
}

issue199_ror_method <- function(a, b, ref, dataset, method) {
    where <- paste(dataset, method)
    fail <- function(reason) stop(paste(where, reason), call. = FALSE)
    ids <- colnames(ref$outList$testData)
    genes <- rownames(ref$outList$testData)
    key_ok <- function(x) is.character(x) && !anyNA(x) && all(nzchar(x)) && !anyDuplicated(x)
    if (!key_ok(ids) || !key_ok(genes)) fail("approved identities")
    for (o in list(a, b)) {
        v <- o$outList$testData; d <- o$outList$distances; s <- o$score.ROR
        if (!is.matrix(v) || !identical(dim(v), dim(ref$outList$testData)) ||
            !identical(rownames(v), genes) || !identical(colnames(v), ids))
            fail("gene/sample identities, dimensions or ordering")
        if (!is.matrix(d) || !identical(dim(d), c(length(ids), 5L)) ||
            !identical(dimnames(d), dimnames(ref$outList$distances)))
            fail("distance identities, dimensions or ordering")
        if (!is.data.frame(s) || !identical(names(s), names(ref$score.ROR)) ||
            !identical(rownames(s), rownames(ref$score.ROR)) ||
            !identical(s$patientID, ref$score.ROR$patientID) ||
            !identical(as.character(s$patientID), ids) ||
            !identical(o$outList$predictions, ref$outList$predictions))
            fail("score/prediction identities or subtype calls")
    }
    x <- a$outList$testData; y <- b$outList$testData
    issue199_numeric_patterns(x, y, paste(where, "parent"))
    ranks <- function(v) apply(v, 2L, rank, ties.method = "min", na.last = "keep")
    if (!identical(ranks(x), ranks(y))) fail("parent ranks/ties")
    if (!identical(a$outList$distances, b$outList$distances,
        num.eq = FALSE, single.NA = FALSE)) fail("five-distance identity")
    if (dataset == "TCGABRCAobj") {
        issue199_portable_compare(x, y, dataset,
            path = paste0("/res/results/", method, "/outList/testData"))
    } else if (!identical(x, y, num.eq = FALSE, single.NA = FALSE)) fail("exact OSLO parent")
    field <- "ROR-P (Subtype + Prolif)"
    sa <- a$score.ROR[[field]]; sb <- b$score.ROR[[field]]
    issue199_numeric_patterns(sa, sb, paste(where, "ROR"))
    if (!all(is.finite(sa)) || !all(is.finite(sb))) fail("finite ROR domain")
    oa <- issue199_ror_domain(x, a$outList$distances, paste(where, "actual"))
    ob <- issue199_ror_domain(y, b$outList$distances, paste(where, "reference"))
    cuts <- c(-.25, .1)
    if (!identical(outer(oa$raw, cuts, "<"), outer(ob$raw, cuts, "<")) ||
        !identical(outer(oa$raw, cuts, ">"), outer(ob$raw, cuts, ">"))) fail("raw cutoff sides")
    group <- function(raw) ifelse(raw < -.25, "low", ifelse(raw > .1, "high", "med"))
    gf <- "ROR-P Group (Subtype + Prolif)"
    if (!identical(a$score.ROR[[gf]], b$score.ROR[[gf]]) ||
        !identical(a$score.ROR[[gf]], group(oa$raw)) ||
        !identical(b$score.ROR[[gf]], group(ob$raw))) fail("ROR risk groups")
    issue199_ror_export_sides(sa, sb, where)
    if (!all(abs(sa - sb) <= 2^-37)) fail("fixture ROR bound exceeded")
    invisible(length(ids))
}

issue199_ror_context <- function(actual, expected, dataset) {
    approved <- issue199_ror_approved(dataset)$res
    routes <- issue199_ror_routes()
    methods <- routes$method[routes$dataset == dataset]
    if (length(methods) != 5L || anyDuplicated(methods))
        stop("ROR route inventory", call. = FALSE)
    # Generic B2 fixtures have only posterior methods and synthetic axes.
    active <- any(c(methods, "cIHC.itr") %in% c(names(actual$results), names(expected$results),
        names(actual$res_subtypes), names(expected$res_subtypes))) ||
        identical(rownames(expected$res_subtypes), rownames(approved$res_subtypes))
    if (!active) return(character())
    if (!identical(expected, approved, num.eq = FALSE, single.NA = FALSE))
        stop(paste(dataset, "approved-example reference provenance"), call. = FALSE)
    if (!all(methods %in% names(actual$results)) || anyDuplicated(names(actual$results)))
        stop(paste(dataset, "missing or duplicated ROR routes"), call. = FALSE)
    if (!identical(actual$res_subtypes, approved$res_subtypes, num.eq = FALSE, single.NA = FALSE))
        stop(paste(dataset, "evaluated/excluded sample sets, subtype calls or entropy"), call. = FALSE)
    counts <- vapply(methods, function(m) issue199_ror_method(actual$results[[m]],
        expected$results[[m]], approved$results[[m]], dataset, m), integer(1))
    required <- if (dataset == "OSLO2EMIT0obj") c(102L, 102L, 94L, 102L, 102L) else
        c(100L, 100L, 96L, 100L, 100L)
    if (!identical(unname(counts), required)) stop("ROR evaluated/excluded counts", call. = FALSE)
    paste0("/res/results/", methods, "/score.ROR/ROR-P (Subtype + Prolif)")
}

issue199_ror_visits <- function(visited, required) {
    if (!identical(sort(visited), sort(required)) || anyDuplicated(visited))
        stop("unused or multiply visited ROR route", call. = FALSE)
    invisible(TRUE)
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
    ror_leaf <- paste0("/res/results/", c("genefu.robust", "ssBC", "ssBC.v2", "cIHC", "PCAPAM50"),
        "/score.ROR/ROR-P (Subtype + Prolif)")
    if (path %in% ror_leaf) stop("ROR scientific comparison requires full parent context", call. = FALSE)
    ror_paths <- character()
    if (scope == "C_vs_approved_C" && path %in% c("", "/res")) {
        ar <- if (path == "") actual$res else actual
        br <- if (path == "") expected$res else expected
        if (!is.null(ar) && !is.null(br)) ror_paths <- issue199_ror_context(ar, br, dataset)
    }
    ror_visited <- character()
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
        issue199_numeric_patterns(a, b, p)
        x <- a[is.finite(a)]; y <- b[is.finite(b)]
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
        if (p %in% ror_paths) {
            ror_visited <<- c(ror_visited, p)
            if (nrow(rule)) visited <<- c(visited, p)
        } else if (nrow(rule)) {
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
    issue199_ror_visits(ror_visited, ror_paths)
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
