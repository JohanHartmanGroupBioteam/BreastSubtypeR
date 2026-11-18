# Shared fixture RNG scope. Restore the full caller configuration even when the
# caller had no seed (withr 3.0.3 only removes the seed in that case).
issue199_rng <- function(seed, code) {
    withr::with_preserve_seed({
        kind <- RNGkind()
        tryCatch({
            RNGkind("Mersenne-Twister", "Inversion", "Rejection")
            set.seed(seed)
            force(code)
        }, finally = do.call(RNGkind, as.list(kind)))
    })
}

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

issue199_posterior_models <- function() {
    hashes <- c(AIMSmodel = "6c541625fe22667540ef7c11d9005290",
                sspbc.models = "1c2ba126cc1ea8ecc15dd6dc4656285d")
    models <- new.env(parent = emptyenv())
    for (name in names(hashes)) {
        file <- system.file("data", paste0(name, ".RData"), package = "BreastSubtypeR")
        if (!nzchar(file) || !identical(unname(tools::md5sum(file)), hashes[[name]]))
            stop(paste(name, "posterior model provenance"), call. = FALSE)
        load(file, envir = models)
    }
    list(AIMS = models$AIMSmodel, sspbc = models$sspbc.models$ssp.pam50)
}

# Fixture proof assumes binary64 round-to-nearest, accumulator precision >= 53
# bits, the reviewed e1071 operation order, normal exp/log relative error <=
# 8u (u = 2^-53), gradual underflow with exp absolute error <= 4*2^-1074,
# exact exp(0)/log(1), and ordinary positive overflow. These probes establish
# observable prerequisites; runtime metadata cannot certify the libm bound.
issue199_posterior_arithmetic <- function() {
    eta <- 2^-1074
    if (!isTRUE(.Machine$double.base == 2L && .Machine$double.digits == 53L &&
        .Machine$double.rounding == 5L && .Machine$double.min.exp == -1022L &&
        .Machine$double.max.exp == 1024L && .Machine$sizeof.pointer == 8L &&
        (is.null(.Machine$longdouble.digits) || .Machine$longdouble.digits >= 53L) &&
        eta > 0 && .Machine$double.xmin / 2^52 == eta && eta/2 == 0 &&
        exp(0) == 1 && log(1) == 0 && exp(1000) == Inf && exp(-1000) == 0 &&
        abs(exp(-744) - 2*eta) <= 4*eta))
        stop("posterior arithmetic model", call. = FALSE)
}

issue199_posterior_x <- function(x) {
    if (!all(is.finite(x) & (x <= 32 | x >= 64)))
        stop("posterior domain guard G-X", call. = FALSE)
}

# Build both possible density columns once per predictor. Use the stored
# predictor order even when the supplied rule matrix has a different order.
issue199_posterior_terms <- function(nb, rules) {
    predictors <- names(nb$tables)
    rule_names <- make.names(rownames(rules))
    index <- match(predictors, rule_names)
    if (!length(index) || anyNA(index) || anyDuplicated(predictors) ||
        anyDuplicated(rule_names) || !all(is.na(rules) | rules == 0 | rules == 1))
        stop("posterior rules: predictor mapping or values", call. = FALSE)
    rules <- rules[index, , drop = FALSE]
    classes <- length(nb$apriori)
    log_density <- lapply(c(0, 1), function(value) {
        vapply(nb$tables, function(tab) {
            sd <- tab[, 2L]
            sd[sd <= 0] <- 0.001
            x <- abs((value - tab[, 1L])/sd)
            issue199_posterior_x(x)
            density <- dnorm(value, tab[, 1L], sd)
            if (!all(is.finite(density)) || any(density[x <= 32] < .Machine$double.xmin) ||
                any(density[x >= 64] != 0))
                stop("posterior domain guard G-X density/floor", call. = FALSE)
            density[density <= 0] <- 0.001
            log(density)
        }, numeric(classes))
    })
    prior <- log(nb$apriori)
    likelihood <- magnitude <- matrix(0, ncol(rules), classes,
        dimnames = list(NULL, names(prior)))
    for (i in seq_len(ncol(rules))) {
        column <- rules[, i]
        ones <- which(column == 1)
        terms <- log_density[[1L]]
        terms[, ones] <- log_density[[2L]][, ones, drop = FALSE]
        terms[, is.na(column)] <- 0 # missing predictor contributes density one
        likelihood[i, ] <- prior + apply(terms, 1L, sum)
        magnitude[i, ] <- abs(prior) + apply(abs(terms), 1L, sum)
    }
    list(L = likelihood, magnitude = magnitude, masked = colSums(!is.na(rules)) == 0L)
}

# Exact same-runtime oracle. Keep the reciprocal form, including overflow.
issue199_posterior_oracle <- function(nb, rules) {
    terms <- issue199_posterior_terms(nb, rules)
    result <- terms$L
    for (i in seq_len(nrow(result))) {
        L <- terms$L[i, ]
        result[i, ] <- if (terms$masked[i]) NA_real_ else
            vapply(L, function(lp) 1/sum(exp(L-lp)), numeric(1L))
    }
    result
}

issue199_posterior_domain <- function(terms, reference, method) {
    sum_cap <- switch(method, AIMS = 520, sspbc = 1520, stop("Unknown posterior method"))
    span_cap <- if (method == "AIMS") 500 else 1520
    L <- terms$L
    if (!all(is.finite(terms$magnitude) & terms$magnitude <= sum_cap))
        stop("posterior domain guard G-S1", call. = FALSE)
    spans <- apply(L, 1L, function(v) max(v) - min(v))
    if (!all(is.finite(spans) & spans <= span_cap))
        stop("posterior domain guard G-S2", call. = FALSE)
    # Subtract the row maximum before exp. The five-term sum is in [1, 5].
    log_sum <- t(vapply(seq_len(nrow(L)), function(i) {
        v <- L[i, ]
        peak <- max(v)
        (peak - v) + log(sum(exp(v - peak)))
    }, numeric(ncol(L))))
    zero <- reference == 0
    tiny <- reference > 0 & reference < .Machine$double.xmin
    normal <- reference >= .Machine$double.xmin
    overflow <- log(.Machine$double.xmax)
    boundary <- 1022*log(2)
    margin <- 2^-4
    if (!isTRUE(all(log_sum[zero] >= overflow + margin,
        log_sum[!zero] <= overflow - margin, log_sum[tiny] >= boundary + margin,
        log_sum[normal] <= boundary - margin)))
        stop("posterior domain guard G-O", call. = FALSE)
    gaps <- apply(reference, 1L, function(v) {
        positive <- sort(v[v > 0])
        if (length(positive) < 2L) return(Inf)
        min(diff(positive)/positive[-1L])
    })
    if (!all(gaps >= 2^-20)) stop("posterior domain guard G-R", call. = FALSE)
    invisible(list(max_sum = max(terms$magnitude), max_span = max(spans),
        min_gap = min(gaps), zero = sum(zero), subnormal = sum(tiny),
        overflow_margin = if (any(zero)) min(log_sum[zero] - overflow) else NA_real_,
        nonzero_margin = min(overflow - log_sum[!zero]),
        subnormal_margin = if (any(tiny)) min(log_sum[tiny] - boundary) else NA_real_,
        normal_margin = min(boundary - log_sum[normal])))
}

issue199_posterior_method <- function(parent, ref, model, method, dataset) {
    k <- switch(method, AIMS = 20L, sspbc = 24L, stop("Unknown posterior method"))
    key <- as.character(k)
    classes <- c("Basal", "Her2", "LumA", "LumB", "Normal")
    nb <- model$one.vs.all.tsp[[k]]
    if (!isTRUE(model$k == k) || length(nb$tables) != 5L*k ||
        !identical(names(nb$apriori), classes) ||
        !identical(names(nb$tables), make.names(model$all.pairs)) ||
        anyDuplicated(names(nb$tables)) ||
        !all(vapply(nb$tables, function(tab) identical(dim(tab), c(5L, 2L)) &&
            identical(rownames(tab), classes) && all(is.finite(tab)), logical(1L))))
        stop(paste(method, "posterior model provenance/structure"), call. = FALSE)
    samples <- switch(dataset, OSLO2EMIT0obj = 102L, TCGABRCAobj = 100L,
        stop("Unknown posterior dataset"))
    p <- parent$all.probs[[key]]
    q <- ref$all.probs[[key]]
    if (!identical(names(parent$all.probs), key) || !identical(dim(p), c(samples, 5L)) ||
        !identical(attributes(p), attributes(q)) || !identical(colnames(p), classes) ||
        !identical(dimnames(parent$prob), dimnames(ref$prob)) ||
        !identical(dimnames(parent$cl), dimnames(ref$cl)))
        stop(paste(method, "posterior structure/axes"), call. = FALSE)
    rules <- parent$rules.matrix
    if (!identical(rules, ref$rules.matrix, num.eq = FALSE, single.NA = FALSE) ||
        !identical(rownames(rules), model$all.pairs) ||
        !identical(dim(rules), c(5L*k, samples)) || !all(rules %in% c(0, 1)))
        stop(paste(method, "posterior rules"), call. = FALSE)
    issue199_posterior_link(parent, paste(method, "actual"))
    issue199_posterior_link(ref, paste(method, "reference"))
    issue199_numeric_patterns(p, q, method)
    ranks <- function(x) apply(x, 1L, rank, ties.method = "min", na.last = "keep")
    if (!identical(ranks(p), ranks(q)) ||
        !identical(max.col(p, ties.method = "first"), max.col(q, ties.method = "first")))
        stop(paste(method, "posterior ranks/ties/argmax"), call. = FALSE)
    terms <- issue199_posterior_terms(nb, rules)
    domain <- issue199_posterior_domain(terms, q, method)
    normal <- is.finite(p) & p >= .Machine$double.xmin
    ceiling <- if (method == "AIMS") 2^-35 else 2^-33
    if (!all(abs(p[normal] - q[normal]) <= ceiling*pmax(p[normal], q[normal])))
        stop(paste(method, "posterior bound exceeded"), call. = FALSE)
    invisible(domain)
}

# Observe the ordered selection at the production calibration boundary. The
# independent replay below uses only the mapped input, model and R primitives.
issue199_itr_protocol <- function(dataset) {
    list(dataset = dataset, mapping_seed = 1992026L, auto_seed = 20260831L,
        iterations = 100L, seed = 118L, ratio = 54/64, er_subsample = "target",
        kind = c("Mersenne-Twister", "Inversion", "Rejection"),
        runtime = R.version.string, platform = R.version$platform,
        scipen = getOption("scipen"), OutDec = getOption("OutDec"),
        collate = Sys.getlocale("LC_COLLATE"))
}

issue199_itr_capture <- function(mapped, dataset, calculate) {
    protocol <- issue199_itr_protocol(dataset)
    # The historical protocol specifies these components only. The audited
    # producer/replay paths use uniform sample() without replacement, not
    # rbinom() or rmultinom(); leave the runtime's binom.kind unchanged.
    if (!identical(RNGkind()[seq_along(protocol$kind)], protocol$kind))
        stop("cIHC.itr producer RNG kinds", call. = FALSE)
    original <- BreastSubtypeR:::docalibration
    draws <- medians <- list()
    observed <- NULL
    observe <- function(y, medians.all, calibration, internal = NA, external = NA) {
        frame <- parent.frame()
        calibrated <- original(y, medians.all, calibration, internal, external)
        if (exists("itr", frame, inherits = FALSE) && exists("keep", frame, inherits = FALSE)) {
            if (!identical(calibration, "Internal") || !identical(internal, "ER.mdns"))
                stop("cIHC.itr producer calibration", call. = FALSE)
            i <- length(draws) + 1L
            draws[[i]] <<- as.character(frame$keep)
            names(draws)[i] <<- as.character(frame$itr)
            if (i == 1L) observed <<- y
            if (!issue199_itr_exact(y, observed))
                stop("cIHC.itr producer assay changed during iterations", call. = FALSE)
            genes <- rownames(calibrated)
            middle <- setNames(medians.all[genes, internal], genes)
            medians[[i]] <<- middle
            if (!issue199_itr_exact(calibrated, sweep(as.matrix(y)[genes, , drop = FALSE],
                    1L, middle, "-")))
                stop("cIHC.itr producer calibrated expression", call. = FALSE)
        }
        calibrated
    }
    testthat::local_mocked_bindings(docalibration = observe, .package = "BreastSubtypeR")
    result <- calculate(mapped)
    if (!identical(names(draws), paste0("itr.", seq_len(protocol$iterations))))
        stop("cIHC.itr producer iteration sequence", call. = FALSE)
    list(data_input = mapped, res = result, protocol = protocol, selections = draws,
        observed_assay = observed, medians = medians)
}

issue199_itr_exact <- function(a, b) {
    identical(a, b, num.eq = FALSE, single.NA = FALSE, attrib.as.set = FALSE)
}

issue199_itr_paths <- function() {
    prefix <- "/res/results/cIHC.itr/"
    paste0(prefix, c("outList/distances", "outList/dist.RORSubtype",
        paste0("score.ROR/", c("Basal", "Her2", "LumA", "LumB", "Normal",
            "ROR-S (Subtype Only)", "ROR-P (Subtype + Prolif)"))))
}

issue199_itr_model <- function() {
    f <- system.file("data", "BreastSubtypeRobj.RData", package = "BreastSubtypeR")
    issue199_itr_model_file(f)
    e <- new.env(parent = emptyenv())
    load(f, envir = e)
    e$BreastSubtypeRobj
}

issue199_itr_model_file <- function(f) {
    if (!nzchar(f) || !file.exists(f) || !identical(unname(tools::md5sum(f)), "d2600fffe228ceec2e1a47b586747543"))
        stop("cIHC.itr centroid/median provenance", call. = FALSE)
    invisible(TRUE)
}

issue199_itr_arithmetic <- function() {
    p <- .Machine$longdouble.digits
    if (.Machine$double.base != 2L || .Machine$double.digits != 53L ||
        .Machine$double.rounding != 5L || .Machine$double.min.exp != -1022L ||
        .Machine$double.max.exp != 1024L || .Machine$sizeof.pointer != 8L ||
        (!is.null(p) && p != 53L && p < 64L) || 2^-1074 == 0)
        stop("cIHC.itr arithmetic model", call. = FALSE)
    if (!identical(getOption("OutDec"), ".") ||
        !identical(as.character(c(1/3, 2/3, -0, .1 + .2)),
            c("0.333333333333333", "0.666666666666667", "0", "0.3")))
        stop("cIHC.itr formatter semantics", call. = FALSE)
}

# Check the actual R token at every distinct contributing value. The integer
# mantissa and 10^k (k <= 21) are exactly representable. Recover the low
# product using split operands before testing the decimal grid.
issue199_itr_roundtrip <- function(x) {
    if (!all(is.finite(x) & abs(x) <= 1)) stop("cIHC.itr correlation domain", call. = FALSE)
    z <- unique(x)
    token <- as.character(z)
    parsed <- issue199_itr_decimal(z, token)
    parsed[match(x, z)]
}

issue199_itr_decimal <- function(z, token) {
    if (length(z) != length(token) || !all(grepl("^-?[0-9]+(\\.[0-9]+)?(e[+-]?[0-9]+)?$", token)))
        stop("cIHC.itr decimal token", call. = FALSE)
    parsed <- as.numeric(token)
    if (!all(is.finite(parsed) & abs(parsed) <= 1)) stop("cIHC.itr decimal parse", call. = FALSE)
    special <- abs(z) %in% c(0, 1)
    if (!all(parsed[special] == z[special])) stop("cIHC.itr decimal endpoint", call. = FALSE)
    v <- abs(z[!special])
    if (length(v)) {
        if (any(v < 2^-20)) stop("cIHC.itr decade-boundary domain", call. = FALSE)
        e <- floor(log10(v))
        # Exact positive integer powers plus an exclusion band much larger
        # than binary64 multiplication error certify floor(log10(v)).
        if (!all(e >= -7 & e <= -1)) stop("cIHC.itr decimal exponent/grid", call. = FALSE)
        normalized <- v * 10^(-e)
        if (!all(normalized > 1 + 2^-40 & normalized < 10*(1 - 2^-40)))
            stop("cIHC.itr decimal exponent/grid", call. = FALSE)
        t <- sub("^-", "", token[!special])
        mant <- sub("e.*$", "", t)
        exponent <- integer(length(t))
        sci <- grepl("e", t, fixed = TRUE)
        exponent[sci] <- as.integer(sub("^.*e", "", t[sci]))
        frac <- nchar(sub("^[^.]*\\.?", "", mant))
        digits <- sub("^0+", "", gsub(".", "", mant, fixed = TRUE))
        places <- exponent - frac
        if (!all(nchar(digits) <= 15L & places >= -21 & places <= 0))
            stop("cIHC.itr decimal significand/grid", call. = FALSE)
        integer_digits <- as.numeric(digits)
        grid <- places - e + 14
        if (!all(grid >= 0 & grid <= 14)) stop("cIHC.itr decimal exponent/grid", call. = FALSE)
        alpha <- integer_digits * 10^grid
        # Reconstruct the 15-digit grid integer. Multiplication residual via
        # Dekker splitting (bounded normal operands) proves Lemma F's limit.
        scale <- 10^(14-e)
        product_error <- function(a, b, integer) {
            # Error-free product expansion on this normal, bounded domain.
            split <- function(x) { t <- (2^27 + 1)*x; t - (t-x) }
            ah <- split(a); al <- a-ah
            bh <- split(b); bl <- b-bh
            product <- a*b
            first <- product - ah*bh
            second <- first - al*bh
            third <- second - ah*bl
            (product-integer) + (al*bl-third)
        }
        formatter_error <- product_error(v, scale, alpha)
        precision <- .Machine$longdouble.digits
        u <- if (is.null(precision)) 2^-53 else 2^-precision
        if (!all(abs(formatter_error) <= .5 + 1e15*(2*u + u*u)))
            stop("cIHC.itr decimal formatter error", call. = FALSE)
        parser_error <- product_error(abs(parsed[!special]), 10^(-places), integer_digits)
        # Compare to the exact rational decimal, including double rounding.
        if (!all(sign(parsed[!special]) == sign(z[!special]) &
                abs(parser_error) <= (2^-53 * (1 + 2^-10)) * integer_digits))
            stop("cIHC.itr decimal parser error", call. = FALSE)
    }
    parsed
}

issue199_itr_selections <- function(se, dataset) {
    cd <- as.data.frame(SummarizedExperiment::colData(se))
    ids <- as.character(cd$PatientID)
    if (!identical(ids, colnames(se)) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
        stop("cIHC.itr input sample identities", call. = FALSE)
    pos <- ids[cd$ER == "ER+"]; neg <- ids[cd$ER == "ER-"]
    oslo <- identical(dataset, "OSLO2EMIT0obj")
    counts <- if (oslo) c(84L, 18L) else c(74L, 26L)
    if (anyNA(c(pos, neg)) || !identical(c(length(pos), length(neg)), counts))
        stop("cIHC.itr ER membership/order", call. = FALSE)
    size <- if (oslo) 16L else 22L
    draws <- issue199_rng(118, lapply(seq_len(100L), function(i)
        c(pos[sample(length(pos), size)], neg)))
    setNames(draws, paste0("itr.", seq_len(100L)))
}

issue199_itr_ranks <- function(y, scale) {
    gap <- apply(y, 2L, function(v) min(diff(sort(v))))
    if (!all(is.finite(gap) & gap > 2^-40 * (1 + scale)))
        stop("cIHC.itr per-iteration rank gap", call. = FALSE)
    apply(y, 2L, rank, ties.method = "min")
}

issue199_itr_replay <- function(se, model, draws) {
    mat <- as.matrix(SummarizedExperiment::assay(se))
    if (!is.double(mat) || !all(is.finite(mat)) || anyDuplicated(rownames(mat)))
        stop("cIHC.itr expression domain", call. = FALSE)
    genes <- sort(intersect(rownames(model$centroid), model$medians$X))
    if (length(genes) != 50L || !identical(sort(rownames(mat)), genes))
        stop("cIHC.itr input gene alignment", call. = FALSE)
    mat <- mat[genes, , drop = FALSE]
    cent <- model$centroid[genes, , drop = FALSE]
    use4 <- !genes %in% c("BIRC5", "CCNB1", "GRB7", "MYBL2")
    if (sum(use4) != 46L || any(!is.finite(as.matrix(cent))) ||
        any(apply(cent, 2L, anyDuplicated) != 0L))
        stop("cIHC.itr centroid rank domain", call. = FALSE)
    ids <- colnames(mat); classes <- colnames(cent)
    if (!identical(classes, c("Basal", "Her2", "LumA", "LumB", "Normal")))
        stop("cIHC.itr centroid classes", call. = FALSE)
    measure <- function(y, centre) {
        d <- vapply(seq_len(ncol(centre)), function(k)
            apply(y, 2L, function(v) -stats::cor(centre[, k], v, method = "spearman",
                use = "pairwise.complete.obs")), numeric(ncol(y)))
        dimnames(d) <- list(ids, colnames(centre))
        d
    }
    n <- ncol(mat)
    iterations <- names(draws)
    if (!identical(iterations, paste0("itr.", seq_len(100L))))
        stop("cIHC.itr ordered iteration sequence", call. = FALSE)
    votes <- matrix("", n, 100L, dimnames = list(ids, iterations))
    five <- four <- ranks5 <- ranks4 <- medians <- vector("list", 100L)
    # A global scale bounds both selected middle operands and each
    # subtraction, on this input, through all iterations and columns.
    bound <- max(abs(mat))
    for (i in seq_along(draws)) {
        selected <- mat[, draws[[i]], drop = FALSE]
        middle <- apply(selected, 1L, sort)
        k <- nrow(middle) %/% 2L
        if (!all(abs(middle[c(k, k+1L), , drop = FALSE]) >= 2^-400) ||
            bound > 2^400) stop("cIHC.itr middle-median domain", call. = FALSE)
        med <- apply(selected, 1L, median)
        medians[[i]] <- med
        y <- sweep(mat, 1L, med, "-")
        ranks5[[i]] <- issue199_itr_ranks(y, 2*bound)
        ranks4[[i]] <- issue199_itr_ranks(y[use4, , drop = FALSE], 2*bound)
        five[[i]] <- measure(y, cent)
        four[[i]] <- measure(y[use4, , drop = FALSE], cent[use4, 1:4, drop = FALSE])
        votes[, i] <- classes[apply(five[[i]], 1L, which.min)]
        td <- if (i == 1L) y else td + y
    }
    consensus <- apply(votes, 1L, function(v) { freq <- table(v); names(freq)[which.max(freq)] })
    keep <- votes == consensus
    counts <- rowSums(keep)
    if (anyNA(keep) || any(counts < 1 | counts > 100))
        stop("cIHC.itr votes/masks/counts", call. = FALSE)
    aggregate <- function(values) {
        contributing <- matrix(0L, nrow(values[[1L]]), ncol(values[[1L]]))
        for (i in seq_along(values)) {
            a <- values[[i]]
            rows <- keep[, i]
            a[!rows, ] <- 0
            a[rows, ] <- issue199_itr_roundtrip(as.vector(a[rows, , drop = FALSE]))
            contributing[rows, ] <- contributing[rows, , drop = FALSE] + is.finite(a[rows, , drop = FALSE])
            total <- if (i == 1L) a else total+a
        }
        issue199_itr_counts(contributing, counts)
        -(total/contributing)
    }
    list(selections = draws, ranks5 = ranks5, ranks4 = ranks4, votes = votes,
        medians = medians, counts = counts, consensus = consensus, D = aggregate(five), D4 = aggregate(four), TD = td/100)
}

issue199_itr_counts <- function(cell_counts, votes) {
    if (length(votes) != nrow(cell_counts) || anyNA(cell_counts) || anyNA(votes) ||
        any(votes < 1L | votes > 100L) || !all(cell_counts == votes))
        stop("cIHC.itr cell contribution counts", call. = FALSE)
    invisible(TRUE)
}

issue199_itr_structure <- function(a, b, dataset) {
    exact <- issue199_itr_exact
    s <- a$score.ROR; r <- b$score.ROR
    same <- function(a, b, fields) all(vapply(fields, function(k) exact(a[[k]], b[[k]]), logical(1)))
    fields <- c("patientID", "Call", "Confidence", "ER", "HER2",
        "ROR-S Group (Subtype Only)", "ROR-P Group (Subtype + Prolif)")
    ids <- as.character(s$patientID)
    n <- if (dataset == "OSLO2EMIT0obj") 102L else 100L
    d <- a$outList$distances; d4 <- a$outList$dist.RORSubtype; td <- a$outList$testData
    if (!identical(dim(d), c(n, 5L)) || !identical(dim(d4), c(n, 4L)) ||
        !identical(dim(td), c(50L, n)) || !identical(ids, rownames(d)) ||
        !identical(ids, rownames(d4)) || !identical(ids, colnames(td)) ||
        !identical(ids, names(a$outList$predictions)) || !exact(names(a), names(b)) || !exact(names(a$outList), names(b$outList)) ||
        !exact(attributes(s), attributes(r)) || !same(s, r, fields) ||
        !exact(a$BS.all, b$BS.all) || !same(a$outList, b$outList, c("predictions", "centroids")) ||
        !all(vapply(c("distances", "dist.RORSubtype", "testData"), function(k)
            exact(attributes(a$outList[[k]]), attributes(b$outList[[k]])), logical(1))) ||
        !exact(s$ER, unname(a$outList$testData["ESR1", ])) ||
        !exact(s$HER2, unname(a$outList$testData["ERBB2", ])))
        stop(paste(dataset, "cIHC.itr structure/identities"), call. = FALSE)
}

issue199_itr_votes <- function(a, b) {
    votes <- a$BS.itr.keep
    ids <- names(a$outList$predictions)
    if (!issue199_itr_exact(votes, b$BS.itr.keep) || anyNA(votes) ||
        !identical(dimnames(votes), list(ids, paste0("itr.", seq_len(100L)))))
        stop("cIHC.itr votes/masks/counts", call. = FALSE)
    calls <- apply(votes, 1L, function(v) { n <- table(v); names(n)[which.max(n)] })
    if (!issue199_itr_exact(calls, a$outList$predictions))
        stop("cIHC.itr votes/masks/counts", call. = FALSE)
}

issue199_itr_links <- function(o) {
    d <- o$outList$distances
    for (k in colnames(d)) {
        if (!issue199_itr_exact(o$score.ROR[[k]], unname(d[, k])))
            stop("cIHC.itr exact score.ROR copies", call. = FALSE)
    }
}

issue199_itr_raw_s <- function(d) {
    unname(.04210193*d[, 1L] + .12466938*d[, 2L] - .35235561*d[, 3L] + .14213283*d[, 4L])
}

issue199_itr_cutoffs <- function(raw, other, score, reference, group, cuts, method) {
    if (!all(is.finite(c(raw, other))) ||
        !all(vapply(cuts, function(k) all(abs(c(raw, other)-k) > 2^-40), logical(1))))
        stop(paste(method, "raw cutoff margin"), call. = FALSE)
    sides <- function(a, b, cuts) identical(outer(a, cuts, "<"), outer(b, cuts, "<")) &&
        identical(outer(a, cuts, ">"), outer(b, cuts, ">"))
    if (!sides(raw, other, cuts)) stop(paste(method, "raw cutoff sides"), call. = FALSE)
    if (!identical(group, ifelse(raw < cuts[1L], "low", ifelse(raw > cuts[2L], "high", "med"))))
        stop(paste(method, "risk groups"), call. = FALSE)
    if (!sides(score, reference, 100*(cuts+.35)/.85))
        stop(paste(method, "export cutoff sides"), call. = FALSE)
}

# Numerical unit controls use this function directly. The scientific route
# additionally requires a matching producer and independent bitwise replay.
issue199_itr_bounds <- function(a, b, dataset) {
    for (field in c("distances", "dist.RORSubtype")) {
        x <- a$outList[[field]]; y <- b$outList[[field]]
        issue199_numeric_patterns(x, y, paste("cIHC.itr", field))
        if (!all(is.finite(c(x,y)) & abs(c(x,y)) <= 1) || !all(abs(x-y) <= 2^-45))
            stop(paste("cIHC.itr", field, "root bound exceeded"), call. = FALSE)
        ranks <- function(m) t(apply(m, 1L, rank, ties.method = "min"))
        if (!identical(ranks(x), ranks(y)) ||
            !identical(max.col(x, ties.method = "first"), max.col(y, ties.method = "first")))
            stop("cIHC.itr distance ordering", call. = FALSE)
    }
    x <- a$outList$testData; y <- b$outList$testData
    issue199_numeric_patterns(x, y, "cIHC.itr parent")
    if (dataset == "TCGABRCAobj")
        issue199_portable_compare(x, y, dataset, path = "/res/results/cIHC.itr/outList/testData")
    else if (!issue199_itr_exact(x, y)) stop("cIHC.itr exact OSLO parent", call. = FALSE)
    if (!identical(apply(x, 2L, rank), apply(y, 2L, rank))) stop("cIHC.itr parent ranks", call. = FALSE)
    da <- a$outList$distances; db <- b$outList$distances
    pa <- issue199_ror_domain(x, da, "cIHC.itr actual")
    pb <- issue199_ror_domain(y, db, "cIHC.itr reference")
    for (method in c("ROR-S", "ROR-P")) {
        subtype <- method == "ROR-S"
        desc <- if (subtype) "(Subtype Only)" else "(Subtype + Prolif)"
        field <- paste(method, desc)
        sa <- a$score.ROR[[field]]; sb <- b$score.ROR[[field]]
        issue199_numeric_patterns(sa, sb, field)
        if (!all(is.finite(c(sa,sb)))) stop(paste(method, "nonfinite score"), call. = FALSE)
        raw <- if (subtype) issue199_itr_raw_s(da) else pa$raw
        ref <- if (subtype) issue199_itr_raw_s(db) else pb$raw
        group <- paste(method, "Group", desc)
        cuts <- c(if (subtype) -.15 else -.25, .1)
        issue199_itr_cutoffs(raw, ref, sa, sb, a$score.ROR[[group]], cuts, method)
        issue199_itr_cutoffs(ref, raw, sb, sa, b$score.ROR[[group]], cuts, method)
        ceiling <- if (subtype) 2^-38 else 2^-37
        if (!all(abs(sa-sb) <= ceiling)) stop(paste(method, "bound exceeded"), call. = FALSE)
    }
    invisible(TRUE)
}

issue199_itr_visits <- function(visited, paths) {
    if (!identical(sort(visited), sort(paths)) || anyDuplicated(visited))
        stop("unused or multiply visited cIHC.itr route", call. = FALSE)
    invisible(TRUE)
}

issue199_itr_trace <- function(a, b) {
    if (!identical(a$selections, b$selections)) stop("cIHC.itr ordered selections", call. = FALSE)
    if (!identical(a$ranks5, b$ranks5)) stop("cIHC.itr 50-gene iteration ranks", call. = FALSE)
    if (!identical(a$ranks4, b$ranks4)) stop("cIHC.itr 46-gene iteration ranks", call. = FALSE)
}

issue199_itr_context <- function(a, fixture, dataset, producer) {
    b <- fixture$res$results$cIHC.itr
    issue199_itr_structure(a, b, dataset)
    issue199_itr_votes(a, b)
    issue199_itr_links(a)
    issue199_itr_links(b)
    paths <- issue199_itr_paths()
    if (is.null(producer)) {
        if (!issue199_itr_exact(a, b)) stop("cIHC.itr missing paired producer", call. = FALSE)
        return(paths)
    }
    issue199_itr_arithmetic()
    if (!identical(names(producer), c("data_input", "res", "protocol", "selections", "observed_assay", "medians")) ||
        !identical(producer$protocol, issue199_itr_protocol(dataset)))
        stop("cIHC.itr unverified producer protocol", call. = FALSE)
    if (!issue199_itr_exact(a, producer$res$results$cIHC.itr))
        stop("cIHC.itr producer result mismatch", call. = FALSE)
    se <- producer$data_input$se_NC
    frozen <- fixture$data_input$se_NC
    if (!identical(SummarizedExperiment::colData(se), SummarizedExperiment::colData(frozen)))
        stop("cIHC.itr input ER membership/order or sample identities", call. = FALSE)
    if (!issue199_itr_exact(producer$observed_assay, SummarizedExperiment::assay(se)))
        stop("cIHC.itr producer assay mismatch", call. = FALSE)
    issue199_portable_compare(producer$data_input, fixture$data_input, dataset, path = "/data_input")
    actual_draws <- issue199_itr_selections(se, dataset)
    reference_draws <- issue199_itr_selections(frozen, dataset)
    if (!identical(producer$selections, actual_draws) || !identical(actual_draws, reference_draws))
        stop("cIHC.itr ordered selections", call. = FALSE)
    model <- issue199_itr_model()
    actual <- issue199_itr_replay(se, model, actual_draws)
    if (!issue199_itr_exact(actual$medians, producer$medians))
        stop("cIHC.itr producer iteration medians", call. = FALSE)
    reference <- issue199_itr_replay(frozen, model, reference_draws)
    issue199_itr_trace(actual, reference)
    if (!identical(actual$votes, a$BS.itr.keep) || !identical(reference$votes, b$BS.itr.keep) ||
        !identical(actual$consensus, a$outList$predictions) ||
        !identical(reference$consensus, b$outList$predictions))
        stop("cIHC.itr replay votes/masks/counts", call. = FALSE)
    for (pair in list(c("D","distances"), c("D4","dist.RORSubtype"), c("TD","testData"))) {
        v <- actual[[pair[1L]]]; field <- pair[2L]
        if (!issue199_itr_exact(v, a$outList[[field]]))
            stop(paste("cIHC.itr same-runtime replay", field), call. = FALSE)
        approved <- b$outList[[field]]
        value <- reference[[pair[1L]]]
        issue199_numeric_patterns(value, approved, paste("cIHC.itr approved replay", field))
        limit <- if (field != "testData") 2^-45 else if (dataset == "TCGABRCAobj") 2^-48 else 0
        if (!all(abs(value-approved) <= limit))
            stop(paste("cIHC.itr approved replay", field), call. = FALSE)
    }
    s <- unname(100*(issue199_itr_raw_s(actual$D)+.35)/.85)
    p <- issue199_ror_formula(actual$TD, actual$D)$score
    if (!issue199_itr_exact(s, a$score.ROR[["ROR-S (Subtype Only)"]]))
        stop("cIHC.itr same-runtime ROR-S formula", call. = FALSE)
    if (!issue199_itr_exact(p, a$score.ROR[["ROR-P (Subtype + Prolif)"]]))
        stop("cIHC.itr same-runtime ROR-P formula", call. = FALSE)
    issue199_itr_bounds(a, b, dataset)
    paths
}

issue199_portable_compare <- function(actual, expected, dataset,
    scope = "C_vs_approved_C", path = "", producer = NULL) {
    if (path %in% issue199_itr_paths())
        stop("cIHC.itr scientific comparison requires full parent context", call. = FALSE)
    keys <- c(AIMS = "20", sspbc = "24")
    posterior_paths <- paste0("/res/results/", names(keys), "/all.probs/", keys)
    if (path %in% posterior_paths)
        stop("Posterior scientific comparison requires full parent context", call. = FALSE)
    ror_paths <- checked_posteriors <- itr_paths <- character()
    if (scope == "C_vs_approved_C" && path %in% c("", "/res")) {
        fixture <- issue199_ror_approved(dataset)
        reference <- if (path == "") fixture else fixture$res
        if (!identical(expected, reference, num.eq = FALSE, single.NA = FALSE))
            stop(paste(dataset, "approved posterior reference provenance"), call. = FALSE)
        res <- if (path == "") actual[["res", exact = TRUE]] else actual
        if (!is.list(res) || !all(names(keys) %in% names(res$results)))
            stop("Missing scientific posterior methods or parent", call. = FALSE)
        ror_paths <- issue199_ror_context(res, fixture$res, dataset)
        if (sum(names(res$results) == "cIHC.itr") != 1L) stop("cIHC.itr route inventory", call. = FALSE)
        if (!is.null(producer) && path == "" &&
            !issue199_itr_exact(actual$data_input, producer$data_input))
            stop("cIHC.itr full input/producer mismatch", call. = FALSE)
        itr_paths <- issue199_itr_context(res$results$cIHC.itr, fixture, dataset, producer)
        issue199_posterior_arithmetic()
        models <- issue199_posterior_models()
        for (method in names(keys))
            issue199_posterior_method(res$results[[method]], fixture$res$results[[method]],
                models[[method]], method, dataset)
        checked_posteriors <- posterior_paths
    }
    issue199_compare_walk(actual, expected, dataset, scope, path, ror_paths, checked_posteriors, itr_paths)
}

# Explicit entry for synthetic tests of the original 19-row policy and B2 links.
issue199_portability_unit <- function(actual, expected, dataset,
    scope = "C_vs_approved_C", path = "") {
    issue199_compare_walk(actual, expected, dataset, scope, path, character(), character())
}

issue199_compare_walk <- function(actual, expected, dataset, scope, path, ror_paths, checked_posteriors, itr_paths = character()) {
    ror_leaf <- paste0("/res/results/", c("genefu.robust", "ssBC", "ssBC.v2", "cIHC", "PCAPAM50"),
        "/score.ROR/ROR-P (Subtype + Prolif)")
    if (path %in% ror_leaf) stop("ROR scientific comparison requires full parent context", call. = FALSE)
    ror_visited <- itr_visited <- character()
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
        if (p %in% itr_paths) {
            itr_visited <<- c(itr_visited, p)
            if (nrow(rule)) visited <<- c(visited, p)
        } else if (p %in% ror_paths) {
            ror_visited <<- c(ror_visited, p)
            if (nrow(rule)) visited <<- c(visited, p)
        } else if (p %in% checked_posteriors) {
            if (nrow(rule) != 1L) fail(p, "missing posterior policy row")
            visited <<- c(visited, p)
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
    if (!identical(sort(visited[visited %in% checked_posteriors]), sort(checked_posteriors)))
        fail(path, "unused or multiply visited posterior route")
    issue199_itr_visits(itr_visited, itr_paths)
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
