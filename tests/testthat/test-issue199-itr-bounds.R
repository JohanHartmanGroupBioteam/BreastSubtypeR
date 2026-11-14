# Use the public AUTO producer on the packaged input. No permission to use a
# bound is cached: each scientific comparison authenticates its own evidence.
itr33_pair <- function(fixture, dataset) {
    calculate <- function(input) suppressMessages(BS_Multi(input, methods = "AUTO",
        Subtype = FALSE, hasClinical = FALSE, er_subsample = "target"))
    withr::with_seed(20260831, issue199_itr_capture(fixture$data_input, dataset, calculate),
        .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
}

itr33_score <- function(result, field, delta) {
    attrs <- attributes(result$score.ROR)
    result$score.ROR[[field]][1L] <- result$score.ROR[[field]][1L] + delta
    attributes(result$score.ROR) <- attrs
    result
}

test_that("scientific iteration routes need independently replayed producer evidence", {
    ds <- "TCGABRCAobj"
    f <- issue199_ror_approved(ds)
    p <- itr33_pair(f, ds)
    o <- p$res$results$cIHC.itr
    b <- f$res$results$cIHC.itr
    check <- function(w, value = o) issue199_itr_context(value, f, ds, w)
    routes <- issue199_itr_paths()
    expect_identical(check(p), routes)
    expect_identical(issue199_itr_context(b, f, ds, NULL), routes)
    expect_length(unique(routes), 9L)
    expect_true(issue199_itr_visits(routes, routes))
    expect_error(issue199_itr_visits(routes[-1L], routes), "^unused or multiply visited cIHC.itr route$")
    expect_error(issue199_itr_visits(c(routes, routes[1L]), routes), "^unused or multiply visited cIHC.itr route$")
    expect_error(issue199_portable_compare(0, 0, ds, path = routes[1L]), "requires full parent context")

    field <- "ROR-P (Subtype + Prolif)"
    ulp <- function(x) 2^(floor(log2(abs(x)))-52)
    a <- itr33_score(o, field, ulp(o$score.ROR[[field]][1L]))
    expect_false(issue199_itr_exact(a, o))
    expect_error(check(NULL, a), "^cIHC.itr missing paired producer$")
    expect_error(check(p, a), "^cIHC.itr producer result mismatch$")
    q <- p
    q$res$results$cIHC.itr <- a
    expect_error(check(q, a), "^cIHC.itr same-runtime ROR-P formula$")
    q <- p
    q$protocol$seed <- 117L
    expect_error(check(q), "^cIHC.itr unverified producer protocol$")
    q <- p
    q$data_input$se_NC <- q$data_input$se_NC[, rev(seq_len(ncol(q$data_input$se_NC)))]
    expect_error(check(q), "^cIHC.itr input ER membership/order or sample identities$")
    q <- p
    cd <- SummarizedExperiment::colData(q$data_input$se_NC)
    j <- c(which(cd$ER == "ER+")[1L], which(cd$ER == "ER-")[1L])
    cd$ER[j] <- rev(cd$ER[j])
    SummarizedExperiment::colData(q$data_input$se_NC) <- cd
    expect_error(check(q), "^cIHC.itr input ER membership/order or sample identities$")
    q <- p
    x <- SummarizedExperiment::assay(q$data_input$se_NC)
    x[1L, 1L] <- x[1L, 1L] + 2^-40
    SummarizedExperiment::assay(q$data_input$se_NC) <- x
    expect_error(check(q), "^cIHC.itr producer assay mismatch$")
    expect_error(issue199_portable_compare(f, f, ds, producer = q), "^cIHC.itr full input/producer mismatch$")
    q <- p
    q$medians[[1L]][1L] <- q$medians[[1L]][1L] + 2^-40
    expect_error(check(q), "^cIHC.itr producer iteration medians$")
    for (order in list(c(10L, 5L), c(5L, 5L))) {
        q <- p
        q$selections[c(5L, 10L)] <- p$selections[order]
        expect_error(check(q), "^cIHC.itr ordered selections$")
    }
    q <- p
    q$selections <- q$selections[-1L]
    expect_error(check(q), "^cIHC.itr ordered selections$")
})

test_that("identical votes and bounded averages do not identify iteration ranks", {
    for (ds in c("OSLO2EMIT0obj", "TCGABRCAobj")) {
        f <- issue199_ror_approved(ds)
        se <- f$data_input$se_NC
        model <- issue199_itr_model()
        draw <- issue199_itr_selections(se, ds)
        a <- issue199_itr_replay(se, model, draw)
        pair <- if (ds == "OSLO2EMIT0obj") c(6L, 9L) else c(5L, 10L)
        expect_false(setequal(draw[[pair[1L]]], draw[[pair[2L]]]))
        draw[pair] <- rev(draw[pair])
        b <- issue199_itr_replay(se, model, draw)
        expect_identical(a$votes, b$votes)
        expect_identical(a$consensus, b$consensus)
        expect_identical(a$counts, b$counts)
        expect_true(max(abs(a$D-b$D), abs(a$D4-b$D4)) <= 2^-45)
        expect_true(max(abs(a$TD-b$TD)) <= 2^-48)
        expect_error(issue199_itr_trace(a, b), "^cIHC.itr ordered selections$")
        # Corrupt ranks alone after checking that the real swap changed both.
        expect_false(identical(a$ranks5, b$ranks5))
        expect_false(identical(a$ranks4, b$ranks4))
        q <- a; q$ranks5 <- b$ranks5
        expect_error(issue199_itr_trace(a, q), "^cIHC.itr 50-gene iteration ranks$")
        q <- a; q$ranks4 <- b$ranks4
        expect_error(issue199_itr_trace(a, q), "^cIHC.itr 46-gene iteration ranks$")
        expect_error(issue199_itr_replay(se, model, rev(draw)), "^cIHC.itr ordered iteration sequence$")
    }
})

test_that("decimal grid, middle operands, gene alignment and masks fail closed", {
    expect_identical(issue199_itr_roundtrip(c(0, 1, -1, .125)), c(0, 1, -1, .125))
    expect_error(issue199_itr_roundtrip(0.1*(1+2^-45)), "^cIHC.itr decimal exponent/grid$")
    expect_error(issue199_itr_decimal(.125, "0.125000000000001"), "^cIHC.itr decimal formatter error$")
    expect_error(issue199_itr_decimal(.125, "-0.125"), "^cIHC.itr decimal parser error$")
    expect_error(issue199_itr_roundtrip(Inf), "^cIHC.itr correlation domain$")
    withr::with_options(list(OutDec = ","), expect_error(issue199_itr_arithmetic(), "^cIHC.itr formatter semantics$"))
    se <- issue199_ror_approved("TCGABRCAobj")$data_input$se_NC
    model <- issue199_itr_model()
    draw <- issue199_itr_selections(se, "TCGABRCAobj")
    expect_error(issue199_itr_replay(se[-1L, ], model, draw), "^cIHC.itr input gene alignment$")
    x <- SummarizedExperiment::assay(se)
    x[1L, ] <- 0
    SummarizedExperiment::assay(se) <- x
    expect_error(issue199_itr_replay(se, model, draw), "^cIHC.itr middle-median domain$")
    expect_error(issue199_itr_ranks(matrix(c(.5,.5), 2L), 1), "^cIHC.itr per-iteration rank gap$")
    expect_true(issue199_itr_counts(matrix(c(75L, 100L), 2L, 4L), c(75L,100L)))
    expect_error(issue199_itr_counts(matrix(74L, 1L, 5L), 75L), "^cIHC.itr cell contribution counts$")
    expect_error(issue199_itr_counts(matrix(NA_integer_), 1L), "^cIHC.itr cell contribution counts$")
    b <- issue199_ror_approved("TCGABRCAobj")$res$results$cIHC.itr
    a <- b; a$BS.itr.keep[1L, 1L] <- NA_character_
    expect_error(issue199_itr_votes(a,b), "^cIHC.itr votes/masks/counts$")
    a <- b; a$outList$distances[1L] <- Inf
    expect_error(issue199_itr_bounds(a,b,"TCGABRCAobj"), "nonfinite masks/signs$")
    expect_error(issue199_numeric_patterns(-0, 0, "root"), "signed zero$")
    expect_error(issue199_numeric_patterns(2^-1074, 2^-1073, "root"), "subnormal support/values$")
})

test_that("bounded numerical units enforce both roots, copies, scores and raw decisions", {
    ds <- "OSLO2EMIT0obj"; b <- issue199_ror_approved(ds)$res$results$cIHC.itr
    for (k in c("distances", "dist.RORSubtype")) {
        a <- b; a$outList[[k]][1L] <- a$outList[[k]][1L] + 2^-44
        expect_error(issue199_itr_bounds(a,b,ds), paste0("^cIHC.itr ", k, " root bound exceeded$"))
    }
    a <- itr33_score(b, "Basal", 2^-50)
    expect_error(issue199_itr_links(a), "^cIHC.itr exact score.ROR copies$")
    s <- "ROR-S (Subtype Only)"; p <- "ROR-P (Subtype + Prolif)"
    value <- b$score.ROR[[p]][1L]
    a <- itr33_score(b, p, 2^(floor(log2(abs(value)))-52))
    expect_false(issue199_itr_exact(a,b))
    expect_true(issue199_itr_bounds(a,b,ds))
    expect_error(issue199_itr_bounds(itr33_score(b,p,2^-36),b,ds), "^ROR-P bound exceeded$")
    expect_error(issue199_itr_bounds(itr33_score(b,s,2^-37),b,ds), "^ROR-S bound exceeded$")
    sides <- function(x,y,group) issue199_itr_cutoffs(x,y,20,20,group,c(-.25,.1),"ROR-P")
    expect_error(sides(.09,.11,"med"), "^ROR-P raw cutoff sides$")
    expect_error(sides(.09,.09,"high"), "^ROR-P risk groups$")
    expect_error(sides(.1,.1,"med"), "^ROR-P raw cutoff margin$")
    a <- b; attrs <- attributes(a$score.ROR)
    a$score.ROR[["ROR-P Group (Subtype + Prolif)"]][1L] <- "changed"
    attributes(a$score.ROR) <- attrs
    expect_error(issue199_itr_structure(a,b,ds), "^OSLO2EMIT0obj cIHC.itr structure/identities$")
})

test_that("iteration capture and replay preserve caller RNG and restore observers", {
    ds <- "OSLO2EMIT0obj"; f <- issue199_ror_approved(ds)
    binding <- BreastSubtypeR:::docalibration
    withr::local_preserve_seed()
    kind <- RNGkind()
    withr::defer(do.call(RNGkind, as.list(kind)))
    RNGkind("L'Ecuyer-CMRG", "Box-Muller", "Rejection")
    set.seed(104)
    seed <- .Random.seed; kinds <- RNGkind()
    p <- itr33_pair(f, ds)
    expect_identical(.Random.seed, seed)
    expect_identical(RNGkind(), kinds)
    expect_identical(BreastSubtypeR:::docalibration, binding)
    expect_identical(issue199_itr_context(p$res$results$cIHC.itr, f, ds, p), issue199_itr_paths())
    expect_identical(.Random.seed, seed)
    expect_identical(RNGkind(), kinds)
    rm(".Random.seed", envir = .GlobalEnv)
    issue199_itr_selections(f$data_input$se_NC, ds)
    expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
    fail <- function(input) stop("producer probe")
    expect_error(withr::with_seed(20260831, issue199_itr_capture(f$data_input, ds, fail),
        .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
        .rng_sample_kind = "Rejection"), "^producer probe$")
    expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
    expect_identical(BreastSubtypeR:::docalibration, binding)
})


test_that("root decisions and the measured parent contract detect independent mutations", {
    for (ds in c("OSLO2EMIT0obj", "TCGABRCAobj")) {
        f <- issue199_ror_approved(ds)
        b <- f$res$results$cIHC.itr
        a <- b
        a$outList$testData[1L] <- a$outList$testData[1L] + 2^-46
        diagnostic <- if (ds == "TCGABRCAobj") "elementwise numeric bound exceeded$" else "^cIHC.itr exact OSLO parent$"
        expect_error(issue199_itr_bounds(a,b,ds), diagnostic)
        a <- b
        a$outList$testData <- a$outList$testData[rev(seq_len(50L)), ]
        expect_error(issue199_itr_structure(a,b,ds), paste0("^", ds, " cIHC.itr structure/identities$"))
        a <- b; a$outList$extra <- 1
        expect_error(issue199_itr_structure(a,b,ds), paste0("^", ds, " cIHC.itr structure/identities$"))
        a <- b; a$outList$distances <- a$outList$distances[rev(seq_len(nrow(a$score.ROR))), ]
        expect_error(issue199_itr_structure(a,b,ds), paste0("^", ds, " cIHC.itr structure/identities$"))
        for (col in colnames(b$outList$distances)) {
            a <- itr33_score(b, col, 2^(floor(log2(abs(b$score.ROR[[col]][1L])))-52))
            expect_error(issue199_itr_links(a), "^cIHC.itr exact score.ROR copies$")
        }
        a <- b
        a$outList$distances[1L, 1L] <- a$outList$distances[1L, 1L] + 2^-46
        attrs <- attributes(a$score.ROR)
        a$score.ROR$Basal <- unname(a$outList$distances[, 1L])
        attributes(a$score.ROR) <- attrs
        expect_no_error(issue199_itr_links(a))
        expect_true(issue199_itr_bounds(a,b,ds))
        a <- b; a$BS.itr.keep[1L, 1L] <- "mutated vote"
        expect_error(issue199_itr_votes(a,b), "^cIHC.itr votes/masks/counts$")
        for (invalid in c(NA_real_, NaN, -Inf)) {
            a <- b; a$outList$dist.RORSubtype[1L] <- invalid
            expect_error(issue199_itr_bounds(a,b,ds), "nonfinite masks/signs$")
        }
    }
    file <- withr::local_tempfile()
    writeLines("changed model", file)
    expect_error(issue199_itr_model_file(file), "^cIHC.itr centroid/median provenance$")
    expect_error(issue199_itr_ranks(matrix(c(1,1+2^-45)),1), "^cIHC.itr per-iteration rank gap$")
})
