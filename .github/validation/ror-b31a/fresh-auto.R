library(BreastSubtypeR)
library(testthat)
script <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1L])
here <- dirname(script)
source(file.path(here, "evidence-utils.R"))
root <- system.file("tests/testthat", package = "BreastSubtypeR")
sys.source(file.path(root, "helper-issue199-portability.R"), globalenv())
sys.source(file.path(root, "helper-auto-ssbc-contract.R"), globalenv())
stopifnot(normalizePath(find.package("BreastSubtypeR")) == normalizePath(file.path(Sys.getenv("R_LIBS"), "BreastSubtypeR")))
out <- file.path(Sys.getenv("EVIDENCE_DIR"), "auto")
dir.create(out, showWarnings = FALSE)
checks <- list(); numeric <- list(); structural <- list(); decision_results <- list()
record <- function(dataset, method, model, check, value) {
    checks[[length(checks) + 1L]] <<- data.frame(dataset, method, model, check, passed = isTRUE(value))
}
seeded <- function(seed, code) withr::with_preserve_seed({
    kind <- RNGkind()
    tryCatch({RNGkind("Mersenne-Twister", "Inversion", "Rejection");set.seed(seed);force(code)},
        finally = do.call(RNGkind, as.list(kind)))
})
first <- function(p) vapply(seq_len(nrow(p)), function(i)
    if (all(is.na(p[i, ]))) NA_integer_ else which.max(p[i, ]), integer(1))
row_ranks <- function(p) t(vapply(seq_len(nrow(p)), function(i)
    rank(-p[i, ], ties.method = "min", na.last = "keep"), numeric(ncol(p))))
tie_pattern <- function(p) lapply(seq_len(nrow(p)), function(i) outer(p[i, ], p[i, ], "=="))
source(file.path(here, "primary-ror.R"))
set.seed(199)
for (dataset in c("OSLO2EMIT0obj", "TCGABRCAobj")) {
    e <- new.env(parent = emptyenv());utils::data(list = dataset, package = "BreastSubtypeR", envir = e)
    obj <- e[[dataset]]
    before <- .Random.seed;kind <- RNGkind()
    mapped <- seeded(1992026, Mapping(obj$se_obj, RawCounts = dataset == "TCGABRCAobj", method = "max", impute = TRUE, verbose = FALSE))
    fresh <- seeded(20260831, BS_Multi(data_input = mapped, methods = "AUTO", Subtype = FALSE, hasClinical = FALSE, er_subsample = "target"))
    record(dataset, "AUTO", "all", "seed preserved by harness", identical(before, .Random.seed))
    record(dataset, "AUTO", "all", "RNG kind preserved", identical(kind, RNGkind()))
    record(dataset, "AUTO", "all", "complete categorical/entropy table exact", identical(fresh$res_subtypes, obj$res$res_subtypes, num.eq = FALSE))
    for (method in c("AIMS", "sspbc")) {
        a <- fresh$results[[method]];b <- obj$res$results[[method]]
        prior <- .Random.seed
        record(dataset, method, "all", "reviewed helper validates actual", tryCatch(issue199_posterior_link(a, method), error = function(e) FALSE))
        record(dataset, method, "all", "reviewed helper validates approved", tryCatch(issue199_posterior_link(b, method), error = function(e) FALSE))
        record(dataset, method, "all", "class matrix exact", identical(a$cl, b$cl))
        record(dataset, method, "all", "rules.matrix bit exact", identical(a$rules.matrix, b$rules.matrix, num.eq = FALSE))
        record(dataset, method, "all", "model names exact", identical(names(a$all.probs), names(b$all.probs)))
        record(dataset, method, "all", "selected probability attributes exact", identical(attributes(a$prob), attributes(b$prob), num.eq = FALSE))
        record(dataset, method, "all", "selected probability missingness exact", identical(is.na(a$prob), is.na(b$prob)))
        for (k in names(a$all.probs)) {
            for (side in c("actual", "approved")) {
                o <- if (side == "actual") a else b
                p <- o$all.probs[[k]];j <- first(p)
                pv <- rep(NA_real_, nrow(p));cl <- rep(NA_character_, nrow(p))
                kept <- which(!is.na(j))
                pv[kept] <- p[cbind(kept, j[kept])];cl[kept] <- colnames(p)[j[kept]]
                masked <- rowSums(is.na(p)) == ncol(p)
                valid <- is.matrix(p) && is.double(p) && nrow(p) == nrow(o$prob) &&
                    is.matrix(o$prob) && is.double(o$prob) && is.matrix(o$cl) && is.character(o$cl) &&
                    identical(dim(o$prob), dim(o$cl)) && identical(dimnames(o$prob), dimnames(o$cl)) &&
                    ncol(o$rules.matrix) == nrow(p) && !anyNA(colnames(p)) && !anyDuplicated(colnames(p)) &&
                    !any(is.nan(p)) && all(rowSums(is.na(p)) == 0L | masked) &&
                    all(is.finite(p[!masked, , drop = FALSE])) && all(p >= 0 & p <= 1, na.rm = TRUE)
                record(dataset, method, k, paste(side, "independent shapes, types and masks"), valid)
                record(dataset, method, k, paste(side, "which.max probability bit exact"), identical(unname(o$prob[, k]), pv, num.eq = FALSE))
                record(dataset, method, k, paste(side, "which.max class exact"), identical(unname(o$cl[, k]), cl))
            }
            p <- a$all.probs[[k]];q <- b$all.probs[[k]]
            record(dataset, method, k, "posterior attributes and axes exact", identical(attributes(p), attributes(q), num.eq = FALSE))
            record(dataset, method, k, "posterior masks exact", identical(is.na(p), is.na(q)) && identical(is.nan(p), is.nan(q)) && identical(is.infinite(p), is.infinite(q)))
            record(dataset, method, k, "independent posterior ranks exact", identical(row_ranks(p), row_ranks(q)))
            record(dataset, method, k, "all posterior ties exact", identical(tie_pattern(p), tie_pattern(q)))
            record(dataset, method, k, "independent first argmax exact", identical(first(p), first(q)))
        }
        record(dataset, method, "all", "posterior checks preserve RNG", identical(prior, .Random.seed) && identical(kind, RNGkind()))
    }
    compatible <- NULL
    reporter <- testthat::ListReporter$new();reporter$start_file(paste0(dataset, "-independent-decisions"))
    testthat::with_reporter(reporter, testthat::test_that("retained axes and all scientific decisions", {
        compatible <<- auto_ssbc_historical_view(fresh, obj$res)
        issue199_expect_decisions(compatible, obj$res)
    }))
    d <- save_test_results(reporter$get_results(), file.path(out, paste0(dataset, "-decisions")))
    decision_results[[dataset]] <- d
    if (is.null(compatible)) stop("Historical-axis projection could not be verified")
    primary_ror(compatible, obj$res, dataset)
    inv <- leaf_inventory(list(data_input = mapped, res = compatible),
        list(data_input = obj$data_input, res = obj$res), dataset, issue199_portability_policy())
    policy <- issue199_portability_policy()
    current <- policy[policy$scope == "C_vs_approved_C" & policy$dataset == dataset, ]
    current$visits <- vapply(current$path, function(p) sum(inv$leaves$path == p), integer(1))
    stopifnot(all(current$visits == 1L))
    write_tsv(current, file.path(out, paste0(dataset, "-policy-visits.tsv")))
    # Preserve old-policy diagnostics, then mark only successfully certified ROR routes.
    ror_paths <- issue199_ror_context(compatible, obj$res, dataset)
    for (p in ror_paths) {
        i <- which(inv$leaves$path == p)
        stopifnot(length(i) == 1L)
        inv$leaves$policy[i] <- "fixture_ROR_2^-37"
        inv$leaves$exceeded[i] <- 0L
        inv$leaves$accepted[i] <- TRUE
    }
    numeric[[dataset]] <- inv$leaves;structural[[dataset]] <- inv$structural
    write_tsv(do.call(rbind, checks), file.path(out, "posterior-checks.tsv"))
    write_tsv(do.call(rbind, numeric), file.path(out, "numeric-leaves.tsv"))
    write_tsv(do.call(rbind, structural), file.path(out, "structural-differences.tsv"))
}
checks <- do.call(rbind, checks);decisions <- do.call(rbind, decision_results)
write_tsv(issue199_portability_policy(), file.path(out, "policy.tsv"))
print(table(checks$passed))
stopifnot(all(checks$passed), !any(decisions$failed | decisions$error | decisions$warning | decisions$skipped),
    !nrow(do.call(rbind, structural)))

stopifnot(sum(do.call(rbind, ror_rows)$n) == 998L, sum(do.call(rbind, ror_rows)$formula_identities) == 998L)
