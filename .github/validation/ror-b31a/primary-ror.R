# Run on freshly computed installed-package outputs before the full AUTO walk.
for (expr in parse(file.path(root, "test-issue199-ror-bounds.R"))) {
    if (!(is.call(expr) && identical(expr[[1L]], as.name("test_that")))) eval(expr, globalenv())
}
base_policy <- new.env(parent = baseenv())
sys.source(file.path(here, "base-policy.R"), base_policy)
stopifnot(identical(issue199_portability_policy(), base_policy$issue199_portability_policy()))
ror_rows <- list(); ror_cells <- list(); ror_excluded <- list(); ror_checks <- list()
primary_ror <- function(actual, approved, dataset) {
    before <- .Random.seed; kind <- RNGkind()
    required <- c("genefu.robust", "ssBC", "ssBC.v2", "cIHC", "PCAPAM50")
    routes <- issue199_ror_routes()
    stopifnot(identical(routes$method[routes$dataset == dataset], required),
        nrow(routes) == 10L, !anyDuplicated(routes))
    paths <- issue199_ror_context(actual, approved, dataset)
    visited <- character()
    for (method in required) {
        a <- actual$results[[method]]; b <- approved$results[[method]]
        n <- issue199_ror_method(a, b, b, dataset, method)
        stopifnot(isTRUE(ror31_formula_check(a)))
        x <- a$outList$testData; d <- a$outList$distances
        formula <- issue199_ror_formula(x, d)
        score <- a$score.ROR[["ROR-P (Subtype + Prolif)"]]
        reference <- b$score.ROR[["ROR-P (Subtype + Prolif)"]]
        same <- vapply(seq_along(score), function(j) identical(score[j], formula$score[j],
            num.eq = FALSE, single.NA = FALSE), logical(1))
        delta <- abs(score - reference)
        stopifnot(n == length(score), all(same), all(delta <= 2^-37))
        ids <- colnames(x)
        excluded <- setdiff(rownames(approved$res_subtypes), ids)
        stopifnot(identical(excluded, rownames(approved$res_subtypes)[is.na(approved$res_subtypes[[method]])]))
        ror_rows[[length(ror_rows) + 1L]] <<- data.frame(dataset, method, n,
            max_abs = max(delta), ceiling = 2^-37, formula_identities = sum(same),
            parents_distances_structure_domain_decisions = TRUE,
            excluded = length(excluded), actual_sha256 = numeric_fingerprint(score),
            reference_sha256 = numeric_fingerprint(reference),
            formula_sha256 = numeric_fingerprint(formula$score))
        ror_cells[[length(ror_cells) + 1L]] <<- data.frame(dataset, method, sample_id = ids,
            score, reference, delta, within_ceiling = delta <= 2^-37, formula_identity = same,
            parent_delta = apply(abs(x - b$outList$testData), 2L, max),
            sd = formula$sd, m1 = formula$m1, raw_cutoff_margin = pmin(abs(formula$raw+.25), abs(formula$raw-.1)))
        if (length(excluded)) ror_excluded[[length(ror_excluded) + 1L]] <<-
            data.frame(dataset, method, sample_id = excluded)
        visited <- c(visited, paste0("/res/results/", method, "/score.ROR/ROR-P (Subtype + Prolif)"))
    }
    stopifnot(issue199_ror_visits(visited, paths), identical(before, .Random.seed), identical(kind, RNGkind()))
    ror_checks[[dataset]] <<- data.frame(dataset, exact_routes = TRUE, no_duplicate_missing_visits = TRUE,
        policy_19_unchanged = TRUE, RNG_preserved = TRUE, eligible = sum(vapply(ror_rows,
            function(x) if (x$dataset == dataset) x$n else 0L, integer(1))))
    write_tsv(do.call(rbind, ror_rows), file.path(out, "ror-primary-routes.tsv"))
    write_tsv(do.call(rbind, ror_cells), file.path(out, "ror-primary-998.tsv"))
    write_tsv(do.call(rbind, ror_excluded), file.path(out, "ror-excluded-identities.tsv"))
    write_tsv(do.call(rbind, ror_checks), file.path(out, "ror-primary-checks.tsv"))
    invisible(TRUE)
}
