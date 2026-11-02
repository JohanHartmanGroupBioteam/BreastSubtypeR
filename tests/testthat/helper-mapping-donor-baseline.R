# Frozen donor-preparation oracle from f878f58c1c12ca6f1cd5bf4489667935bf25969b.
# Keep this implementation independent of the optimized production helper.
mapping219_baseline <- function(x, y, method, required) {
    y <- y[!is.na(y$ENTREZID), , drop = FALSE]
    groups <- split(seq_len(nrow(y)),
        factor(y$ENTREZID, levels = unique(y$ENTREZID)))
    values <- list()
    excluded <- data.frame(gene = character(), reason = character(),
        stringsAsFactors = FALSE)
    for (id in names(groups)) {
        annotation <- y[groups[[id]], , drop = FALSE]
        collapse <- function() duplicate_genes(
            x[annotation$probe, , drop = FALSE], annotation, method)
        if (id %in% required) {
            # Never drop a protected row or conceal its duplicate-resolution error.
            values[[id]] <- collapse()
            next
        }
        gene <- tryCatch(collapse(), error = function(e) {
            if (method %in% c("stdev", "iqr") &&
                inherits(e, "BreastSubtypeR_duplicate_undefined")) return(e)
            # Preserve recognition of the existing max error by its full
            # message. All unrelated failures must still propagate.
            expected <- sprintf(paste0(
                "Cannot resolve Entrez ID %s with method = \"max\": ",
                "none of its %d probes has a defined total across samples ",
                "(all row sums are NA/NaN). Review missing or nonfinite ",
                "measurements and the preprocessing strategy. ",
                "No probe was selected for this Entrez ID."), id, nrow(annotation))
            if (identical(method, "max") &&
                identical(conditionMessage(e), expected)) return(e)
            stop(e)
        })
        reason <- if (inherits(gene, "error")) "duplicate_resolution_failed" else
            if (any(!is.finite(gene) & !is.na(gene))) "input_infinity" else
            if (all(is.na(gene))) "entirely_missing" else NULL
        if (is.null(reason)) {
            # No missing-fraction cutoff: otherwise usable partial rows remain.
            values[[id]] <- gene
        } else {
            excluded[nrow(excluded) + 1L, ] <- list(id, reason)
        }
    }
    if (nrow(excluded)) warning(structure(list(
        message = paste0("Optional KNN excluded unusable background genes: ",
            paste0(excluded$gene, " (", excluded$reason, ")", collapse = "; "),
            ". Protected genes were retained."),
        call = NULL, exclusions = excluded),
        class = c("BreastSubtypeR_knn_donor_exclusion", "warning", "condition")))
    matrix <- if (length(values)) do.call(rbind, unname(values)) else
        matrix(numeric(), nrow = 0L, ncol = ncol(x),
            dimnames = list(NULL, colnames(x)))
    list(matrix = matrix, exclusions = excluded)
}

environment(mapping219_baseline) <- asNamespace("BreastSubtypeR")
