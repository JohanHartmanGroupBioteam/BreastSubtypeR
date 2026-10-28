# Private Mapping preparation and acceptance around the frozen B3 estimator.
# Donor policy, output projection and Mapping errors belong here, not in B3.

.mapping_required_ids <- function(reference) {
    models <- new.env(parent = emptyenv())
    data(list = c("AIMSmodel", "sspbc.models"), envir = models,
        package = "BreastSubtypeR")
    pairs <- c(models$AIMSmodel$all.pairs,
        models$sspbc.models$ssp.pam50$all.pairs,
        models$sspbc.models$ssp.subtype$all.pairs)
    unique(c(as.character(reference$genes.signature$EntrezGene.ID),
        as.character(reference$genes.sig50$EntrezGene.ID),
        sub("^e", "", unlist(strsplit(pairs, "<", fixed = TRUE),
            use.names = FALSE))))
}

.mapping_output_targets <- function(x, reference) {
    genes <- reference$genes.signature
    output_ids <- unique(c(as.character(reference$genes.sig50$EntrezGene.ID),
        as.character(genes$EntrezGene.ID[genes$SSP_based == "Yes"])))
    rownames(x)[rownames(x) %in% output_ids]
}

.mapping_knn_stop <- function(message, targets = character(),
    reasons = character(), diagnostics = NULL, exclusions = NULL,
    class = "BreastSubtypeR_knn_unsupported") {
    stop(structure(list(message = message, call = NULL,
        targets = unname(unique(targets)), reasons = reasons,
        diagnostics = diagnostics, exclusions = exclusions),
        class = c(class, "error", "condition")))
}

.mapping_donor_matrix <- function(x, y, method, required) {
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
            # duplicate_genes has one intentional unsupported-collapse error.
            # Its full message is matched so unrelated failures are propagated.
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

# Arithmetic-only guard, exposed privately so integer boundaries need no large
# allocations. All validation and multiplication here use double arithmetic.
.mapping_knn_capacity <- function(n_genes, n_samples, n_targets, n_computed,
    input_counts) {
    limit <- as.double(.Machine$integer.max)
    fail <- function(reason) .mapping_knn_stop(
        paste0("Optional KNN integer capacity: ", reason),
        reasons = c(capacity = reason), class = "BreastSubtypeR_knn_capacity")
    scalar <- function(value, name) {
        if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
            !is.finite(value) || value < 0 || value != floor(value))
            fail(paste0(name, " must be one finite nonnegative whole number."))
        value <- as.double(value)
        if (value > limit) fail(paste0(name, " exceeds the existing integer limit."))
        value
    }
    n_genes <- scalar(n_genes, "gene count")
    n_samples <- scalar(n_samples, "sample/overlap count")
    n_targets <- scalar(n_targets, "target count")
    n_computed <- scalar(n_computed, "computed target count")
    if (n_targets > n_genes || n_computed > n_targets)
        fail("represented/computed target counts are inconsistent.")
    count_names <- c("NA_only", "NaN", "posInf", "negInf")
    if (!is.numeric(input_counts) || length(input_counts) != 4L ||
        !identical(names(input_counts), count_names))
        fail("input counts must be named NA_only, NaN, posInf, negInf in order.")
    counts <- vapply(seq_along(input_counts), function(i)
        scalar(input_counts[[i]], count_names[[i]]), numeric(1))
    names(counts) <- count_names
    infinity_count <- counts[["posInf"]] + counts[["negInf"]]
    if (infinity_count > limit)
        fail("combined positive/negative infinity count exceeds the integer limit.")
    donor_pairs <- n_computed * max(n_genes - 1, 0)
    if (!is.finite(donor_pairs) || donor_pairs > limit)
        fail("computed targets times other donors exceeds the integer limit.")
    invisible(c(gene_count = n_genes, sample_count = n_samples,
        target_count = n_targets, computed_targets = n_computed,
        donor_pairs = donor_pairs, infinity_count = infinity_count))
}

.mapping_impute_log <- function(x, targets, ssp_targets = character(),
    exclusions = NULL, verbose = FALSE) {
    withr::local_preserve_seed()
    boundary <- function(reason, ids = targets, diagnostics = NULL)
        .mapping_knn_stop(paste0("Optional KNN boundary: ", reason),
            targets = ids, reasons = stats::setNames(rep(reason, length(ids)), ids),
            diagnostics = diagnostics, exclusions = exclusions,
            class = "BreastSubtypeR_knn_boundary")
    if (!is.matrix(x) || typeof(x) != "double")
        boundary("prepared expression must be a double matrix.")
    if (nrow(x) > 0L && (is.null(rownames(x)) || anyNA(rownames(x)) ||
        any(!nzchar(rownames(x))) || anyDuplicated(rownames(x))))
        boundary("prepared gene IDs must be unique, nonmissing and nonempty.")
    if (!is.character(targets) || anyNA(targets) || any(!nzchar(targets)) ||
        anyDuplicated(targets))
        boundary("target IDs must be distinct, nonmissing character IDs.", character())
    targets <- unname(targets)
    target_rows <- match(targets, rownames(x))
    if (anyNA(target_rows))
        boundary("a frozen represented target is absent from the donor matrix.",
            targets[is.na(target_rows)])
    target_holes <- is.na(x[target_rows, , drop = FALSE])
    if (!any(target_holes))
        return(list(matrix = x, diagnostics = NULL, exclusions = exclusions))

    # Guard every count/position category before B3 casts counts to integer.
    counts <- c(NA_only = as.double(sum(is.na(x) & !is.nan(x))),
        `NaN` = as.double(sum(is.nan(x))),
        posInf = as.double(sum(x == Inf, na.rm = TRUE)),
        negInf = as.double(sum(x == -Inf, na.rm = TRUE)))
    missing_per_target <- rowSums(target_holes)
    n_computed <- if (counts[["posInf"]] + counts[["negInf"]] > 0) 0 else
        as.double(sum(missing_per_target > 0 & missing_per_target < ncol(x)))
    .mapping_knn_capacity(as.double(nrow(x)), as.double(ncol(x)),
        as.double(length(targets)), n_computed, counts)

    diagnostic <- knn_vectorized_reference(x, targets)
    # Validate the frozen return schema before indexing fields or testing flags.
    fields <- function(value, expected) is.list(value) &&
        identical(names(value), expected)
    vector_type <- function(value, type, size, missing = FALSE) {
        typeof(value) == type && is.null(dim(value)) && length(value) == size &&
            (missing || !anyNA(value))
    }
    table_type <- function(value, types, size = NULL) {
        is.data.frame(value) && identical(names(value), names(types)) &&
            (is.null(size) || nrow(value) == size) &&
            all(vapply(names(types), function(name)
                vector_type(value[[name]], types[[name]], nrow(value), TRUE),
                logical(1)))
    }
    invalid <- function(reason, ids = targets)
        boundary(reason, ids, diagnostics = diagnostic)
    count_names <- c("NA_only", "NaN", "posInf", "negInf")
    if (!fields(diagnostic, c("matrix", "input_counts", "input_positions", "targets", "work")) ||
        !is.list(diagnostic$targets) || is.data.frame(diagnostic$targets) ||
        !vector_type(diagnostic$input_counts, "integer", 4L) ||
        !identical(names(diagnostic$input_counts), count_names) ||
        any(diagnostic$input_counts < 0L) ||
        !fields(diagnostic$input_positions, count_names) ||
        !vector_type(diagnostic$work, "integer", 2L) ||
        !identical(names(diagnostic$work), c("target_rows", "donor_pairs")) ||
        any(diagnostic$work < 0L))
        invalid("estimator returned an invalid result structure.")
    for (name in count_names) {
        positions <- diagnostic$input_positions[[name]]
        if (!is.matrix(positions) || typeof(positions) != "integer" ||
            ncol(positions) != 2L || anyNA(positions) ||
            nrow(positions) != diagnostic$input_counts[[name]] ||
            any(positions[, 1] < 1L | positions[, 1] > nrow(x)) ||
            any(positions[, 2] < 1L | positions[, 2] > ncol(x)))
            invalid("estimator returned invalid input-position diagnostics.")
    }
    completed <- diagnostic$matrix
    if (!is.matrix(completed) || typeof(completed) != "double" ||
        !identical(dim(completed), dim(x)) ||
        !identical(dimnames(completed), dimnames(x)))
        boundary("estimator dimensions, names/order or storage changed.",
            diagnostics = diagnostic)
    writable <- matrix(FALSE, nrow(x), ncol(x))
    writable[target_rows, ] <- target_holes
    if (!identical(completed[!writable], x[!writable]))
        boundary("estimator changed observed or nontarget values.",
            diagnostics = diagnostic)
    if (!identical(names(diagnostic$targets), targets))
        boundary("estimator changed target identities/order.", diagnostics = diagnostic)

    failure <- character()
    for (i in seq_along(targets)) {
        id <- targets[[i]]
        record <- diagnostic$targets[[i]]
        missing <- unname(as.integer(which(target_holes[i, ])))
        statuses <- c("complete", "imputed", "invalid_input_infinity",
            "entirely_missing", "insufficient_neighbors", "partial",
            "unresolved_coordinates")
        if (!fields(record, c("gene", "index", "status", "missing_columns",
            "neighbors", "candidates", "cells")) ||
            !vector_type(record$gene, "character", 1L) ||
            !vector_type(record$index, "integer", 1L) ||
            !vector_type(record$status, "character", 1L) ||
            !record$status %in% statuses ||
            !vector_type(record$missing_columns, "integer", length(missing)) ||
            !table_type(record$neighbors, c(index = "integer", gene = "character",
                distance = "double", overlap = "integer")) ||
            !table_type(record$candidates, c(index = "integer", gene = "character",
                overlap = "integer", distance = "double", reason = "character")) ||
            !table_type(record$cells, c(column = "integer", contributors = "integer",
                estimate = "double", resolved = "logical", reason = "character"),
                length(missing)))
            invalid("estimator returned malformed target diagnostics.", id)
        computed <- record$status %in% c("imputed", "partial",
            "unresolved_coordinates", "insufficient_neighbors")
        selected <- record$status %in% c("imputed", "partial", "unresolved_coordinates")
        if (nrow(record$candidates) != (if (computed) nrow(x) else 0L) ||
            nrow(record$neighbors) != (if (selected) 10L else 0L))
            invalid("estimator returned invalid donor diagnostic lengths.", id)
        if (computed && (!identical(record$candidates$index, seq_len(nrow(x))) ||
            !identical(record$candidates$gene, rownames(x)) ||
            anyNA(record$candidates$overlap) || anyNA(record$candidates$reason) ||
            any(record$candidates$overlap < 0L | record$candidates$overlap > ncol(x))))
            invalid("estimator returned invalid candidate coordinates.", id)
        if (selected && (anyNA(record$neighbors$index) ||
            any(record$neighbors$index < 1L | record$neighbors$index > nrow(x)) ||
            anyDuplicated(record$neighbors$index) ||
            any(record$neighbors$index == target_rows[[i]]) ||
            !identical(record$neighbors$gene, rownames(x)[record$neighbors$index]) ||
            anyNA(record$neighbors$overlap) ||
            any(record$neighbors$overlap < 1L | record$neighbors$overlap > ncol(x)) ||
            any(!is.finite(record$neighbors$distance))))
            invalid("estimator returned invalid selected-neighbor diagnostics.", id)
        cells <- record$cells
        if (anyNA(cells$column) || anyNA(cells$resolved) || anyNA(cells$reason) ||
            any(cells$contributors < 0L | cells$contributors > 10L, na.rm = TRUE) ||
            any(!cells$reason %in% c("estimated", "invalid_input_infinity",
                "entirely_missing", "insufficient_neighbors", "no_contributors",
                "nonfinite_estimate")))
            invalid("estimator returned malformed cell diagnostics.", id)
        if (!identical(record$gene, id) ||
            !identical(record$index, as.integer(target_rows[[i]])) ||
            !identical(record$missing_columns, missing) ||
            !identical(record$cells$column, missing))
            boundary("estimator target/cell coordinates changed.", id, diagnostic)
        resolved <- cells$resolved
        if (any(resolved)) {
            if (anyNA(cells$contributors[resolved]) ||
                any(cells$contributors[resolved] < 1L) ||
                any(cells$reason[resolved] != "estimated") ||
                any(!is.finite(cells$estimate[resolved])) ||
                !identical(unname(completed[target_rows[[i]], missing[resolved]]),
                    unname(cells$estimate[resolved])))
                invalid("resolved cell estimates disagree with matrix values or diagnostics.", id)
        }
        if (any(!resolved & cells$reason == "estimated") ||
            (record$status == "complete" && length(missing) != 0L) ||
            (record$status == "imputed" && (!length(missing) || !all(resolved))) ||
            (record$status == "partial" && (!any(resolved) || all(resolved))))
            invalid("estimator target status disagrees with cell diagnostics.", id)
        if (!record$status %in% c("complete", "imputed") ||
            !all(record$cells$resolved) ||
            any(!is.finite(record$cells$estimate)) ||
            any(!is.finite(completed[target_rows[[i]], missing]))) {
            reasons <- unique(c(record$status,
                record$cells$reason[!record$cells$resolved]))
            failure[id] <- paste(reasons, collapse = "/")
        }
    }
    if (length(failure)) .mapping_knn_stop(paste0(
        "Optional KNN cannot resolve required represented targets: ",
        paste0(names(failure), " (", failure, ")", collapse = "; "),
        ". No Mapping output was accepted."), targets = names(failure),
        reasons = failure, diagnostics = diagnostic, exclusions = exclusions)

    # The output formulas remain unchanged. Only newly estimated coordinates
    # are checked here, so complete-input behavior is not tightened incidentally.
    projected_ids <- intersect(targets, as.character(ssp_targets))
    projection_failure <- projected_ids[vapply(projected_ids, function(id) {
        row <- match(id, rownames(x))
        holes <- is.na(x[row, ])
        any(!is.finite(2^completed[row, holes]))
    }, logical(1))]
    if (length(projection_failure)) .mapping_knn_stop(paste0(
        "Optional KNN estimates have a nonfinite SSP 2^x projection for targets: ",
        paste(projection_failure, collapse = ", "),
        ". No Mapping output was accepted."), targets = projection_failure,
        reasons = stats::setNames(rep("nonfinite_ssp_projection",
            length(projection_failure)), projection_failure),
        diagnostics = diagnostic, exclusions = exclusions)
    if (isTRUE(verbose)) message("Optional KNN resolved ", sum(target_holes),
        " originally missing represented output measurements.")
    list(matrix = completed, diagnostics = diagnostic, exclusions = exclusions)
}
