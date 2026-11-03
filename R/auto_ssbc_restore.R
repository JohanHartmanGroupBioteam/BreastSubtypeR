# Output assembly only: never pass padded values back to a classifier.
.auto_ssbc_ids <- function(ids, n, field, universe = ids) {
    if (n == 0L && is.null(ids)) return(invisible(ids))
    if (!is.character(ids) || length(ids) != n || anyNA(ids) ||
        any(!nzchar(ids)) || anyDuplicated(ids) || any(!ids %in% universe)) {
        stop("Cannot restore AUTO ssBC results: invalid sample IDs in ", field,
             ". IDs must be complete, unique and drawn from the input assay.",
             call. = FALSE)
    }
    invisible(ids)
}

.restore_auto_ssbc <- function(result, ids) {
    .auto_ssbc_ids(ids, length(ids), "input assay")
    if (anyDuplicated(names(result)) || anyDuplicated(names(result$outList)))
        stop("Cannot restore AUTO ssBC results: duplicated field names.", call. = FALSE)

    restore_frame <- function(x, key, field, use_key = FALSE) {
        if (!is.data.frame(x) || !key %in% names(x) || anyDuplicated(names(x)) ||
            !(is.character(x[[key]]) || is.factor(x[[key]])))
            stop("Cannot restore AUTO ssBC results: invalid ", field, " schema.", call. = FALSE)
        row_ids <- if (use_key) as.character(x[[key]]) else rownames(x)
        .auto_ssbc_ids(row_ids, nrow(x), field, ids)
        labels <- as.character(x[[key]])
        if (!identical(rownames(x), row_ids) ||
            any(!is.na(labels) & labels != row_ids))
            stop("Cannot restore AUTO ssBC results: conflicting IDs in ", field, ".", call. = FALSE)
        index <- match(ids, row_ids)
        restored <- x[index, , drop = FALSE]
        rownames(restored) <- ids
        # Existing NA identifier payloads (native missing-receptor rows) stay NA.
        # Only newly added rows acquire their actual identifier label.
        added <- is.na(index)
        if (any(added)) {
            if (is.factor(restored[[key]]))
                levels(restored[[key]]) <- union(levels(restored[[key]]), ids[added])
            restored[[key]][added] <- ids[added]
        }
        restored
    }
    for (field in c("BS.all", "score.ROR")) {
        if (field %in% names(result))
            result[[field]] <- restore_frame(result[[field]],
                if (field == "BS.all") "PatientID" else "patientID", field,
                use_key = field == "BS.all")
    }
    for (field in c("predictions", "predictions.FourSubtype")) {
        if (field %in% names(result$outList)) {
            x <- result$outList[[field]]
            if (!is.atomic(x) || !is.null(dim(x)))
                stop("Cannot restore AUTO ssBC results: invalid ", field, " vector.", call. = FALSE)
            .auto_ssbc_ids(names(x), length(x), field, ids)
            result$outList[[field]] <- stats::setNames(x[match(ids, names(x))], ids)
        }
    }
    for (field in c("distances", "dist.RORSubtype", "testData")) {
        if (field %in% names(result$outList)) {
            x <- result$outList[[field]]
            if (!is.matrix(x))
                stop("Cannot restore AUTO ssBC results: ", field,
                     " must retain its native matrix dimensions (see issue #234).", call. = FALSE)
            columns <- field == "testData"
            axis_ids <- if (columns) colnames(x) else rownames(x)
            .auto_ssbc_ids(axis_ids, if (columns) ncol(x) else nrow(x), field, ids)
            index <- match(ids, axis_ids)
            if (columns) {
                x <- x[, index, drop = FALSE]
                colnames(x) <- ids
            } else {
                x <- x[index, , drop = FALSE]
                rownames(x) <- ids
            }
            result$outList[[field]] <- x
        }
    }
    result
}
