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

# With no recognized pair, native ER.v2 has already restored every sample.
# Its score table has no patientID column and keeps logical NA in every
# column. Validate this complete return before passing it through verbatim.
# The four-class return (NCP-02) adds only BS.Subtype, Call.Subtype after Call
# and predictions.FourSubtype; the five-class schema is unchanged.
.validate_empty_auto_ssbc <- function(result, ids, clinical, subtype = FALSE) {
    classes <- c("Basal", "Her2", "LumA", "LumB", "Normal")
    fields <- c(classes, "Call", if (subtype) "Call.Subtype", "Confidence",
        "ROR-S (Subtype Only)", "ROR-S Group (Subtype Only)",
        "ROR-P (Subtype + Prolif)", "ROR-P Group (Subtype + Prolif)",
        if (clinical) c("ROR-C (Subtype + Clinic)", "ROR-C Group (Subtype + Clinic)",
            "ROR-PC (Subtype + Clinic + Prolif)", "ROR-PC Group (Subtype + Clinic + Prolif)",
            "ROR-PC (Subtype + Clinic + Prolif.Subtype)",
            "ROR-PC Group (Subtype + Clinic + Prolif.Subtype)"), "ER", "HER2")
    n <- length(ids)
    missing <- rep(NA, n)
    calls <- if (subtype) {
        data.frame(PatientID = ids, BS = missing, BS.Subtype = missing, row.names = ids)
    } else {
        data.frame(PatientID = ids, BS = missing, row.names = ids)
    }
    scores <- as.data.frame(stats::setNames(rep(list(missing), length(fields)), fields),
        row.names = ids, check.names = FALSE)
    if (!identical(attributes(result), list(names = c("BS.all", "score.ROR", "mdns", "outList"))) ||
        !identical(result$BS.all, calls) || !identical(result$score.ROR, scores))
        stop("Invalid zero-pair AUTO ssBC.v2 result frames.", call. = FALSE)

    correlations <- function(labels) {
        matrix(NA_real_, n, length(labels), dimnames = list(ids, labels))
    }
    models <- new.env(parent = emptyenv())
    utils::data("BreastSubtypeRobj", package = "BreastSubtypeR", envir = models)
    model <- models$BreastSubtypeRobj
    no_call <- stats::setNames(rep(NA_character_, n), ids)
    expected <- c(list(predictions = no_call),
        if (subtype) list(predictions.FourSubtype = no_call),
        list(testData = matrix(numeric(), 0L, n, dimnames = list(NULL, ids)),
        distances = correlations(classes), dist.RORSubtype = correlations(classes[1:4]),
        centroids = model$centroid))
    if (!identical(result$outList, expected) || !identical(result$mdns, model$ssBC.subgroupQuantile))
        stop("Invalid zero-pair AUTO ssBC.v2 prediction or model fields.", call. = FALSE)
    invisible(NULL)
}

.restore_auto_ssbc <- function(result, ids, zero_pair_clinical = NULL,
                               zero_pair_subtype = FALSE) {
    .auto_ssbc_ids(ids, length(ids), "input assay")
    if (!is.null(zero_pair_clinical)) {
        .validate_empty_auto_ssbc(result, ids, zero_pair_clinical, zero_pair_subtype)
        return(result)
    }
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
