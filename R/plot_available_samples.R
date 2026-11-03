# Only wholly unavailable sample measurements are removed, from plot copies.
# Complete inputs keep their existing calculation and ordering paths.
.plot_available_samples <- function(x, out, rows = FALSE, plot) {
    unavailable <- if (rows) rowSums(is.na(x)) == ncol(x) else
        colSums(is.na(x)) == nrow(x)
    if (!any(unavailable)) return(list(x = x, out = out))
    ids <- if (rows) rownames(x) else colnames(x)
    labels <- as.character(out$PatientID)
    valid <- function(z) !is.null(z) && !anyNA(z) && all(nzchar(z)) && !anyDuplicated(z)
    if (!valid(ids) || !valid(labels) || (rows && any(!ids %in% labels)))
        stop(plot, ": filtering unavailable measurements requires unique, matching sample IDs.",
             call. = FALSE)
    message(plot, ": excluded ", sum(unavailable),
            " sample(s) with entirely unavailable measurements: ",
            paste(ids[unavailable], collapse = ", "))
    x <- if (rows) x[!unavailable, , drop = FALSE] else x[, !unavailable, drop = FALSE]
    # Boxes pair rows positionally; expression plots retain their established
    # annotation order / named join. Never filter on the subtype call.
    out <- if (rows) out[match(ids[!unavailable], labels), , drop = FALSE] else
        out[!labels %in% ids[unavailable], , drop = FALSE]
    out$PatientID <- as.character(out$PatientID)
    list(x = x, out = out)
}
