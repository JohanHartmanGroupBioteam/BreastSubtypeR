# B3 experimental optimization of frozen B2. No Mapping integration.
# Column-wise donor vectorization preserves each pair's B2 summation order.
knn_vectorized_reference <- function(x, targets) {
    if (!is.matrix(x) || typeof(x) != "double")
        stop("x must be a double-precision numeric matrix")
    genes <- rownames(x)
    if (nrow(x) > 0L && (is.null(genes) || anyNA(genes) ||
        any(!nzchar(genes)) || anyDuplicated(genes)))
        stop("x must have unique, nonmissing, nonempty gene row names")
    if (!is.character(targets) || anyNA(targets) ||
        any(!nzchar(targets)) || anyDuplicated(targets))
        stop("targets must be distinct, nonmissing, nonempty character IDs")

    # Classify original values once; never read from completed when finding donors.
    masks <- list(NA_only = is.na(x) & !is.nan(x), `NaN` = is.nan(x),
                  posInf = x == Inf, negInf = x == -Inf)
    input_counts <- vapply(masks, function(z) as.integer(sum(z, na.rm = TRUE)), integer(1))
    input_positions <- lapply(masks, which, arr.ind = TRUE)
    finite <- is.finite(x)
    invalid <- input_counts[["posInf"]] + input_counts[["negInf"]] > 0L
    completed <- x
    records <- setNames(vector("list", length(targets)), unname(targets))
    work <- c(target_rows = 0L, donor_pairs = 0L)
    empty_neighbors <- function() data.frame(index = integer(), gene = character(),
        distance = double(), overlap = integer(), stringsAsFactors = FALSE)
    empty_candidates <- function() data.frame(index = integer(), gene = character(),
        overlap = integer(), distance = double(), reason = character(),
        stringsAsFactors = FALSE)
    # Explicit left-to-right ordinary double summation, also used for means.
    ordered_sum <- function(values) {
        total <- 0
        for (value in values) total <- total + value
        total
    }
    for (q in seq_along(targets)) {
        gene <- unname(targets[[q]])
        ti <- match(gene, genes)
        missing <- if (is.na(ti)) integer() else which(is.na(x[ti, ]))
        missing <- unname(as.integer(missing))
        record <- list(gene = gene, index = as.integer(ti), status = "absent",
            missing_columns = missing, neighbors = empty_neighbors(),
            candidates = empty_candidates(),
            cells = data.frame(column = missing,
                contributors = rep(NA_integer_, length(missing)),
                estimate = rep(NA_real_, length(missing)),
                resolved = rep(FALSE, length(missing)),
                reason = rep("", length(missing)), stringsAsFactors = FALSE))
        if (is.na(ti)) {
            records[[q]] <- record
            next
        }
        if (invalid) {
            record$status <- "invalid_input_infinity"
            record$cells$reason[] <- "invalid_input_infinity"
            records[[q]] <- record
            next
        }
        if (!length(missing)) {
            record$status <- "complete"
            records[[q]] <- record
            next
        }
        if (length(missing) == ncol(x)) {
            record$status <- "entirely_missing"
            record$cells$reason[] <- "entirely_missing"
            records[[q]] <- record
            next
        }

        # Same ordered arithmetic per pair as B2, vectorized ACROSS donors.
        # Iterate columns in input order. Add only shared observed contributions.
        # Every other row remains in the donor universe; no shortlist or pruning.
        distance_sums <- double(nrow(x))
        overlaps <- integer(nrow(x))
        other <- seq_len(nrow(x)) != ti
        for (column in which(finite[ti, ])) {
            shared <- finite[, column] & other
            differences <- x[ti, column] - x[shared, column]
            squared <- differences * differences
            distance_sums[shared] <- distance_sums[shared] + squared
            overlaps[shared] <- overlaps[shared] + 1L
        }
        distances <- rep(NA_real_, nrow(x))
        supported_overlap <- overlaps > 0L
        distances[supported_overlap] <- distance_sums[supported_overlap] /
                                         overlaps[supported_overlap]
        reasons <- rep("no_overlap", nrow(x))
        reasons[supported_overlap] <- ifelse(is.finite(distances[supported_overlap]),
                                            "supported", "nonfinite_distance")
        reasons[ti] <- "self"
        # Construct diagnostics once, instead of copying fields for each donor.
        candidates <- data.frame(index = seq_len(nrow(x)), gene = genes,
            overlap = overlaps, distance = distances, reason = reasons,
            stringsAsFactors = FALSE)
        work[["target_rows"]] <- work[["target_rows"]] + 1L
        work[["donor_pairs"]] <- work[["donor_pairs"]] + as.integer(nrow(x) - 1L)
        record$candidates <- candidates
        supported <- which(candidates$reason == "supported")
        if (length(supported) < 10L) {
            record$status <- "insufficient_neighbors"
            record$cells$reason[] <- "insufficient_neighbors"
            records[[q]] <- record
            next
        }
        ranked <- supported[order(candidates$distance[supported],
                                  candidates$index[supported], method = "radix")]
        chosen <- ranked[seq_len(10L)]
        record$neighbors <- candidates[chosen, c("index", "gene", "distance", "overlap")]
        rownames(record$neighbors) <- NULL
        for (ci in seq_along(missing)) {
            column <- missing[[ci]]
            available <- x[chosen, column]
            available <- available[is.finite(available)]
            contributors <- length(available)
            record$cells$contributors[ci] <- as.integer(contributors)
            if (contributors == 0L) {
                record$cells$reason[ci] <- "no_contributors"
                next
            }
            estimate <- ordered_sum(available) / contributors
            record$cells$estimate[ci] <- estimate
            if (!is.finite(estimate)) {
                record$cells$reason[ci] <- "nonfinite_estimate"
                next
            }
            record$cells$reason[ci] <- "estimated"
            record$cells$resolved[ci] <- TRUE
            completed[ti, column] <- estimate
        }
        record$status <- if (all(record$cells$resolved)) "imputed" else
            if (any(record$cells$resolved)) "partial" else "unresolved_coordinates"
        records[[q]] <- record
    }
    list(matrix = completed, input_counts = input_counts,
         input_positions = input_positions, targets = records, work = work)
}
