write_tsv <- function(x, path) utils::write.table(x, path, sep = "\t", row.names = FALSE, na = "NA")
save_test_results <- function(result, prefix) {
    d <- as.data.frame(result)
    write_tsv(d[c("file", "test", "nb", "passed", "failed", "error", "warning", "skipped")], paste0(prefix, "-blocks.tsv"))
    conditions <- list()
    for (b in result) for (j in seq_along(b$results)) {
        cnd <- b$results[[j]]
        if (inherits(cnd, "expectation_success")) next
        conditions[[length(conditions) + 1L]] <- data.frame(file = b$file,
            test = b$test, index = j, class = paste(class(cnd), collapse = "/"),
            message = conditionMessage(cnd), stringsAsFactors = FALSE)
    }
    cond <- if (length(conditions)) do.call(rbind, conditions) else
        data.frame(file = character(), test = character(), index = integer(), class = character(), message = character())
    write_tsv(cond, paste0(prefix, "-conditions.tsv"))
    summary <- aggregate(d[c("passed", "failed", "error", "warning", "skipped")], list(file = d$file), sum)
    summary$blocks <- as.integer(table(d$file)[summary$file])
    write_tsv(summary, paste0(prefix, "-summary.tsv"))
    print(summary)
    d
}
numeric_fingerprint <- function(x) digest::digest(
    writeBin(as.double(x), raw(), size = 8L, endian = "little"), algo = "sha256", serialize = FALSE)

# A read-only, full-leaf diagnostic. This does not replace any package assertion.
leaf_inventory <- function(actual, reference, dataset, policy) {
    leaves <- list(); structural <- list()
    structure_error <- function(path, reason) {
        structural[[length(structural) + 1L]] <<- data.frame(dataset, path, reason)
    }
    numeric_leaf <- function(a, b, path) {
        if (!is.double(a) || !is.double(b)) {
            if (!identical(a, b, num.eq = FALSE)) structure_error(path, "non-double value changed")
            return(invisible(NULL))
        }
        finite <- is.finite(a) & is.finite(b)
        masks <- identical(is.na(a), is.na(b)) && identical(is.nan(a), is.nan(b)) &&
            identical(is.infinite(a), is.infinite(b)) && identical(a[is.infinite(a)], b[is.infinite(b)])
        x <- a[finite]; y <- b[finite]
        zeros <- identical(x == 0, y == 0) && identical(sign(x), sign(y)) &&
            identical(1 / x[x == 0], 1 / y[y == 0])
        tiny_x <- x != 0 & abs(x) < .Machine$double.xmin
        tiny_y <- y != 0 & abs(y) < .Machine$double.xmin
        tiny <- identical(tiny_x, tiny_y) && identical(x[tiny_x], y[tiny_y], num.eq = FALSE)
        rule <- policy[policy$scope == "C_vs_approved_C" & policy$dataset == dataset & policy$path == path, , drop = FALSE]
        derived <- path %in% c("/res/results/AIMS/prob", "/res/results/sspbc/prob")
        status <- if (nrow(rule)) "bounded" else if (derived) "derived" else "exact"
        delta <- abs(x - y)
        exceeded <- if (nrow(rule)) sum(delta > rule$atol + rule$rtol * pmax(abs(x), abs(y))) else
            if (derived) 0L else sum(x != y)
        exact <- identical(a, b, num.eq = FALSE)
        accepted <- masks && zeros && tiny && if (derived) TRUE else if (nrow(rule)) exceeded == 0L else exact
        leaves[[length(leaves) + 1L]] <<- data.frame(dataset, path, n = length(a),
            changed = sum(x != y), max_abs = if (length(delta)) max(delta) else 0,
            policy = status, exceeded, accepted, exact, masks, zeros, tiny,
            actual_sha256 = numeric_fingerprint(a), reference_sha256 = numeric_fingerprint(b), stringsAsFactors = FALSE)
    }
    walk <- function(a, b, path) {
        if (!identical(typeof(a), typeof(b)) || !identical(class(a), class(b)) || length(a) != length(b)) {
            structure_error(path, "type/class/length"); return(invisible(NULL))
        }
        if (isS4(a)) {
            slots <- methods::slotNames(a)
            if (!identical(slots, methods::slotNames(b))) {structure_error(path, "S4 slots");return(invisible(NULL))}
            ax <- attributes(a); bx <- attributes(b); ax[slots] <- NULL; bx[slots] <- NULL
            if (!identical(ax, bx, num.eq = FALSE)) structure_error(path, "S4 attributes")
            for (n in slots) walk(methods::slot(a, n), methods::slot(b, n), paste0(path, "/", n))
            return(invisible(NULL))
        }
        if (!identical(attributes(a), attributes(b), num.eq = FALSE)) structure_error(path, "attributes/names/dimensions")
        paths <- c("/data_input/se_NC/assays/data/listData/counts", "/data_input/se_SSP/assays/data/listData/counts")
        i <- match(path, paths)
        if (!is.na(i)) path <- c("/data_input/se_NC/assay", "/data_input/se_SSP/assay")[i]
        rule <- policy[policy$scope == "C_vs_approved_C" & policy$dataset == dataset & policy$path == path, , drop = FALSE]
        if (is.data.frame(a) && nrow(rule) && all(vapply(a, is.double, logical(1))) && all(vapply(b, is.double, logical(1)))) {
            numeric_leaf(as.matrix(a), as.matrix(b), path)
        } else if (is.list(a)) {
            for (j in seq_along(a)) walk(a[[j]], b[[j]], paste0(path, "/", if (length(names(a))) names(a)[j] else j))
        } else if (is.double(a)) numeric_leaf(a, b, path)
        else if (!identical(a, b, num.eq = FALSE)) structure_error(path, "non-double value changed")
        invisible(NULL)
    }
    walk(actual, reference, "")
    list(leaves = do.call(rbind, leaves), structural = if (length(structural)) do.call(rbind, structural) else
        data.frame(dataset = character(), path = character(), reason = character()))
}
