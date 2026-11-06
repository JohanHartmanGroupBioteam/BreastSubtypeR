.shiny_input_lines <- function(path) {
    lines <- tryCatch(readLines(path), warning = function(w) {
        # A missing final newline is valid. readLines(warn = FALSE) would
        # also silently truncate a record at an embedded NUL.
        con <- file(path, "rb")
        on.exit(close(con))
        repeat {
            bytes <- readBin(con, "raw", 65536L)
            if (any(bytes == as.raw(0L)))
                stop("Input file contains a NUL byte. Export the table as plain text.", call. = FALSE)
            if (length(bytes) < 65536L) break
        }
        readLines(path, warn = FALSE)
    })
    # Mark UTF-8 records explicitly, including on hosts using the C locale.
    # Other records keep the input connection's native encoding.
    utf8 <- validUTF8(lines)
    Encoding(lines[utf8]) <- "UTF-8"
    if (length(lines) && utf8[1L]) lines[1L] <- sub("^\ufeff", "", lines[1L])
    lines
}

# Read records before associating fields with headers. Appending the separator
# preserves an explicit empty last field in strsplit().
.shiny_delimited_records <- function(lines, sep, kind = "Expression") {
    line <- seq_along(lines)
    quoted <- grepl('"', lines, fixed = TRUE)
    split_at <- sep
    if (any(quoted)) {
        # A double-quoted field can contain delimiters, doubled quotes and
        # newlines. Only quotes enclosing a complete field have this meaning.
        quoted_field <- '"(?:[^"]++|"")*+"'
        field <- paste0("(?:", quoted_field, '|(?!")[^', sep, "\r\n]*+)")
        record <- paste0("\\A", field, "(?:", sep, field, ")*+\\z")
        invalid <- which(quoted & !grepl(record, lines, perl = TRUE))
        if (length(invalid)) {
            con <- textConnection(lines, encoding = "UTF-8")
            on.exit(close(con))
            widths <- utils::count.fields(con, sep = sep, quote = '"',
                comment.char = "", blank.lines.skip = FALSE)
            ends <- which(!is.na(widths))
            if (length(widths) == length(lines) && length(ends) &&
                ends[length(ends)] == length(lines)) {
                line <- c(1L, utils::head(ends, -1L) + 1L)
                for (i in which(ends > line)) {
                    lines[line[i]] <- paste(lines[line[i]:ends[i]], collapse = "\n")
                }
                lines <- lines[line]
                invalid <- line[!grepl(record, lines, perl = TRUE)]
            }
            if (length(invalid)) {
                stop(sprintf("%s file line %d has an invalid quoted field. Enclose the whole field in double quotes and double embedded quotes.",
                             kind, invalid[1L]), call. = FALSE)
            }
        }
        split_at <- paste0("(?<![^", sep, "])", quoted_field, "(*SKIP)(*F)|", sep)
    }
    keep <- nzchar(trimws(lines)) | grepl(sep, lines, fixed = TRUE)
    lines <- lines[keep]
    fields <- strsplit(paste0(lines, sep), split_at, fixed = !any(quoted),
                       perl = any(quoted))
    for (i in which(grepl('"', lines, fixed = TRUE))) {
        x <- fields[[i]]
        q <- startsWith(x, '"')
        x[q] <- gsub('""', '"', substring(x[q], 2L, nchar(x[q]) - 1L), fixed = TRUE)
        fields[[i]] <- x
    }
    list(fields = fields, line = line[keep])
}

.read_shiny_metadata <- function(path, sep, key) {
    lines <- .shiny_input_lines(path)
    if (any(grepl('"', lines, fixed = TRUE))) {
        # Give read.table unambiguous quoting, including when a legacy bare
        # ID contains a literal quote alongside conventionally quoted fields.
        records <- .shiny_delimited_records(lines, sep, "Metadata")
        lines <- vapply(records$fields, function(x) {
            paste0('"', gsub('"', '""', x, fixed = TRUE), '"', collapse = sep)
        }, character(1L))
    }
    df <- utils::read.table(text = lines, header = TRUE, sep = sep, row.names = NULL,
                            check.names = FALSE, quote = '"', comment.char = "", fill = TRUE,
                            colClasses = "character", na.strings = character())
    # Only the matching key bypasses read.table's usual type conversion.
    other <- names(df) != key
    df[other] <- lapply(df[other], utils::type.convert, as.is = TRUE)
    df
}

.read_shiny_expression <- function(path, sep, sample_ids, probe_ids, raw_counts) {
    records <- .shiny_delimited_records(.shiny_input_lines(path), sep)
    if (length(records$line) < 2L)
        stop("Expression file needs a header and at least one data row.", call. = FALSE)
    header <- records$fields[[1L]]
    rows <- records$fields[-1L]
    line <- records$line[-1L]
    widths <- lengths(rows)
    width <- length(header)
    probe <- which(tolower(header) == "probe")
    if (length(probe) > 1L)
        stop("Use exactly one probe header for the gene IDs.", call. = FALSE)

    numbered <- as.character(seq_along(rows))
    first <- vapply(rows, '[[', character(1L), 1L)
    annotated <- first %in% probe_ids
    leading_ids <- all(nzchar(first)) && !anyDuplicated(first) && any(annotated)
    numbered_ids <- any(numbered %in% probe_ids)
    extra <- which(widths > width)
    width_error <- function(i) {
        stop(sprintf("Expression file line %d has %d fields, expected %d. Check the delimiter and use a decimal point. Label gene IDs with an explicit probe header.",
                     line[i], widths[i], width), call. = FALSE)
    }
    # Used only to veto an inferred ID layout. The original tokens are kept
    # for the measurement parser, which rejects decimal commas.
    measurement_like <- function(x) {
        token <- trimws(x)
        number <- suppressWarnings(as.double(chartr(",", ".", token)))
        !is.na(number) | is.nan(number) |
            tolower(token) %in% c("", "na", "n/a", "#n/a", "-")
    }
    short <- which(widths < width)
    if (length(short)) {
        i <- short[1L]
        stop(sprintf("Expression file line %d has %d fields, expected %d. Preserve each sample's position with its delimiter. Use NA or an empty field for missing normalized expression.",
                     line[i], widths[i], width), call. = FALSE)
    }
    if (length(extra)) {
        # R row-name exports may omit the corner header. An extra field alone
        # is insufficient evidence for shifting every sample one place.
        too_wide <- which(widths > width + 1L)
        if (length(too_wide)) width_error(too_wide[1L])
        if (any(widths != width + 1L) || length(probe) || !nzchar(header[1L]) ||
            !any(header %in% sample_ids) ||
            !leading_ids || numbered_ids ||
            (!all(annotated) && all(measurement_like(first[annotated])))) {
            width_error(extra[1L])
        }
        header <- c("", header)
    }
    cells <- do.call(rbind, rows)
    last_named <- max(c(0L, which(nzchar(header))))
    # A data-only trailing delimiter can mimic an omitted corner header.
    # In this ambiguous case require an explicit probe header or final NA.
    if (length(extra) && all(cells[, last_named:ncol(cells), drop = FALSE] == ""))
        width_error(extra[1L])
    if (!length(probe)) {
        if (!nzchar(header[1L])) {
            probe <- 1L
        } else if (!header[1L] %in% sample_ids) {
            leading_ids <- leading_ids && any(header[-1L] %in% sample_ids)
            # An unmatched sample's measurements can coincide with annotation
            # IDs at every row. Only infer a named column from textual IDs.
            if (leading_ids == numbered_ids || (leading_ids && any(measurement_like(first))))
                stop("Ambiguous first expression column. Label gene IDs with an explicit probe header.", call. = FALSE)
            if (leading_ids) probe <- 1L
        }
    }
    samples <- setdiff(seq_along(header), probe)
    trailing <- which(seq_along(header) > last_named)
    if (length(trailing) && all(cells[, trailing, drop = FALSE] == ""))
        samples <- setdiff(samples, trailing)
    if (!length(samples) || any(!nzchar(header[samples])) || anyDuplicated(header[samples]))
        stop("Expression file needs a unique, nonempty header for each sample.", call. = FALSE)

    ids <- if (length(probe)) cells[, probe] else numbered
    if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
        stop("Expression probe IDs must be nonempty and unique.", call. = FALSE)
    result <- matrix(NA_real_, length(ids), length(samples),
                     dimnames = list(ids, header[samples]))
    errors <- character()
    count <- 0L
    unmatched_error <- FALSE
    for (j in seq_along(samples)) {
        token <- trimws(cells[, samples[j]])
        missing <- is.na(token) | token %in% c("", "NA")
        value <- suppressWarnings(as.double(token))
        # A numeric literal contains a digit. R's explicit infinity spellings
        # do not, regardless of case, sign or accepted surrounding whitespace.
        overflow <- is.infinite(value) & grepl("[0-9]", token)
        bad <- which(overflow | (is.na(value) & !is.nan(value) & !missing) |
                     (raw_counts & is.na(value)))
        count <- count + length(bad)
        if (length(bad) && !header[samples[j]] %in% sample_ids) unmatched_error <- TRUE
        for (i in utils::head(bad, 3L - length(errors))) {
            shown <- if (is.na(token[i]) || !nzchar(token[i])) "<blank>" else
                encodeString(substr(token[i], 1L, 64L), quote = '"')
            reason <- if (overflow[i]) " (numeric overflow)" else ""
            errors <- c(errors, sprintf("line %d, probe %s, sample %s = %s%s",
                line[i], ids[i], header[samples[j]], shown, reason))
        }
        result[, j] <- value
    }
    if (count) {
        hint <- if (raw_counts) "Raw counts require a number at every measurement." else
            "Use NA, NaN or an empty field for missing normalized expression."
        if (unmatched_error)
            hint <- paste(hint, "Check unmatched column names against clinical PatientID. Put feature metadata in the annotation table.")
        noun <- if (count == 1L) "measurement" else "measurements"
        stop(sprintf("Expression file has %d invalid %s (%s). Use a decimal point. %s",
                     count, noun, paste(errors, collapse = " / "), hint), call. = FALSE)
    }
    result
}
