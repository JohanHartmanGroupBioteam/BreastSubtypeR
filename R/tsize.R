# Names and display labels do not change a stored clinical code. Other
# attributes must survive until the decoder has checked their semantics.
.tsize_display_attrs <- c("names", "comment", "label", "labels", "value.labels", "units",
    "format.spss", "format.sas", "format.stata", "display_width")

.needs_tsize_decode <- function(x) {
    (is.object(x) && !identical(class(x), "factor") &&
        !identical(class(x), c("ordered", "factor"))) ||
        length(setdiff(names(attributes(x)), c("class", if (is.factor(x)) "levels",
                                                .tsize_display_attrs))) > 0L
}

.tsize_metadata_error <- function(detail) {
    stop("TSIZE has unverified metadata (", detail,
         "). Resolve special missing codes to NA and supply ordinary binary codes.",
         call. = FALSE)
}

.check_tsize_attrs <- function(keys, allowed) {
    extra <- setdiff(keys, allowed)
    if (length(extra)) .tsize_metadata_error(paste(extra, collapse = ", "))
}

# Read only transparent declaration storage. Check the declaration's own
# attributes first: a labelled range can carry a second missing-value rule.
.tsize_declared_codes <- function(x, key) {
    if (is.null(x)) return(NULL)
    extra <- setdiff(names(attributes(x)), c("class", .tsize_display_attrs))
    if (length(extra))
        .tsize_metadata_error(paste(key, "contains", paste(extra, collapse = ", ")))
    storage <- typeof(x)
    wrappers <- c("AsIs", "labelled", "haven_labelled", "haven_labelled_spss",
                  "vctrs_vctr", storage, if (storage == "double") "numeric")
    classes <- oldClass(x)
    haven <- any(classes %in% c("haven_labelled", "haven_labelled_spss"))
    if (haven != ("vctrs_vctr" %in% classes))
        .tsize_metadata_error(paste("haven cannot validate", key, "class stack"))
    if (isS4(x) || !storage %in% c("integer", "double", "character") ||
        anyDuplicated(classes) || !all(classes %in% wrappers))
        .tsize_metadata_error(paste("haven cannot validate", key, "class or storage"))
    as.vector(unclass(x))
}

# Inspect the original column before row selection or class removal. Only
# haven's SPSS family can attach a missingness rule to stored binary codes.
.preflight_tsize <- function(x) {
    classes <- oldClass(x)
    spss <- "haven_labelled_spss" %in% classes
    .check_tsize_attrs(names(attributes(x)), c("class", if (is.factor(x)) "levels",
        .tsize_display_attrs, if (spss) c("na_values", "na_range")))
    if (spss || "haven_labelled" %in% classes) {
        storage <- typeof(x)
        core <- setdiff(classes, c("AsIs", "labelled"))
        if (is.factor(x)) {
            if (spss || storage != "integer" ||
                !all(core %in% c("ordered", "factor", "haven_labelled")))
                .tsize_metadata_error("haven cannot validate factor storage in this class stack")
            oldClass(x) <- c(if ("ordered" %in% core) "ordered", "factor")
            return(x)
        }
        allowed <- c("haven_labelled", "vctrs_vctr", storage,
                     if (storage == "double") "numeric", if (spss) "haven_labelled_spss")
        if (anyDuplicated(classes) || !storage %in% c("integer", "double", "character") ||
            (spss && !"vctrs_vctr" %in% core) || !all(core %in% allowed))
            .tsize_metadata_error("unsupported haven class stack or storage")
        if (spss) {
            values <- .tsize_declared_codes(attr(x, "na_values", exact = TRUE), "na_values")
            range <- .tsize_declared_codes(attr(x, "na_range", exact = TRUE), "na_range")
            if (!requireNamespace("haven", quietly = TRUE))
                stop("An SPSS-labelled TSIZE column requires the optional haven package. ",
                     "Install haven or convert SPSS user-missing values to NA before supplying ",
                     "ordinary binary TSIZE codes.", call. = FALSE)
            # Let haven validate and cast the declarations together. Its native
            # is.na method then evaluates this local copy in both cold and warm sessions.
            codes <- as.vector(unclass(x))
            names(codes) <- names(x)
            return(tryCatch(haven::labelled_spss(codes, na_values = values, na_range = range),
                error = function(e) .tsize_metadata_error(paste(
                    "haven cannot validate na_values or na_range:", conditionMessage(e)))))
        }
        oldClass(x) <- c("haven_labelled", "vctrs_vctr", storage)
    } else if ("AsIs" %in% classes) {
        oldClass(x) <- setdiff(classes, "AsIs")
    }
    x
}


# A plain column survives data.frame coercion and assay subsetting. Keep
# invalid distinct from missing until the evaluated PatientIDs are known.
.scoring_tsize <- function(x) {
    size <- .decode_tsize(x)
    codes <- as.character(size$value)
    codes[size$invalid] <- "invalid"
    codes
}

.nc_pheno <- function(se_obj, hasClinical) {
    cd <- colData(se_obj)
    if (.needs_tsize_decode(cd$TSIZE)) {
        cd$TSIZE <- if (hasClinical) .scoring_tsize(cd$TSIZE) else NULL
    }
    data.frame(cd)
}

# Decode a clinical covariate without changing the supplied annotations.
.decode_tsize <- function(x) {
    if (isS4(x)) {
        if (!methods::is(x, "Rle"))
            stop("Unsupported S4 TSIZE class. Supply binary codes in a supported vector or Rle.",
                 call. = FALSE)
        if (!identical(as.character(class(x)), "Rle"))
            .tsize_metadata_error("unsupported Rle subclass")
        .check_tsize_attrs(names(attributes(x)), c("class", "values", "lengths",
            "metadata", "elementMetadata", .tsize_display_attrs))
        .check_tsize_attrs(c(names(S4Vectors::metadata(x)), names(S4Vectors::mcols(x))),
                           .tsize_display_attrs)
        runs <- S4Vectors::runValue(x)
        .check_tsize_attrs(names(attributes(runs)), c("class", if (is.factor(runs)) "levels",
                                                     .tsize_display_attrs))
        if (!(identical(class(runs), "factor") || identical(class(runs), c("ordered", "factor")) ||
              (!is.object(runs) && typeof(runs) %in% c("double", "integer", "logical", "character"))))
            stop("Unsupported TSIZE Rle run values. Supply numeric, logical, factor or character codes.",
                 call. = FALSE)
        x <- as.vector(x)
    } else {
        x <- .preflight_tsize(x)
    }
    if (inherits(x, "labelled") && !inherits(x, "haven_labelled")) {
        classes <- oldClass(x)
        if (typeof(x) %in% c("integer", "double") && classes[1L] == "labelled" &&
            all(classes[-1L] %in% c("numeric", "integer")))
            x <- unclass(x)
    }
    if (inherits(x, "integer64")) {
        # The double storage contains integer bits. Decode through the class
        # method, including for a serialized column read in a fresh session.
        if (!requireNamespace("bit64", quietly = TRUE))
            stop("An integer64 TSIZE column requires the optional bit64 package. ",
                 "Install bit64 or supply ordinary binary TSIZE codes.", call. = FALSE)
        x <- as.character(x)
    }
    if (is.factor(x)) x <- as.character(x)
    missing <- is.na(x)
    # Preserve class-specific missingness, then parse stored character labels
    # or haven's known numeric storage without a class-specific double cast.
    if (is.character(x) ||
        (inherits(x, "haven_labelled") && typeof(x) %in% c("integer", "double")))
        x <- unclass(x)
    if (is.character(x)) {
        text <- trimws(x, whitespace = "[ \t\r\n\v\f]")
        value <- suppressWarnings(as.numeric(text))
        # A parsed double can round to 1 or underflow to 0. An exact decimal
        # one has just one nonzero mantissa digit (1); for hex it is 1/2/4/8.
        hex <- grepl("^[+-]?0[xX]", text)
        mantissa <- sub("[eE].*$", "", text)
        mantissa[hex] <- sub("[pP].*$", "", sub("^[+-]?0[xX]", "", text[hex]))
        digits <- gsub("[.+-]", "", mantissa)
        one <- grepl("^0*10*$", digits)
        one[hex] <- grepl("^0*[1248]0*$", digits[hex])
        valid <- (value %in% 0 & grepl("^0+$", digits)) | (value %in% 1 & one)
    } else if ((!is.object(x) && is.numeric(x)) || is.logical(x) || is.raw(x)) {
        value <- tryCatch(as.numeric(x), error = function(e) rep(NA_real_, length(x)))
        valid <- value %in% c(0, 1)
    } else {
        value <- rep(NA_real_, length(x))
        valid <- rep(FALSE, length(x))
    }
    value[!valid | missing] <- NA_real_
    list(value = value, invalid = !missing & !valid)
}

.warn_invalid_tsize <- function(ids) {
    ids <- unique(as.character(ids))
    if (length(ids)) {
        msg <- paste0("Invalid binary TSIZE for ", length(ids),
            " evaluated patient(s): ", .short_ids(ids),
            ". Use 0 for tumor size <= 2 cm or 1 for > 2 cm. ",
            "If the stored codes are already binary, convert the unsupported column class ",
            "to a supported representation, preserving missing values. ",
            "Only TSIZE-dependent clinical research ROR scores and groups ",
            "are unavailable for these patients.")
        warning(structure(list(message = msg, call = NULL, patient_ids = ids),
            class = c("bsr_invalid_tsize_warning", "warning", "condition")))
    }
    invisible(NULL)
}

# Keep separate public invocations independent. A multi-stage or multi-method
# caller collects affected IDs, then issues the warning once before returning.
.tsize_warning_collector <- function() {
    ids <- character()
    list(
        collect = function(w) {
            ids <<- union(ids, w$patient_ids)
            invokeRestart("muffleWarning")
        },
        warn = function() .warn_invalid_tsize(ids)
    )
}
