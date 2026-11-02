#!/usr/bin/env Rscript
# Issue #199: export the six Shiny tables from the approved complete objects.
# Requires SummarizedExperiment and digest; no classifier or Mapping is called.
# Usage: Rscript export-shiny-examples.R --package-root=PATH
#        --historical-oslo-clinical=PATH [--output-dir=PATH] [--manifest=PATH]
args <- commandArgs(TRUE)
get_arg <- function(key, default = NULL) {
    value <- args[startsWith(args, paste0("--", key, "="))]
    if (length(value) > 1L) stop("Duplicate argument: ", key)
    if (!length(value)) return(default)
    substring(value, nchar(key) + 4L)
}
root <- get_arg("package-root", ".")
legacy_path <- get_arg("historical-oslo-clinical")
if (is.null(legacy_path)) stop("--historical-oslo-clinical is required")
out <- get_arg("output-dir", file.path(root, "inst", "RshinyTest"))
manifest <- get_arg("manifest")
stopifnot(requireNamespace("SummarizedExperiment", quietly = TRUE),
          requireNamespace("digest", quietly = TRUE), dir.exists(root))
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
expected <- c(
    OSLO2EMIT0obj = "1b7c1a6a2b1fb21f4968f185e6a64e90c0282ccdd7205aaf78a04cbdf75d6439",
    TCGABRCAobj = "b4ac2d8c898557453d69f5c7a3da81a1fbce4db4486bf0623adffcb1597f5d7d")
legacy_sha <- "00a24cfde1bc7837a5f5cc46d9c2e4ce55a60a624b99caf5ee276005d89b8996"
stopifnot(identical(sha(legacy_path), legacy_sha))
objects <- lapply(names(expected), function(n) {
    path <- file.path(root, "data", paste0(n, ".RData"))
    stopifnot(identical(sha(path), unname(expected[n])))
    e <- new.env(parent = emptyenv())
    stopifnot(identical(load(path, envir = e), n))
    e[[n]]
})
names(objects) <- names(expected)
legacy <- read.delim(legacy_path, check.names = FALSE, colClasses = "character")
legacy_cols <- c("PatientID", "ER", "PR", "HER2", "TSIZE", "NODE", "IHC",
                 "Prosigna.Subtype", "Prosigna.ROR", "MKI67.corrected",
                 "Ki67.estimated", "IHC.estimated")
stopifnot(identical(names(legacy), legacy_cols),
          !anyNA(legacy$PatientID), !anyDuplicated(legacy$PatientID))
# 17 significant decimal digits round-trip IEEE doubles on supported R readers.
# Integer coverage sums are emitted without changing their numerical values.
numeric_text <- function(x) {
    result <- sprintf("%.17g", x)
    result[is.na(x)] <- NA_character_
    result
}
as_text <- function(df) {
    df[] <- lapply(df, function(x) if (is.numeric(x)) numeric_text(x) else as.character(x))
    df
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
records <- list()
for (name in names(objects)) {
    se <- objects[[name]]$se_obj
    raw <- identical(name, "TCGABRCAobj")
    prefix <- if (raw) "TCGABRCAsub" else "OSLO2EMIT0"
    x <- as.matrix(SummarizedExperiment::assay(se))
    ann <- as.data.frame(SummarizedExperiment::rowData(se))
    clinical <- as.data.frame(SummarizedExperiment::colData(se))
    stopifnot(identical(ann$probe, rownames(x)),
              identical(clinical$PatientID, colnames(x)),
              !anyDuplicated(rownames(x)), !anyDuplicated(colnames(x)),
              all(is.finite(x)), if (raw) all(x >= 0 & x == floor(x)) else TRUE)
    if (!raw) {
        stopifnot(setequal(legacy$PatientID, clinical$PatientID))
        legacy_aligned <- legacy[match(clinical$PatientID, legacy$PatientID), , drop = FALSE]
        stopifnot(!anyNA(legacy_aligned$PatientID),
                  identical(legacy_aligned$PatientID, clinical$PatientID))
        # Six candidate fields are authoritative; retain six ancillary fields
        # from the pinned historical public table to preserve its 12-column schema.
        clinical <- cbind(as_text(clinical),
                          legacy_aligned[, setdiff(legacy_cols, names(clinical)), drop = FALSE])
        stopifnot(identical(names(clinical), legacy_cols))
    }
    x_text <- matrix(numeric_text(as.numeric(x)), nrow = nrow(x), dimnames = dimnames(x))
    tables <- list(
        GEX = list(value = x_text, suffix = if (raw) "_GEX_RawCounts.txt" else "_GEX_log2.FPKM.txt",
                   row_names = TRUE, col_names = NA),
        annotation = list(value = as_text(ann), suffix = "_anno.txt",
                          row_names = FALSE, col_names = TRUE),
        clinical = list(value = as_text(clinical), suffix = "_clinical.txt",
                        row_names = FALSE, col_names = TRUE))
    for (kind in names(tables)) {
        tab <- tables[[kind]]
        path <- file.path(out, paste0(prefix, tab$suffix))
        write.table(tab$value, path, sep = "\t", quote = FALSE,
                    row.names = tab$row_names, col.names = tab$col_names,
                    na = "NA", eol = "\n", fileEncoding = "UTF-8")
        records[[length(records) + 1L]] <- data.frame(
            file = basename(path), object = name, source_sha256 = unname(expected[name]),
            kind = kind, rows = nrow(tab$value), columns = ncol(tab$value),
            units = if (kind != "GEX") "annotation" else if (raw)
                "recount3 gene-level base-pair coverage sums" else
                "log2(protocol-adjusted FPKM + 1)",
            sha256 = sha(path), bytes = file.info(path)$size)
    }
    reread <- as.matrix(read.delim(file.path(out, paste0(prefix, tables$GEX$suffix)),
                                  row.names = 1L, check.names = FALSE))
    storage.mode(reread) <- "double"
    storage.mode(x) <- "double"
    stopifnot(identical(reread, x))
}
if (!is.null(manifest)) {
    dir.create(dirname(manifest), recursive = TRUE, showWarnings = FALSE)
    write.table(do.call(rbind, records), manifest, sep = "\t", quote = FALSE,
                row.names = FALSE, eol = "\n")
}
print(do.call(rbind, records), row.names = FALSE)

