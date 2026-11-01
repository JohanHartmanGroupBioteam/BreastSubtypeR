#!/usr/bin/env Rscript
# Build compact historical fixtures from immutable originals, never classifiers.
# Usage: Rscript compact-historical-fixtures.R --historical-dir=/originals
#        --package-dir=/checkout --out=/new/fixture-build
options(stringsAsFactors = FALSE)
args <- commandArgs(TRUE)
opt <- setNames(sub("^[^=]+=", "", args), sub("=.*$", "", sub("^--", "", args)))
if (!all(c("historical-dir", "package-dir", "out") %in% names(opt))) {
    stop("Required: --historical-dir=DIR --package-dir=DIR --out=NEW_DIR")
}
original_dir <- normalizePath(opt[["historical-dir"]], winslash = "/", mustWork = TRUE)
package_dir <- normalizePath(opt[["package-dir"]], winslash = "/", mustWork = TRUE)
out <- opt[["out"]]
if (file.exists(out) || dir.exists(out)) stop("Output already exists.")
if (!dir.create(out, recursive = TRUE)) stop("Could not create output.")
suppressPackageStartupMessages(library(SummarizedExperiment))
if (!requireNamespace("digest", quietly = TRUE)) stop("Builder requires digest.")
source(file.path(package_dir, "tests", "testthat", "helper-historical-examples.R"))

source_sha <- c(OSLO2EMIT0obj =
    "43393c2bf06f0ea44b262d5525b75154ab99162fac5de898433556b6805450ec",
    TCGABRCAobj =
    "b94309c6f915d5f5fcb18368f3876f281b2111acf50eeb25b167f201f499a45d")
public_sha <- c(OSLO2EMIT0obj =
    "1b7c1a6a2b1fb21f4968f185e6a64e90c0282ccdd7205aaf78a04cbdf75d6439",
    TCGABRCAobj =
    "b4ac2d8c898557453d69f5c7a3da81a1fbce4db4486bf0623adffcb1597f5d7d")
sha <- function(path) digest::digest(file = path, algo = "sha256")
raw_sha <- function(bytes) digest::digest(bytes, serialize = FALSE, algo = "sha256")
checks <- list()
check <- function(object, requirement, value) {
    checks[[length(checks) + 1L]] <<- data.frame(object, requirement,
        status = if (isTRUE(value)) "PASS" else "FAIL")
    if (!isTRUE(value)) stop(object, ": ", requirement)
}
read_object <- function(path, name) {
    e <- new.env(parent = emptyenv())
    stopifnot(identical(load(path, envir = e), name))
    e[[name]]
}
read_payload <- function(path) {
    con <- gzfile(path, "rb")
    on.exit(close(con))
    readBin(con, "raw", n = .Machine$integer.max)
}
manifest <- list()
for (name in names(source_sha)) {
    original_path <- file.path(original_dir, paste0(name, ".RData"))
    public_path <- file.path(package_dir, "data", paste0(name, ".RData"))
    check(name, "original_container_sha256", identical(sha(original_path), unname(source_sha[name])))
    check(name, "approved_public_container_sha256", identical(sha(public_path), unname(public_sha[name])))
    old <- read_object(original_path, name)
    current <- read_object(public_path, name)
    reference <- historical_fixture_reference(name)
    old_assays <- assays(old$se_obj, withDimnames = FALSE)
    check(name, "one_original_assay", length(old_assays) == 1L)
    check(name, "original_complete_object_digest", identical(
        historical_canonical_md5(old), reference$object_md5))
    check(name, "original_assay_digest", identical(
        historical_canonical_md5(old_assays[[1L]]), reference$assay_md5))
    first <- seq_len(nrow(old$se_obj))
    check(name, "approved_public_prefix_is_original_assay", identical(
        assays(current$se_obj, withDimnames = FALSE)[[1L]][first, , drop = FALSE],
        old_assays[[1L]]))
    fixture <- list(format = "issue199-shared-assay-v1", object_name = name,
        original_nrow = nrow(old$se_obj), assay_names = names(old_assays),
        historical = old)
    assays(fixture$historical$se_obj, withDimnames = FALSE) <- S4Vectors::SimpleList()
    check(name, "assayless_template_is_valid", isTRUE(methods::validObject(fixture$historical$se_obj)))
    target <- file.path(out, paste0(name, ".rds"))
    con <- xzfile(target, "wb", compression = 9L)
    tryCatch(invisible(serialize(fixture, con, version = 3L)), finally = close(con))
    saved <- readRDS(target)
    restored <- reconstruct_historical_example(saved, current)
    check(name, "entire_historical_object_identical", identical(restored, old))
    check(name, "historical_se_obj_identical", identical(restored$se_obj, old$se_obj))
    check(name, "historical_mapped_inputs_identical", identical(restored$data_input, old$data_input))
    check(name, "historical_frozen_results_identical", identical(restored$res, old$res))
    check(name, "historical_annotations_identical", identical(rowData(restored$se_obj), rowData(old$se_obj)))
    check(name, "historical_sample_metadata_identical", identical(colData(restored$se_obj), colData(old$se_obj)))
    check(name, "historical_assay_storage_type_identical", identical(
        typeof(assay(restored$se_obj)), typeof(assay(old$se_obj))))
    manifest[[name]] <- data.frame(file = basename(target),
        sha256 = sha(target), md5 = unname(tools::md5sum(target)),
        bytes = file.info(target)$size, object = name,
        original_file = basename(original_path), original_sha256 = sha(original_path),
        original_md5 = unname(tools::md5sum(original_path)),
        original_bytes = file.info(original_path)$size,
        original_payload_sha256 = raw_sha(read_payload(original_path)),
        canonical_object_sha256 = raw_sha(historical_canonical_bytes(old)),
        canonical_object_md5 = reference$object_md5,
        canonical_assay_sha256 = raw_sha(historical_canonical_bytes(old_assays[[1L]])),
        canonical_assay_md5 = reference$assay_md5,
        rows = nrow(old$se_obj), columns = ncol(old$se_obj))
}
write.table(do.call(rbind, manifest), file.path(out, "manifest.tsv"),
    sep = "\t", quote = TRUE, row.names = FALSE)
write.table(do.call(rbind, checks), file.path(out, "construction-checks.tsv"),
    sep = "\t", quote = TRUE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
versions <- as.data.frame(installed.packages()[intersect(loadedNamespaces(),
    rownames(installed.packages())), c("Package", "Version", "Built"), drop = FALSE])
write.table(versions, file.path(out, "software.tsv"), sep = "\t",
    quote = TRUE, row.names = FALSE)
cat(sum(vapply(checks, function(x) x$status == "PASS", logical(1))), "checks PASS\n")
