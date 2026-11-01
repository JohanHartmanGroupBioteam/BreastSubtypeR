# Historical regression objects keep independent annotations, metadata,
# mapped inputs and frozen results. Only unchanged assay prefixes are shared
# with the public examples; independent frozen digests reject any drift.
historical_fixture_reference <- function(name) {
    switch(name,
        OSLO2EMIT0obj = list(nrow = 356L, ncol = 102L,
            assay_md5 = "d7e5e33355b130dd552f989643a4ce5e",
            object_md5 = "b074b7fc768e44df13f84c0beabaeff3"),
        TCGABRCAobj = list(nrow = 5000L, ncol = 100L,
            assay_md5 = "2b0e969408e8983d3291f3ebc49ae82c",
            object_md5 = "87926c829eb05788d9f060716f79b133"),
        stop("Unknown historical fixture name.", call. = FALSE))
}

historical_canonical_bytes <- function(object) {
    bytes <- serialize(object, NULL, version = 2L, xdr = TRUE)
    stopifnot(identical(bytes[1:2], charToRaw("X\n")),
        identical(bytes[3:6], as.raw(c(0, 0, 0, 2))))
    # Normalize only the writer/minimum-reader version integers in the XDR
    # header. Values, types, attributes, ordering and metadata are untouched.
    bytes[7:14] <- as.raw(0)
    bytes
}

historical_canonical_md5 <- function(object) {
    path <- tempfile("issue199-fixture-digest-")
    on.exit(unlink(path), add = TRUE)
    writeBin(historical_canonical_bytes(object), path)
    unname(tools::md5sum(path))
}

reconstruct_historical_example <- function(fixture, current) {
    if (!identical(fixture$format, "issue199-shared-assay-v1")) {
        stop("Unsupported historical fixture format.", call. = FALSE)
    }
    reference <- historical_fixture_reference(fixture$object_name)
    old <- fixture$historical
    if (!identical(names(old), c("se_obj", "data_input", "res")) ||
        !identical(dim(old$se_obj), c(reference$nrow, reference$ncol)) ||
        length(SummarizedExperiment::assays(old$se_obj)) != 0L ||
        !identical(fixture$original_nrow, reference$nrow)) {
        stop("Historical fixture schema integrity failure.", call. = FALSE)
    }
    first <- seq_len(reference$nrow)
    if (nrow(current$se_obj) < reference$nrow ||
        !identical(colnames(current$se_obj), colnames(old$se_obj)) ||
        !identical(rownames(current$se_obj)[first], rownames(old$se_obj))) {
        stop("Historical assay row/sample order integrity failure.", call. = FALSE)
    }
    supplied <- SummarizedExperiment::assays(current$se_obj,
        withDimnames = FALSE)
    if (length(supplied) != 1L ||
        !identical(names(supplied), fixture$assay_names)) {
        stop("Historical assay schema integrity failure.", call. = FALSE)
    }
    values <- supplied[[1L]][first, , drop = FALSE]
    if (!identical(historical_canonical_md5(values), reference$assay_md5)) {
        stop("Historical assay content/type integrity failure.", call. = FALSE)
    }
    restored <- S4Vectors::SimpleList(values)
    names(restored) <- fixture$assay_names
    SummarizedExperiment::assays(old$se_obj, withDimnames = FALSE) <- restored
    if (!identical(historical_canonical_md5(old), reference$object_md5)) {
        stop("Historical complete-object integrity failure.", call. = FALSE)
    }
    old
}

current_example <- function(name) {
    env <- new.env(parent = emptyenv())
    utils::data(list = name, package = "BreastSubtypeR", envir = env)
    env[[name]]
}

load_historical_example <- function(name, envir = parent.frame()) {
    stopifnot(length(name) == 1L,
        name %in% c("OSLO2EMIT0obj", "TCGABRCAobj"), is.environment(envir))
    path <- testthat::test_path("fixtures", "historical-cb7a7c70",
        paste0(name, ".rds"))
    fixture <- readRDS(path)
    stopifnot(identical(fixture$object_name, name))
    object <- reconstruct_historical_example(fixture, current_example(name))
    assign(name, object, envir = envir)
    invisible(name)
}

historical_example <- function(name) {
    env <- new.env(parent = emptyenv())
    load_historical_example(name, envir = env)
    env[[name]]
}
