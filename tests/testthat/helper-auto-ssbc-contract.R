# Independent contract assertions, also used to project fresh AUTO results onto
# immutable historical axes. No production restoration code is used as oracle.
expect_auto_ssbc_axes <- function(actual, native, ids, payload = TRUE) {
    testthat::expect_identical(names(actual), names(native))
    testthat::expect_identical(names(actual$outList), names(native$outList))
    testthat::expect_identical(actual$mdns, native$mdns)
    testthat::expect_identical(actual$outList$centroids, native$outList$centroids)
    projected <- actual
    for (field in c("BS.all", "score.ROR")) {
        a <- actual[[field]]; b <- native[[field]]
        key <- if (field == "BS.all") "PatientID" else "patientID"
        testthat::expect_identical(rownames(a), ids)
        testthat::expect_identical(names(a), names(b))
        testthat::expect_identical(lapply(a, class), lapply(b, class))
        testthat::expect_identical(vapply(a, typeof, ""), vapply(b, typeof, ""))
        extra <- setdiff(ids, rownames(b))
        testthat::expect_identical(as.character(a[extra, key]), extra)
        testthat::expect_true(all(is.na(a[extra, setdiff(names(a), key), drop = FALSE])))
        projected[[field]] <- a[rownames(b), , drop = FALSE]
        if (payload) testthat::expect_identical(projected[[field]], b)
    }
    for (field in c("predictions", "predictions.FourSubtype")) {
        if (!field %in% names(native$outList)) next
        a <- actual$outList[[field]]; b <- native$outList[[field]]
        testthat::expect_identical(names(a), ids)
        testthat::expect_identical(typeof(a), typeof(b))
        testthat::expect_true(all(is.na(a[setdiff(ids, names(b))])))
        projected$outList[[field]] <- a[names(b)]
        if (payload) testthat::expect_identical(projected$outList[[field]], b)
    }
    for (field in c("distances", "dist.RORSubtype", "testData")) {
        a <- actual$outList[[field]]; b <- native$outList[[field]]
        columns <- field == "testData"
        testthat::expect_true(is.matrix(a))
        testthat::expect_identical(typeof(a), typeof(b))
        testthat::expect_identical(if (columns) rownames(a) else colnames(a),
                                  if (columns) rownames(b) else colnames(b))
        testthat::expect_identical(if (columns) colnames(a) else rownames(a), ids)
        old_ids <- if (columns) colnames(b) else rownames(b)
        extra <- setdiff(ids, old_ids)
        added <- if (columns) a[, extra, drop = FALSE] else a[extra, , drop = FALSE]
        testthat::expect_true(all(is.na(added)))
        projected$outList[[field]] <- if (columns) a[, old_ids, drop = FALSE] else
            a[old_ids, , drop = FALSE]
        if (payload) testthat::expect_identical(projected$outList[[field]], b)
    }
    if (payload) testthat::expect_identical(projected, native)
    projected
}

auto_ssbc_historical_view <- function(actual, frozen) {
    # Only #233's explicitly approved axes are projected. All old payloads,
    # other fields/methods and aggregates still go through the existing oracle.
    for (method in intersect(c("ssBC", "ssBC.v2"), names(frozen$results))) {
        actual$results[[method]] <- expect_auto_ssbc_axes(actual$results[[method]],
            frozen$results[[method]], rownames(frozen$res_subtypes), payload = FALSE)
    }
    actual
}
