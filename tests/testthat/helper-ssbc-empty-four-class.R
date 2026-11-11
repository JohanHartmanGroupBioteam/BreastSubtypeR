## NCP-02 empty four-class ssBC return. Two independent constructions: one
## written out field by field, and one derived from the native five-class
## return. Every sample is unevaluated, so every payload is missing.
ncp02_model <- function() {
    e <- new.env(parent = emptyenv())
    utils::data("BreastSubtypeRobj", package = "BreastSubtypeR", envir = e)
    e$BreastSubtypeRobj
}

ncp02_explicit_empty <- function(ids, clinical, subtype = TRUE) {
    model <- ncp02_model()
    classes <- c("Basal", "Her2", "LumA", "LumB", "Normal")
    fields <- c(classes, "Call", if (subtype) "Call.Subtype", "Confidence",
        "ROR-S (Subtype Only)", "ROR-S Group (Subtype Only)",
        "ROR-P (Subtype + Prolif)", "ROR-P Group (Subtype + Prolif)",
        if (clinical) c("ROR-C (Subtype + Clinic)", "ROR-C Group (Subtype + Clinic)",
            "ROR-PC (Subtype + Clinic + Prolif)", "ROR-PC Group (Subtype + Clinic + Prolif)",
            "ROR-PC (Subtype + Clinic + Prolif.Subtype)",
            "ROR-PC Group (Subtype + Clinic + Prolif.Subtype)"),
        "ER", "HER2")
    n <- length(ids)
    na <- rep(NA, n)
    calls <- if (subtype) data.frame(PatientID = ids, BS = na, BS.Subtype = na, row.names = ids) else
        data.frame(PatientID = ids, BS = na, row.names = ids)
    scores <- as.data.frame(stats::setNames(rep(list(na), length(fields)), fields),
        row.names = ids, check.names = FALSE)
    no_call <- stats::setNames(rep(NA_character_, n), ids)
    out <- c(list(predictions = no_call), if (subtype) list(predictions.FourSubtype = no_call),
        list(testData = matrix(numeric(), 0L, n, dimnames = list(NULL, ids)),
            distances = matrix(NA_real_, n, 5L, dimnames = list(ids, classes)),
            dist.RORSubtype = matrix(NA_real_, n, 4L, dimnames = list(ids, classes[1:4])),
            centroids = model$centroid))
    list(BS.all = calls, score.ROR = scores, mdns = model$ssBC.subgroupQuantile, outList = out)
}

ncp02_four_from_five <- function(five) {
    ids <- rownames(five$BS.all)
    four <- five
    four$BS.all$BS.Subtype <- rep(NA, length(ids))
    at <- match("Call", names(five$score.ROR))
    four$score.ROR <- cbind(five$score.ROR[seq_len(at)], Call.Subtype = rep(NA, length(ids)),
        five$score.ROR[-seq_len(at)])
    four$outList <- append(five$outList,
        list(predictions.FourSubtype = five$outList$predictions), after = 1L)
    four
}

## OSLO2EMIT0 public expression with synthetic annotations. `ix` permutes the
## assay; PatientID labels always follow the assay order.
ncp02_input <- function(ix = c(seq(2, 36, 2), seq(1, 35, 2)), factor_ids = FALSE) {
    e <- new.env(parent = emptyenv())
    utils::data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = e)
    x <- lapply(e$OSLO2EMIT0obj$data_input, function(z) {
        z <- z[, ix, drop = FALSE]
        ids <- colnames(z)
        z$PatientID <- if (factor_ids) factor(ids, levels = c("unused", rev(ids))) else ids
        z$ER <- rep(c("ER+", "ER-"), length.out = ncol(z))
        z$HER2 <- rep(c("HER2-", "HER2-", "HER2+"), length.out = ncol(z))
        z$TN <- ifelse(z$ER == "ER-" & z$HER2 == "HER2-", "TN", "nonTN")
        z$TSIZE <- rep(c(0, 1), length.out = ncol(z))
        z$NODE <- rep(c(4, 0, 2, 7), length.out = ncol(z))
        z
    })
    x
}

## Run a public call; return the value and the warning classes and messages.
ncp02_capture <- function(expr) {
    warnings <- list()
    value <- withCallingHandlers(suppressMessages(force(expr)), warning = function(w) {
        warnings[[length(warnings) + 1L]] <<- w
        invokeRestart("muffleWarning")
    })
    list(value = value, classes = vapply(warnings, function(w) class(w)[1L], ""),
        messages = vapply(warnings, conditionMessage, ""), warnings = warnings)
}

## A call expected to raise no warning at all; returns its value.
ncp02_no_warning <- function(expr) {
    captured <- ncp02_capture(expr)
    testthat::expect_identical(captured$classes, character())
    captured$value
}
