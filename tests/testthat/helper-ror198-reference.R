## Independent test oracle for #198. Formula/coefficients transcribed from the
## original Parker scripts bundled in official CRAN PCAPAM50 1.0.4:
## https://cran.r-project.org/src/contrib/PCAPAM50_1.0.4.tar.gz
## SHA256 2ea54117a6f6557396607fb619280e7f1c9a4e56a3b240379c3d18b8ae6fc002.
## Use mean()/sd() explicitly so a wrong scale() axis or panel-only scaling in
## production cannot also change the expected result. No classifier is replaced.
ror198_panel <- c("CCNB1", "UBE2C", "BIRC5", "KNTC2", "CDC20", "PTTG1",
    "RRM2", "MKI67", "TYMS", "CEP55", "CDCA1")
ror198_changed <- c("ROR-P (Subtype + Prolif)", "ROR-P Group (Subtype + Prolif)",
    "ROR-PC (Subtype + Clinic + Prolif)", "ROR-PC Group (Subtype + Clinic + Prolif)")

ror198_proliferation <- function(x) {
    x <- as.matrix(x)
    present <- rownames(x) %in% ror198_panel
    vapply(seq_len(ncol(x)), function(j) {
        sample <- x[, j]
        observed <- sample[!is.na(sample)]
        if (length(observed) < 2L || any(!is.finite(observed))) return(NA_real_)
        spread <- stats::sd(observed)
        if (!is.finite(spread) || spread <= 0) return(NA_real_)
        p <- mean((sample[present] - mean(observed)) / spread, na.rm = TRUE)
        if (is.finite(p)) p else NA_real_
    }, numeric(1), USE.NAMES = FALSE)
}

ror198_groups <- function(raw, low, high) {
    group <- raw
    group[raw < low] <- "low"
    group[raw > low & raw < high] <- "med"
    group[raw > high] <- "high"
    group[!is.finite(raw)] <- NA_character_
    as.character(group) # retain the deferred exact-boundary text (#172)
}

ror198_reference <- function(out, clinical = NULL) {
    d <- out$distances
    p <- setNames(ror198_proliferation(out$testData), colnames(out$testData))
    p <- p[rownames(d)]
    raw <- -0.0009299747 * d[, 1] + 0.0692289192 * d[, 2] -
        0.0951505484 * d[, 3] + 0.0493487685 * d[, 4] + 0.3385116381 * p
    scaled <- function(x) {
        x <- 100 * (x + 0.35) / 0.85
        x[!is.finite(x)] <- NA_real_
        x
    }
    result <- data.frame(patientID = rownames(d),
        "ROR-P (Subtype + Prolif)" = scaled(raw),
        "ROR-P Group (Subtype + Prolif)" = ror198_groups(raw, -0.25, 0.1),
        check.names = FALSE)
    if (!is.null(clinical)) {
        tsize <- as.numeric(as.character(clinical$TSIZE[
            match(rownames(d), as.character(clinical$PatientID))]))
        rawpc <- -0.009383416 * d[, 1] + 0.073725503 * d[, 2] -
            0.090436516 * d[, 3] + 0.053013865 * d[, 4] +
            0.131605960 * tsize + 0.327259375 * p
        result[[ror198_changed[3]]] <- scaled(rawpc)
        result[[ror198_changed[4]]] <- ror198_groups(rawpc, -0.2, 0.2)
    }
    for (pair in list(ror198_changed[1:2], ror198_changed[3:4])) {
        if (pair[1] %in% names(result)) result[[pair[2]]][is.na(result[[pair[1]]])] <- NA_character_
    }
    result
}

ror198_expect_reference <- function(result, clinical = NULL) {
    expected <- ror198_reference(result$outList, clinical)
    ids <- as.character(result$score.ROR$patientID)
    expect_false(anyNA(ids) || anyDuplicated(ids) > 0L)
    expect_setequal(ids, expected$patientID)
    expected <- expected[match(ids, expected$patientID), , drop = FALSE]
    for (nm in intersect(ror198_changed, names(expected))) {
        if (is.numeric(expected[[nm]])) {
            expect_equal(unname(result$score.ROR[[nm]]), unname(expected[[nm]]),
                tolerance = 1e-12, info = nm)
        } else {
            expect_identical(as.character(result$score.ROR[[nm]]),
                as.character(expected[[nm]]), info = nm)
        }
    }
}

ror198_run <- function(method, se, clinical = TRUE, subtype = TRUE, iterations = 3L) {
    common <- list(se_obj = se, hasClinical = clinical, Subtype = subtype)
    switch(method,
        parker.original = do.call(BS_parker, c(common, list(calibration = "Internal", internal = "medianCtr"))),
        genefu.scale = do.call(BS_parker, c(common, list(calibration = "Internal", internal = "meanCtr"))),
        genefu.robust = do.call(BS_parker, c(common, list(calibration = "Internal", internal = "qCtr"))),
        cIHC = do.call(BS_cIHC, c(common, list(seed = 118))),
        cIHC.itr = do.call(BS_cIHC.itr, c(common, list(seed = 118, iteration = iterations))),
        PCAPAM50 = do.call(BS_PCAPAM50, c(common, list(seed = 118))),
        ssBC = do.call(BS_ssBC, c(common, list(s = "ER"))),
        ssBC.v2 = do.call(BS_ssBC, c(common, list(s = "ER.v2"))))
}
