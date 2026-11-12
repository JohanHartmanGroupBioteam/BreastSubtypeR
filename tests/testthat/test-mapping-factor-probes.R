map01_probe_input <- function(x, probe, entrez) {
    SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = x),
        rowData = S4Vectors::DataFrame(probe = probe, ENTREZID = entrez,
            Length = rep(1000, nrow(x))),
        colData = S4Vectors::DataFrame(PatientID = colnames(x),
            series = seq_len(ncol(x)), row.names = colnames(x)))
}

test_that("Mapping resolves factor probe labels in annotation order for all methods", {
    x <- rbind(z = c(8, 2, 2), b = c(2, 8, 2), a = c(0, 0, 0),
        q = c(3, 4, 5), background = c(17, 16, 15))
    colnames(x) <- c("S3", "S1", "S2")
    entrez <- c(z = 596L, b = 596L, a = 596L, q = 332L, background = 800001L)

    for (permuted in c(FALSE, TRUE)) {
        values <- x
        labels <- rownames(x)
        if (permuted) {
            values <- x[c("q", "a", "background", "z", "b"), c("S2", "S3", "S1")]
            labels <- c("b", "background", "q", "z", "a")
        }
        se <- map01_probe_input(values, labels, unname(entrez[labels]))
        encodings <- list(ordinary = factor(labels),
            reversed = factor(labels, levels = rev(levels(factor(labels)))),
            ordered = ordered(labels),
            unused = factor(labels, levels = c("unused", labels)))

        # Match the assay by name. Both z and b tie for greatest total, SD and
        # IQR; the first in the ordered target annotation wins on the log scale.
        duplicates <- values[labels[entrez[labels] == 596L], , drop = FALSE]
        first_tied <- rownames(duplicates)[rownames(duplicates) %in% c("z", "b")][1L]
        for (method in c("max", "mean", "median", "stdev", "iqr")) {
            bcl2 <- switch(method, mean = apply(duplicates, 2L, mean, na.rm = TRUE),
                median = apply(duplicates, 2L, stats::median), values[first_tied, ])
            nc <- rbind(BCL2 = bcl2, BIRC5 = values["q", ])
            ssp <- 2^nc[c("BIRC5", "BCL2"), , drop = FALSE]
            rownames(ssp) <- c("332", "596")
            control <- Mapping(se, method = method, verbose = FALSE)
            expect_identical(as.matrix(SummarizedExperiment::assay(control$se_NC)), nc)
            expect_identical(SummarizedExperiment::assay(control$se_SSP), ssp)
            for (name in names(encodings)) {
                input <- se
                SummarizedExperiment::rowData(input)$probe <- encodings[[name]]
                before <- serialize(input, NULL)
                mapped <- Mapping(input, method = method, verbose = FALSE)
                expect_identical(mapped, control, info = paste(name, method, permuted))
                expect_identical(serialize(input, NULL), before)
            }
        }
    }
})

test_that("the factor correction retains current numeric probe and raw NA behavior", {
    x <- rbind(c(8, 2), c(3, 3), c(5, 7), c(2, 2))
    dimnames(x) <- list(c("2", "1", "3", "4"), c("S2", "S1"))
    ids <- c(596L, 800001L, 332L, 800002L)
    se <- map01_probe_input(x, rownames(x), ids)
    expect_identical(rownames(Mapping(se, verbose = FALSE)$se_NC), c("BCL2", "BIRC5"))

    # Numeric probes still use the established positional subscript. A separate
    # API decision is needed before changing or rejecting them.
    for (probe in list(as.integer(rownames(x)), as.double(rownames(x)))) {
        input <- map01_probe_input(x, probe, ids)
        mapped <- Mapping(input, verbose = FALSE)
        expected <- x["3", , drop = FALSE]
        rownames(expected) <- "BIRC5"
        expect_identical(as.matrix(SummarizedExperiment::assay(mapped$se_NC)), expected)
    }

    x[1L, 1L] <- NA_real_
    for (probe in list(rownames(x), factor(rownames(x)))) {
        input <- map01_probe_input(x, probe, ids)
        for (impute in c(FALSE, TRUE)) {
            expect_error(Mapping(input, RawCounts = TRUE, impute = impute, verbose = FALSE),
                "NA counts not allowed: RawCounts=TRUE requires complete counts", fixed = TRUE)
        }
    }
})
