# All donors tie on S2:S5, so KNN takes the first ten Entrez groups.
# The expected probe names below are independent of the factor level codes.
test_that("Mapping uses probe labels for KNN targets and background donors", {
    donor_ids <- setNames(800001:800012, paste0("donor", 1:12))
    for (case in c("background duplicates", "observed target", "missing target")) {
        ids <- c(target = 596L, donor_ids)
        x <- cbind(S1 = seq_along(ids), S2 = 2, S3 = 3, S4 = 4, S5 = 5)
        rownames(x) <- names(ids)
        x["target", "S1"] <- NA_real_
        targets <- c(BCL2 = "target")
        neighbors <- names(donor_ids)[1:10]

        if (case == "background duplicates") {
            # Max selects the new probes for donors 1 and 2. mean(c(82, 93, 4:11)) is 23.5.
            x <- rbind(x, first = c(82, 2:5), second = c(93, 2:5))
            ids <- c(ids, first = 800001L, second = 800002L)
            neighbors[1:2] <- c("first", "second")
            swap <- c("first", "second")
        } else if (case == "observed target") {
            # Max selects the complete BCL2 probe. S1 must remain 60 with KNN on or off.
            x <- rbind(x, first = c(60, 2:5), second = c(93, 2:5))
            ids <- c(ids, first = 596L, second = 800001L)
            targets["BCL2"] <- "first"
            swap <- c("first", "second")
        } else {
            # Keep a second represented target when a positional lookup drops BCL2.
            ids["donor1"] <- 332L
            targets <- c(targets, BIRC5 = "donor1")
            swap <- c("target", "donor2")
        }

        keys <- rownames(x)
        levels <- keys
        levels[match(swap, keys)] <- rev(swap)
        se <- SummarizedExperiment::SummarizedExperiment(assays = list(expr = x),
            rowData = S4Vectors::DataFrame(probe = keys, ENTREZID = unname(ids)),
            colData = S4Vectors::DataFrame(PatientID = colnames(x)))
        factor_se <- se
        SummarizedExperiment::rowData(factor_se)$probe <- factor(keys, levels = levels)
        before <- serialize(factor_se, NULL)

        for (impute in c(FALSE, TRUE)) {
            expected <- x[unname(targets), , drop = FALSE]
            rownames(expected) <- names(targets)
            if (impute && is.na(expected["BCL2", "S1"])) {
                expected["BCL2", "S1"] <- mean(x[neighbors, "S1"])
            }
            # The SSP panel orders BIRC5 (332) before BCL2 (596).
            ssp_ids <- c(BIRC5 = "332", BCL2 = "596")
            included <- intersect(names(ssp_ids), rownames(expected))
            expected_ssp <- 2^expected[included, , drop = FALSE]
            rownames(expected_ssp) <- unname(ssp_ids[included])

            reference <- Mapping(se, impute = impute, verbose = FALSE)
            result <- Mapping(factor_se, impute = impute, verbose = FALSE)
            info <- paste(case, "impute =", impute)
            expect_identical(as.matrix(SummarizedExperiment::assay(reference$se_NC)),
                expected, info = info)
            expect_identical(as.matrix(SummarizedExperiment::assay(reference$se_SSP)),
                expected_ssp, info = info)
            expect_identical(result, reference, info = info)
            expect_identical(serialize(factor_se, NULL), before, info = info)
        }
    }
})
