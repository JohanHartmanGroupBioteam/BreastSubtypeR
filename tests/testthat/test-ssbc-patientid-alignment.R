# Ordered sample-key regressions; public expression and synthetic metadata.
local({
fx <- local({
    e <- new.env(); data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = e)
    se <- e$OSLO2EMIT0obj$data_input$se_NC[, 1:24]
    ids <- sprintf("k%02d", 1:24); colnames(se) <- ids
    SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(PatientID = ids,
        ER = rep(c("ER+", "ER+", "ER-", "ER-", "ER+", "ER-"), 4),
        HER2 = rep(c("HER2-", "HER2+", "HER2-", "HER2+"), 6),
        TSIZE = rep(c(0, 1, 1), 8), NODE = rep(c(0, 2, 5, 1), 6), row.names = ids)
    se$TN <- ifelse(se$ER == "ER-" & se$HER2 == "HER2-", "TN", "nonTN")
    se
})
ids <- colnames(fx)
cell <- paste(fx$ER, fx$HER2, fx$TN)
within <- integer(length(ids))
for (i in split(seq_along(ids), cell)) within[i] <- i[c(2:length(i), 1)]
P_CROSS <- ids[c(2:24, 1)]          # crosses receptor strata
P_WITHIN <- ids[within]             # keeps every ER x HER2 x TN membership
stopifnot(setequal(P_CROSS, ids), setequal(P_WITHIN, ids), all(P_WITHIN != ids),
          identical(cell[match(P_WITHIN, ids)], cell), !identical(cell[match(P_CROSS, ids)], cell))
swap <- function(i, j) { p <- ids; p[c(i, j)] <- ids[c(j, i)]; p }
PARTIAL <- list(swap(23, 24), swap(11, 12), swap(1, 24), swap(1, 4))   # few evaluated rows; (1,4) crosses ER inside nonTN
grid <- expand.grid(s = c("ER", "ER.v2", "TN", "TN.v2"), st = c(FALSE, TRUE), cl = c(FALSE, TRUE), stringsAsFactors = FALSE)
canon <- function(r, keys) {
    r$BS.all <- r$BS.all[keys, , drop = FALSE]; r$score.ROR <- r$score.ROR[keys, , drop = FALSE]
    for (k in intersect(c("predictions", "predictions.FourSubtype"), names(r$outList))) r$outList[[k]] <- r$outList[[k]][keys]
    for (k in c("distances", "dist.RORSubtype")) r$outList[[k]] <- r$outList[[k]][keys, , drop = FALSE]
    r$outList$testData <- r$outList$testData[, keys, drop = FALSE]; r
}
quiet <- function(expr) {
    withCallingHandlers(force(expr),
        bsr_experimental_ror_pc_warning = function(w) {
            expect_identical(conditionMessage(w), ror221_notice)
            invokeRestart("muffleWarning")
        },
        message = function(m) {
            expect_match(conditionMessage(m), "ssBC.*running|Running AUTO mode for subtyping")
            invokeRestart("muffleMessage")
        })
}
pid_expect_error <- function(expr, regexp = "PatientID", info = NULL) {
    seen <- list()
    error <- tryCatch(withCallingHandlers(expr, warning = function(w) {
        seen[[length(seen) + 1L]] <<- w
        invokeRestart("muffleWarning")
    }), error = identity)
    expect_s3_class(error, "bsr_patientid_mismatch_error")
    expect_length(seen, 0L)
    if (inherits(error, "error")) {
        expect_match(conditionMessage(error), regexp)
        expect_match(conditionMessage(error), "assay")
        expect_match(conditionMessage(error), "complete, unique and in the same order", fixed = TRUE)
        expect_null(conditionCall(error))
    }
    invisible(error)
}
id_error <- "PatientID"

test_that("order-only PatientID permutations are refused for every strategy and flag", {
    # Partial swaps are used only where the swapped rows are evaluated (all PARTIAL rows
    # are nonTN, i.e. excluded under TN.v2); excluded-row cases belong to block [D1].
    for (i in seq_len(nrow(grid))) for (p in c(list(P_CROSS, P_WITHIN), if (grid$s[i] != "TN.v2") PARTIAL)) {
        bad <- fx; bad$PatientID <- p
        pid_expect_error(quiet(BS_ssBC(bad, grid$s[i], Subtype = grid$st[i], hasClinical = grid$cl[i])), id_error,
                     info = paste(grid[i, ], collapse = "/"))
    }
})

test_that("a swapped singleton stratum is refused instead of moving the no-call", {
    q <- fx; q$ER <- rep("ER+", 24); q$ER[3] <- "ER-"; q$PatientID <- swap(3, 5)
    pid_expect_error(quiet(BS_ssBC(q, "ER", Subtype = TRUE, hasClinical = TRUE)), id_error)
})

test_that("factor PatientID with contradictory labels is refused (labels, not codes, matter)", {
    bad <- fx; bad$PatientID <- factor(P_CROSS, levels = c("unused", rev(ids)))
    pid_expect_error(quiet(BS_ssBC(bad, "ER", Subtype = TRUE, hasClinical = TRUE)), id_error)
})

test_that("refusal happens before ssBC centering and clinical scoring (public behaviour)", {
    # Sentinels: a singleton in the LAST ER stratum warns only after the first stratum
    # was centred; an invalid TSIZE warns only during clinical scoring. A refusal before
    # calculation signals neither.
    q <- fx; q$ER <- rep("ER-", 24); q$ER[24] <- "ER+"; q$TSIZE[5] <- 2
    for (p in list(P_CROSS, swap(1, 2))) {
        q$PatientID <- p; seen <- character()
        err <- tryCatch(withCallingHandlers(suppressMessages(BS_ssBC(q, "ER", Subtype = TRUE, hasClinical = TRUE)),
            warning = function(w) { seen <<- c(seen, paste(class(w), collapse = "/"), conditionMessage(w)); invokeRestart("muffleWarning") }),
            error = identity)
        expect_s3_class(err, "bsr_patientid_mismatch_error")
        expect_length(seen, 0L)
    }
})

test_that("caller object and RNG are unchanged when the identity check refuses", {
    bad <- fx; bad$PatientID <- P_WITHIN
    before <- serialize(bad, NULL); set.seed(7); rng <- .Random.seed
    pid_expect_error(quiet(BS_ssBC(bad, "ER.v2", Subtype = TRUE, hasClinical = TRUE)), id_error)
    expect_identical(serialize(bad, NULL), before); expect_identical(.Random.seed, rng)
    d <- list(se_NC = bad, se_SSP = bad); before <- serialize(d, NULL); set.seed(8); rng <- .Random.seed
    pid_expect_error(quiet(BS_Multi(d, methods = c("ssBC", "AIMS"), Subtype = TRUE, hasClinical = TRUE)), id_error)
    expect_identical(serialize(d, NULL), before); expect_identical(.Random.seed, rng)
})

test_that("manual BS_Multi propagates the refusal; AUTO keeps its early error", {
    bad <- fx; bad$PatientID <- P_CROSS
    d <- list(se_NC = bad, se_SSP = bad)
    pid_expect_error(quiet(BS_Multi(d, methods = c("ssBC", "ssBC.v2"), Subtype = TRUE, hasClinical = TRUE)), id_error)
    expect_error(quiet(BS_Multi(d, methods = "AUTO", Subtype = TRUE, hasClinical = TRUE)),
                 "AUTO requires PatientID labels to match the input assay sample order")
})

test_that("valid identities are preserved exactly: factor labels, joint reordering, consistent renaming, numeric labels", {
    for (i in seq_len(nrow(grid))) {
        f <- function(x) quiet(BS_ssBC(x, grid$s[i], Subtype = grid$st[i], hasClinical = grid$cl[i]))
        ref <- f(fx)
        fac <- fx; fac$PatientID <- factor(ids, levels = c("unused", rev(ids)))
        expect_identical(f(fac), ref)
        o <- c(13:24, 1:12); expect_identical(canon(f(fx[, o]), ids), ref)
        sr_same <- function(a, b) { k <- setdiff(names(b), "patientID"); expect_identical(unname(as.list(a[k])), unname(as.list(b[k]))) }
        ren <- fx; new <- paste0("x_", ids); colnames(ren) <- new; ren$PatientID <- new
        r2 <- f(ren); expect_identical(unname(r2$outList$distances), unname(ref$outList$distances))
        expect_identical(unname(r2$BS.all$BS), unname(ref$BS.all$BS)); sr_same(r2$score.ROR, ref$score.ROR)
        num <- fx; colnames(num) <- as.character(1:24); num$PatientID <- 1:24
        r3 <- f(num); expect_identical(unname(r3$outList$distances), unname(ref$outList$distances)); sr_same(r3$score.ROR, ref$score.ROR)
        for (store in list(matrix(ids, ncol = 1), as.list(ids), S4Vectors::Rle(ids), I(ids))) {
            st <- fx; cd <- SummarizedExperiment::colData(st); cd$PatientID <- store; SummarizedExperiment::colData(st) <- cd
            expect_identical(f(st), ref)
        }
    }
})

test_that("ROR-01 provenance: invalid TSIZE masks and names the physical sample on aligned input", {
    q <- fx; q$TSIZE[5] <- 2
    w <- NULL
    r <- withCallingHandlers(suppressMessages(BS_ssBC(q, "ER", Subtype = TRUE, hasClinical = TRUE)),
        bsr_invalid_tsize_warning = function(c) { w <<- c; invokeRestart("muffleWarning") },
        bsr_experimental_ror_pc_warning = function(c) {
            expect_identical(conditionMessage(c), ror221_notice); invokeRestart("muffleWarning") })
    clin <- grep("Clinic", names(r$score.ROR), value = TRUE)
    expect_length(clin, 6L); expect_true(all(is.na(r$score.ROR[ids[5], clin])))
    expect_identical(as.character(w$patient_ids), ids[5])
    q2 <- fx; q2$NODE[5] <- 9   # valid NODE >= 4 stays high on the physical sample
    r2 <- quiet(BS_ssBC(q2, "ER", Subtype = TRUE, hasClinical = TRUE))
    g <- grep("Group.*Prolif.Subtype", names(r2$score.ROR), value = TRUE)
    expect_identical(r2$score.ROR[ids[5], g], "high")
})

test_that("valid singleton and NCP-02 empty outputs keep their identities", {
    q <- fx; q$ER <- rep("ER+", 24); q$ER[3] <- "ER-"
    expect_warning(r <- suppressMessages(BS_ssBC(q, "ER", Subtype = TRUE, hasClinical = FALSE)), "sample k03")
    expect_true(is.na(r$BS.all["k03", "BS"]))
    e <- fx; e$ER <- NA_character_
    r <- quiet(BS_ssBC(e, "ER", Subtype = TRUE, hasClinical = TRUE))
    expect_identical(rownames(r$BS.all), ids); expect_true(all(is.na(r$BS.all$BS.Subtype)))
    expect_true("Call.Subtype" %in% names(r$score.ROR))
})

test_that("excluded and wholly empty receptor paths also refuse inconsistent identities", {
    q <- fx; q$ER[c(23, 24)] <- NA_character_
    b1 <- q; b1$PatientID[23] <- ""
    b2 <- q; b2$PatientID[23] <- "not_in_assay"
    b3 <- q; b3$PatientID[c(23, 24)] <- ids[c(24, 23)]
    for (b in list(b1, b2, b3)) pid_expect_error(quiet(BS_ssBC(b, "ER", Subtype = FALSE, hasClinical = FALSE)), id_error)
    t2 <- fx; nt <- which(t2$TN == "nonTN")[1:2]; t2$PatientID[nt] <- ids[rev(nt)]
    pid_expect_error(quiet(BS_ssBC(t2, "TN.v2", Subtype = TRUE, hasClinical = TRUE)), id_error)
    e <- fx; e$ER <- NA_character_; e$PatientID <- P_CROSS
    pid_expect_error(quiet(BS_ssBC(e, "ER", Subtype = FALSE, hasClinical = TRUE)), id_error)
})

test_that("a second PatientID column is refused rather than silently ignored", {
    q <- fx; cd <- SummarizedExperiment::colData(q)
    SummarizedExperiment::colData(q) <- cbind(cd, S4Vectors::DataFrame(PatientID = P_CROSS, check.names = FALSE))
    pid_expect_error(quiet(BS_ssBC(q, "ER", Subtype = TRUE, hasClinical = FALSE)), id_error)
    cd2 <- cd; cd2$PatientID <- P_CROSS
    SummarizedExperiment::colData(q) <- cbind(cd2, S4Vectors::DataFrame(PatientID = ids, check.names = FALSE))
    pid_expect_error(quiet(BS_ssBC(q, "ER", Subtype = TRUE, hasClinical = FALSE)), id_error)
})


test_that("all strategies reject incomplete keys and original duplicate key columns", {
    for (i in seq_len(nrow(grid))) {
        inputs <- list()
        for (kind in c("missing", "blank", "duplicate", "disjoint", "missing-column",
                       "double-aligned", "double-permuted")) {
            z <- fx
            if (kind == "missing") z$PatientID[24] <- NA_character_
            if (kind == "blank") z$PatientID[24] <- ""
            if (kind == "duplicate") z$PatientID[24] <- z$PatientID[23]
            if (kind == "disjoint") z$PatientID <- paste0("other_", ids)
            if (kind == "missing-column") z$PatientID <- NULL
            if (startsWith(kind, "double")) {
                cd <- SummarizedExperiment::colData(z)
                if (kind == "double-permuted") cd$PatientID <- P_CROSS
                SummarizedExperiment::colData(z) <- cbind(cd, S4Vectors::DataFrame(PatientID = ids, check.names = FALSE))
            }
            pid_expect_error(BS_ssBC(z, grid$s[i], Subtype = grid$st[i], hasClinical = grid$cl[i]))
        }
    }
})

test_that("ordered keys are required on every empty and excluded path", {
    for (i in seq_len(nrow(grid))) for (kind in c("empty", "excluded", "singleton")) {
        z <- fx
        if (kind == "empty") { z$ER <- NA_character_; z$HER2 <- NA_character_; z$TN <- NA_character_ }
        if (kind == "excluded") { z$ER[c(23, 24)] <- NA_character_; z$HER2[c(23, 24)] <- NA_character_; z$TN[c(23, 24)] <- NA_character_ }
        if (kind == "singleton") { z$ER <- "ER+"; z$ER[23] <- "ER-"; z$HER2 <- "HER2-"; z$TN <- "nonTN"; z$TN[23] <- "TN" }
        z$PatientID <- swap(23, 24)
        pid_expect_error(BS_ssBC(z, grid$s[i], Subtype = grid$st[i], hasClinical = grid$cl[i]))
    }
})

test_that("clinical TSIZE metadata refusal retains precedence over identity", {
    for (kind in c("permuted", "duplicate-column")) {
        z <- fx; z$PatientID <- P_CROSS
        if (kind == "duplicate-column") {
            cd <- SummarizedExperiment::colData(z)
            SummarizedExperiment::colData(z) <- cbind(cd, S4Vectors::DataFrame(PatientID = ids, check.names = FALSE))
        }
        attr(z$TSIZE, "na_values") <- 1
        for (st in c(FALSE, TRUE)) for (s in c("ER", "ER.v2", "TN", "TN.v2")) {
            err <- tryCatch(BS_ssBC(z, s, Subtype = st, hasClinical = TRUE), error = identity)
            expect_s3_class(err, "error")
            expect_false(inherits(err, "bsr_patientid_mismatch_error"))
            expect_match(conditionMessage(err), "TSIZE.*na_values")
            pid_expect_error(BS_ssBC(z, s, Subtype = st, hasClinical = FALSE))
        }
    }
})

test_that("assay keys are validated without manufacturing labels or rewriting metadata", {
    for (kind in c("duplicate", "blank")) {
        z <- fx; key <- ids
        if (kind == "duplicate") key[24] <- key[23]
        if (kind == "blank") key[24] <- ""
        colnames(z) <- key
        pid_expect_error(BS_ssBC(z, "ER"))
    }
    # Invalid container axes can fail in construction or .nc_pheno(), before
    # the approved after-.nc_pheno guard is reached.
    z <- fx; colnames(z) <- NULL
    expect_error(BS_ssBC(z, "ER"), "dimnames")
    z <- fx; key <- ids; key[24] <- NA_character_
    expect_error(colnames(z) <- key, "missing values")
    z <- fx; colnames(z) <- as.character(100001:100024)
    z$PatientID <- as.numeric(100001:100024)
    reference <- z; reference$PatientID <- colnames(z)
    expect_identical(quiet(BS_ssBC(z, "ER")), quiet(BS_ssBC(reference, "ER")))
    z$PatientID <- paste0(" ", z$PatientID)
    pid_expect_error(BS_ssBC(z, "ER"))
    z <- fx; z$PatientID <- toupper(ids)
    pid_expect_error(BS_ssBC(z, "ER"))
    z <- fx; z$OtherPatientID <- rev(ids)
    expect_identical(quiet(BS_ssBC(z, "ER")), quiet(BS_ssBC(fx, "ER")))
    before <- serialize(z, NULL); set.seed(934); rng <- .Random.seed
    quiet(BS_ssBC(z, "ER"))
    expect_identical(serialize(z, NULL), before); expect_identical(.Random.seed, rng)
})


# Conversion fixtures deliberately retain one genuine aligned original key.
converted_fixture <- function(kind) {
    z <- fx; cd <- SummarizedExperiment::colData(z)
    if (startsWith(kind, "C1")) {
        key <- cd[["PatientID"]]; cd$PatientID <- NULL
        cd$meta <- if (kind == "C1-matrix")
            matrix(P_CROSS, dimnames = list(NULL, "PatientID")) else
            S4Vectors::DataFrame(PatientID = if (kind == "C1-factor")
                factor(P_CROSS, levels = rev(ids)) else P_CROSS)
        cd$PatientID <- key
    } else {
        cd$PatientID <- matrix(ids, dimnames = list(NULL, "id"))
        if (kind != "B1-none") {
            prefix <- sub("^B1-", "", kind)
            cd[[prefix]] <- if (prefix == "matching") ids else P_CROSS
            if (prefix == "matching") names(cd)[names(cd) == "matching"] <- "PatientID_old"
        }
    }
    SummarizedExperiment::colData(z) <- cd
    z
}
conversion_refusal <- function(z, s, subtype, clinical, manual = FALSE) {
    before <- serialize(z, NULL); withr::local_seed(322); rng <- .Random.seed
    expression <- function() {
        if (!manual) BS_ssBC(z, s, Subtype = subtype, hasClinical = clinical) else
            BS_Multi(list(se_NC = z, se_SSP = z), methods = c(s, "AIMS"),
                     Subtype = subtype, hasClinical = clinical)
    }
    if (manual && !"PatientID" %in% names(data.frame(SummarizedExperiment::colData(z)))) {
        # BS_Multi refuses a missing converted key before wrapper dispatch.
        warnings <- list()
        err <- tryCatch(withCallingHandlers(expression(), warning = function(w) {
            warnings[[length(warnings) + 1L]] <<- w; invokeRestart("muffleWarning")
        }), error = identity)
        expect_s3_class(err, "error")
        expect_false(inherits(err, "bsr_patientid_mismatch_error"))
        expect_match(conditionMessage(err), "colData(se_NC) must have a 'PatientID' column", fixed = TRUE)
        expect_length(warnings, 0L)
    } else {
        pid_expect_error(withCallingHandlers(expression(), message = function(m) {
            expect_match(conditionMessage(m), "ssBC.*running")
            invokeRestart("muffleMessage")
        }), "[Mm]etadata conversion")
    }
    expect_identical(serialize(z, NULL), before)
    expect_identical(.Random.seed, rng)
}

test_that("B1 requires an exact converted key, never a PatientID prefix", {
    for (kind in c("B1-none", "B1-PatientID_old", "B1-PatientIDx",
                   "B1-PatientID.1", "B1-matching")) {
        z <- converted_fixture(kind)
        cd <- SummarizedExperiment::colData(z); pheno <- data.frame(cd)
        expect_identical(sum(names(cd) == "PatientID"), 1L)
        expect_identical(as.character(cd[["PatientID"]]), ids)
        expect_false("PatientID" %in% names(pheno))
        for (i in seq_len(nrow(grid)))
            conversion_refusal(z, grid$s[i], grid$st[i], grid$cl[i])
        for (s in c("ssBC", "ssBC.v2")) for (st in c(FALSE, TRUE))
            for (cl in c(FALSE, TRUE)) conversion_refusal(z, s, st, cl, manual = TRUE)
    }
})

test_that("C1 refuses a flattened decoy that replaces the exact converted key", {
    for (kind in c("C1-matrix", "C1-nested", "C1-factor")) {
        z <- converted_fixture(kind)
        cd <- SummarizedExperiment::colData(z); pheno <- data.frame(cd)
        expect_identical(sum(names(cd) == "PatientID"), 1L)
        expect_identical(as.character(cd[["PatientID"]]), ids)
        expect_identical(sum(names(pheno) == "PatientID"), 1L)
        expect_identical(as.character(pheno[["PatientID"]]), P_CROSS)
        if (kind != "C1-factor") expect_identical(pheno[["PatientID"]], P_CROSS)
        if (kind == "C1-factor") {
            expect_s3_class(pheno[["PatientID"]], "factor")
            expect_s3_class(cd$meta[["PatientID"]], "factor")
        }
        expect_identical(pheno[["PatientID.1"]], ids)
        for (i in seq_len(nrow(grid)))
            conversion_refusal(z, grid$s[i], grid$st[i], grid$cl[i])
        for (s in c("ssBC", "ssBC.v2")) for (st in c(FALSE, TRUE))
            for (cl in c(FALSE, TRUE)) conversion_refusal(z, s, st, cl, manual = TRUE)
    }
})

test_that("M9 rejects identical malformed assay and PatientID labels", {
    for (kind in c("duplicate", "blank")) {
        z <- fx; keys <- ids
        keys[24] <- if (kind == "duplicate") keys[23] else ""
        colnames(z) <- keys; z$PatientID <- keys
        # Construction succeeds for these cases; both malformed keys agree.
        expect_identical(as.character(SummarizedExperiment::colData(z)[["PatientID"]]), colnames(z))
        expect_identical(data.frame(SummarizedExperiment::colData(z))[["PatientID"]], keys)
        before <- serialize(z, NULL); withr::local_seed(611); rng <- .Random.seed
        for (i in seq_len(nrow(grid)))
            pid_expect_error(BS_ssBC(z, grid$s[i], Subtype = grid$st[i], hasClinical = grid$cl[i]))
        expect_identical(serialize(z, NULL), before)
        expect_identical(.Random.seed, rng)
    }
})

test_that("exact valid keys ignore unrelated PatientID prefixes", {
    capture_valid <- function(z, row) {
        before <- serialize(z, NULL); withr::local_seed(512); rng <- .Random.seed; warnings <- list()
        result <- withCallingHandlers(
            BS_ssBC(z, row$s, Subtype = row$st, hasClinical = row$cl),
            warning = function(w) { warnings[[length(warnings) + 1L]] <<- w; invokeRestart("muffleWarning") })
        expect_identical(vapply(warnings, function(w) class(w)[1], ""),
                         if (row$cl) "bsr_experimental_ror_pc_warning" else character())
        if (row$cl) expect_identical(conditionMessage(warnings[[1]]), ror221_notice)
        expect_identical(serialize(z, NULL), before); expect_identical(.Random.seed, rng)
        list(value = result, warnings = warnings)
    }
    for (i in seq_len(nrow(grid))) {
        ref <- capture_valid(fx, grid[i, ])
        for (name in c("PatientID_old", "PatientIDx")) for (values in list(ids, P_CROSS)) {
            z <- fx; SummarizedExperiment::colData(z)[[name]] <- values
            expect_identical(capture_valid(z, grid[i, ]), ref)
        }
        z <- fx
        SummarizedExperiment::colData(z)$PatientID <- matrix(ids, dimnames = list(NULL, "PatientID"))
        expect_identical(capture_valid(z, grid[i, ]), ref)
    }
})

test_that("TSIZE metadata preparation still precedes converted-key refusal", {
    for (kind in c("B1-none", "B1-PatientID_old", "C1-matrix", "C1-nested", "C1-factor")) {
        z <- converted_fixture(kind)
        if (kind == "C1-factor") {
            cd <- SummarizedExperiment::colData(z); pheno <- data.frame(cd)
            expect_identical(sum(names(cd) == "PatientID"), 1L)
            expect_identical(as.character(cd[["PatientID"]]), ids)
            expect_identical(tail(names(pheno), 2), c("PatientID", "PatientID.1"))
            expect_s3_class(pheno[["PatientID"]], "factor")
            expect_identical(as.character(pheno[["PatientID"]]), P_CROSS)
            expect_identical(pheno[["PatientID.1"]], ids)
        }
        attr(z$TSIZE, "na_values") <- 1
        before <- serialize(z, NULL); withr::local_seed(911); rng <- .Random.seed
        for (s in c("ER", "ER.v2", "TN", "TN.v2")) for (st in c(FALSE, TRUE)) {
            seen <- list()
            err <- tryCatch(withCallingHandlers(BS_ssBC(z, s, Subtype = st, hasClinical = TRUE),
                warning = function(w) { seen[[length(seen) + 1L]] <<- w; invokeRestart("muffleWarning") }), error = identity)
            expect_s3_class(err, "error")
            expect_false(inherits(err, "bsr_patientid_mismatch_error"))
            expect_match(conditionMessage(err), "TSIZE.*na_values")
            expect_length(seen, 0L)
            conversion_refusal(z, s, st, FALSE)
        }
        expect_identical(serialize(z, NULL), before); expect_identical(.Random.seed, rng)
    }
})

test_that("manual preconditions and AUTO wrapper boundaries remain distinct", {
    z <- fx; z$PatientID <- NULL
    err <- tryCatch(BS_Multi(list(se_NC = z, se_SSP = z), methods = c("ssBC", "ssBC.v2")), error = identity)
    expect_s3_class(err, "error")
    expect_false(inherits(err, "bsr_patientid_mismatch_error"))
    expect_match(conditionMessage(err), "colData(se_NC) must have a 'PatientID' column", fixed = TRUE)
    e <- new.env(); data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = e)
    auto <- e$OSLO2EMIT0obj$data_input
    auto_ids <- sprintf("auto%03d", seq_len(ncol(auto$se_NC)))
    n <- length(auto_ids)
    for (k in names(auto)) {
        colnames(auto[[k]]) <- auto_ids
        cd <- S4Vectors::DataFrame(PatientID = auto_ids,
            ER = rep("ER-", n),
            HER2 = rep("HER2-", n),
            TSIZE = rep(0:1, length.out = n), NODE = rep(c(0, 2, 5), length.out = n),
            row.names = auto_ids)
        cd$TN <- ifelse(cd$ER == "ER-" & cd$HER2 == "HER2-", "TN", "nonTN")
        SummarizedExperiment::colData(auto[[k]]) <- cd
    }
    expect_true(all(auto$se_NC$TN == "TN") && n >= 18L)
    auto_call <- function(d) {
        messages <- character()
        result <- tryCatch(withCallingHandlers(BS_Multi(d, methods = "AUTO"),
            message = function(m) {
                messages <<- c(messages, conditionMessage(m)); invokeRestart("muffleMessage")
            }), error = identity)
        expect_true(all(grepl("Running AUTO mode|\\[AUTO\\]|ssBC.*running", messages)))
        if (inherits(result, "error")) stop(result)
        result
    }
    for (kind in c("duplicate-original", "C1-matrix")) {
        d <- auto; cd <- SummarizedExperiment::colData(d$se_NC)
        if (kind == "duplicate-original") {
            cd <- cbind(cd, S4Vectors::DataFrame(PatientID = auto_ids, check.names = FALSE))
        } else {
            key <- cd$PatientID; cd$PatientID <- NULL
            cd$meta <- matrix(rev(auto_ids), dimnames = list(NULL, "PatientID"))
            cd$PatientID <- key
        }
        SummarizedExperiment::colData(d$se_NC) <- cd
        # Eligible strata select ssBC when the existing AUTO precheck passes.
        expect_identical(as.character(cd$PatientID), auto_ids)
        before <- serialize(d, NULL); withr::local_seed(337); rng <- .Random.seed
        if (kind == "duplicate-original") {
            expect_identical(data.frame(cd)[["PatientID"]], auto_ids)
            pid_expect_error(auto_call(d))
        } else {
            expect_false(identical(data.frame(cd)[["PatientID"]], auto_ids))
            err <- tryCatch(auto_call(d), error = identity)
            expect_s3_class(err, "error")
            expect_false(inherits(err, "bsr_patientid_mismatch_error"))
            expect_identical(conditionMessage(err),
                "AUTO requires PatientID labels to match the input assay sample order.")
        }
        expect_identical(serialize(d, NULL), before); expect_identical(.Random.seed, rng)
    }
})


test_that("positive converted keys preserve physical-sample clinical results", {
    capture <- function(z, s, st) {
        before <- serialize(z, NULL); withr::local_seed(761); rng <- .Random.seed
        warnings <- messages <- list()
        value <- withCallingHandlers(BS_ssBC(z, s, Subtype = st, hasClinical = TRUE),
            warning = function(w) {
                warnings[[length(warnings) + 1L]] <<- list(class = class(w),
                    message = conditionMessage(w), patient_ids = w$patient_ids)
                invokeRestart("muffleWarning")
            }, message = function(m) {
                messages[[length(messages) + 1L]] <<- list(class = class(m), message = conditionMessage(m))
                invokeRestart("muffleMessage")
            })
        expect_identical(serialize(z, NULL), before); expect_identical(.Random.seed, rng)
        expect_setequal(vapply(warnings, function(w) w$class[1], ""),
            c("bsr_invalid_tsize_warning", "bsr_experimental_ror_pc_warning"))
        invalid <- Filter(function(w) "bsr_invalid_tsize_warning" %in% w$class, warnings)
        expect_length(invalid, 1L); expect_identical(invalid[[1]]$patient_ids, ids[3])
        clinical <- grep("Clinic", names(value$score.ROR), value = TRUE)
        expect_length(clinical, 6L)
        expect_true(all(is.na(value$score.ROR[ids[3], clinical])))
        expect_false(anyNA(value$score.ROR[ids[9], clinical]))
        list(value = value, warnings = warnings, messages = messages)
    }
    z <- fx; z$TSIZE[3] <- 2; z$NODE[9] <- 7
    # k03 and k09 both belong to TN, so all four strategies evaluate them.
    expect_identical(z$TN[c(3, 9)], rep("TN", 2))
    for (s in c("ER", "ER.v2", "TN", "TN.v2")) for (st in c(FALSE, TRUE)) {
        reference <- capture(z, s, st)
        for (kind in c("unused-dot-one", "matching-flattened")) {
            q <- z; cd <- SummarizedExperiment::colData(q)
            if (kind == "unused-dot-one") cd$PatientID.1 <- P_CROSS else {
                key <- cd$PatientID; cd$PatientID <- NULL
                cd$meta <- S4Vectors::DataFrame(PatientID = key)
                cd$PatientID <- key
            }
            SummarizedExperiment::colData(q) <- cd
            pheno <- data.frame(cd)
            expect_identical(sum(names(cd) == "PatientID"), 1L)
            expect_identical(as.character(cd[["PatientID"]]), ids)
            expect_identical(sum(names(pheno) == "PatientID"), 1L)
            expect_identical(pheno[["PatientID"]], ids)
            expect_identical(pheno[["PatientID.1"]], if (kind == "unused-dot-one") P_CROSS else ids)
            if (kind == "matching-flattened") expect_s4_class(cd$meta, "DataFrame")
            actual <- capture(q, s, st)
            expect_identical(actual, reference)
        }
    }
})

})
