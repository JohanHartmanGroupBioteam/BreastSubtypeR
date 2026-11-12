# NCP-02: four-class ER, ER.v2 and TN ssBC with no sample in any selected
# receptor stratum. Public OSLO2EMIT0 expression with synthetic annotations.
# Expectations are full returns, warning sets and the preconditions that keep
# other inputs on their existing paths.

ncp02_strategies <- list(
    ER = function(z) { z$ER <- rep(NA_character_, ncol(z)); z },
    ER.v2 = function(z) { z$HER2 <- rep(NA_character_, ncol(z)); z },
    ER.v2 = function(z) { z$HER2 <- rep("2+", ncol(z)); z },
    TN = function(z) { z$TN <- rep(NA_character_, ncol(z)); z })

test_that("empty four-class strata return the five-class empty object plus four-class fields", {
    withr::local_seed(202611)
    for (ix in list(seq_len(36), c(seq(2, 36, 2), seq(1, 35, 2)))) for (factor_ids in c(FALSE, TRUE))
        for (k in seq_along(ncp02_strategies)) for (clinical in c(FALSE, TRUE)) {
            s <- names(ncp02_strategies)[k]
            z <- ncp02_strategies[[k]](ncp02_input(ix, factor_ids)$se_NC)
            ids <- colnames(z)
            before <- z
            seed <- .Random.seed
            five <- ncp02_capture(BS_ssBC(z, s, Subtype = FALSE, hasClinical = clinical))
            four <- ncp02_capture(BS_ssBC(z, s, Subtype = TRUE, hasClinical = clinical))
            expect_identical(.Random.seed, seed)
            expect_identical(z, before)
            expected_classes <- if (clinical) "bsr_experimental_ror_pc_warning" else character()
            expect_identical(four$classes, expected_classes)
            expect_identical(five$classes, expected_classes)
            a <- four$value
            expect_identical(a, ncp02_explicit_empty(ids, clinical))
            expect_identical(a, ncp02_four_from_five(five$value))
            expect_identical(five$value, ncp02_explicit_empty(ids, clinical, subtype = FALSE))
            expect_identical(ncol(a$score.ROR), if (clinical) 20L else 14L)
            expect_identical(match("Call.Subtype", names(a$score.ROR)), match("Call", names(a$score.ROR)) + 1L)
            expect_false("patientID" %in% names(a$score.ROR))
            expect_true(all(vapply(a$score.ROR, is.logical, TRUE)))
            if (clinical)
                expect_true(all(is.na(a$score.ROR$`ROR-PC Group (Subtype + Clinic + Prolif.Subtype)`)))
        }
})

test_that("unrecognized structurally valid receptor labels use literal membership", {
    z0 <- ncp02_input()$se_NC
    cases <- list(
        ER = function(z) { z$ER <- rep(c("er+", "er-"), 18); z },
        ER = function(z) { z$ER <- rep(c(1, 0), 18); z },
        ER = function(z) { z$ER <- factor(rep(c("positive", "negative"), 18)); z },
        ER.v2 = function(z) { z$HER2 <- rep(c("negative", "unknown"), 18); z },
        TN = function(z) { z$TN <- rep(c(TRUE, FALSE), 18); z })
    for (k in seq_along(cases)) for (clinical in c(FALSE, TRUE)) {
        z <- cases[[k]](z0)
        s <- names(cases)[k]
        a <- ncp02_capture(BS_ssBC(z, s, Subtype = TRUE, hasClinical = clinical))
        expect_identical(a$value, ncp02_explicit_empty(colnames(z), clinical))
        expect_identical(a$classes, if (clinical) "bsr_experimental_ror_pc_warning" else character())
    }
})

test_that("one sample and one represented PAM50 gene qualify only when no sample is a member", {
    x <- ncp02_input()$se_NC
    one <- x[, 1]
    one$ER <- NA_character_
    expect_identical(ncp02_no_warning(BS_ssBC(one, "ER", Subtype = TRUE)),
        ncp02_explicit_empty(colnames(one), FALSE))
    gene <- intersect(rownames(x), rownames(ncp02_model()$centroid))[1]
    single <- x[gene, ]
    empty <- single
    empty$ER <- rep(NA_character_, ncol(empty))
    expect_identical(ncp02_no_warning(BS_ssBC(empty, "ER", Subtype = TRUE)),
        ncp02_explicit_empty(colnames(empty), FALSE))
    # A populated one-gene input keeps its existing behavior in both settings.
    for (subtype in c(FALSE, TRUE))
        expect_error(BS_ssBC(single, "ER", Subtype = subtype), "'names' attribute", fixed = TRUE)
    # A recognized singleton keeps its own no-call path and warning.
    single_member <- x
    single_member$ER <- c("ER+", rep(NA_character_, ncol(x) - 1L))
    w <- ncp02_capture(BS_ssBC(single_member, "ER", Subtype = TRUE))
    expect_identical(w$messages, sprintf(
        "ssBC/ssBC.v2: subgroup 'ER_pos' has n=1 (sample %s); no subtype or risk estimate.",
        colnames(x)[1]))
    expect_identical(as.character(w$value$score.ROR$patientID[1]), colnames(x)[1])
    expect_true(all(is.na(w$value$BS.all$BS.Subtype)))
    # One recognized pair is evaluated by the ordinary path.
    pair <- x
    pair$ER <- c("ER+", "ER+", rep(NA_character_, ncol(x) - 2L))
    p <- ncp02_no_warning(BS_ssBC(pair, "ER", Subtype = TRUE))
    expect_identical(sum(!is.na(p$BS.all$BS.Subtype)), 2L)
    expect_true("patientID" %in% names(p$score.ROR))
})

test_that("inputs outside the empty four-class preconditions do not become successful", {
    x <- ncp02_input()$se_NC
    empty <- x
    empty$ER <- rep(NA_character_, ncol(x))
    four <- function(z) BS_ssBC(z, "ER", Subtype = TRUE)
    # Zero samples and zero PAM50 overlap, with empty or populated strata.
    expect_error(four(empty[, 0]), "undefined columns selected", fixed = TRUE)
    absent <- empty
    rownames(absent) <- paste0("absent_", seq_len(nrow(absent)))
    expect_error(four(absent), "undefined columns selected", fixed = TRUE)
    populated_absent <- x
    rownames(populated_absent) <- rownames(absent)
    expect_error(four(populated_absent), "undefined columns selected", fixed = TRUE)
    # Sample identities: PatientID permuted alone, disjoint, missing, empty or duplicated.
    ids <- colnames(x)
    for (labels in list(rev(ids), paste0("X", ids), c("", ids[-1]))) {
        z <- empty
        z$PatientID <- labels
        expect_error(four(z), class = "bsr_patientid_mismatch_error")
    }
    z <- empty
    z$PatientID <- c(NA, ids[-1])
    expect_error(four(z), class = "bsr_patientid_mismatch_error")
    z$PatientID <- c(ids[2], ids[-1])
    expect_error(four(z), class = "bsr_patientid_mismatch_error")
    # A consistently permuted assay with its own metadata remains valid.
    perm <- empty[, rev(seq_len(ncol(empty)))]
    expect_identical(ncp02_no_warning(four(perm)), ncp02_explicit_empty(rev(ids), FALSE))
    # A receptor column that is not a plain vector is outside the branch.
    listed <- empty
    listed$ER <- I(as.list(rep(NA_character_, ncol(x))))
    expect_error(four(listed), "undefined columns selected", fixed = TRUE)
    # Structure is checked on the supplied colData, before data.frame() coercion.
    for (column in list(matrix(NA_character_, ncol(x), 1L), S4Vectors::Rle(rep(NA_character_, ncol(x))))) {
        shaped <- empty
        SummarizedExperiment::colData(shaped)$ER <- column
        expect_error(four(shaped), "undefined columns selected", fixed = TRUE)
    }
    # PatientID labels must be character or factor values identical to the assay names.
    numbered <- empty
    colnames(numbered) <- as.character(seq_len(ncol(x)))
    numbered$PatientID <- seq_len(ncol(x))
    expect_error(four(numbered), "undefined columns selected", fixed = TRUE)
    numbered$PatientID <- as.character(seq_len(ncol(x)))
    expect_identical(ncp02_no_warning(four(numbered)),
        ncp02_explicit_empty(as.character(seq_len(ncol(x))), FALSE))
    # PatientID structure is checked on the supplied colData too.
    for (column in list(S4Vectors::Rle(ids), matrix(ids, ncol(x), 1L))) {
        shaped <- empty
        SummarizedExperiment::colData(shaped)$PatientID <- column
        expect_error(four(shaped), "undefined columns selected", fixed = TRUE)
    }
    second <- empty
    SummarizedExperiment::colData(second) <- cbind(SummarizedExperiment::colData(second),
        S4Vectors::DataFrame(PatientID = paste0("X", ids), check.names = FALSE))
    expect_error(four(second), class = "bsr_patientid_mismatch_error")
    # TN.v2 and five-class returns are not routed through the new branch.
    tn <- x
    tn$TN <- rep(NA_character_, ncol(x))
    a <- BS_ssBC(tn, "TN.v2", Subtype = TRUE)
    expect_true("patientID" %in% names(a$score.ROR))
    expect_identical(dim(a$outList$testData), c(50L, ncol(x)))
})

test_that("TSIZE handling, warnings and NODE grouping are unchanged for empty strata", {
    d <- ncp02_input()
    z <- d$se_NC
    z$HER2 <- rep(NA_character_, ncol(z))
    ids <- colnames(z)
    # Ambiguous metadata stops the clinical call before receptor handling.
    for (subtype in c(FALSE, TRUE)) {
        bad <- z
        size <- bad$TSIZE
        attr(size, "na_values") <- 1
        bad$TSIZE <- size
        expect_error(BS_ssBC(bad, "ER.v2", Subtype = subtype, hasClinical = TRUE), "na_values")
    }
    # Invalid codes on unevaluated rows raise no TSIZE warning from ssBC.
    invalid <- z
    invalid$TSIZE[1:3] <- c(2, 5, 2)
    a <- ncp02_capture(BS_ssBC(invalid, "ER.v2", Subtype = TRUE, hasClinical = TRUE))
    expect_identical(a$classes, "bsr_experimental_ror_pc_warning")
    expect_identical(a$value, ncp02_explicit_empty(ids, TRUE))
    # A peer that evaluates those rows reports them once, with exact IDs.
    d$se_NC <- invalid
    m <- ncp02_capture(BS_Multi(d, methods = c("ssBC.v2", "parker.original"),
        Subtype = TRUE, hasClinical = TRUE))
    expect_identical(sort(m$classes), c("bsr_experimental_ror_pc_warning", "bsr_invalid_tsize_warning"))
    tsize <- m$warnings[[match("bsr_invalid_tsize_warning", m$classes)]]
    expect_identical(sort(tsize$patient_ids), sort(ids[1:3]))
    expect_identical(m$value$results$ssBC.v2, ncp02_explicit_empty(ids, TRUE))
    grp <- "ROR-PC Group (Subtype + Clinic + Prolif.Subtype)"
    # The evaluating peer is unchanged in every NODE category (0, 1-3, >= 4).
    peer <- ncp02_capture(BS_Multi(d, methods = c("parker.original", "AIMS"),
        Subtype = TRUE, hasClinical = TRUE))
    expect_identical(sort(peer$classes), c("bsr_experimental_ror_pc_warning", "bsr_invalid_tsize_warning"))
    expect_identical(m$value$results$parker.original, peer$value$results$parker.original)
    expect_true(all(c(0, 2, 4, 7) %in% invalid$NODE))
    parker <- m$value$results$parker.original$score.ROR
    evaluated_high <- ids[-(1:3)][invalid$NODE[-(1:3)] >= 4]
    expect_true(all(parker[evaluated_high, grp] == "high"))
    expect_true(all(is.na(parker[ids[1:3], grp])))
    # Warning-as-error behavior for the experimental notice is preserved.
    withr::local_options(warn = 2)
    expect_error(BS_ssBC(z, "ER.v2", Subtype = TRUE, hasClinical = TRUE), "Experimental 18-gene ROR-PC")
    expect_identical(BS_ssBC(z, "ER.v2", Subtype = TRUE, hasClinical = FALSE), ncp02_explicit_empty(ids, FALSE))
})

test_that("manual BS_Multi keeps peers unchanged around an empty four-class method", {
    d <- ncp02_input()
    for (n in names(d)) d[[n]]$ER <- rep(NA_character_, 36)
    ids <- colnames(d$se_NC)
    before <- d
    # With no available vote for any patient, both summaries have NA entropy.
    none <- ncp02_no_warning(BS_Multi(d, methods = c("ssBC", "ssBC.v2"), Subtype = TRUE))
    for (field in c("res_subtypes", "res_subtypes.Subtype")) {
        expect_identical(none[[field]]$ssBC, rep(NA, 36))
        expect_identical(none[[field]]$entropy, rep(NA_real_, 36))
    }
    expect_identical(d, before)
    # A duplicated receptor column is outside the branch here as well.
    dup <- d
    SummarizedExperiment::colData(dup$se_NC) <- cbind(SummarizedExperiment::colData(dup$se_NC),
        S4Vectors::DataFrame(ER = rep(c("ER+", "ER-"), 18), check.names = FALSE))
    expect_error(suppressMessages(BS_Multi(dup, methods = c("ssBC", "AIMS"), Subtype = TRUE)),
        "undefined columns selected", fixed = TRUE)
    control <- ncp02_capture(BS_Multi(d, methods = c("AIMS", "sspbc"), Subtype = TRUE))
    expect_identical(control$classes, character())
    for (method in c("ssBC", "ssBC.v2")) {
        a <- ncp02_capture(BS_Multi(d, methods = c(method, "AIMS", "sspbc"), Subtype = TRUE))
        expect_identical(a$classes, control$classes)
        expect_identical(a$value$results[[method]], ncp02_explicit_empty(ids, FALSE))
        expect_identical(a$value$results[c("AIMS", "sspbc")], control$value$results)
        for (field in c("res_subtypes", "res_subtypes.Subtype")) {
            summary <- a$value[[field]]
            expect_identical(summary[[method]], rep(NA, 36))
            peers <- summary
            peers[[method]] <- NULL
            expect_identical(peers, control$value[[field]])
            calls <- summary[setdiff(names(summary), "entropy")]
            entropy <- apply(calls, 1, function(row) {
                row <- row[!is.na(row)]
                if (!length(row)) return(NA_real_)
                p <- prop.table(table(row))
                -sum(p * log2(p))
            })
            expect_equal(summary$entropy, unname(entropy), tolerance = 1e-14)
        }
    }
})

test_that("the AUTO zero-pair exception admits only the exact four-class schema", {
    ids <- c("a", "b", "c")
    for (clinical in c(FALSE, TRUE)) {
        valid <- ncp02_explicit_empty(ids, clinical)
        expect_identical(BreastSubtypeR:::.restore_auto_ssbc(valid, ids, clinical, zero_pair_subtype = TRUE), valid)
        expect_error(BreastSubtypeR:::.restore_auto_ssbc(valid, ids, clinical), "result frames")
        expect_error(BreastSubtypeR:::.restore_auto_ssbc(ncp02_explicit_empty(ids, clinical, FALSE),
            ids, clinical, zero_pair_subtype = TRUE), "result frames")
        expect_error(BreastSubtypeR:::.restore_auto_ssbc(valid, ids), "invalid score.ROR schema")
        broken <- valid; broken$outList$predictions.FourSubtype <- NULL
        expect_error(BreastSubtypeR:::.restore_auto_ssbc(broken, ids, clinical, zero_pair_subtype = TRUE), "prediction or model")
        broken <- valid; broken$score.ROR <- broken$score.ROR[c(1:5, 7, 6, 8:ncol(broken$score.ROR))]
        expect_error(BreastSubtypeR:::.restore_auto_ssbc(broken, ids, clinical, zero_pair_subtype = TRUE), "result frames")
        broken <- valid; broken$score.ROR$patientID <- NA
        expect_error(BreastSubtypeR:::.restore_auto_ssbc(broken, ids, clinical, zero_pair_subtype = TRUE), "result frames")
        broken <- valid; broken$outList$testData <- matrix(NA_real_, 50L, 3L, dimnames = list(rownames(valid$outList$centroids), ids))
        expect_error(BreastSubtypeR:::.restore_auto_ssbc(broken, ids, clinical, zero_pair_subtype = TRUE), "prediction or model")
    }
})

test_that("alias-coded AUTO routes outside the zero-pair exception still fail", {
    base <- ncp02_input(seq_len(36))
    relabel <- function(f) { x <- base; for (n in names(x)) x[[n]] <- f(x[[n]]); x }
    partial <- relabel(function(z) { h <- rep(NA_character_, 36); h[c(1, 3)] <- "negative"; z$HER2 <- h; z })
    tn_yes <- relabel(function(z) { z$TN <- rep("yes", 36); z })
    worded <- relabel(function(z) {
        z$ER <- ifelse(z$ER == "ER+", "positive", "negative")
        z$HER2 <- ifelse(z$HER2 == "HER2+", "positive", "negative")
        z$TN <- rep(NA_character_, 36)
        z
    })
    # Routing coerces these labels with its own warning; nothing else may warn.
    coerced_only <- function(expr, pattern) {
        seen <- character()
        expect_error(withCallingHandlers(suppressMessages(expr), warning = function(w) {
            seen <<- c(seen, conditionMessage(w))
            invokeRestart("muffleWarning")
        }), pattern)
        expect_true(length(seen) > 0L && all(startsWith(seen, "Phenodata: coerced")))
    }
    for (subtype in c(FALSE, TRUE)) {
        coerced_only(BS_Multi(worded, methods = "AUTO", Subtype = subtype), "ER balancing")
        coerced_only(BS_Multi(partial, methods = "AUTO", Subtype = subtype), "invalid score.ROR schema")
        coerced_only(BS_Multi(tn_yes, methods = "AUTO", Subtype = subtype), "invalid score.ROR schema")
    }
})

test_that("the Subtype argument keeps its established if() coercion", {
    x <- ncp02_input()$se_NC
    for (s in c("ER", "TN.v2")) {
        logical_four <- BS_ssBC(x, s, Subtype = TRUE)
        logical_five <- BS_ssBC(x, s, Subtype = FALSE)
        expect_identical(BS_ssBC(x, s, Subtype = "TRUE"), logical_four)
        expect_identical(BS_ssBC(x, s, Subtype = "T"), logical_four)
        expect_identical(BS_ssBC(x, s, Subtype = "FALSE"), logical_five)
    }
    empty <- x
    empty$ER <- rep(NA_character_, ncol(x))
    expect_identical(ncp02_no_warning(BS_ssBC(empty, "ER", Subtype = "TRUE")),
        ncp02_explicit_empty(colnames(x), FALSE))
})

test_that("a duplicated receptor column is outside the new branch", {
    x <- ncp02_input()$se_NC
    x$ER <- rep(NA_character_, ncol(x))
    SummarizedExperiment::colData(x) <- cbind(SummarizedExperiment::colData(x),
        S4Vectors::DataFrame(ER = rep(c("ER+", "ER-"), 18), check.names = FALSE))
    expect_identical(sum(names(SummarizedExperiment::colData(x)) == "ER"), 2L)
    expect_error(BS_ssBC(x, "ER", Subtype = TRUE), "undefined columns selected", fixed = TRUE)
})

test_that("AUTO without ER annotations still routes only to the single-sample predictors", {
    d <- ncp02_input()
    for (n in names(d)) d[[n]]$ER <- rep(NA_character_, 36)
    for (subtype in c(FALSE, TRUE)) {
        a <- ncp02_capture(BS_Multi(d, methods = "AUTO", Subtype = subtype))
        expect_identical(names(a$value$results), c("AIMS", "sspbc"))
        expect_identical(a$classes, character())
    }
})
