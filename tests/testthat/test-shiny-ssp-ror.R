# Independent projection of the public API for TSV comparisons. Compare the
# serialized decimal values exactly, including tiny posteriors and NA fields.
ssp_reference_export <- function(subtype, ror, four) {
    ids <- rownames(subtype$cl)
    stopifnot(identical(rownames(ror$cl), ids))
    k <- if (four) "21" else "24"
    labels <- c("Basal", "Her2", "LumA", "LumB", if (!four) "Normal")
    bins <- sprintf("c%03d", seq(5, 90, 5))
    tab <- data.frame(ids, subtype$cl[, k], subtype$prob[, k],
        subtype$all.probs[[k]][, labels, drop = FALSE],
        colSums(!is.na(subtype$rules.matrix)), nrow(subtype$rules.matrix),
        ror$cl[, "21"], ror$prob[, "21"], ror$all.probs[["21"]][, bins, drop = FALSE],
        colSums(!is.na(ror$rules.matrix)), nrow(ror$rules.matrix),
        row.names = NULL, check.names = FALSE)
    names(tab) <- c("PatientID", if (four) "Call_4class" else "Call_5class",
        "SSP_Call_Prob", paste0("SSP_Prob_", labels), "SSP_Rules_Evaluable",
        "SSP_Rules_Total", "SSP_ROR_Bin", "SSP_ROR_Bin_Prob", paste0("SSP_ROR_Prob_", bins),
        "SSP_ROR_Rules_Evaluable", "SSP_ROR_Rules_Total")
    tab
}

ssp_read_export <- function(path) {
    utils::read.delim(path, check.names = FALSE, colClasses = "character",
                      comment.char = "#", na.strings = character())
}

ssp_reference_tsv <- function(tab) {
    path <- tempfile()
    on.exit(unlink(path))
    utils::write.table(tab, path, sep = "\t", row.names = FALSE, qmethod = "double")
    ssp_read_export(path)
}

ssp_small_upload <- function(se) {
    table_file <- function(df) {
        path <- tempfile(fileext = ".txt")
        utils::write.table(df, path, sep = "\t", row.names = FALSE, qmethod = "double")
        data.frame(name = basename(path), datapath = path, size = file.size(path), type = "text/plain")
    }
    list(GEX = table_file(data.frame(probe = rownames(se), SummarizedExperiment::assay(se),
                                     check.names = FALSE)),
         clinic = table_file(as.data.frame(SummarizedExperiment::colData(se))),
         anno = table_file(as.data.frame(SummarizedExperiment::rowData(se))))
}

test_that("both SSP exports match the public API on the normalized example and raw uploads", {
    env <- shiny_app_env()
    map <- BreastSubtypeR::Mapping
    predict <- BreastSubtypeR::BS_sspbc
    mapped_ror <- logical()
    models <- character()
    local_mocked_bindings(Mapping = function(se_obj, ..., include_ssp_ror = FALSE) {
        mapped_ror <<- c(mapped_ror, include_ssp_ror)
        map(se_obj, ..., include_ssp_ror = include_ssp_ror)
    }, BS_sspbc = function(se_obj, ssp.name) {
        models <<- c(models, ssp.name)
        predict(se_obj, ssp.name)
    }, .package = "BreastSubtypeR")

    for (raw in c(FALSE, TRUE)) {
        se <- current_example(if (raw) "TCGABRCAobj" else "OSLO2EMIT0obj")$se_obj
        expect_identical(dim(se), if (raw) c(5350L, 100L) else c(569L, 102L))
        expr <- shiny_upload(if (raw) "TCGABRCAsub_GEX_RawCounts.txt" else "OSLO2EMIT0_GEX_log2.FPKM.txt")
        # Use the uploaded decimal values on this R platform.
        x <- as.matrix(utils::read.delim(expr$datapath, row.names = 1L, check.names = FALSE))
        expect_identical(dimnames(x), dimnames(se))
        SummarizedExperiment::assay(se) <- x
        reference <- map(se, RawCounts = raw, include_ssp_ror = TRUE, impute = FALSE, verbose = FALSE)
        ror <- predict(reference$se_SSP_ROR, "ssp.ror")
        expect_identical(nrow(reference$se_SSP_ROR), 296L)
        expect_true(all(colSums(!is.na(ror$rules.matrix)) == 361L))
        models <- character()
        mapped_ror <- logical()
        shiny::testServer(env$server, {
            session$setInputs(BSmethod = "sspbc", k_subtypes = "5",
                              is_raw_counts = if (raw) "raw" else "norm")
            if (raw) {
                session$setInputs(GEX = expr, clinic = shiny_upload("TCGABRCAsub_clinical.txt"),
                                  anno = shiny_upload("TCGABRCAsub_anno.txt"))
            } else session$setInputs(load_example_step1 = 1)
            session$setInputs(map = 1)
            expect_identical(mapped_ror, FALSE)
            for (i in 1:2) {
                four <- i == 2L
                model <- if (four) "ssp.subtype" else "ssp.pam50"
                expected <- ssp_reference_export(predict(reference$se_SSP, model), ror, four)
                expected_tsv <- ssp_reference_tsv(expected)
                session$setInputs(k_subtypes = if (four) "4" else "5", run = i)
                p <- current_prediction()
                html <- output$export_selector$html
                expect_match(html, 'value="calls" checked', fixed = TRUE)
                expect_match(html, "Full metrics (incl. ROR)", fixed = TRUE)
                expect_equal(sum(mapped_ror), i - 1L)
                expect_equal(sum(models == "ssp.ror"), i - 1L)
                session$setInputs(export_kind = "calls")
                expect_identical(ssp_read_export(output$download), expected_tsv[, 1:2])
                session$setInputs(export_kind = p$ssp$full_choice)
                path <- output$download
                expect_identical(ssp_read_export(path), expected_tsv)
                expect_identical(ncol(expected_tsv), if (four) 31L else 32L)
                header <- readLines(path, n = 3)
                expect_match(header[1], "10.1093/nargab/lqaf131", fixed = TRUE)
                expect_match(header[2], model, fixed = TRUE)
                expect_match(header[2], "categorical ssp.ror", fixed = TRUE)
                expect_false(grepl("ROR: OFF", header[2], fixed = TRUE))
                expect_match(header[3], "Export: full", fixed = TRUE)
                expect_null(output$ssp_ror_notice)
                session$setInputs(export_kind = "calls")
                expect_identical(ssp_read_export(output$download), expected_tsv[, 1:2])
                session$setInputs(export_kind = p$ssp$full_choice)
                expect_identical(ssp_read_export(output$download), expected_tsv)
                expect_identical(current_prediction(), p)
                expect_equal(sum(mapped_ror), i)
                expect_equal(sum(models == "ssp.ror"), i)
                expect_equal(sum(models == model), 1L)
            }
        })
    }
})

test_that("partial and zero evidence preserve native cells, sample IDs, and singleton dimensions", {
    env <- shiny_app_env()
    se <- current_example("OSLO2EMIT0obj")$se_obj[, c(7L, 2L, 5L, 1L)]
    ids <- c('z#01/"A"', "007", "a.b / 1", "patient-B")
    colnames(se) <- ids
    SummarizedExperiment::colData(se)$PatientID <- ids
    ror_ids <- rownames(Mapping(se, include_ssp_ror = TRUE, verbose = FALSE)$se_SSP_ROR)
    genes <- SummarizedExperiment::rowData(se)$ENTREZID
    x <- round(SummarizedExperiment::assay(se), 6)
    x[genes %in% ror_ids, 1] <- NA_real_
    x[genes %in% head(ror_ids, 10), 2] <- NA_real_
    x[, 3] <- NA_real_
    SummarizedExperiment::assay(se) <- x
    for (cols in list(seq_along(ids), 2L)) {
        input_se <- se[, cols, drop = FALSE]
        mapped <- Mapping(input_se, include_ssp_ror = TRUE, verbose = FALSE)
        ror <- BS_sspbc(mapped$se_SSP_ROR, "ssp.ror")
        subtype <- BS_sspbc(mapped$se_SSP, "ssp.subtype")
        expected <- ssp_reference_export(subtype, ror, TRUE)
        # Serialize once through the upload table as well as through the TSV.
        uploads <- ssp_small_upload(input_se)
        shiny::testServer(env$server, {
            session$setInputs(BSmethod = "sspbc", k_subtypes = "4", is_raw_counts = "norm")
            do.call(session$setInputs, uploads)
            session$setInputs(map = 1)
            session$setInputs(run = 1)
            expect_identical(ssp_read_export(output$download), ssp_reference_tsv(expected[, 1:2]))
            session$setInputs(export_kind = current_prediction()$ssp$full_choice)
            path <- output$download
            tab <- ssp_read_export(path)
            expect_identical(tab, ssp_reference_tsv(expected))
            expect_identical(tab$PatientID, ids[cols])
            note <- readLines(path, n = 4)[4]
            expect_match(note, "^# Incomplete SSP.ROR rule coverage")
            expect_match(output$ssp_ror_notice$html, "Incomplete SSP.ROR rule coverage", fixed = TRUE)
            if (length(cols) > 1L) {
                expect_identical(as.integer(tab$SSP_ROR_Rules_Evaluable[c(1, 3, 4)]), c(0L, 0L, 361L))
                expect_true(as.integer(tab$SSP_ROR_Rules_Evaluable[2]) %in% 1:360)
                expect_true(all(!is.na(ror$cl[c(1, 3), 1])))
                expect_true(is.na(subtype$cl[3, 1]))
                expect_match(note, "SSP_ROR_Rules_Evaluable = 0 (2)", fixed = TRUE)
                expect_match(output$ssp_ror_notice$html, "no expression support", fixed = TRUE)
                expect_match(note, "Check probe-to-Entrez annotation", fixed = TRUE)
            } else expect_false(grepl("zero evaluable rules", note, fixed = TRUE))
            session$setInputs(export_kind = "calls")
            expect_null(output$ssp_ror_notice)
            expect_false(any(grepl("SSP.ROR rule coverage", readLines(output$download), fixed = TRUE)))
        })
    }
})

test_that("optional mapping, prediction and alignment errors are cached and leave calls usable", {
    env <- shiny_app_env()
    se <- current_example("OSLO2EMIT0obj")$se_obj[, 1:3]
    clean <- ssp_small_upload(se)
    reference <- Mapping(se, include_ssp_ror = TRUE, verbose = FALSE)
    ror_only <- setdiff(rownames(reference$se_SSP_ROR), rownames(reference$se_SSP))
    ids <- as.character(SummarizedExperiment::rowData(se)$ENTREZID)
    # Two all-missing duplicate probes for an Entrez ID required only by ROR.
    rows <- match(head(intersect(ror_only[!ror_only %in% rownames(reference$se_NC)],
                                 ids[!duplicated(ids) & !duplicated(ids, fromLast = TRUE)]), 2L), ids)
    SummarizedExperiment::rowData(se)$ENTREZID[rows[2]] <- SummarizedExperiment::rowData(se)$ENTREZID[rows[1]]
    SummarizedExperiment::assay(se)[rows, ] <- NA_real_
    bad <- ssp_small_upload(se)

    state <- new.env()
    state$fault <- ""
    state$maps <- state$rors <- 0L
    map <- BreastSubtypeR::Mapping
    predict <- BreastSubtypeR::BS_sspbc
    local_mocked_bindings(Mapping = function(se_obj, ..., include_ssp_ror = FALSE) {
        state$maps <- state$maps + as.integer(include_ssp_ror)
        map(se_obj, ..., include_ssp_ror = include_ssp_ror)
    }, BS_sspbc = function(se_obj, ssp.name) {
        if (ssp.name != "ssp.ror") {
            if (state$fault == "subtype") stop("injected subtype error")
            return(predict(se_obj, ssp.name))
        }
        state$rors <- state$rors + 1L
        if (state$fault == "predict") stop("injected predictor error")
        result <- predict(se_obj, ssp.name)
        if (state$fault == "align") colnames(result$data.used) <- rev(colnames(result$data.used))
        result
    }, .package = "BreastSubtypeR")

    shiny::testServer(env$server, {
        session$setInputs(BSmethod = "sspbc", k_subtypes = "5", is_raw_counts = "norm")
        do.call(session$setInputs, bad)
        session$setInputs(map = 1)
        faults <- c("map", "predict", "align")
        for (i in seq_along(faults)) {
            state$fault <- faults[i]
            if (i == 2L) {
                do.call(session$setInputs, clean)
                expect_null(current_prediction())
                session$setInputs(map = 2)
            }
            session$setInputs(run = i)
            p <- current_prediction()
            calls <- ssp_read_export(output$download)
            choice <- p$ssp$full_choice
            session$setInputs(export_kind = choice)
            error <- c("Entrez ID", "injected predictor error", "Cannot align ssp.ror")[i]
            expect_match(output$ssp_ror_notice$html, error, fixed = TRUE)
            expect_null(output$download_control)
            expect_error(download_filename(), error, fixed = TRUE)
            path <- tempfile()
            expect_error(write_download(path), error, fixed = TRUE)
            expect_false(file.exists(path))
            expect_error(output$download, error, fixed = TRUE)
            expect_identical(current_prediction(), p)
            expect_match(output$plotSection$html, "Subtype calls", fixed = TRUE)
            expect_null(reactive_files$error)
            session$setInputs(export_kind = "calls")
            expect_null(output$ssp_ror_notice)
            expect_identical(ssp_read_export(output$download), calls)
            session$setInputs(export_kind = choice)
            expect_error(write_download(path), error, fixed = TRUE)
            expect_equal(state$maps, i)
            expect_equal(state$rors, i - 1L)
            session$setInputs(export_kind = "calls")
        }
        state$fault <- ""
        session$setInputs(run = 4)
        session$setInputs(export_kind = current_prediction()$ssp$full_choice)
        expect_null(output$ssp_ror_notice)
        expect_identical(download_filename(), "results-sspbc-5class-full.txt")
        expect_equal(ncol(ssp_read_export(output$download)), 32L)
        expect_equal(state$maps, 4L)
        expect_equal(state$rors, 3L)
        state$fault <- "subtype"
        session$setInputs(run = 5)
        expect_true(mapping_current())
        expect_null(current_prediction())
        expect_error(download_filename(), "No current prediction")
        expect_error(write_download(tempfile()), "No current prediction")
    })
})

test_that("both callbacks revoke SSP exports across input and prediction changes", {
    env <- shiny_app_env()
    se <- current_example("OSLO2EMIT0obj")$se_obj[, 1:2]
    uploads <- ssp_small_upload(se)
    for (change in c("GEX", "clinic", "anno", "example", "raw", "map", "classes", "method", "rerun")) {
        shiny::testServer(env$server, {
            session$setInputs(BSmethod = "sspbc", k_subtypes = "5", is_raw_counts = "norm")
            do.call(session$setInputs, uploads)
            session$setInputs(map = 1)
            session$setInputs(run = 1)
            old <- current_prediction()
            session$setInputs(export_kind = old$ssp$full_choice)
            expect_identical(download_filename(), "results-sspbc-5class-full.txt")
            path <- tempfile()
            write_download(path)
            expect_equal(ncol(ssp_read_export(path)), 32L)
            if (change %in% names(uploads)) {
                replacement <- uploads[change]
                replacement[[1]]$name <- "replacement.txt"
                do.call(session$setInputs, replacement)
            } else switch(change,
                example = session$setInputs(load_example_step1 = 1),
                raw = session$setInputs(is_raw_counts = "raw"),
                map = session$setInputs(map = 2),
                classes = session$setInputs(k_subtypes = "4"),
                method = session$setInputs(BSmethod = "AIMS"),
                rerun = session$setInputs(run = 2))
            if (change == "rerun") {
                expect_false(is.null(current_prediction()))
                expect_match(output$export_selector$html, 'value="calls" checked', fixed = TRUE)
                msg <- "does not provide this export"
            } else {
                expect_null(current_prediction())
                msg <- "No current prediction"
            }
            expect_null(output$download_control)
            expect_error(download_filename(), msg, fixed = TRUE)
            dest <- tempfile()
            expect_error(write_download(dest), msg, fixed = TRUE)
            expect_false(file.exists(dest))
            if (change == "classes") {
                # Bypass the clearing observer. The provenance check still wins.
                reactive_files$prediction <- old
                expect_error(write_download(dest), msg, fixed = TRUE)
            }
            if (change == "method") {
                session$setInputs(run = 2, export_kind = "full")
                expect_identical(current_prediction()$method, "AIMS")
                expect_false(grepl("Full metrics (incl. ROR)", output$export_selector$html, fixed = TRUE))
                expect_error(write_download(dest), "does not provide this export")
            }
            if (change == "rerun") {
                session$setInputs(export_kind = "calls")
                expect_identical(download_filename(), "results-sspbc-5class-calls.txt")
                session$setInputs(export_kind = current_prediction()$ssp$full_choice)
                expect_identical(download_filename(), "results-sspbc-5class-full.txt")
            }
        })
    }
})

test_that("SSP positional metrics reject inconsistent axes and align verified IDs to calls", {
    env <- shiny_app_env()
    se <- current_example("OSLO2EMIT0obj")$data_input$se_SSP[, c(6L, 2L, 1L)]
    native <- BS_sspbc(se, "ssp.pam50")
    ids <- colnames(se)
    project <- function(result, calls = ids) env$.ssp_export_metrics(result, "ssp.pam50", ids, calls)
    reordered <- project(native, rev(ids))
    expect_identical(rownames(reordered), rev(ids))
    expect_identical(reordered$SSP_Call_Prob, unname(rev(native$prob[, 1])))
    expect_error(project(native, c(ids[1], ids[1], ids[3])), "Cannot align")
    bad <- native
    rownames(bad$prob) <- rev(ids)
    expect_error(project(bad), "Cannot align")
    bad <- native
    bad$all.probs[["24"]] <- bad$all.probs[["24"]][-1L, , drop = FALSE]
    expect_error(project(bad), "Cannot align")
    bad <- native
    colnames(bad$rules.matrix) <- rev(ids)
    expect_error(project(bad), "Cannot align")
})
