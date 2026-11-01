test_that("compact historical fixtures restore independently frozen objects", {
    withr::local_seed(199)
    for (name in c("OSLO2EMIT0obj", "TCGABRCAobj")) {
        fixture <- readRDS(test_path("fixtures", "historical-cb7a7c70",
            paste0(name, ".rds")))
        current <- current_example(name)
        before <- current
        rng <- .Random.seed
        reference <- historical_fixture_reference(name)
        expect_identical(length(SummarizedExperiment::assays(
            fixture$historical$se_obj)), 0L, info = name)
        restored <- reconstruct_historical_example(fixture, current)
        expect_identical(names(restored), c("se_obj", "data_input", "res"), info = name)
        expect_identical(dim(restored$se_obj),
            c(reference$nrow, reference$ncol), info = name)
        expect_identical(historical_canonical_md5(restored),
            reference$object_md5, info = name)
        expect_identical(historical_canonical_md5(
            SummarizedExperiment::assays(restored$se_obj,
                withDimnames = FALSE)[[1L]]), reference$assay_md5, info = name)
        expect_identical(current, before, info = name)
        expect_identical(.Random.seed, rng, info = name)
        expect_true(methods::validObject(restored$se_obj), info = name)
    }
})

test_that("historical fixture sharing fails closed on public assay drift", {
    for (name in c("OSLO2EMIT0obj", "TCGABRCAobj")) {
        fixture <- readRDS(test_path("fixtures", "historical-cb7a7c70",
            paste0(name, ".rds")))
        current <- current_example(name)

        changed <- current
        values <- SummarizedExperiment::assay(changed$se_obj)
        values[1L, 1L] <- values[1L, 1L] + 1
        SummarizedExperiment::assay(changed$se_obj) <- values
        expect_error(reconstruct_historical_example(fixture, changed),
            "assay content/type integrity", info = name)

        changed <- current
        order <- seq_len(nrow(changed$se_obj))
        order[1:2] <- c(2L, 1L)
        changed$se_obj <- changed$se_obj[order, ]
        expect_error(reconstruct_historical_example(fixture, changed),
            "row/sample order integrity", info = name)

        changed <- current
        changed$se_obj <- changed$se_obj[, rev(seq_len(ncol(changed$se_obj)))]
        expect_error(reconstruct_historical_example(fixture, changed),
            "row/sample order integrity", info = name)

        changed <- current
        values <- SummarizedExperiment::assay(changed$se_obj)
        storage.mode(values) <- if (typeof(values) == "integer") "double" else "character"
        SummarizedExperiment::assay(changed$se_obj) <- values
        expect_error(reconstruct_historical_example(fixture, changed),
            "assay content/type integrity", info = name)
    }
})

test_that("all independently stored historical metadata and outputs remain pinned", {
    for (name in c("OSLO2EMIT0obj", "TCGABRCAobj")) {
        fixture <- readRDS(test_path("fixtures", "historical-cb7a7c70",
            paste0(name, ".rds")))
        current <- current_example(name)
        changed <- fixture
        S4Vectors::metadata(changed$historical$se_obj)$issue199_corruption <- TRUE
        expect_error(reconstruct_historical_example(changed, current),
            "complete-object integrity", info = name)

        changed <- fixture
        SummarizedExperiment::rowData(changed$historical$se_obj)$ENTREZID[1L] <- "changed"
        expect_error(reconstruct_historical_example(changed, current),
            "complete-object integrity", info = name)

        changed <- fixture
        SummarizedExperiment::colData(changed$historical$se_obj)$issue199_corruption <- TRUE
        expect_error(reconstruct_historical_example(changed, current),
            "complete-object integrity", info = name)

        changed <- fixture
        values <- SummarizedExperiment::assay(changed$historical$data_input$se_NC)
        values[1L, 1L] <- values[1L, 1L] + 1
        SummarizedExperiment::assay(changed$historical$data_input$se_NC) <- values
        expect_error(reconstruct_historical_example(changed, current),
            "complete-object integrity", info = name)

        changed <- fixture
        changed$historical$res$res_subtypes$entropy[1L] <-
            changed$historical$res$res_subtypes$entropy[1L] + 1
        expect_error(reconstruct_historical_example(changed, current),
            "complete-object integrity", info = name)
    }
})
