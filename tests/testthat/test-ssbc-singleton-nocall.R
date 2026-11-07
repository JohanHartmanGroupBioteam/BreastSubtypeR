test_that("the high-node exception applies only to the deliberate singleton no-call", {
    env <- new.env(parent = emptyenv())
    utils::data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = env)
    se <- env$OSLO2EMIT0obj$data_input$se_NC[, c(1, 2, 5, 7)]
    se$ER <- c("ER+", "ER-", "ER+", "ER+")
    se$TSIZE <- rep(1, 4)
    se$NODE <- rep(4, 4)
    singleton <- colnames(se)[2]
    populated <- which(se$ER == "ER+")
    for (four in c(FALSE, TRUE)) {
        control <- expect_experimental_ror_pc(
            BS_ssBC(se[, populated], "ER", Subtype = four, hasClinical = TRUE))
        expect_warning(mixed <- expect_experimental_ror_pc(
            BS_ssBC(se, "ER", Subtype = four, hasClinical = TRUE)), "has n=1")
        scores <- mixed$score.ROR
        expect_identical(scores[singleton, ror221_fields[2]], NA_character_)
        expect_identical(scores[singleton, ror221_fields[1]], NA_real_)
        expect_true(all(is.na(scores[singleton, setdiff(names(scores), "patientID")])))
        expect_identical(scores[populated, ror221_fields[2]], rep("high", length(populated)))
        expect_false(anyNA(scores$Call[populated]))
        expect_identical(scores[populated, , drop = FALSE], control$score.ROR)
        expect_identical(mixed$BS.all[populated, , drop = FALSE], control$BS.all)
    }
})

test_that("manual and AUTO summaries align singleton no-calls and use available votes", {
    env <- new.env(parent = emptyenv())
    utils::data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = env)
    di <- env$OSLO2EMIT0obj$data_input
    se <- di$se_NC
    pos <- colnames(se)[se$ER == "ER+"]
    neg <- colnames(se)[se$ER == "ER-"]
    id <- neg[which(se$HER2[match(neg, colnames(se))] == "HER2+")[1]]
    withr::local_seed(73)
    for (auto in c(FALSE, TRUE)) {
        ids <- if (auto) rev(c(pos[1:20], neg)) else
            c(pos[se$HER2[match(pos, colnames(se))] == "HER2-"][1:11], id)
        input <- lapply(di, function(x) x[, ids, drop = FALSE])
        # Seven known ER+/HER2- and one ER-/HER2+, below AUTO subset minima.
        if (auto) {
            h <- rep(NA_character_, length(ids))
            h[match(pos[1:7], ids)] <- "HER2-"
            h[match(id, ids)] <- "HER2+"
            input$se_NC$HER2 <- h
        }
        input$se_NC$NODE <- rep(4, length(ids))
        input$se_NC$TSIZE <- rep(1, length(ids))
        methods <- if (auto) "AUTO" else c("ssBC", "ssBC.v2")
        before <- .Random.seed
        output <- suppressMessages(suppressWarnings(BS_Multi(input, methods = methods,
            hasClinical = TRUE, Subtype = TRUE)))
        expect_identical(.Random.seed, before)
        for (m in if (auto) "ssBC.v2" else methods) {
            result <- output$results[[m]]
            expect_identical(rownames(result$score.ROR), ids)
            expect_identical(result$score.ROR[id, ror221_fields[2]], NA_character_)
            expect_true(all(is.na(result$outList$testData[, id])))
            for (table in c("res_subtypes", "res_subtypes.Subtype"))
                expect_identical(output[[table]][id, m], NA_character_)
        }
        for (table in c("res_subtypes", "res_subtypes.Subtype")) {
            tab <- output[[table]]
            expect_identical(rownames(tab), ids)
            votes <- unlist(tab[id, setdiff(names(tab), "entropy"), drop = FALSE])
            counts <- as.numeric(base::table(votes))
            if (auto) expect_gt(sum(counts), 0) else expect_length(counts, 0)
            p <- counts / sum(counts)
            expected <- if (length(p)) -sum(p * log2(p)) else NA_real_
            expect_equal(tab[id, "entropy"], expected)
        }
    }
})

test_that("Shiny displays and exports actual all-no-call ssBC results", {
    skip_if_not_installed("shiny")
    skip_if_not_installed("bslib")
    env <- new.env(parent = asNamespace("shiny"))
    eval(parse(system.file("ShinyBreastSubtypeR/server.R", package = "BreastSubtypeR"),
               encoding = "UTF-8", keep.source = FALSE), env)
    utils::data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = env)
    input <- lapply(env$OSLO2EMIT0obj$data_input, function(x) x[, 6:1, drop = FALSE])
    ids <- colnames(input$se_NC)
    input$se_NC$TN <- c("nonTN", "TN", rep("nonTN", 4))
    input$se_NC$NODE <- rep(4, 6)
    input$se_NC$TSIZE <- rep(1, 6)
    input_data_for_test <- input
    shiny::testServer(env$server, {
        session$setInputs(BSmethod = "ssBC", s = "TN.v2", k_subtypes = "4",
                          hasClinical = TRUE, is_raw_counts = "norm", export_kind = "calls")
        reactive_files$data_input <- input_data_for_test
        reactive_files$mapped_key <- input_key()
        session$flushReact()
        session$setInputs(run = 1)
        expect_null(reactive_files$error)
        expect_false(is.null(output$pie1))
        expect_false(is.null(output$heat2))
        pie <- Vis_pie(current_prediction()$plot)
        expect_identical(pie$labels$title, "No assigned subtype calls")
        for (kind in c("calls", "full")) {
            session$setInputs(export_kind = kind)
            tab <- utils::read.delim(output$download, comment.char = "#", check.names = FALSE)
            expect_identical(tab$PatientID, ids)
            expect_true(all(is.na(tab$Call_4class)))
            if (kind == "full") {
                expect_true(all(is.na(tab[[ror221_fields[2]]])))
                expect_true(all(is.na(tab$Confidence)))
            }
        }
    })
})
