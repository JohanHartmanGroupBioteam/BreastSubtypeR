test_that("Spearman batching accepts only complete nonconstant plain columns", {
    gate <- BreastSubtypeR:::.can_batch_spearman
    samples <- matrix(seq_len(12L), 4L, 3L)
    centres <- as.data.frame(samples)
    expect_true(gate(centres, samples))
    expect_true(gate(as.matrix(centres), samples + 0.25))

    # A bare levels attribute is enough to change as.matrix.data.frame coercion.
    levels_column <- centres
    attr(levels_column[[1L]], "levels") <- letters[1:4]
    classed_column <- centres
    class(classed_column[[1L]]) <- "ssp_test_numeric"
    matrix_column <- centres
    matrix_column[[1L]] <- I(cbind(1:4, 4:1))
    frame_column <- centres
    frame_column[[1L]] <- centres[, 1:2]
    bad_centres <- list(
        NULL, 1:4, centres[1:3, ], centres[, FALSE, drop = FALSE],
        structure(centres, class = c("ssp_test_frame", "data.frame")),
        data.table::as.data.table(centres), S4Vectors::DataFrame(centres),
        levels_column, classed_column, matrix_column, frame_column,
        data.frame(label = letters[1:4]), data.frame(group = factor(1:4))
    )
    # Unsupported shapes, types and classes can occur on either side.
    bad_matrices <- list(
        as.logical(samples), matrix(TRUE, 4L, 2L), matrix(1+2i, 4L, 2L),
        matrix("1", 4L, 2L), array(1:12, c(4L, 3L, 1L)),
        structure(samples, class = "ssp_test_matrix"), samples[, integer(), drop = FALSE],
        samples[-1L, ], replace(samples, seq_len(4L), rep(c(0, -0), 2L))
    )
    for (value in c(NA_real_, NaN, Inf, -Inf)) {
        bad_matrices <- c(bad_matrices, list(replace(samples, 1L, value)))
    }
    for (x in c(bad_centres, bad_matrices)) expect_false(gate(x, samples))
    for (y in c(list(NULL, centres), bad_matrices)) expect_false(gate(centres, y))
    for (n in 0:1) {
        rows <- seq_len(n)
        expect_false(gate(centres[rows, , drop = FALSE], samples[rows, , drop = FALSE]))
    }
})

test_that("sspPredict executes the expected correlation at each reached site", {
    predict <- BreastSubtypeR:::sspPredict
    gate <- BreastSubtypeR:::.can_batch_spearman
    correlate <- stats::cor
    distance <- stats::dist

    # Probe the frozen first-centroid loop before mocking. Across R runtimes,
    # apply() can pass difftime columns to cor() with or without their class.
    cases <- ssp_predict_cases()
    classed <- cases$classed_matrix$args
    aligned <- environment(ssp_scalar_reference)$overlapSets(classed$x, classed$y)
    probe_calls <- 0L
    probe_failed <- withRestarts({
        apply(aligned$y, 2L, function(xx) {
            probe_calls <<- probe_calls + 1L
            -tryCatch(correlate(aligned$x[, 1L], xx, method = "spearman",
                                use = "pairwise.complete.obs"),
                      error = function(e) invokeRestart("scalar_cor_error"))
        })
        FALSE
    }, scalar_cor_error = function() TRUE)
    if (probe_failed) {
        expect_identical(probe_calls, 1L)
        cases$classed_matrix$sites <- c(main = "S")
    } else {
        expect_identical(probe_calls, ncol(aligned$y))
    }

    # Locate the active predictor frame once per call. Its output allocations
    # distinguish sites even when all three have the same operand dimensions.
    record_site <- function() {
        if (is.null(predict_frame)) {
            for (i in rev(seq_along(sys.frames()))) {
                if (identical(sys.function(i), predict)) {
                    predict_frame <<- sys.frame(i)
                    break
                }
            }
        }
        stopifnot(!is.null(predict_frame))
        if (exists("dist.RORSubtype", predict_frame, inherits = FALSE)) return("ror")
        if (exists("dist4", predict_frame, inherits = FALSE)) return("four")
        "main"
    }
    log_gate <- function(centres, samples) {
        eligible <- gate(centres, samples)
        key <- record_site()
        gate_events[[key]] <<- c(gate_events[[key]], eligible)
        gate_shapes[[key]] <<- c(dim(centres), dim(samples))
        eligible
    }
    log_cor <- function(x, y, use, method) {
        key <- record_site()
        # Record the call before delegation, including scalar calls that fail.
        if (identical(use, "all.obs") && identical(method, "spearman") &&
            is.matrix(x) && is.matrix(y)) {
            route <- "B"
            batch_shapes[[key]] <<- c(dim(x), dim(y))
        } else if (identical(use, "pairwise.complete.obs") &&
                   identical(method, expected_method[[key]])) {
            route <- "S"
        } else {
            route <- "unexpected"
        }
        events[[key]] <<- c(events[[key]], route)
        correlate(x, y, use = use, method = method)
    }
    log_dist <- function(x) {
        key <- record_site()
        events[[key]] <<- c(events[[key]], "E")
        distance(x)
    }
    # Intercept both the qualified matrix call and the imported scalar call.
    # All replacements expire with this test; delegate each call to base stats.
    local_mocked_bindings(cor = log_cor, .package = "stats")
    local_mocked_bindings(.can_batch_spearman = log_gate, cor = log_cor,
                          dist = log_dist, .package = "BreastSubtypeR")

    for (label in names(cases)) {
        item <- cases[[label]]
        args <- item$args
        for (subtype in c(TRUE, FALSE)) {
            predict_frame <- NULL
            events <- gate_events <- gate_shapes <- batch_shapes <- list()
            expected <- item$sites
            if (!subtype) expected <- expected[names(expected) != "four"]
            method <- if (args$distm == "spearman") "spearman" else "pearson"
            expected_method <- c(main = method, four = method, ror = "spearman")
            # Reference-equivalence tests compare all conditions. Here, keep
            # correlation calls up to and including the call that fails.
            failed <- tryCatch({
                suppressWarnings(do.call(predict, c(args, list(Subtype = subtype))))
                FALSE
            }, error = function(e) TRUE)
            info <- paste(label, "Subtype", subtype)
            if (label == "classed_matrix") {
                expect_identical(failed, probe_failed, info = info)
                if (probe_failed) expect_identical(events$main, "S", info = info)
            }
            expect_identical(as.character(names(events)), as.character(names(expected)),
                             info = info)
            gated <- if (args$distm == "spearman") expected else expected[names(expected) == "ror"]
            expect_identical(as.character(names(gate_events)), as.character(names(gated)),
                             info = info)
            for (site in names(gated)) {
                expect_identical(gate_events[[site]], unname(gated[[site]] == "B"), info = info)
            }

            retained <- intersect(rownames(args$x), rownames(args$y))
            for (site in names(expected)) {
                path <- expected[[site]]
                expect_identical(unique(events[[site]]), path, info = paste(info, site))
                class_count <- if (site == "main") ncol(args$x) else min(4L, ncol(args$x))
                if (path == "B") {
                    genes <- if (site == "ror") {
                        setdiff(retained, c("BIRC5", "CCNB1", "GRB7", "MYBL2"))
                    } else retained
                    shape <- c(length(genes), class_count, length(genes), ncol(args$y))
                    expect_identical(gate_shapes[[site]], shape, info = info)
                    expect_identical(batch_shapes[[site]], shape, info = info)
                    expect_length(events[[site]], 1L)
                } else if (!failed) {
                    calls <- if (path == "E") class_count else class_count * ncol(args$y)
                    expect_length(events[[site]], calls)
                }
            }
        }
    }
})
