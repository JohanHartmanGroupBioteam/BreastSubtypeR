# Manual performance experiment; never run by testthat or R CMD check.
# Run from the repository root with Rscript --vanilla tools/issue-219/benchmark.R
# Optional arguments: baseline|candidate|both output-directory repetitions
# First run baseline to record pre-change profiles/allocations. Then run both
# for alternating paired timings and candidate profiles/allocations. The both
# mode does not repeat baseline profiling (which is deliberately untimed).
# Baseline is read from the pinned upstream commit, not a locally installed build.
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1L]] else "both"
out <- if (length(args) > 1L) args[[2L]] else "tools/issue-219/results"
repetitions <- if (length(args) > 2L) as.integer(args[[3L]]) else 3L
stopifnot(mode %in% c("baseline", "candidate", "both"), repetitions >= 3L)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
pkgload::load_all(quiet = TRUE, export_all = FALSE)
ns <- asNamespace("BreastSubtypeR")
base_ref <- "f878f58c1c12ca6f1cd5bf4489667935bf25969b"
base_source <- system2("git", c("show", paste0(base_ref, ":R/mapping_imputation.R")),
    stdout = TRUE)
stopifnot(is.null(attr(base_source, "status")))
base_env <- new.env(parent = ns)
eval(parse(text = base_source), envir = base_env)
baseline <- base_env$.mapping_donor_matrix
candidate <- get(".mapping_donor_matrix", ns)
reference <- new.env(parent = emptyenv())
data("BreastSubtypeRobj", envir = reference, package = "BreastSubtypeR")
reference <- reference$BreastSubtypeRobj
required <- get(".mapping_required_ids", ns)(reference)
target <- intersect(as.character(reference$genes.sig50$EntrezGene.ID),
    as.character(reference$genes.signature$EntrezGene.ID[
        reference$genes.signature$SSP_based == "Yes"]))[[1L]]

fixture <- function(n, shape) {
    set.seed(219L + n)
    x <- matrix(rnorm(n * 100L, mean = 5, sd = 1), nrow = n,
        dimnames = list(paste0("probe", seq_len(n)), paste0("sample", seq_len(100L))))
    # Keep the sole represented target a singleton in both shapes.
    ids <- as.character(9000000L + seq_len(n))
    if (shape == "duplicate-heavy") {
        start <- as.integer(n / 10L) + 1L
        ids[start:n] <- rep(as.character(10000000L + seq_len((n - start + 1L) / 2L)),
            each = 2L)
    }
    ids[[1L]] <- target
    stopifnot(sum(ids %in% required) == 1L)
    x[1L, 50L] <- NA_real_
    y <- data.frame(probe = rownames(x), ENTREZID = as.integer(ids))
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = x),
        rowData = S4Vectors::DataFrame(y),
        colData = S4Vectors::DataFrame(PatientID = colnames(x)))
    list(x = x, y = y, se = se, unique_genes = length(unique(ids)))
}

with_donor <- function(fun, code) {
    old <- get(".mapping_donor_matrix", ns)
    was_locked <- bindingIsLocked(".mapping_donor_matrix", ns)
    if (was_locked) unlockBinding(".mapping_donor_matrix", ns)
    assign(".mapping_donor_matrix", fun, ns)
    if (was_locked) lockBinding(".mapping_donor_matrix", ns)
    on.exit({
        if (was_locked) unlockBinding(".mapping_donor_matrix", ns)
        assign(".mapping_donor_matrix", old, ns)
        if (was_locked) lockBinding(".mapping_donor_matrix", ns)
    })
    force(code)
}

profile <- function(fun, prefix) {
    path <- tempfile(fileext = ".Rprof")
    on.exit(unlink(path))
    gc()
    Rprof(path, interval = 0.01, memory.profiling = TRUE)
    started <- proc.time()[["elapsed"]]
    iterations <- 0L
    repeat {
        value <- fun()
        iterations <- iterations + 1L
        if (proc.time()[["elapsed"]] - started >= 2) break
    }
    Rprof(NULL)
    p <- summaryRprof(path, memory = "both")
    write.csv(p$by.total, paste0(prefix, "-profile-total.csv"))
    write.csv(p$by.self, paste0(prefix, "-profile-self.csv"))
    writeLines(paste("Sampling time:", p$sampling.time, "seconds; interval:",
        p$sample.interval, "; iterations:", iterations),
        paste0(prefix, "-profile.txt"))
    invisible(value)
}

allocation <- function(fun) {
    # A measurable lower bound on cumulative allocations, not peak resident RAM.
    # Threshold avoids enormous logs from small, per-group temporary objects.
    path <- tempfile(fileext = ".Rprofmem")
    on.exit(unlink(path))
    gc()
    Rprofmem(path, threshold = 10000)
    value <- fun()
    Rprofmem(NULL)
    lines <- readLines(path, warn = FALSE)
    bytes <- suppressWarnings(as.numeric(sub(" .*", "", lines)))
    c(large_allocation_bytes = sum(bytes, na.rm = TRUE),
        large_allocation_count = sum(!is.na(bytes)),
        result_bytes = as.numeric(object.size(value)))
}

writeLines(c(paste("baseline:", base_ref), paste("mode:", mode),
    paste("date:", Sys.time()),
    paste("method: max; represented target:", target, "; samples: 100"),
    paste("candidate source MD5:", tools::md5sum("R/mapping_imputation.R")),
    capture.output(sessionInfo()),
    paste("OS:", paste(Sys.info(), collapse = " | "))),
    file.path(out, paste0(mode, "-session.txt")))
implementations <- switch(mode, baseline = list(baseline = baseline),
    candidate = list(candidate = candidate), both = list(baseline = baseline, candidate = candidate))
timings <- list()
allocations <- list()
for (n in c(20000L, 5000L)) for (shape in c("singleton-heavy", "duplicate-heavy")) {
    f <- fixture(n, shape)
    # Always gate candidate performance on exact donor and public-output equality.
    expected <- baseline(f$x, f$y, "max", required)
    if (mode != "baseline") {
        stopifnot(identical(expected, candidate(f$x, f$y, "max", required)))
        old_output <- with_donor(baseline, BreastSubtypeR::Mapping(f$se,
            impute = TRUE, verbose = FALSE))
        new_output <- with_donor(candidate, BreastSubtypeR::Mapping(f$se,
            impute = TRUE, verbose = FALSE))
        stopifnot(identical(old_output, new_output))
    }
    estimator <- function() get("knn_vectorized_reference", ns)(expected$matrix, target)
    stage_function <- function(implementation, stage) {
        donor <- implementations[[implementation]]
        switch(stage,
            preparation = function() donor(f$x, f$y, "max", required),
            estimator = estimator,
            Mapping = function() with_donor(donor, BreastSubtypeR::Mapping(f$se,
                impute = TRUE, verbose = FALSE)))
    }
    for (stage in c("preparation", "estimator", "Mapping")) {
        # The equality calls above warm both preparation and Mapping. Warm the
        # timing control separately; baseline-only Mapping has no equality call.
        # Alternate paired execution order to reduce systematic drift.
        if (stage == "estimator" || (mode == "baseline" && stage == "Mapping"))
            for (implementation in names(implementations))
                invisible(stage_function(implementation, stage)())
        for (iteration in seq_len(repetitions)) {
            execution_order <- names(implementations)
            if (iteration %% 2L == 0L) execution_order <- rev(execution_order)
            for (implementation in execution_order) {
                fun <- stage_function(implementation, stage)
                gc()
                elapsed <- unname(system.time(value <- fun())[["elapsed"]])
                record <- data.frame(implementation, rows = n, samples = 100L, shape,
                    unique_genes = f$unique_genes, stage, iteration, elapsed)
                timings[[length(timings) + 1L]] <- record
                write.csv(do.call(rbind, timings), file.path(out, paste0(mode, "-timings.csv")),
                    row.names = FALSE)
                cat(implementation, n, shape, stage, iteration, elapsed, "seconds\n")
                flush.console()
            }
        }
        for (implementation in names(implementations)) {
            if (n == 20000L &&
                (mode != "both" || implementation == "candidate")) {
                fun <- stage_function(implementation, stage)
                prefix <- file.path(out, paste(implementation, n, shape, stage, sep = "-"))
                profile(fun, prefix)
                cat("PROFILE", prefix, "\n")
                flush.console()
                if (stage != "estimator") {
                    memory <- allocation(fun)
                    allocations[[length(allocations) + 1L]] <- data.frame(implementation,
                        rows = n, shape, stage, as.list(memory))
                    write.csv(do.call(rbind, allocations),
                        file.path(out, paste0(mode, "-allocations.csv")), row.names = FALSE)
                }
            }
        }
    }
}
medians <- aggregate(elapsed ~ implementation + rows + samples + shape + unique_genes + stage,
    data = do.call(rbind, timings), FUN = median)
write.csv(medians, file.path(out, paste0(mode, "-medians.csv")), row.names = FALSE)
print(medians, row.names = FALSE)

if (mode == "both") {
    all_timings <- do.call(rbind, timings)
    baseline_memory_path <- file.path(out, "baseline-allocations.csv")
    baseline_memory <- if (file.exists(baseline_memory_path))
        read.csv(baseline_memory_path) else NULL
    memory <- rbind(baseline_memory, do.call(rbind, allocations))
    cases <- unique(all_timings[c("rows", "samples", "shape", "unique_genes", "stage")])
    measurements <- lapply(seq_len(nrow(cases)), function(i) {
        case <- cases[i, ]
        selected <- all_timings[all_timings$rows == case$rows &
            all_timings$shape == case$shape & all_timings$stage == case$stage, ]
        result <- case
        for (implementation in names(implementations)) {
            times <- selected[selected$implementation == implementation, ]
            times <- times$elapsed[order(times$iteration)]
            for (j in seq_along(times))
                result[[paste0(implementation, "_elapsed_", j, "_s")]] <- times[[j]]
            result[[paste0(implementation, "_median_s")]] <- median(times)
            allocation_row <- memory[memory$implementation == implementation &
                memory$rows == case$rows & memory$shape == case$shape &
                memory$stage == case$stage, ]
            for (field in c("large_allocation_bytes", "large_allocation_count", "result_bytes"))
                result[[paste0(implementation, "_", field)]] <-
                    if (nrow(allocation_row)) allocation_row[[field]][[1L]] else NA_real_
        }
        # The estimator is an unchanged timing control, not an optimization.
        result$speedup <- if (case$stage == "estimator") NA_real_ else
            result$baseline_median_s / result$candidate_median_s
        result
    })
    write.csv(do.call(rbind, measurements), file.path(out, "measurements.csv"),
        row.names = FALSE)
}
