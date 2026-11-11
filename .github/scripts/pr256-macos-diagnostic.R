args <- commandArgs(TRUE)
mode <- args[1L]
out <- normalizePath(args[2L], winslash = "/", mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE,
    null = "null", na = "null", digits = NA)
strict <- function(a, b) identical(a, b, num.eq = FALSE)
if (mode == "compare") {
    roots <- args[3:4]
    files <- lapply(roots, function(p) list.files(file.path(p, "collect/pairs"),
        pattern = "[.]rds$", recursive = TRUE))
    stopifnot(identical(files[[1]], files[[2]]))
    equal <- lapply(files[[1]], function(f) {
        a <- readRDS(file.path(roots[1], "collect/pairs", f))
        b <- readRDS(file.path(roots[2], "collect/pairs", f))
        data.frame(pair = f, actual_identical = strict(a$actual, b$actual),
            reference_identical = strict(a$reference, b$reference))
    })
    write.csv(do.call(rbind, equal), file.path(out, "baseline-candidate-pairs.csv"), row.names = FALSE)
    checks <- lapply(roots, function(p) jsonlite::read_json(file.path(p, "collect/checks.json"), simplifyVector = TRUE))
    json(list(pair_count = length(equal), all_pair_payloads_identical =
        all(vapply(equal, function(x) x$actual_identical && x$reference_identical, logical(1))),
        diagnostic_checks_identical = strict(checks[[1]], checks[[2]])),
        file.path(out, "baseline-candidate-summary.json"))
    quit(status = 0L)
}
lib <- normalizePath(args[3L], winslash = "/")
source <- normalizePath(args[4L], winslash = "/")
python_script <- normalizePath(args[5L], winslash = "/", mustWork = FALSE)
.libPaths(c(lib, .libPaths()))
library(BreastSubtypeR)
library(testthat)
stopifnot(normalizePath(find.package("BreastSubtypeR"), winslash = "/") ==
    paste0(lib, "/BreastSubtypeR"))
tests <- system.file("tests/testthat", package = "BreastSubtypeR")
stopifnot(nzchar(tests))
selected <- c("test-expanded-examples.R", "test-mapping-factor-probes.R",
    "test-shiny-example-parity.R")
stopifnot(all(file.exists(file.path(tests, selected))))
ip <- installed.packages()
deps <- ip[ip[, "Package"] != "BreastSubtypeR", c("Package", "Version", "LibPath", "Built"), drop = FALSE]
deps <- deps[order(deps[, "Package"], deps[, "LibPath"]), , drop = FALSE]
saveRDS(deps, file.path(out, "dependency-inventory.rds"))
write.csv(deps, file.path(out, "dependency-inventory.csv"))
capture.output(sessionInfo(), file = file.path(out, "session-info.txt"))
capture.output(.Machine, capabilities(), extSoftVersion(),
    file = file.path(out, "machine-capabilities.txt"))
saveRDS(list(machine = .Machine, capabilities = capabilities(), R = R.version,
    locale = Sys.getlocale(), rng_kind = RNGkind(), lib_paths = .libPaths(),
    package_path = find.package("BreastSubtypeR"), pid = Sys.getpid()),
    file.path(out, "environment.rds"))
json(list(source = source, library = lib, package_path = find.package("BreastSubtypeR"),
    pid = Sys.getpid(), R = R.version.string, arch = R.version$arch,
    Bioconductor = as.character(BiocManager::version()),
    test_file_md5 = as.list(tools::md5sum(file.path(tests, selected)))),
    file.path(out, "process.json"))
if (mode == "focused") {
    results <- testthat::test_dir(tests,
        filter = "^(expanded-examples|mapping-factor-probes|shiny-example-parity)$",
        reporter = "silent", package = "BreastSubtypeR", load_package = "installed",
        load_helpers = TRUE, stop_on_failure = FALSE, stop_on_warning = FALSE)
    saveRDS(results, file.path(out, "focused-results.rds"))
    d <- as.data.frame(results)
    write.csv(d[, !vapply(d, is.list, logical(1)), drop = FALSE],
        file.path(out, "focused-results.csv"), row.names = FALSE)
    json(list(files = unique(d$file), passes = sum(d$passed),
        failed_expectations = sum(d$failed), errors = sum(d$error),
        warnings = sum(d$warning), skips = sum(d$skipped),
        test_blocks = nrow(d)), file.path(out, "focused-summary.json"))
    # Expected test failures are data; collection failure still raises an error.
    quit(status = 0L)
}
stopifnot(mode == "collect")
setwd(tests)
Sys.setenv(TESTTHAT = "true")
env <- new.env(parent = globalenv())
testthat::source_test_helpers(tests, env)
for (f in selected) {
    for (expr in parse(file.path(tests, f))) {
        if (is.call(expr) && identical(expr[[1L]], as.name("test_that"))) break
        eval(expr, env)
    }
}
options(digits = 22)
checks <- list()
stages <- list()
check <- function(label, value, detail = NULL) {
    checks[[length(checks) + 1L]] <<- list(check = label, passed = isTRUE(value), detail = detail)
    invisible(value)
}
stage <- function(name, expr) {
    tryCatch({
        value <- force(expr)
        stages[[name]] <<- list(completed = TRUE)
        json(stages, file.path(out, "stages.json"))
        value
    }, error = function(e) {
        stages[[name]] <<- list(completed = FALSE, error = conditionMessage(e),
            calls = vapply(sys.calls(), function(x) paste(deparse(x), collapse = " "), ""))
        json(stages, file.path(out, "stages.json"))
        message("Diagnostic stage failed: ", name, ": ", conditionMessage(e))
        NULL
    })
}
bits <- function(x) matrix(as.integer(writeBin(as.double(x), raw(), size = 8L,
    endian = "big")), nrow = 8L)
hex <- function(x) apply(bits(x), 2L, function(z) paste(sprintf("%02x", z), collapse = ""))
ulp <- function(a, b) {
    ordered <- function(x) {
        z <- bits(x)
        hi <- as.vector(c(256^3, 256^2, 256, 1) %*% z[1:4, , drop = FALSE])
        lo <- as.vector(c(256^3, 256^2, 256, 1) %*% z[5:8, , drop = FALSE])
        neg <- hi >= 2^31
        hi[!neg] <- hi[!neg] + 2^31
        hi[neg] <- 2^32 - 1 - hi[neg]
        lo[neg] <- 2^32 - 1 - lo[neg]
        list(hi = hi, lo = lo)
    }
    x <- ordered(a); y <- ordered(b)
    abs((x$hi - y$hi) * 2^32 + (x$lo - y$lo))
}
stopifnot(ulp(1, 1 + .Machine$double.eps) == 1,
    ulp(-1, -1 - .Machine$double.eps) == 1,
    ulp(2 - .Machine$double.eps, 2) == 1, ulp(-0, 0) == 1)
kind <- function(x) ifelse(is.nan(x), "NaN", ifelse(is.na(x), "NA",
    ifelse(is.infinite(x), ifelse(x > 0, "+Inf", "-Inf"), "finite")))
pair_number <- 0L
pair_index <- list()
pair <- function(label, actual, reference) {
    pair_number <<- pair_number + 1L
    id <- sprintf("%03d-%s", pair_number, gsub("[^A-Za-z0-9._-]", "_", label))
    d <- file.path(out, "pairs", id)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    saveRDS(list(actual = actual, reference = reference), file.path(d, "complete.rds"))
    leaves <- list()
    differences <- list()
    walk <- function(a, b, path) {
        if (strict(a, b)) return(invisible(NULL))
        if (isS4(a) && isS4(b)) {
            slots <- union(methods::slotNames(a), methods::slotNames(b))
            ax <- attributes(a); bx <- attributes(b)
            ax[methods::slotNames(a)] <- NULL; bx[methods::slotNames(b)] <- NULL
            walk(ax, bx, paste0(path, "/@attributes"))
            for (n in slots) {
                if (n %in% methods::slotNames(a) && n %in% methods::slotNames(b))
                    walk(methods::slot(a, n), methods::slot(b, n), paste0(path, "/", n))
                else leaves[[length(leaves) + 1L]] <<- list(path = paste0(path, "/", n), missing_slot = TRUE)
            }
            return(invisible(NULL))
        }
        meta <- list(path = path, type_actual = typeof(a), type_reference = typeof(b),
            class_exact = strict(class(a), class(b)), type_exact = identical(typeof(a), typeof(b)),
            length_actual = length(a), length_reference = length(b),
            attributes_exact = strict(attributes(a), attributes(b)),
            names_exact = strict(names(a), names(b)), dimensions_exact = strict(dim(a), dim(b)),
            dimnames_exact = strict(dimnames(a), dimnames(b)))
        if (is.list(a) && is.list(b)) {
            if (!meta$attributes_exact || !meta$type_exact || length(a) != length(b))
                leaves[[length(leaves) + 1L]] <<- meta
            for (j in seq_len(min(length(a), length(b)))) {
                n <- if (length(names(a))) names(a)[j] else as.character(j)
                walk(a[[j]], b[[j]], paste0(path, "/", n))
            }
            return(invisible(NULL))
        }
        if (is.numeric(a) && is.numeric(b) && length(a) == length(b)) {
            x <- as.vector(a); y <- as.vector(b)
            ka <- kind(x); kb <- kind(y)
            finite <- is.finite(x) & is.finite(y)
            changed <- ka != kb
            changed[finite] <- x[finite] != y[finite]
            zeros <- finite & x == 0 & y == 0
            changed[zeros] <- 1 / x[zeros] != 1 / y[zeros]
            ix <- which(changed)
            fi <- which(finite & changed)
            delta <- abs(x[fi] - y[fi])
            rel <- delta / pmax(abs(x[fi]), abs(y[fi]))
            units <- if (length(fi)) ulp(x[fi], y[fi]) else numeric()
            meta$nonfinite_patterns_exact <- identical(ka, kb)
            meta$nonfinite_values_exact <- strict(x[!finite], y[!finite])
            meta$zero_signs_exact <- strict(1 / x[zeros], 1 / y[zeros])
            meta$differing_values <- length(ix)
            meta$max_absolute <- if (length(delta)) max(delta) else 0
            meta$max_relative <- if (length(rel)) max(rel) else 0
            meta$max_ulp_distance <- if (length(units)) max(units) else 0
            meta$ulp_integer_exact <- all(units <= 2^53)
            meta$global_ranks_ties_exact <- strict(rank(x, ties.method = "min", na.last = "keep"),
                rank(y, ties.method = "min", na.last = "keep"))
            if (is.matrix(a) && identical(dim(a), dim(b))) {
                ranks <- function(z) t(apply(z, 1L, rank, ties.method = "min", na.last = "keep"))
                meta$row_ranks_ties_exact <- strict(ranks(a), ranks(b))
            }
            if (length(ix)) {
                nr <- if (is.matrix(a)) nrow(a) else length(a)
                row <- (ix - 1L) %% nr + 1L
                column <- if (is.matrix(a)) (ix - 1L) %/% nr + 1L else rep(NA_integer_, length(ix))
                absd <- reld <- ud <- rep(NA_real_, length(ix))
                ok <- finite[ix]
                absd[ok] <- abs(x[ix[ok]] - y[ix[ok]])
                reld[ok] <- absd[ok] / pmax(abs(x[ix[ok]]), abs(y[ix[ok]]))
                ud[ok] <- ulp(x[ix[ok]], y[ix[ok]])
                differences[[length(differences) + 1L]] <<- data.frame(path = path, index = ix,
                    row = row, column = column,
                    row_id = if (is.matrix(a) && !is.null(rownames(a))) rownames(a)[row] else
                        if (!is.null(names(a))) names(a)[ix] else NA_character_,
                    column_id = if (is.matrix(a) && !is.null(colnames(a))) colnames(a)[column] else NA_character_,
                    actual = sprintf("%.17g", x[ix]), reference = sprintf("%.17g", y[ix]),
                    actual_hex = hex(x[ix]), reference_hex = hex(y[ix]),
                    actual_kind = ka[ix], reference_kind = kb[ix],
                    absolute = absd, relative = reld, ulp_distance = ud, stringsAsFactors = FALSE)
            }
        } else {
            meta$value_exact <- FALSE
        }
        leaves[[length(leaves) + 1L]] <<- meta
        invisible(NULL)
    }
    walk(actual, reference, "")
    if (length(differences)) write.csv(do.call(rbind, differences),
        file.path(d, "numeric-differences.csv"), row.names = FALSE)
    json(list(label = label, identical = strict(actual, reference), unequal_components = leaves),
        file.path(d, "summary.json"))
    pair_index[[id]] <<- list(label = label, identical = strict(actual, reference))
    json(pair_index, file.path(out, "pair-index.json"))
    invisible(strict(actual, reference))
}
decisions <- function(label, a, b) {
    rows <- function(x) t(apply(as.matrix(x), 1L, function(v) rank(-v, ties.method = "min", na.last = "keep")))
    ranks <- function(tag, x, y) {
        check(paste(label, tag, "ranks/ties"), strict(rows(x), rows(y)))
        check(paste(label, tag, "argmax"), strict(max.col(as.matrix(x), ties.method = "first"),
            max.col(as.matrix(y), ties.method = "first")))
    }
    check(paste(label, "res_subtypes"), strict(a$res_subtypes, b$res_subtypes))
    pair(paste(label, "res_subtypes"), a$res_subtypes, b$res_subtypes)
    for (method in intersect(names(a$results), names(b$results))) {
        x <- a$results[[method]]; y <- b$results[[method]]
        if (method %in% c("AIMS", "sspbc")) {
            for (n in c("cl", "rules.matrix", "prob"))
                check(paste(label, method, n), strict(x[[n]], y[[n]]))
            for (k in names(x$all.probs)) ranks(paste(method, k), x$all.probs[[k]], y$all.probs[[k]])
        } else {
            for (n in c("distances", "dist.RORSubtype")) ranks(paste(method, n), x$outList[[n]], y$outList[[n]])
            check(paste(label, method, "Confidence"), strict(x$score.ROR$Confidence, y$score.ROR$Confidence))
            for (n in grep("^ROR-[SP] ", names(x$score.ROR), value = TRUE)) {
                if (!is.numeric(x$score.ROR[[n]])) next
                cuts <- 100 * ((if (startsWith(n, "ROR-S ")) c(-.15, .1) else c(-.25, .1)) + .35) / .85
                check(paste(label, method, n, "cutoff sides"), strict(outer(x$score.ROR[[n]], cuts, "<"), outer(y$score.ROR[[n]], cuts, "<")) &&
                    strict(outer(x$score.ROR[[n]], cuts, ">"), outer(y$score.ROR[[n]], cuts, ">")))
            }
            groups <- grep("^ROR-[SP] Group ", names(x$score.ROR), value = TRUE)
            check(paste(label, method, "risk groups"), strict(x$score.ROR[, groups, drop = FALSE], y$score.ROR[, groups, drop = FALSE]))
        }
    }
}
stage("mapping", {
    capture.output(BreastSubtypeR:::duplicate_genes, mean.default, colMeans,
        file = file.path(out, "collapse-functions.txt"))
    x <- rbind(z = c(8, 2, 2), b = c(2, 8, 2), a = c(0, 0, 0),
        q = c(3, 4, 5), background = c(17, 16, 15))
    colnames(x) <- c("S3", "S1", "S2")
    entrez <- c(z = 596L, b = 596L, a = 596L, q = 332L, background = 800001L)
    cases <- list()
    for (permuted in c(FALSE, TRUE)) {
        values <- x; labels <- rownames(x)
        if (permuted) {
            values <- x[c("q", "a", "background", "z", "b"), c("S2", "S3", "S1")]
            labels <- c("b", "background", "q", "z", "a")
        }
        se <- env$map01_probe_input(values, labels, unname(entrez[labels]))
        encodings <- list(ordinary = factor(labels),
            reversed = factor(labels, levels = rev(levels(factor(labels)))),
            ordered = ordered(labels), unused = factor(labels, levels = c("unused", labels)))
        duplicates <- values[labels[entrez[labels] == 596L], , drop = FALSE]
        cm <- colMeans(duplicates)
        mn <- apply(duplicates, 2L, mean, na.rm = TRUE)
        controls <- list()
        for (method in c("max", "mean", "median", "stdev", "iqr")) {
            control <- Mapping(se, method = method, verbose = FALSE)
            controls[[method]] <- control
            for (encoding in names(encodings)) {
                input <- se
                SummarizedExperiment::rowData(input)$probe <- encodings[[encoding]]
                before <- serialize(input, NULL)
                mapped <- Mapping(input, method = method, verbose = FALSE)
                check(paste("mapping factor exact", permuted, method, encoding), strict(mapped, control))
                check(paste("mapping input preserved", permuted, method, encoding), identical(serialize(input, NULL), before))
                if (!strict(mapped, control)) pair(paste("factor", permuted, method, encoding), mapped, control)
            }
        }
        prod <- as.matrix(SummarizedExperiment::assay(controls$mean$se_NC))["BCL2", ]
        prod_ssp <- SummarizedExperiment::assay(controls$mean$se_SSP)["596", ]
        pair(paste("mapping", permuted, "colMeans-vs-mean"), cm, mn)
        pair(paste("mapping", permuted, "production-vs-mean"), prod, mn)
        pair(paste("mapping", permuted, "SSP-vs-exp-mean"), prod_ssp, 2^mn)
        pair(paste("mapping", permuted, "SSP-vs-exp-colMeans"), prod_ssp, 2^cm)
        check(paste("production collapse equals columnwise mean", permuted), strict(prod, mn))
        check(paste("SSP differences fully explained by exponentiation", permuted), strict(prod_ssp, 2^mn))
        cases[[as.character(permuted)]] <- list(values = values, labels = labels, entrez = entrez,
            input = se, encodings = encodings, duplicates = duplicates, colMeans = cm,
            columnwise_mean = mn, production = controls)
    }
    saveRDS(cases, file.path(out, "mapping-complete.rds"))
})
for (name in c("OSLO2EMIT0obj", "TCGABRCAobj")) {
    stage(paste0("expanded-", name), {
        obj <- env$current_example(name)
        fresh <- env$issue199_fresh(name)
        saveRDS(list(actual = fresh, reference = obj),
            file.path(out, paste0(name, "-expanded-unprojected-complete.rds")))
        compatible <- env$auto_ssbc_historical_view(fresh$res, obj$res)
        pair(paste(name, "expanded-complete"),
            list(data_input = fresh$data_input, res = compatible),
            list(data_input = obj$data_input, res = obj$res))
        saveRDS(fresh, file.path(out, paste0(name, "-fresh-unprojected.rds")))
        decisions(paste("expanded", name), compatible, obj$res)
        calls <- env$issue199_calls()[[name]]$C
        pair(paste(name, "expanded-reference-calls"), fresh$res$res_subtypes, calls)
        check(paste(name, "stored calls equal fixture"), strict(obj$res$res_subtypes, calls))
        methods <- setdiff(names(obj$res$res_subtypes), "entropy")
        check(paste(name, "stored method calls equal fixture"),
            strict(obj$res$res_subtypes[, methods], calls[, methods]))
        entropy <- apply(obj$res$res_subtypes[, methods, drop = FALSE], 1L, function(x) {
            x <- x[!is.na(x)]
            if (!length(x)) return(NA_real_)
            p <- as.numeric(table(x)) / length(x)
            -sum(p * log2(p))
        })
        check(paste(name, "original entropy reconstruction tolerance"),
            isTRUE(all.equal(unname(entropy), obj$res$res_subtypes$entropy, tolerance = 1e-12)))
        check(paste(name, "eight methods"), length(obj$res$results) == 8L)
        check(paste(name, "ssBC.v2 missing count"), sum(is.na(obj$res$res_subtypes$ssBC.v2)) ==
            if (name == "OSLO2EMIT0obj") 8L else 4L)
        vv <- cbind(actual = compatible$results$genefu.robust$score.ROR[["ROR-P (Subtype + Prolif)"]],
            reference = obj$res$results$genefu.robust$score.ROR[["ROR-P (Subtype + Prolif)"]])
        rownames(vv) <- obj$res$results$genefu.robust$score.ROR$patientID
        saveRDS(vv, file.path(out, paste0(name, "-genefu-robust-ROR-P.rds")))
        write.csv(vv, file.path(out, paste0(name, "-genefu-robust-ROR-P.csv")))
        cuts <- 100 * (c(-.25, .1) + .35) / .85
        json(list(cuts = cuts, minimum_reference_distance = min(abs(outer(vv[, "reference"], cuts, "-")), na.rm = TRUE)),
            file.path(out, paste0(name, "-ROR-cutoffs.json")))
    })
    stage(paste0("shiny-", name), {
        obj <- env$shiny199_object(name)
        raw <- identical(name, "TCGABRCAobj")
        prefix <- if (raw) "TCGABRCAsub" else "OSLO2EMIT0"
        path <- system.file("RshinyTest", paste0(prefix, if (raw) "_GEX_RawCounts.txt" else "_GEX_log2.FPKM.txt"),
            package = "BreastSubtypeR")
        reference <- as.matrix(SummarizedExperiment::assay(obj$se_obj))
        storage.mode(reference) <- "double"
        read <- as.matrix(utils::read.delim(path, row.names = 1L, check.names = FALSE))
        storage.mode(read) <- "double"
        parsed <- BreastSubtypeR:::.read_shiny_expression(path, "\t",
            colnames(reference), rownames(reference), raw)
        tokens <- as.matrix(utils::read.delim(path, row.names = 1L, check.names = FALSE, colClasses = "character"))
        decimal_dir <- file.path(out, paste0("decimal-", name))
        dir.create(decimal_dir)
        for (tag in c("reference", "read", "parsed")) writeBin(as.double(get(tag)),
            file.path(decimal_dir, paste0(tag, ".bin")), size = 8L, endian = "little")
        saveRDS(list(tokens = tokens, read = read, parsed = parsed, reference = reference),
            file.path(decimal_dir, "complete-arrays.rds"))
        py <- Sys.getenv("DIAG_PYTHON", Sys.which("python3"))
        stopifnot(nzchar(py))
        status <- system2(py, c(shQuote(python_script), shQuote(path), shQuote(decimal_dir)))
        stopifnot(status == 0L)
        nearest <- matrix(readBin(file.path(decimal_dir, "nearest.bin"), "double",
            n = length(reference), size = 8L, endian = "little"), nrow(reference), ncol(reference), byrow = TRUE,
            dimnames = dimnames(reference))
        pair(paste(name, "read.delim-vs-binary"), read, reference)
        pair(paste(name, "parser-vs-binary"), parsed, reference)
        pair(paste(name, "parser-vs-read.delim"), parsed, read)
        pair(paste(name, "independent-decimal-vs-binary"), nearest, reference)
        pair(paste(name, "parser-vs-independent-decimal"), parsed, nearest)
        parsed_se <- obj$se_obj
        SummarizedExperiment::assay(parsed_se) <- parsed
        binary_mapping <- env$issue199_rng(1992026, suppressMessages(Mapping(obj$se_obj,
            RawCounts = raw, method = "max", impute = TRUE, verbose = FALSE)))
        text_mapping <- env$issue199_rng(1992026, suppressMessages(Mapping(parsed_se,
            RawCounts = raw, method = "max", impute = TRUE, verbose = FALSE)))
        pair(paste(name, "Mapping-text-vs-binary-input"), text_mapping, binary_mapping)
        pair(paste(name, "Mapping-binary-vs-bundled"), binary_mapping, obj$data_input)
        for (component in c("se_NC", "se_SSP")) {
            pair(paste(name, component, "text-assay-vs-bundled"),
                SummarizedExperiment::assay(text_mapping[[component]]),
                SummarizedExperiment::assay(obj$data_input[[component]]))
        }
        inspect_server <- function() {
          original_mapping <- BreastSubtypeR::Mapping
          observed <- new.env(parent = emptyenv())
          testthat::local_mocked_bindings(Mapping = function(se_obj, ...) {
              observed$input <- SummarizedExperiment::assay(se_obj)
              original_mapping(se_obj = se_obj, ...)
          }, .package = "BreastSubtypeR")
          shiny::testServer(env$shiny199_server(), {
            session$setInputs(BSmethod = "AUTO Mode", k_subtypes = "5",
                is_raw_counts = if (raw) "raw" else "norm",
                hasClinical = FALSE, export_kind = "calls")
            if (raw) session$setInputs(
                GEX = env$shiny199_upload("TCGABRCAsub_GEX_RawCounts.txt"),
                clinic = env$shiny199_upload("TCGABRCAsub_clinical.txt"),
                anno = env$shiny199_upload("TCGABRCAsub_anno.txt"))
            else session$setInputs(load_example_step1 = 1)
            RNGkind("Mersenne-Twister", "Inversion", "Rejection")
            set.seed(1992026)
            session$setInputs(map = 1)
            pair(paste(name, "server-ingestion-vs-parser"), observed$input, parsed)
            pair(paste(name, "server-ingestion-vs-binary"), observed$input, reference)
            check(paste(name, "server map error absent"), is.null(reactive_files$error))
            check(paste(name, "server mapping current"), mapping_current())
            check(paste(name, "server raw flag"), identical(reactive_files$RawCounts, raw))
            pair(paste(name, "server-Mapping-vs-direct-text"), reactive_files$data_input, text_mapping)
            check(paste(name, "server AUTO ready"), auto_requirements()$ready)
            set.seed(20260831)
            session$setInputs(run = 1)
            check(paste(name, "server run error absent"), is.null(reactive_files$error))
            prediction <- current_prediction()
            check(paste(name, "prediction present"), !is.null(prediction))
            expected_calls <- data.frame(PatientID = rownames(obj$res$res_subtypes),
                obj$res$res_subtypes, check.names = FALSE)
            pair(paste(name, "prediction-calls"), prediction$calls, expected_calls)
            pair(paste(name, "prediction-plot"), prediction$plot, obj$res$res_subtypes)
            check(paste(name, "prediction five classes"), identical(prediction$want_4, FALSE))
            check(paste(name, "prediction full null"), is.null(prediction$full))
            export_path <- output$download
            check(paste(name, "download exists"), file.exists(export_path))
            file.copy(export_path, file.path(out, paste0(name, "-download.tsv")), overwrite = FALSE)
            export <- utils::read.delim(export_path, comment.char = "#", check.names = FALSE)
            pair(paste(name, "download-values"), export, expected_calls)
            fields <- setdiff(names(export), "entropy")
            check(paste(name, "download nonentropy fields as character"),
                all(vapply(fields, function(n) identical(as.character(export[[n]]),
                    as.character(expected_calls[[n]])), logical(1))))
            check(paste(name, "download entropy original tolerance"),
                isTRUE(all.equal(export$entropy, expected_calls$entropy, tolerance = 1e-14)))
            check(paste(name, "download columns"), identical(names(export), names(expected_calls)))
            check(paste(name, "download rows"), nrow(export) == ncol(obj$se_obj))
            headers <- readLines(export_path, n = 2L)
            check(paste(name, "download package header"), grepl("BreastSubtypeR", headers[1], fixed = TRUE))
            check(paste(name, "download class header"), grepl("Classes: 5", headers[2], fixed = TRUE))
            saveRDS(prediction, file.path(out, paste0(name, "-prediction-complete.rds")))
          })
        }
        inspect_server()
    })
}
json(checks, file.path(out, "checks.json"))
saveRDS(checks, file.path(out, "checks.rds"))
capture.output(sessionInfo(), file = file.path(out, "session-info-after.txt"))
json(list(stages_completed = all(vapply(stages, function(x) isTRUE(x$completed), logical(1))),
    stage_count = length(stages), checks = length(checks),
    false_checks = sum(!vapply(checks, function(x) x$passed, logical(1)))),
    file.path(out, "collection-summary.json"))
if (!all(vapply(stages, function(x) isTRUE(x$completed), logical(1)))) quit(status = 1L)
