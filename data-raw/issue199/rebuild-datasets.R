# Reconstruct #199 examples in a NEW directory; never writes data/ or installs software.
# Rscript --vanilla rebuild-datasets.R --archive=/retained/2026-10-26_expanded-inputs --historical-dir=/original/complete-objects --out=/new/run
# Optional --package-source=/checkout (defaults to this script's repository).
options(stringsAsFactors = FALSE, digits = 17)
arg <- commandArgs(TRUE)
if (!all(grepl("^--[^=]+=.+$", arg))) stop("Use --name=value arguments.")
key <- sub("^--([^=]+)=.*$", "\\1", arg)
opt <- setNames(sub("^--[^=]+=", "", arg), key)
if (anyDuplicated(key) || !all(c("archive", "historical-dir", "out") %in% key) ||
    any(!key %in% c("archive", "historical-dir", "out", "package-source", "verify-only"))) stop("Required: --archive, --historical-dir and --out.")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
here <- dirname(normalizePath(script, winslash = "/", mustWork = TRUE))
repo <- if ("package-source" %in% key) opt[["package-source"]] else file.path(here, "../..")
repo <- normalizePath(repo, winslash = "/", mustWork = TRUE)
archive <- normalizePath(opt[["archive"]], winslash = "/", mustWork = TRUE)
historical_dir <- normalizePath(opt[["historical-dir"]], winslash = "/", mustWork = TRUE)
out <- opt[["out"]]
if (file.exists(out) || !dir.exists(dirname(out))) stop("Use a NEW output directory with an existing parent.")
for (p in c("BreastSubtypeR", "SummarizedExperiment", "digest", "jsonlite"))
    if (!requireNamespace(p, quietly = TRUE)) stop("Missing dependency: ", p)
sha <- function(p) {
    con <- file(p, "rb"); on.exit(close(con))
    bytes <- readBin(con, "raw", n = file.info(p)$size)
    digest::digest(bytes, algo = "sha256", serialize = FALSE)
}
# Verify the retained archive manifest itself, every external plan/source/model,
# and both approved candidates/historical fixtures before any reconstruction.
resolve <- function(p) {
    if (startsWith(p, "archive/")) return(file.path(archive, substring(p, 9L)))
    if (startsWith(p, "package/")) return(file.path(repo, substring(p, 9L)))
    if (startsWith(p, "historical/")) return(file.path(historical_dir, substring(p, 12L)))
    stop("Unexpected manifest path: ", p)
}
checked <- 0L
for (manifest_name in c("source-inputs.tsv", "approved-artifacts.tsv")) {
    m <- read.delim(file.path(here, "manifests", manifest_name), check.names = FALSE)
    for (i in seq_len(nrow(m))) {
        p <- resolve(m$path[i])
        if (!file.exists(p) || sha(p) != m$sha256[i])
            stop("Pinned source/artifact hash mismatch: ", m$path[i])
        if ("bytes" %in% names(m) && file.info(p)$size != m$bytes[i])
            stop("Pinned source size mismatch: ", m$path[i])
        checked <- checked + 1L
    }
}
cat("PASS: ", checked, " pinned source/artifact manifest entries.\n", sep = "")
if ("verify-only" %in% key) {
    if (!identical(opt[["verify-only"]], "true")) stop("--verify-only must be true.")
    quit(status = 0)
}
# Pin supported implementations; keep the approved generation baseline explicit.
implementations <- c(
    cb7a7c70e0f94e1690d15ac2ec2b464a40934f36 = "implementation-function-hashes.tsv",
    f878f58c1c12ca6f1cd5bf4489667935bf25969b = "implementation-function-hashes-f878f58c.tsv")
matches <- vapply(implementations, function(manifest) {
    baseline <- read.delim(file.path(here, "manifests", manifest))
    all(vapply(seq_len(nrow(baseline)), function(i) {
        f <- tryCatch(get(baseline$name[i], asNamespace("BreastSubtypeR")),
            error = function(e) NULL)
        is.function(f) &&
            identical(digest::digest(body(f), algo = "sha256", serializeVersion = 3),
                baseline$body_sha256[i]) &&
            identical(digest::digest(formals(f), algo = "sha256", serializeVersion = 3),
                baseline$formals_sha256[i])
    }, logical(1)))
}, logical(1))
if (sum(matches) != 1L) stop("Installed implementation does not match a supported pinned source.")
implementation_commit <- names(implementations)[matches]
stopifnot(dir.create(out))
out <- normalizePath(out, winslash = "/", mustWork = TRUE)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
run <- function(name, args) {
    log <- file.path(out, paste0(name, ".log"))
    status <- system2(rscript, c("--vanilla", shQuote(file.path(here, paste0(name, ".R"))),
        vapply(args, shQuote, "")), stdout = log, stderr = log)
    if (status != 0L) stop(name, " failed; see ", log)
}
historical <- historical_dir
run("reconstruct_oslo", c(
    paste0("--baseline=", file.path(historical, "OSLO2EMIT0obj.RData")),
    paste0("--expression=", file.path(archive, "source-data/oslo/OSLO2EMIT0.103.genematrix_noNeg.Rdata")),
    paste0("--annotation=", file.path(archive, "source-data/oslo/Gene.ID.ann.Rdata")),
    paste0("--model-dir=", file.path(repo, "data")),
    paste0("--discrepancies=", file.path(archive, "provenance/reference-tables/OSLO_exact_discrepancies.tsv")),
    paste0("--feature-plan=", file.path(archive, "provenance/reference-tables/OSLO_planned_features.tsv")),
    paste0("--out=", file.path(out, "oslo-input"))))
run("reconstruct_tcga", c(
    paste0("--baseline=", file.path(historical, "TCGABRCAobj.RData")),
    paste0("--archive=", archive), paste0("--model-dir=", file.path(repo, "data")),
    paste0("--out=", file.path(out, "tcga-input"))))
suppressPackageStartupMessages(library(BreastSubtypeR))
settings <- list(source_commit = implementation_commit,
    approved_generation_commit = "cb7a7c70e0f94e1690d15ac2ec2b464a40934f36",
    input_model_commit = "6b2a4ab16798493f87713a91cee4674b2ae4df8e",
    RNGkind = c("Mersenne-Twister", "Inversion", "Rejection"),
    mapping_seed = 1992026L, prediction_seed = 20260831L,
    Mapping = list(method = "max", impute = TRUE, verbose = FALSE),
    BS_Multi = list(methods = "AUTO", Subtype = FALSE, hasClinical = FALSE, er_subsample = "target"))
spec <- list(OSLO2EMIT0obj = c("oslo-input", "OSLO2EMIT0_stage1_input.rds"),
    TCGABRCAobj = c("tcga-input", "TCGABRCA_stage2_input.rds"))
checks <- data.frame(dataset=character(), check=character(), status=character())
outputs <- list()
for (nm in names(spec)) {
    se <- readRDS(do.call(file.path, as.list(c(out, spec[[nm]]))))
    RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(settings$mapping_seed)
    mapped <- Mapping(se, RawCounts = nm == "TCGABRCAobj", method = "max",
        impute = TRUE, verbose = FALSE)
    RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(settings$prediction_seed)
    result <- BS_Multi(mapped, methods = "AUTO", Subtype = FALSE,
        hasClinical = FALSE, er_subsample = "target")
    obj <- list(se_obj = se, data_input = mapped, res = result)
    e <- new.env(); load(file.path(repo, "data", paste0(nm, ".RData")), envir = e)
    for (component in names(obj)) {
        ok <- identical(obj[[component]], e[[nm]][[component]])
        checks[nrow(checks)+1L,] <- list(nm, paste0(component, "_canonical_identity"), if(ok) "PASS" else "FAIL")
    }
    if (any(checks$status == "FAIL")) stop("Reproduction differs from approved artifact: ", nm)
    e <- new.env(parent = emptyenv()); e[[nm]] <- obj
    target <- file.path(out, paste0(nm, ".RData"))
    save(list = nm, envir = e, file = target, compress = "xz", version = 3)
    outputs[[nm]] <- data.frame(dataset = nm, file = basename(target),
        sha256 = sha(target), bytes = file.info(target)$size)
}
write.table(checks, file.path(out, "canonical-identity.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(do.call(rbind, outputs), file.path(out, "outputs-sha256.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
jsonlite::write_json(settings, file.path(out, "settings.json"), pretty = TRUE, auto_unbox = TRUE)
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
cat("PASS: six complete component identity checks. File hashes identify this run; byte identity across R environments is not promised.\n")
