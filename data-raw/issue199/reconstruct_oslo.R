# Issue #199 Stage 1: preserve historical OSLO rows; append adjusted-source genes.
# This script creates an INPUT-ONLY SummarizedExperiment in a NEW output directory.
# No Mapping, classifier, package mutation, downloads, or historical result regeneration.
# Usage (quote each --name=value argument when paths contain spaces):
# Rscript --vanilla reconstruct_oslo.R --baseline=FILE --expression=FILE \
#   --annotation=FILE --model-dir=DIR --discrepancies=FILE --feature-plan=FILE --out=NEW_DIR
# --library=DIR may select an isolated R library; runtime libs can also be set by R_LIBS.
options(stringsAsFactors = FALSE, digits = 17)
args <- commandArgs(trailingOnly = TRUE)
if (any(!grepl("^--[^=]+=.+$", args))) stop("Arguments must use --name=value.")
keys <- sub("^--([^=]+)=.*$", "\\1", args)
if (anyDuplicated(keys)) stop("Duplicate argument.")
opt <- setNames(sub("^--[^=]+=", "", args), keys)
required <- c("baseline", "expression", "annotation", "model-dir",
              "discrepancies", "feature-plan", "out")
if (!all(required %in% keys) || any(!keys %in% c(required, "library")))
    stop("Required: ", paste(required, collapse = ", "))
if ("library" %in% keys) .libPaths(c(opt[["library"]], .libPaths()))
for (p in c("SummarizedExperiment", "S4Vectors", "digest", "jsonlite"))
    if (!requireNamespace(p, quietly = TRUE)) stop("Required package unavailable: ", p)
suppressPackageStartupMessages(library(SummarizedExperiment))
sha <- function(p) {
    # Read bytes directly: Windows file.access() can report a false negative
    # for readable OneDrive sources. The bytes are still SHA256-verified.
    con <- file(p, "rb"); on.exit(close(con))
    n <- file.info(p)$size
    bytes <- readBin(con, "raw", n = n)
    if (length(bytes) != n) stop("Incomplete source read: ", p)
    digest::digest(bytes, algo = "sha256", serialize = FALSE)
}
load_one <- function(p) {
    e <- new.env(parent = emptyenv()); nm <- load(p, envir = e)
    if (length(nm) != 1L) stop("Expected one object in ", p)
    e[[nm]]
}
commit <- "6b2a4ab16798493f87713a91cee4674b2ae4df8e"
expected <- c(
    baseline = "43393c2bf06f0ea44b262d5525b75154ab99162fac5de898433556b6805450ec",
    expression = "4311bf87dc81c01ae37069a6acac0d25d1e1a087350c324ee2a3f68286459b94",
    annotation = "6f10528db362160313c2f1345770c1333f5bab30a19e4f68b0ceaf0ab595d678",
    discrepancies = "13930cad46ab40b5c57b8475bf63cbeef2b5c0a890a1f1135cc95f43fb47e65f",
    feature_plan = "fb4ad21a557ef96b6ccb52b3a0755f7529efe177675178e7ca3b1132e1f95ebd",
    BreastSubtypeRobj = "dac97591d52816121c6d06ab2c7cbca16385f205616c9fc8fa73ae9793c03096",
    AIMSmodel = "52b1190a09535d1df0ff09c06a4e0a180374da81649faaaf00194fd178b07b92",
    sspbc.models = "9cf0c03578c4ed7fc5a0b8b2c8a1eb683778df039e6b02d53ddc5d65d500634c"
)
paths <- c(baseline = opt[["baseline"]], expression = opt[["expression"]],
    annotation = opt[["annotation"]], discrepancies = opt[["discrepancies"]],
    feature_plan = opt[["feature-plan"]])
for (nm in c("BreastSubtypeRobj", "AIMSmodel", "sspbc.models"))
    paths[nm] <- file.path(opt[["model-dir"]], paste0(nm, ".RData"))
paths <- vapply(paths, normalizePath, "", winslash = "/", mustWork = TRUE)
actual <- vapply(paths, sha, "")
if (!identical(unname(actual), unname(expected[names(paths)])))
    stop("An input SHA256 does not match the pinned manifest.")
if (file.exists(opt[["out"]]) || dir.exists(opt[["out"]]))
    stop("Output path exists; refusing to overwrite anything.")
if (!dir.exists(dirname(opt[["out"]]))) stop("Output parent must already exist.")
out_parent <- normalizePath(dirname(opt[["out"]]), winslash = "/", mustWork = TRUE)
out <- file.path(out_parent, basename(opt[["out"]]))
# A new output directory may not contain any source input.
if (any(startsWith(tolower(paths), paste0(tolower(out), "/"))))
    stop("Output directory contains a scientific input.")
checks <- data.frame(check = character(), status = character(), detail = character())
check <- function(name, condition, detail) {
    ok <- isTRUE(condition)
    checks[nrow(checks) + 1L, ] <<- list(name, if (ok) "PASS" else "FAIL", detail)
    if (!ok) stop(name, ": ", detail)
}
old_object <- load_one(paths[["baseline"]])
baseline <- old_object$se_obj
check("baseline_schema", inherits(baseline, "SummarizedExperiment") &&
      identical(dim(baseline), c(356L, 102L)) &&
      identical(assayNames(baseline), "counts"),
      "Original 356 x 102 SummarizedExperiment with counts assay.")
old_values <- assay(baseline)
old_rd <- rowData(baseline)
old_cd <- colData(baseline)
check("baseline_ids", !anyDuplicated(rownames(baseline)) &&
      !anyDuplicated(colnames(baseline)) &&
      !anyDuplicated(as.character(old_rd$ENTREZID)),
      "Original feature, sample and Entrez identifiers are unique.")
source_values <- load_one(paths[["expression"]])
ann <- load_one(paths[["annotation"]])
check("adjusted_source_schema", identical(dim(source_values), c(19675L, 103L)) &&
      !anyDuplicated(rownames(source_values)) &&
      !anyDuplicated(colnames(source_values)) &&
      all(is.finite(source_values)) && all(source_values >= 0),
      "Verified adjusted linear FPKM source; no unadjusted source is read.")
check("sample_identity", identical(
    colnames(source_values)[colnames(source_values) != "OSLO2EMIT0.073"],
    colnames(baseline)), "Historical 102 samples in exact source-relative order.")
myo <- which(ann$HGNC == "MYO10")
check("MYO10_correction", length(myo) == 1L &&
      identical(ann$Gene.ID[myo], "ENSG00000145555.14") &&
      identical(ann$EntrezGene[myo], "NA") &&
      identical(as.character(old_rd$ENTREZID[
          old_rd$probe == "ENSG00000145555.14"]), "4651"),
      "Public annotation literal NA is the documented MYO10 correction to 4651.")
ann$EntrezGene[myo] <- "4651"
check("annotation_source_alignment", !anyDuplicated(ann$Gene.ID) &&
      all(rownames(source_values) %in% ann$Gene.ID),
      "Source features have unambiguous annotation rows.")
ann <- ann[match(rownames(source_values), ann$Gene.ID), , drop = FALSE]
ref <- load_one(paths[["BreastSubtypeRobj"]])
aims <- load_one(paths[["AIMSmodel"]])
ssp <- load_one(paths[["sspbc.models"]])
models <- c(list(AIMS = aims), ssp[c("ssp.subtype", "ssp.pam50", "ssp.ror")])
rule_tables <- lapply(names(models), function(nm) {
    pairs <- strsplit(models[[nm]]$all.pairs, "<", fixed = TRUE)
    if (!all(lengths(pairs) == 2L)) stop("Invalid pair specification.")
    data.frame(model = nm, rule_index = seq_along(pairs),
        rule = models[[nm]]$all.pairs,
        gene_a = sub("^e", "", vapply(pairs, "[[", "", 1L)),
        gene_b = sub("^e", "", vapply(pairs, "[[", "", 2L)))
})
rules <- do.call(rbind, rule_tables)
req <- c(list(PAM50 = as.character(ref$genes.sig50$EntrezGene.ID)),
    setNames(lapply(rule_tables, function(t) unique(c(t$gene_a, t$gene_b))),
             names(models)))
union_ids <- sort(unique(unlist(req, use.names = FALSE)))
subtype_ids <- unique(unlist(req[setdiff(names(req), "ssp.ror")], use.names = FALSE))
check("model_union", length(union_ids) == 569L && length(subtype_ids) == 356L &&
      setequal(subtype_ids, as.character(ref$genes.signature$EntrezGene.ID)) &&
      length(req[["ssp.ror"]]) == 296L && sum(rules$model == "ssp.ror") == 361L,
      "PAM50 + AIMS + SSP subtype/PAM50/ROR: 569 genes; ROR 296 genes/361 rules.")
check("historical_feature_union", setequal(as.character(old_rd$ENTREZID), subtype_ids),
      "Historical 356 genes equal the verified subtype union.")
new_ids <- setdiff(union_ids, as.character(old_rd$ENTREZID))
extra <- ann[which(ann$EntrezGene %in% new_ids), , drop = FALSE]
check("additional_gene_identity", length(new_ids) == 213L && nrow(extra) == 213L &&
      !anyDuplicated(extra$EntrezGene) && setequal(extra$EntrezGene, new_ids),
      "Exactly one adjusted-source feature for each of 213 additional Entrez IDs.")
source_old_ann <- ann[match(rownames(baseline), ann$Gene.ID), , drop = FALSE]
check("original_annotation_source_check",
      identical(as.character(old_rd$probe), source_old_ann$Gene.ID) &&
      identical(as.character(old_rd$SYMBOL), source_old_ann$HGNC) &&
      identical(as.character(old_rd$ENTREZID), as.character(source_old_ann$EntrezGene)),
      "Original feature annotations agree with source after documented MYO10 correction.")
features <- c(rownames(baseline), extra$Gene.ID)
entrez <- c(as.character(old_rd$ENTREZID), as.character(extra$EntrezGene))
plan <- read.delim(paths[["feature_plan"]], check.names = FALSE, colClasses = "character")
check("retained_feature_plan", identical(features, plan$probe) &&
      identical(entrez, plan$ENTREZID),
      "Original 356 rows followed by 213 additions in adjusted-source order.")

# Only APPENDED rows are transformed for the candidate. Historical assay values
# below are copied from the verified baseline, never recomputed or replaced.
extra_values <- log2(source_values[extra$Gene.ID, colnames(baseline), drop = FALSE] + 1)
new_rd <- S4Vectors::DataFrame(probe = extra$Gene.ID, SYMBOL = extra$HGNC,
    ENTREZID = as.character(extra$EntrezGene), row.names = extra$Gene.ID)
extra_se <- SummarizedExperiment(assays = setNames(list(extra_values), assayNames(baseline)),
    rowData = new_rd, colData = old_cd, metadata = metadata(baseline))
candidate <- rbind(baseline, extra_se)
metadata(candidate) <- metadata(baseline)
check("candidate_dimensions_and_ids", identical(dim(candidate), c(569L, 102L)) &&
      !anyDuplicated(as.character(rowData(candidate)$ENTREZID)) &&
      setequal(as.character(rowData(candidate)$ENTREZID), union_ids),
      "Input-only candidate: exactly 569 unique model genes x 102 samples.")
check("historical_values_exact", identical(assay(candidate)[seq_len(356L), , drop = FALSE],
      old_values), "All 36,312 historical expression values and dimnames are identical.")
check("sample_metadata_exact", identical(colnames(candidate), colnames(baseline)) &&
      identical(colData(candidate), old_cd), "Sample identifiers, order and colData unchanged.")
check("feature_metadata_exact", identical(rowData(candidate)[seq_len(356L), , drop = FALSE],
      old_rd), "All original rowData rows and column types unchanged.")
check("finite_expression", all(is.finite(assay(candidate))) &&
      !anyNA(assay(candidate)), "All 58,038 expression values are finite and nonmissing.")
check("new_annotation_complete", !anyNA(new_rd) &&
      all(nzchar(as.character(new_rd$SYMBOL))) &&
      all(grepl("^[0-9]+$", as.character(new_rd$ENTREZID))),
      "Added probe, symbol and Entrez annotations are complete.")

# Recompute original-source values ONLY for a validation comparison.
# This matrix never supplies candidate historical expression values.
recomputed <- log2(source_values[rownames(baseline), colnames(baseline), drop = FALSE] + 1)
idx <- which(recomputed != old_values, arr.ind = TRUE)
old_discrepancies <- read.delim(paths[["discrepancies"]], check.names = FALSE,
    colClasses = "character")
observed <- data.frame(probe = rownames(baseline)[idx[, 1L]],
    sample_id = colnames(baseline)[idx[, 2L]],
    packaged = sprintf("%.17g", old_values[idx]),
    recomputed = sprintf("%.17g", recomputed[idx]),
    difference = sprintf("%.17g", recomputed[idx] - old_values[idx]),
    packaged_hex = sprintf("%a", old_values[idx]),
    recomputed_hex = sprintf("%a", recomputed[idx]))
key <- function(d) paste(d$probe, d$sample_id, sep = "\t")
check("historical_source_discrepancy_coordinates", nrow(observed) == 51L &&
      identical(key(observed), key(old_discrepancies)),
      "The same 51 documented original-source discrepancies, in matrix order.")
check("historical_source_discrepancy_values",
      identical(observed$packaged_hex, old_discrepancies$packaged_hex) &&
      identical(observed$recomputed_hex, old_discrepancies$recomputed_hex),
      "Both hexadecimal doubles match the retained discrepancy evidence.")
# For positive finite doubles, exponent/mantissa unsigned bit patterns are ordered.
next_positive_double_raw <- function(v) {
    b <- as.integer(writeBin(v, raw(), size = 8L, endian = "little"))
    for (i in seq_along(b)) {
        if (b[i] < 255L) { b[i] <- b[i] + 1L; break }
        b[i] <- 0L
    }
    as.raw(b)
}
one_ulp <- vapply(seq_len(nrow(idx)), function(i) {
    a <- old_values[idx[i, 1L], idx[i, 2L]]
    b <- recomputed[idx[i, 1L], idx[i, 2L]]
    a > 0 && b > 0 && identical(next_positive_double_raw(min(a, b)),
        writeBin(max(a, b), raw(), size = 8L, endian = "little"))
}, logical(1))
observed$ulp_distance <- ifelse(one_ulp, 1L, NA_integer_)
check("historical_source_one_ULP", all(one_ulp),
      "Each recorded discrepancy is exactly one adjacent positive-double step.")
provenance <- list(
    schema_version = 1L, task = "Issue #199 Stage 1",
    object_kind = "input-only SummarizedExperiment",
    expression_units = "log2(protocol-adjusted StringTie FPKM + 1)",
    historical_rows = 356L, appended_rows = 213L, samples = 102L,
    preservation = "Original 356 assay rows, rowData and colData retained exactly.",
    appended_transformation = "log2(adjusted source FPKM + 1), applied only to appended rows",
    feature_order = "Original 356 feature order, then 213 additions in source row order",
    sample_order = "Original 102 sample IDs and order; source OSLO2EMIT0.073 excluded",
    annotation_correction = list(probe = "ENSG00000145555.14", symbol = "MYO10",
        source_EntrezGene = "NA (literal token)", corrected_EntrezGene = "4651"),
    model_source_commit = commit, package_source_version = "1.5.12",
    source_dataset = "Mendeley Data 10.17632/yzxtxn4nmd.4",
    source_file_id = "216a4cc7-e688-4c72-9513-2ad081e2c801",
    source_folder = "StringTie FPKM Gene Data LibProtocol adjusted",
    source_sha256 = as.list(actual),
    historical_source_discrepancies = list(count = 51L, each_ULP = 1L,
        historical_values_replaced = FALSE),
    limitations = c(
        "Complete source-gene coverage does not establish native SSP-ROR compatibility.",
        "SSPBC development used unadjusted StringTie FPKM; adjusted-input prediction equivalence is unestablished.",
        "Model-aware Mapping remains #131; normalization sensitivity remains #208.",
        "No predictions, mapped inputs or historical result lists are attached.")
)
metadata(candidate)$issue199_reconstruction <- provenance
check("SE_structure", identical(class(candidate), class(baseline)) &&
      identical(assayNames(candidate), assayNames(baseline)) &&
      identical(typeof(assay(candidate)), typeof(old_values)) &&
      isTRUE(validObject(candidate, test = TRUE)),
      "Original SummarizedExperiment class, assay name/storage preserved; object valid.")

if (!dir.create(out, recursive = FALSE)) stop("Could not create a new output directory.")
tab <- function(x, name) write.table(x, file.path(out, name), sep = "\t",
    quote = FALSE, row.names = FALSE, na = "NA")
json <- function(x, name) jsonlite::write_json(x, file.path(out, name), pretty = TRUE,
    auto_unbox = TRUE, null = "null", na = "null", digits = NA)
output_file <- file.path(out, "OSLO2EMIT0_stage1_input.rds")
saveRDS(candidate, output_file, compress = "xz", version = 3)
check("serialization_roundtrip", identical(readRDS(output_file), candidate),
      "Saved RDS reloads identically; it contains only the candidate SE.")
feature_manifest <- data.frame(output_position = seq_along(features), probe = features,
    SYMBOL = as.character(rowData(candidate)$SYMBOL), ENTREZID = entrez,
    source_position = match(features, rownames(source_values)),
    historical_feature = seq_along(features) <= 356L,
    transformation = ifelse(seq_along(features) <= 356L,
        "preserved_historical_log2_adjusted_FPKM_plus_1", "log2_adjusted_FPKM_plus_1"),
    annotation_correction = ifelse(features == "ENSG00000145555.14", "MYO10_to_4651", "none"),
    models = vapply(entrez, function(id) paste(names(req)[vapply(req,
        function(ids) id %in% ids, logical(1))], collapse = ";"), ""))
sample_manifest <- cbind(data.frame(output_position = seq_len(ncol(candidate)),
    sample_id = colnames(candidate), source_position = match(colnames(candidate),
        colnames(source_values))), as.data.frame(colData(candidate)))
req_manifest <- do.call(rbind, lapply(names(req), function(nm)
    data.frame(model = nm, ENTREZID = req[[nm]])))
coverage <- do.call(rbind, lapply(names(req), function(nm) {
    r <- rules[rules$model == nm, , drop = FALSE]
    data.frame(model = nm, genes_required = length(req[[nm]]),
        genes_present = sum(req[[nm]] %in% entrez),
        rules_required = if (nrow(r)) nrow(r) else NA_integer_,
        rules_structurally_present = if (nrow(r))
            sum(r$gene_a %in% entrez & r$gene_b %in% entrez) else NA_integer_)
}))
tab(feature_manifest, "feature_manifest.tsv")
tab(sample_manifest, "sample_manifest.tsv")
tab(req_manifest, "model_requirements.tsv")
tab(rules, "model_rules.tsv")
tab(coverage, "source_coverage.tsv")
tab(observed, "historical_source_discrepancies.tsv")
json(provenance, "provenance.json")
source_manifest <- data.frame(role = names(paths), path = unname(paths),
    bytes = unname(file.info(paths)$size), sha256 = unname(actual))
tab(source_manifest, "input_sha256.tsv")
json(list(input = source_manifest, candidate = list(path = basename(output_file),
    sha256 = sha(output_file), bytes = file.info(output_file)$size,
    dimensions = dim(candidate), in_memory_bytes = as.numeric(object.size(candidate)))),
    "manifest.json")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
ip <- as.data.frame(installed.packages()[, c("Package", "Version", "LibPath", "Built"),
    drop = FALSE], stringsAsFactors = FALSE)
tab(ip, "package_versions.tsv")
# Recheck all external inputs after generation.
check("inputs_unchanged", identical(vapply(paths, sha, ""), actual),
      "All consumed external scientific/model inputs retain their original hashes.")
tab(checks, "construction_checks.tsv")
cat("PASS: ", nrow(checks), " construction checks.\n", sep = "")
cat("Candidate: ", output_file, "\nSHA256: ", sha(output_file), "\n", sep = "")

