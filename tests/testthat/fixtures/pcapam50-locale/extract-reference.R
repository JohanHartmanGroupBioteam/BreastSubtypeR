# Extract frozen expectations from the completed #232 investigation, not from
# current BreastSubtypeR. Not sourced by tests; takes an archive and output path.
# Rscript extract-reference.R /path/to/retained-investigation /new/reference.rds
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, !file.exists(args[2]))
archive <- args[1]
checksums <- c(`trace-reference-results.rds` = "297b5f56ee9fb07bc578cbfdae959e43",
    `permutation.rds` = "3bba1998e456dca6427db6df970ebf42",
    `sample-id-map.csv` = "355d12d0e42bac5b67a9e23c7ea0dca5")
for (f in names(checksums)) stopifnot(identical(
    unname(tools::md5sum(file.path(archive, f))), unname(checksums[f])))
records <- readRDS(file.path(archive, "trace-reference-results.rds"))
ids <- read.csv(file.path(archive, "sample-id-map.csv"), stringsAsFactors = FALSE)
extract_stage <- function(record, stage, canonical_ids) {
    first <- stage == 1L
    q <- record$trace[[if (first) "REF.makeCalls.PC1ihc" else "REF.makeCalls.v1PAM"]]
    result <- record[[paste0("ref", stage)]]
    correlation <- result$outList$distances
    # The original returned matrix has no dimnames. These labels come from
    # that same original output; distances are already positive correlations.
    dimnames(correlation) <- list(names(result$outList$predictions),
        colnames(result$outList$centroids))
    list(ordered_ids = if (first) q$df.pca1$PatientID else q$df.pam$PatientID,
        positive = if (first) q$ERP.pc1ihc$PatientID else q$ERP.pam$PatientID,
        negative = if (first) q$ERN.pc1ihc$PatientID else q$ERN.pam$PatientID,
        draw = q$i,
        selected = colnames(if (first) q$mbal.pc1ihc else q$mbal.pam),
        medians = q$mdns,
        calls = result$outList$predictions[canonical_ids],
        correlations = correlation[canonical_ids, , drop = FALSE])
}
fixture <- list(
    reference_commit = "8920f64bb7cc669c1cd1fe09370aaa6029877118",
    package_commit = "0d3e64d68a83319bb8bf711958304e855b371f05",
    seed = 118L, rng_kind = c("Mersenne-Twister", "Inversion", "Rejection"),
    ids = ids, permutation = readRDS(file.path(archive, "permutation.rds")),
    cases = list())
for (idkind in c("original", "mixed")) for (loc in c("C", "en_US.UTF-8")) {
    key <- paste(idkind, loc, sep = "__")
    cases <- lapply(c("original", "reverse", "permuted"), function(ordering) {
        record <- records[[paste(idkind, ordering, loc, sep = "__")]]
        list(first = extract_stage(record, 1L, ids[[idkind]]),
             second = extract_stage(record, 2L, ids[[idkind]]))
    })
    # All 12 original-reference captures independently support these four
    # compact expectations. No BSR result or current classifier is consulted.
    stopifnot(identical(cases[[1]], cases[[2]]), identical(cases[[1]], cases[[3]]))
    fixture$cases[[key]] <- cases[[1]]
}
q <- records[["original__original__C"]]$trace[["REF.makeCalls.PC1ihc"]]
fixture$clinical <- q$df.cln[match(ids$original, q$df.cln$PatientID),
    c("PatientID", "ER", "PR", "HER2", "TSIZE", "NODE"), drop = FALSE]
rownames(fixture$clinical) <- NULL
fixture$genes <- names(q$mdns)
fixture$er_order <- unique(q$df.cln$ER_status)
saveRDS(fixture, args[2], version = 2, compress = "xz")
cat("Extracted four expectations from 12 original-reference captures\n")
