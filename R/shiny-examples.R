# The UI, expression reader and Mapping use the same packaged example record.
.shiny_examples <- list(
    oslo = list(name = "OSLO2-EMIT0", prefix = "OSLO2EMIT0",
                expression = "OSLO2EMIT0_GEX_log2.FPKM.txt",
                units = "log2(protocol-adjusted FPKM + 1)", raw = FALSE,
                dimensions = c(569L, 102L)),
    tcga = list(name = "TCGA-BRCA", prefix = "TCGABRCAsub",
                expression = "TCGABRCAsub_GEX_RawCounts.txt",
                units = "recount3 G026 gene-level base-pair coverage sums", raw = TRUE,
                dimensions = c(5350L, 100L))
)

.read_shiny_example <- function(id, directory = system.file("RshinyTest", package = "BreastSubtypeR")) {
    if (!is.character(id) || length(id) != 1L || !id %in% names(.shiny_examples))
        stop("Select a supported packaged example.", call. = FALSE)
    spec <- .shiny_examples[[id]]
    files <- c(GEX = spec$expression, anno = paste0(spec$prefix, "_anno.txt"),
               clinic = paste0(spec$prefix, "_clinical.txt"))
    if (!nzchar(directory) || !all(file.exists(file.path(directory, files))))
        stop("Packaged example files are missing. Reinstall BreastSubtypeR or upload all three tables.",
             call. = FALSE)
    anno <- .read_shiny_metadata(file.path(directory, files["anno"]), "\t", "probe")
    clinic <- .read_shiny_metadata(file.path(directory, files["clinic"]), "\t", "PatientID")
    if (!all(c("probe", "ENTREZID") %in% names(anno)) || !"PatientID" %in% names(clinic))
        stop("Example tables require probe, ENTREZID and PatientID columns.", call. = FALSE)
    if (spec$raw && (!is.numeric(anno$Length) || any(!is.finite(anno$Length) | anno$Length <= 0)))
        stop("The raw example requires a positive, finite Length for every feature.", call. = FALSE)
    gex <- .read_shiny_expression(file.path(directory, files["GEX"]), "\t",
                                  clinic$PatientID, anno$probe, spec$raw)
    if (!identical(dim(gex), spec$dimensions) ||
        !identical(rownames(gex), anno$probe) || !identical(colnames(gex), clinic$PatientID))
        stop("Example dimensions or identifiers differ from the packaged cohort.", call. = FALSE)
    rownames(anno) <- anno$probe
    rownames(clinic) <- clinic$PatientID
    list(GEX = gex, anno = anno, clinic = clinic)
}
