#' Describe model-specific coverage of processed expression inputs
#'
#' @param data_input The named list returned by [Mapping()], containing
#'   `se_NC` and `se_SSP` `SummarizedExperiment` objects.
#'
#' @return A data.frame with rows in this order: `NC-based (PAM50)`,
#'   `AIMS`, `SSPBC (ssp.subtype)`, and `SSPBC (ssp.pam50)`. The `model`
#'   column is character; all count columns are integer:
#'   \describe{
#'     \item{genes_present, genes_total, genes_absent}{Unique model genes
#'       represented in the processed input, included in the bundled model,
#'       and absent from the input, respectively.}
#'     \item{rules_evaluable_min, rules_evaluable_max}{Smallest and largest
#'       number of evaluable gene-pair rules across samples.}
#'     \item{rules_total}{Number of gene-pair rules in the bundled model.}
#'   }
#'   All three rule columns are `NA` for `NC-based (PAM50)`. Malformed or
#'   unavailable processed inputs produce an error; an input with no gene rows
#'   reports zero represented genes and zero evaluable rules.
#'
#' @details
#' This is descriptive reporting for four input models, independently of
#' which methods AUTO selects or a user runs. It does not run a classifier,
#' change predictions, assess calibration or prediction reliability, or
#' impose coverage thresholds or method-exclusion policies.
#'
#' NC-based coverage concerns the shared 50-gene PAM50 subtype reference.
#' It does not assess calibration feasibility or coverage of every ROR component.
#'
#' Gene representation and rule evaluability are different properties.
#' Gene and rule denominators come from the bundled references and fitted
#' model rule sets. Gene presence is structural: a represented gene with
#' zero or missing expression is distinct from an absent gene. Genes are
#' matched by their classifier-specific identifiers, including the native
#' SSPBC Entrez-prefix translation. Duplicate identifiers are counted once;
#' SSP duplicate expression rows follow the native per-sample selection.
#'
#' A rule is evaluable when the native strict `gene A < gene B` comparison
#' yields either `TRUE` or `FALSE`. Equal values, including zero versus zero,
#' give an evaluable `FALSE`. Missing genes and `NA`/`NaN` comparisons are
#' unavailable; infinities follow the native comparison semantics. Counts
#' may therefore differ across samples. No low-expression cutoff is used.
#'
#' These diagnostics describe the processed matrices in `data_input`.
#' [Mapping()] may already have performed imputation, so the report cannot
#' distinguish originally observed values from imputed values. It neither
#' preprocesses the expression again nor changes inputs, models or RNG state.
#'
#' @examples
#' data("OSLO2EMIT0obj")
#' data_input <- Mapping(OSLO2EMIT0obj$se_obj[, 1:3], verbose = FALSE)
#' InputCoverage(data_input)
#'
#' @md
#' @export
InputCoverage <- function(data_input) {
    required <- c("se_NC", "se_SSP")
    if (!is.list(data_input) ||
        any(vapply(required, function(nm) {
            sum(names(data_input) == nm, na.rm = TRUE) != 1L
        }, logical(1)))) {
        stop("InputCoverage requires a named list containing se_NC and se_SSP.",
            call. = FALSE)
    }
    nc <- .input_coverage_matrix(data_input$se_NC, "se_NC")
    ssp <- .input_coverage_matrix(data_input$se_SSP, "se_SSP")
    if (ncol(nc) != ncol(ssp) || !identical(colnames(nc), colnames(ssp))) {
        stop("InputCoverage requires matching sample identifiers and order in se_NC and se_SSP.",
            call. = FALSE)
    }

    models <- new.env(parent = emptyenv())
    utils::data(list = c("BreastSubtypeRobj", "AIMSmodel", "sspbc.models",
        "Gene.ID.ann"), envir = models, package = "BreastSubtypeR")
    reference <- models$BreastSubtypeRobj
    nc_genes <- unique(as.character(reference$genes.sig50$Symbol))
    nc_present <- as.integer(sum(nc_genes %in% rownames(nc)))
    coverage <- data.frame(
        model = "NC-based (PAM50)",
        genes_present = nc_present,
        genes_total = length(nc_genes),
        genes_absent = length(nc_genes) - nc_present,
        rules_evaluable_min = NA_integer_,
        rules_evaluable_max = NA_integer_,
        rules_total = NA_integer_,
        stringsAsFactors = FALSE
    )

    # BS_AIMS first restricts rows to its signature's unprefixed Entrez IDs.
    aims_genes <- as.character(reference$genes.signature$EntrezGene.ID[
        which(reference$genes.signature$AIMS == "Yes")])
    aims_rows <- rownames(ssp) %in% aims_genes
    coverage <- rbind(coverage, .input_coverage_pairs(
        ssp[aims_rows, , drop = FALSE], rownames(ssp)[aims_rows],
        models$AIMSmodel, "AIMS", .smaller_AIMS))

    ssp_ids <- translate_id2entrez(rownames(ssp), models$Gene.ID.ann,
        id.type = "EntrezGene", e = TRUE)
    for (key in c("ssp.subtype", "ssp.pam50")) {
        coverage <- rbind(coverage, .input_coverage_pairs(
            ssp, ssp_ids, models$sspbc.models[[key]],
            paste0("SSPBC (", key, ")"), .smaller))
    }
    rownames(coverage) <- NULL
    coverage
}

.input_coverage_matrix <- function(se, label) {
    fail <- function(reason) {
        stop("InputCoverage: ", label, " ", reason, call. = FALSE)
    }
    if (!inherits(se, "SummarizedExperiment")) {
        fail("must be a SummarizedExperiment object.")
    }
    if (length(SummarizedExperiment::assays(se)) == 0L) {
        fail("must contain a numeric expression assay.")
    }
    x <- as.matrix(SummarizedExperiment::assay(se))
    if (length(dim(x)) != 2L || (length(x) > 0L && !is.numeric(x))) {
        fail("must contain a numeric expression matrix.")
    }
    if (ncol(x) == 0L) fail("must contain at least one sample.")
    ids <- rownames(x)
    if (nrow(x) > 0L && (is.null(ids) || anyNA(ids) || any(!nzchar(ids)))) {
        fail("must have nonmissing, nonempty gene identifiers as row names.")
    }
    x
}

.input_coverage_pairs <- function(x, ids, model, label, compare) {
    pairs <- strsplit(model$all.pairs, "<", fixed = TRUE)
    genes <- unique(unlist(pairs, use.names = FALSE))
    keep <- ids %in% genes
    native <- .removeDuplicatedEntrezPerPatients(
        x[keep, , drop = FALSE], as.character(ids[keep]))
    present <- as.integer(sum(genes %in% native$EntrezID))
    evaluable <- integer(ncol(x))
    for (pair in pairs) {
        rows <- match(pair, native$EntrezID)
        if (!anyNA(rows)) {
            evaluable <- evaluable + as.integer(!is.na(compare(
                native$dataset[rows[1], ], native$dataset[rows[2], ])))
        }
    }
    data.frame(
        model = label,
        genes_present = present,
        genes_total = length(genes),
        genes_absent = length(genes) - present,
        rules_evaluable_min = min(evaluable),
        rules_evaluable_max = max(evaluable),
        rules_total = length(pairs),
        stringsAsFactors = FALSE
    )
}

.message_input_coverage <- function(coverage) {
    lines <- character(nrow(coverage))
    for (i in seq_len(nrow(coverage))) {
        row <- coverage[i, ]
        if (is.na(row$rules_total)) {
            lines[i] <- sprintf(
                "%s: %d/%d reference genes represented; %d absent.",
                row$model, row$genes_present, row$genes_total, row$genes_absent)
        } else {
            varies <- row$rules_evaluable_min != row$rules_evaluable_max
            count <- if (varies) {
                paste(row$rules_evaluable_min, row$rules_evaluable_max, sep = "-")
            } else as.character(row$rules_evaluable_min)
            lines[i] <- sprintf(
                "%s: %d/%d model genes represented; %s/%d rules evaluable%s.",
                row$model, row$genes_present, row$genes_total,
                count, row$rules_total, if (varies) " per sample" else "")
        }
    }
    message(paste(lines, collapse = "\n"))
}
