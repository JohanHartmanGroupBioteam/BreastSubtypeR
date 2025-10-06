#' @import stringr
#' @import Biobase
#' @importFrom e1071 naiveBayes
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment assay
#' @importFrom data.table data.table
#' @importFrom data.table set
#' @importFrom rlang call2
#' @importFrom rlang dots_list
#' @importFrom utils data installed.packages read.delim write.table
#' @importFrom graphics barplot mtext par
#' @importFrom grDevices dev.off pdf
#' @importFrom stats prcomp cor cor.test dist quantile median na.omit
#' @importFrom methods is
#' 
NULL

#' BreastSubtypeR: A Unified R/Bioconductor Package for Intrinsic Molecular Subtyping in Breast Cancer Research
#'
#'
#' @name BreastSubtypeR
#' @aliases NULL
#' @docType package
#'
#' @description
#' **BreastSubtypeR** is an R/Bioconductor package that unifies multiple
#' published intrinsic subtyping (IS) methods for breast cancer into a single,
#' reproducible framework. It supports both nearest-centroid (NC-based) and
#' single-sample predictor (SSP-based) classifiers and introduces an
#' assumption-aware **AUTO mode** that dynamically selects methods compatible
#' with the input cohort.
#'
#' By standardising input handling, applying method-specific normalisation,
#' and providing optimised probe-to-gene mapping, BreastSubtypeR reduces
#' inconsistencies across platforms and improves reproducibility in translational
#' research. A companion Shiny app (**iBreastSubtypeR**) offers an intuitive GUI
#' for non-programmers while preserving data privacy.
#'
#' ## Workflow
#' 1. **Data Input**: Supply a gene expression dataset as a `SummarizedExperiment`.
#'    Supported inputs include raw RNA-seq counts (with gene lengths),
#'    log2(FPKM+1) RNA-seq, or log2-normalised microarray/nCounter data.
#' 2. **Gene Mapping**: Prepare expression data with \code{\link{Mapping}},
#'    including Entrez ID-based resolution of duplicates.
#' 3. **Subtyping**: Apply multiple classifiers simultaneously using
#'    \code{\link{BS_Multi}}, or enable **AUTO mode** for
#'    cohort-aware method selection.
#' 4. **Visualisation**: Summarise and interpret subtyping results with
#'    \code{\link{Vis_Multi}}.
#'
#'
#' ## Key Features
#' - **Multi-method framework**: Ten published NC- and SSP-based classifiers,
#'   harmonised under one interface.
#' - **AUTO mode**: Evaluates cohort composition (e.g., ER/HER2 prevalence,
#'   subtype purity, subgroup sizes) and disables classifiers with violated
#'   assumptions; improves accuracy, Cohen’s kappa, and IHC concordance.
#' - **Standardised normalisation**: Upper-quartile log2-CPM for NC-based
#'   methods; FPKM for SSP-based methods.
#' - **Optimised gene mapping**: Entrez ID-based mapping with conflict resolution.
#' - **Dual accessibility**: A Bioconductor-compliant R API and a local Shiny app
#'   (iBreastSubtypeR).
#'
#' @seealso \code{\link{Mapping}}, \code{\link{BS_Multi}}, \code{\link{Vis_Multi}}
#'
#'
"_PACKAGE"


#' Gene ID Mapping
#'
#' @name Mapping
#' @description Preprocesses and maps gene expression input to prepare for
#'   intrinsic subtyping workflows (NC- and SSP-based).
#'
#' @param se_obj A `SummarizedExperiment` object containing:
#'   - **Assay data**:
#'     - If `RawCounts = FALSE`: `assay()` must contain log2-normalized expression (e.g., pre-normalized microarray/nCounter, or log2(FPKM+1) RNAseq).
#'     - If `RawCounts = TRUE`: `assay()` contains raw RNA-seq counts (see `RawCounts`).
#'   - **Row metadata** (required):
#'     - `"probe"`: feature identifiers (e.g., gene symbols or probe IDs)
#'     - `"ENTREZID"`: corresponding Entrez Gene IDs.
#'     - If row names are gene symbols, provide an additional `SYMBOL` column,
#'       renamed as `probe`.
#'   - **Column metadata** (optional): sample-level metadata in `colData()`.
#'
#' @param RawCounts Logical. If `TRUE`, indicates that `assay()` holds raw RNA-seq counts.
#'   In this case, `rowData()` must also provide gene lengths (column `"Length"`, in base pairs), used for:
#'   - NC-based methods: log2-CPM (upper-quartile normalization).
#'   - SSP-based methods: linear FPKM (not log-transformed).
#'
#' @param method Strategy for resolving duplicate probes/genes. Options:
#'   - `"iqr"`: probe with highest interquartile range (short-oligo arrays, e.g., Affymetrix).
#'   - `"mean"`: probe with highest mean expression (long-oligo arrays, e.g., Agilent/Illumina).
#'   - `"max"`: probe with highest expression value (often used for RNA-seq).
#'   - `"stdev"`: probe with highest standard deviation.
#'   - `"median"`: probe with highest median expression.
#'
#' @param impute Logical. If `TRUE`, applies KNN-based imputation to missing values.
#'
#' @param verbose Logical. If `TRUE`, prints progress messages during execution.
#'
#' @return A named list with:
#' \describe{
#'   \item{se_NC}{`SummarizedExperiment` holding log2-transformed data prepared for NC-based methods (assay name: `counts`).}
#'   \item{se_SSP}{`SummarizedExperiment` holding linear-scale data prepared for SSP-based methods (assay name: `counts`).}
#' }
#'
#' @references
#' Yang Q, Hartman J, Sifakis EG.
#' *BreastSubtypeR: A Unified R/Bioconductor Package for Intrinsic Molecular Subtyping in Breast Cancer Research.*
#' NAR Genomics and Bioinformatics. 2025. https://doi.org/10.1093/nargab/lqaf131. Selected as Editor’s Choice.
#'
#' @details
#' `Mapping()` supports multiple input types:
#' - **Raw RNA-seq counts** (with gene lengths): normalized to CPM (NC) or FPKM (SSP).
#' - **Precomputed log2(FPKM+1)**: used directly for NC; back-transformed for SSP.
#' - **log2-normalized microarray/nCounter data**: used directly for NC; back-transformed for SSP.
#'
#' This design allows users to supply a single expression format, while
#' BreastSubtypeR automatically applies method-specific preprocessing.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'
#' # Using example raw RNA-seq counts (with gene lengths)
#' data("TCGABRCAobj")
#' se_obj_counts <- TCGABRCAobj$se_obj
#' res <- Mapping(se_obj_counts, RawCounts = TRUE)
#'
#' # Using example pre-normalized log2(FPKM+0.1)
#' data("OSLO2EMIT0obj")
#' se_obj_fpkm <- OSLO2EMIT0obj$se_obj
#' res <- Mapping(se_obj_fpkm, RawCounts = FALSE)
#' }
#' }
#' 
#' @export

Mapping <- function(
        se_obj,
        RawCounts = FALSE,
        method = c("max", "mean", "median", "iqr", "stdev"),
        impute = TRUE,
        verbose = TRUE) {
    method <- match.arg(method)

    arguments <- rlang::dots_list(
        se_obj = se_obj,
        RawCounts = RawCounts,
        method = method,
        impute = impute,
        verbose = verbose,
        .homonyms = "last"
    )

    call <- rlang::call2(domapping, !!!arguments)
    res <- eval(call)

    # Extract results for NC-based and SSP-based methods
    x_NC <- res$x_NC
    x_SSP <- res$x_SSP

    ## create se_obj for NC_based methods
    se_NC <- SummarizedExperiment(
        assays = list(counts = x_NC),
        colData = colData(se_obj)
    )
    ## create se_obj for SSP_based methods
    se_SSP <- SummarizedExperiment(
        assays = list(counts = x_SSP),
        colData = colData(se_obj)
    )
    # Return both objects as a list
    return(list(se_NC = se_NC, se_SSP = se_SSP))
}


#' Original Parker Intrinsic Subtyping (BS_parker)
#'
#' @name BS_parker
#' @description
#' Implements the original PAM50 nearest-centroid classifier as described by
#' Parker et al. (2009), along with supported calibration strategies and
#' variations. This function assigns intrinsic breast cancer subtypes
#' (Luminal A, Luminal B, HER2-enriched, Basal-like, and optionally Normal-like).
#'
#' @param se_obj A `SummarizedExperiment` object containing:
#'   - **Assay data**: A log-transformed, normalized gene expression matrix with genes (Gene Symbols) as rows
#'   and samples as columns.
#'   - **Column metadata** (`colData`): Optional sample- or patient-level
#'     information.
#'
#' @param calibration Character. One of:
#'   - `"None"`: no centering/scaling.
#'   - `"Internal"`: center by a method derived from the current cohort (see `internal`).
#'   - `"External"`: center by medians from an external cohort (see `external`).
#'
#' @param internal Internal calibration method used when `calibration = "Internal"`.
#'   Accepts:
#'   - `NA` or `"medianCtr"` (identical): gene-wise median centering (as in Parker et al.).
#'   - `"meanCtr"`: gene-wise z-scoring (mean 0, sd 1; as implemented in `genefu.scale`).
#'   - `"qCtr"`: robust centering (quantile rescale with mq = 0.05; as in `genefu.robust`).
#'   Defaults to `NA` (median centering).
#'
#' @param external Character string specifying the external calibration source.
#'   - To use training cohort medians, provide the platform/column name.
#'   - To supply user-defined medians, set `external = "Given.mdns"` and pass
#'     values via `medians`.
#'
#' @param medians A matrix or data.frame of user-provided medians (required if
#'   `external = "Given.mdns"`).
#'   - First column: 50 PAM50 genes.
#'   - Second column: Corresponding median expression values.
#'
#' @param Subtype Logical. If `TRUE`, assigns only the four main intrinsic
#'   subtypes (Luminal A, Luminal B, HER2-enriched, Basal-like),
#'   excluding Normal-like.
#'
#' @param hasClinical Logical. If `TRUE`, incorporates clinical variables from
#'   `colData(se_obj)`. Required columns:
#' - "TSIZE": Tumor size (0 = \eqn{\le 2}{<= 2} cm; 1 = \eqn{> 2}{> 2} cm).
#' - "NODE": Lymph node status (0 = negative; \eqn{\ge 1}{>= 1} = positive). Must be numeric.
#'
#' @return A list containing PAM50 intrinsic subtype calls using the Parker
#'   classifier and selected calibration strategy.
#'
#'
#' @references
#' - Parker JS, Mullins M, Cheung MCU, Leung S, Voduc D, et al.
#' *Supervised risk predictor of breast cancer based on intrinsic subtypes*.
#' Journal of Clinical Oncology. 2009;27(8).
#' https://doi.org/10.1200/JCO.2008.18.1370
#'
#' - Gendoo DMA, Ratanasirigulchai N, Schröder MS, Paré L, Parker JS, Prat A, et al.
#' *Genefu: An R/Bioconductor package for computation of gene expression-based signatures in breast cancer*.
#' Bioinformatics. 2016;32(7).
#' https://doi.org/10.1093/bioinformatics/btv693
#'
#' @examples
#' data("OSLO2EMIT0obj")
#' res <- BS_parker(
#'     se_obj = OSLO2EMIT0obj$data_input$se_NC,
#'     calibration = "Internal",
#'     internal = NA, # NA is equal to "medianCtr"
#'     Subtype = FALSE,
#'     hasClinical = FALSE
#' )
#'
#' @export

BS_parker <- function(
        se_obj,
        calibration = "None",
        internal = NA,
        external = NA,
        medians = NA,
        Subtype = FALSE,
        hasClinical = FALSE) {
    # Check if input is a SummarizedExperiment object
    if (!inherits(se_obj, "SummarizedExperiment")) {
        stop("Input must be a SummarizedExperiment object.")
    }

    ## --- robust calibration handling & legacy -1 guard ---
    # tolerate NULL/NA/missing and unnamed vectors
    if (missing(calibration) || is.null(calibration) || is.na(calibration)) {
        calibration <- "None"
    } else if (length(calibration) > 1L) {
        calibration <- calibration[1L]
    }

    # legacy semantics: internal == -1 or "-1" means *no* adjustment
    if (!missing(internal) && !is.null(internal) && !all(is.na(internal))) {
        if (any(internal %in% c(-1, "-1"))) {
            calibration <- "None"
            internal <- NULL
        }
    }

    # Default Parker behavior: internal = NA → median centering
    if (identical(calibration, "Internal") &&
        (length(internal) == 0L || all(is.na(internal)))) {
        internal <- "medianCtr"
    }

    ## Extract data from SummarizedExperiment
    gene_expr <- assay(se_obj)
    pheno <- colData(se_obj) %>% data.frame()

    # Handle clinical metadata if required
    if (hasClinical) {
        req <- c("PatientID", "TSIZE", "NODE")
        miss <- setdiff(req, colnames(pheno))
        if (length(miss)) {
            stop(
                "When hasClinical = TRUE, colData(se_obj) must include: ",
                paste(req, collapse = ", "),
                ". Missing: ", paste(miss, collapse = ", "), "."
            )
        }
        if (!is.numeric(pheno$NODE)) {
            stop("colData(se_obj)$NODE must be numeric (0 for negative, >=1 for positive).")
        }
        rownames(pheno) <- pheno$PatientID
    } else {
        pheno <- NULL
    }

    arguments <- rlang::dots_list(
        mat = gene_expr,
        df.cln = pheno,
        calibration = calibration,
        internal = internal,
        external = external,
        medians = medians,
        Subtype = Subtype,
        hasClinical = hasClinical,
        .homonyms = "last"
    )

    call <- rlang::call2(makeCalls.parker, !!!arguments)
    res_parker <- eval(call)

    return(res_parker)
}

#' Conventional IHC Intrinsic Subtyping (BS_cIHC)
#'
#' @name BS_cIHC
#' @description
#' Implements the conventional immunohistochemistry-based (cIHC) intrinsic
#' subtyping approach, which balances cohorts by estrogen receptor (ER) status
#' before applying gene-expression–based subtyping. This method is useful for
#' ER-skewed cohorts where assumptions of nearest-centroid classifiers are
#' violated.
#'
#' @param se_obj A `SummarizedExperiment` object containing:
#'   - **Assay data**: A log2-transformed, normalised expression matrix with
#'     genes (Gene Symbols) as rows and samples as columns.
#'   - **Column metadata** (`colData`): Must include:
#'     - `"PatientID"`: Unique sample or patient identifier.
#'     - `"ER"`: Estrogen receptor status, coded as `"ER+"` or `"ER-"`.
#'
#' @param Subtype Logical. If `TRUE`, returns only the four main subtypes
#'   (Luminal A, Luminal B, HER2-enriched, Basal-like), excluding Normal-like.
#'
#' @param hasClinical Logical. If `TRUE`, incorporates additional clinical
#'   variables from `colData(se_obj)`. Required columns:
#' - "TSIZE": Tumor size (0 = \eqn{\le 2}{<= 2} cm; 1 = \eqn{> 2}{> 2} cm).
#' - "NODE": Lymph node status (0 = negative; \eqn{\ge 1}{>= 1} = positive). Must be numeric.
#'
#' @param seed Integer. Random seed for reproducibility of ER-balancing.
#'
#' @return A `data.frame` containing intrinsic subtype assignments estimated
#'   using the conventional IHC (cIHC) approach.
#'
#' @references
#' Ciriello G, Gatza ML, Beck AH, Wilkerson MD, Rhie SK, Pastore A, et al.
#' *Comprehensive Molecular Portraits of Invasive Lobular Breast Cancer.*
#' Cell. 2015;163(2):506–519.
#' https://doi.org/10.1016/j.cell.2015.09.033
#'
#' @examples
#' data("OSLO2EMIT0obj")
#' res <- BS_cIHC(
#'     se_obj = OSLO2EMIT0obj$data_input$se_NC,
#'     Subtype = FALSE,
#'     hasClinical = FALSE
#' )
#'
#' @export


BS_cIHC <- function(se_obj,
    Subtype = FALSE,
    hasClinical = FALSE,
    seed = 118) {
    # Check if input is a SummarizedExperiment object
    if (!inherits(se_obj, "SummarizedExperiment")) {
        stop("Input must be a SummarizedExperiment object.")
    }

    ## Extract data from SummarizedExperiment
    gene_expr <- assay(se_obj)
    # Extract clinical metadata if required
    pheno <- colData(se_obj) %>% data.frame()

    if (!all(c("PatientID", "ER") %in% colnames(pheno))) {
        stop("The 'colData' of 'se_obj' must include 'PatientID' and 'ER' columns.")
    }
    rownames(pheno) <- pheno$PatientID

    if (hasClinical) {
        req <- c("TSIZE", "NODE")
        miss <- setdiff(req, colnames(pheno))
        if (length(miss)) stop("When hasClinical = TRUE, colData(se_obj) must include: TSIZE and NODE. Missing: ", paste(miss, collapse = ", "), ".")
        if (!is.numeric(pheno$NODE)) stop("colData(se_obj)$NODE must be numeric.")
    }

    arguments <- rlang::dots_list(
        mat = gene_expr,
        df.cln = pheno,
        Subtype = Subtype,
        hasClinical = hasClinical,
        seed = seed,
        .homonyms = "last"
    )

    call <- rlang::call2(makeCalls_ihc, !!!arguments)
    res_IHC <- eval(call)

    return(res_IHC)
}


#' Iterative Conventional IHC Intrinsic Subtyping (BS_cIHC.itr)
#'
#' @name BS_cIHC.itr
#' @description
#' Implements an **iterative** version of the conventional IHC-based intrinsic
#' subtyping approach. This method repeatedly balances samples by estrogen
#' receptor (ER) status across multiple iterations, allowing refinement of
#' subtype calls in ER-skewed cohorts. Users can customise the ER+/ER– ratio to
#' match specific cohort assumptions (e.g., training distribution).
#'
#' @param se_obj A `SummarizedExperiment` object containing:
#'   - **Assay data**: A log2-transformed, normalised expression matrix with
#'     genes (Gene Symbols) as rows and samples as columns.
#'   - **Column metadata** (`colData`): Must include:
#'     - `"PatientID"`: Unique sample or patient identifier.
#'     - `"ER"`: Estrogen receptor status, coded as `"ER+"` or `"ER-"`.
#'
#' @param iteration Integer. Number of iterations for the ER-balancing procedure.
#'   Default: 100.
#'
#' @param ratio Numeric. Target ER+/ER– ratio for balancing. Options:
#'   - `1:1`: Equal balancing.
#'   - `54:64`: Default; reflects the ER+/ER– ratio in the UNC232 training cohort.
#'
#' @param Subtype Logical. If `TRUE`, returns only the four main subtypes
#'   (Luminal A, Luminal B, HER2-enriched, Basal-like), excluding Normal-like.
#'
#' @param hasClinical Logical. If `TRUE`, incorporates additional clinical
#'   variables from `colData(se_obj)`. Required columns:
#' - "TSIZE": Tumor size (0 = \eqn{\le 2}{<= 2} cm; 1 = \eqn{> 2}{> 2} cm).
#' - "NODE": Lymph node status (0 = negative; \eqn{\ge 1}{>= 1} = positive). Must be numeric.
#'
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return A list containing:
#'   - `subtypes`: Intrinsic subtype predictions across iterations.
#'   - `confidence`: Confidence estimates for each assigned subtype.
#'   - `ER_balance`: Proportions of ER+ and ER– subsets observed across iterations.
#'
#' @references
#' Curtis C, Shah SP, Chin SF, Turashvili G, Rueda OM, Dunning MJ, et al.
#' *The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups.*
#' Nature. 2012;486(7403):346–352.
#' https://doi.org/10.1038/nature10983
#'
#' @examples
#' data("OSLO2EMIT0obj")
#' res <- BS_cIHC.itr(
#'     se_obj = OSLO2EMIT0obj$data_input$se_NC,
#'     iteration = 10, ## for final analysis, use iteration = 100
#'     Subtype = FALSE,
#'     hasClinical = FALSE
#' )
#'
#' @export

BS_cIHC.itr <- function(se_obj,
    iteration = 100,
    ratio = 54 / 64,
    Subtype = FALSE,
    hasClinical = FALSE,
    seed = 118) {
    # Check if input is a SummarizedExperiment object
    if (!inherits(se_obj, "SummarizedExperiment")) {
        stop("Input must be a SummarizedExperiment object.")
    }

    ## Extract data from SummarizedExperiment
    gene_expr <- assay(se_obj)
    # Extract clinical metadata if required
    pheno <- colData(se_obj) %>% data.frame()

    if (!all(c("PatientID", "ER") %in% colnames(pheno))) {
        stop("The 'colData' of 'se_obj' must include 'PatientID' and 'ER' columns.")
    }
    rownames(pheno) <- pheno$PatientID

    if (hasClinical) {
        req <- c("TSIZE", "NODE")
        miss <- setdiff(req, colnames(pheno))
        if (length(miss)) stop("When hasClinical = TRUE, colData(se_obj) must include: TSIZE and NODE. Missing: ", paste(miss, collapse = ", "), ".")
        if (!is.numeric(pheno$NODE)) stop("colData(se_obj)$NODE must be numeric.")
    }

    arguments <- rlang::dots_list(
        mat = gene_expr,
        df.cln = pheno,
        iteration = iteration,
        ratio = ratio,
        Subtype = Subtype,
        hasClinical = hasClinical,
        seed = seed,
        .homonyms = "last"
    )

    call <- rlang::call2(makeCalls_ihc.iterative, !!!arguments)
    res_IHC.itr <- eval(call)

    return(res_IHC.itr)
}


#' PCA-PAM50 Intrinsic Subtyping (BS_PCAPAM50)
#'
#' @name BS_PCAPAM50
#' @description
#' Implements the PCA-PAM50 method, which integrates **Principal Component
#' Analysis (PCA)** of ESR1 expression to adjust for estrogen receptor (ER)
#' imbalance prior to applying the PAM50 nearest-centroid classifier. This
#' approach improves subtype consistency, particularly in ER-skewed cohorts.
#'
#' @param se_obj A `SummarizedExperiment` object containing:
#'   - **Assay data**: A log2-transformed, normalised expression matrix with
#'     genes (Gene Symbols) as rows and samples as columns.
#'   - **Column metadata** (`colData`): Must include:
#'     - `"PatientID"`: Unique sample or patient identifier.
#'     - `"ER"`: Estrogen receptor status, coded as `"ER+"` or `"ER-"`.
#'
#' @param Subtype Logical. If `TRUE`, returns only the four main subtypes
#'   (Luminal A, Luminal B, HER2-enriched, Basal-like), excluding Normal-like.
#'
#' @param hasClinical Logical. If `TRUE`, incorporates additional clinical
#'   variables from `colData(se_obj)`. Required columns:
#' - "TSIZE": Tumor size (0 = \eqn{\le 2}{<= 2} cm; 1 = \eqn{> 2}{> 2} cm).
#' - "NODE": Lymph node status (0 = negative; \eqn{\ge 1}{>= 1} = positive). Must be numeric.
#'
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return A character vector of intrinsic subtype predictions assigned to each
#'   sample using the PCA-PAM50 method.
#'
#' @references
#' Raj-Kumar PK, Liu J, Hooke JA, Kovatich AJ, Kvecher L, Shriver CD, et al.
#' *PCA-PAM50 improves consistency between breast cancer intrinsic and clinical subtyping, reclassifying a subset of luminal A tumors as luminal B.*
#' Scientific Reports. 2019;9(1):1–12.
#' https://doi.org/10.1038/s41598-019-44339-4
#'
#' @examples
#' data("OSLO2EMIT0obj")
#' res <- BS_PCAPAM50(
#'     se_obj = OSLO2EMIT0obj$data_input$se_NC,
#'     Subtype = FALSE,
#'     hasClinical = FALSE
#' )
#'
#' @export

BS_PCAPAM50 <- function(se_obj,
    Subtype = FALSE,
    hasClinical = FALSE,
    seed = 118) {
    # Check if input is a SummarizedExperiment object
    if (!inherits(se_obj, "SummarizedExperiment")) {
        stop("Input must be a SummarizedExperiment object.")
    }

    ## Extract data from SummarizedExperiment
    gene_expr <- assay(se_obj)
    # Extract clinical metadata if required
    pheno <- colData(se_obj) %>% data.frame()

    if (!all(c("PatientID", "ER") %in% colnames(pheno))) {
        stop("The 'colData' of 'se_obj' must include 'PatientID' and 'ER' columns.")
    }
    rownames(pheno) <- pheno$PatientID

    if (hasClinical) {
        req <- c("TSIZE", "NODE")
        miss <- setdiff(req, colnames(pheno))
        if (length(miss)) stop("When hasClinical = TRUE, colData(se_obj) must include: TSIZE and NODE. Missing: ", paste(miss, collapse = ", "), ".")
        if (!is.numeric(pheno$NODE)) stop("colData(se_obj)$NODE must be numeric.")
    }

    samples <- pheno$PatientID

    if ("ER" %in% colnames(pheno)) {
        ## create IHC column for PCAPAM50
        pheno$IHC <- dplyr::case_when(pheno$ER == "ER+" ~ "Luminal",
            pheno$ER == "ER-" ~ "non-Luminal",
            .default = NA
        )

        pheno$ER_status <- NA

        pheno$ER_status[which(pheno$ER == "ER+")] <- "pos"
        pheno$ER_status[which(pheno$ER == "ER-")] <- "neg"
        pheno <- pheno[order(pheno$ER_status, decreasing = TRUE), ]
    }

    gene_expr <- as.matrix(gene_expr[, pheno$PatientID])

    ## first step
    arguments <- rlang::dots_list(
        mat = gene_expr,
        df.cln = pheno,
        Subtype = Subtype,
        hasClinical = hasClinical,
        seed = seed,
        .homonyms = "last"
    )

    call <- rlang::call2(makeCalls.PC1ihc, !!!arguments)
    res_PC1IHC <- eval(call)

    ## second step
    if (hasClinical) {
        df.pc1pam <- data.frame(
            PatientID = res_PC1IHC$BS.all$PatientID,
            PAM50 = res_PC1IHC$BS.all$BS,
            TSIZE = pheno[res_PC1IHC$BS.all$PatientID, ]$TSIZE,
            NODE = pheno[res_PC1IHC$BS.all$PatientID, ]$NODE,
            stringsAsFactors = FALSE
        )
    } else {
        df.pc1pam <- data.frame(
            PatientID = res_PC1IHC$BS.all$PatientID,
            PAM50 = res_PC1IHC$BS.all$BS,
            stringsAsFactors = FALSE
        )
    }

    arguments2 <- rlang::dots_list(
        mat = gene_expr,
        df.pam = df.pc1pam,
        Subtype = Subtype,
        hasClinical = hasClinical,
        seed = seed,
        .homonyms = "last"
    )

    call <- rlang::call2(makeCalls.v1PAM, !!!arguments2)
    res_PCAPAM50 <- eval(call)

    ## reorder
    idx <- na.omit(match(samples, res_PCAPAM50$BS.all$PatientID))
    res_PCAPAM50$BS.all <- res_PCAPAM50$BS.all[idx, ]

    return(res_PCAPAM50)
}


#' Subgroup-Specific Gene-Centering Intrinsic Subtyping (BS_ssBC)
#'
#' @name BS_ssBC
#' @description
#' Implements the **subgroup-specific gene-centering (ssBC)** method for breast
#' cancer intrinsic subtyping. The ssBC approach applies precomputed,
#' subgroup-specific centering values to adjust PAM50 nearest-centroid
#' classification when the study cohort is skewed relative to the original
#' training cohort (e.g., ER-selected, HER2-enriched, or triple-negative cohorts).
#'
#' @param se_obj A `SummarizedExperiment` object containing:
#'   - **Assay data**: A log2-transformed, normalised expression matrix with
#'     genes (Gene Symbols) as rows and samples as columns.
#'   - **Column metadata** (`colData`): If `hasClinical = TRUE`, must include:
#'     - `"PatientID"`: Unique patient/sample identifier.
#'     - Depending on the chosen `s` parameter:
#'       - `"ER"`: Estrogen receptor status (`"ER+"` or `"ER-"`) if `s = "ER"`.
#'       - `"HER2"`: HER2 status (`"HER2+"` or `"HER2-"`) if `s = "ER.v2"`.
#'       - `"TN"`: Triple-negative status (`"TN"` or `"nonTN"`) if `s = "TN"` or `"TN.v2"`.
#'
#' @param s Character. Specifies which subgroup-specific quantiles to use:
#'   - `"ER"`, `"TN"`: Original subgroup-specific quantiles (*Breast Cancer Research*, 2015).
#'   - `"ER.v2"`, `"TN.v2"`: Updated subgroup-specific quantiles (*Journal of Clinical Oncology*, 2020).
#'
#' @param Subtype Logical. If `TRUE`, returns only the four main subtypes
#'   (Luminal A, Luminal B, HER2-enriched, Basal-like), excluding Normal-like.
#'
#' @param hasClinical Logical. If `TRUE`, incorporates additional clinical
#'   variables from `colData(se_obj)`. Required columns:
#' - "TSIZE": Tumor size (0 = \eqn{\le 2}{<= 2} cm; 1 = \eqn{> 2}{> 2} cm).
#' - "NODE": Lymph node status (0 = negative; \eqn{\ge 1}{>= 1} = positive). Must be numeric.
#'
#' @return A character vector of intrinsic subtype predictions assigned to each
#'   sample using the ssBC method.
#'
#' @references
#' Zhao X, Rodland EA, Tibshirani R, Plevritis S.
#' *Molecular subtyping for clinically defined breast cancer subgroups.*
#' Breast Cancer Research. 2015;17(1):29.
#' https://doi.org/10.1186/s13058-015-0520-4
#'
#' Fernandez-Martinez A, Krop IE, Hillman DW, Polley MY, Parker JS, Huebner L, et al.
#' *Survival, pathologic response, and genomics in CALGB 40601 (Alliance), a neoadjuvant Phase III trial of paclitaxel–trastuzumab with or without lapatinib in HER2-positive breast cancer.*
#' Journal of Clinical Oncology. 2020;38(36):4184–4197.
#' https://doi.org/10.1200/JCO.20.01276
#'
#' @examples
#' ## Example: Updated subgroup-specific quantiles (ER.v2)
#' data("OSLO2EMIT0obj")
#' res <- BS_ssBC(
#'     se_obj = OSLO2EMIT0obj$data_input$se_NC,
#'     s = "ER.v2",
#'     Subtype = FALSE,
#'     hasClinical = FALSE
#' )
#'
#' @export

BS_ssBC <- function(se_obj,
    s,
    Subtype = FALSE,
    hasClinical = FALSE) {
    # Check that input is a SummarizedExperiment object
    if (!inherits(se_obj, "SummarizedExperiment")) {
        stop("Input must be a SummarizedExperiment object.")
    }

    # Extract gene expression matrix
    gene_expr <- assay(se_obj)

    # Extract clinical metadata if hasClinical is TRUE
    pheno <- colData(se_obj) %>% data.frame()
    if (!"PatientID" %in% colnames(pheno)) {
        stop("The 'colData' of 'se_obj' must include a 'PatientID' column.")
    }
    rownames(pheno) <- pheno$PatientID

    if (hasClinical) {
        req <- c("TSIZE", "NODE")
        miss <- setdiff(req, colnames(pheno))
        if (length(miss)) stop("When hasClinical = TRUE, colData(se_obj) must include: TSIZE and NODE. Missing: ", paste(miss, collapse = ", "), ".")
        if (!is.numeric(pheno$NODE)) stop("colData(se_obj)$NODE must be numeric.")
    }

    # Additional checks based on `s`
    required_columns <- switch(s,
        "ER" = c("ER"),
        "ER.v2" = c("ER", "HER2"),
        "TN" = c("TN"),
        "TN.v2" = c("TN"),
        stop("Invalid value for 's'. Must be one of 'ER', 'ER.v2', 'TN', or 'TN.v2'.")
    )

    missing_columns <- setdiff(required_columns, colnames(pheno))
    if (length(missing_columns) > 0) {
        stop(paste("The following required columns are missing from 'colData':", paste(missing_columns, collapse = ", ")))
    }

    arguments <- rlang::dots_list(
        mat = gene_expr,
        df.cln = pheno,
        s = s,
        Subtype = Subtype,
        hasClinical = hasClinical,
        .homonyms = "last"
    )

    call <- rlang::call2(makeCalls.ssBC, !!!arguments)
    res_ssBC <- eval(call)

    return(res_ssBC)
}


#' AIMS Intrinsic Subtyping (BS_AIMS)
#'
#' @name BS_AIMS
#' @description
#' Implements the **AIMS (Absolute Assignment of Intrinsic Molecular Subtype)**
#' method for breast cancer intrinsic subtyping. Unlike nearest-centroid (NC)
#' approaches, AIMS is a single-sample predictor (SSP): it assigns subtypes
#' independently for each sample using within-sample, pairwise gene expression
#' rules. This makes it robust to cohort composition and scaling.
#'
#' @param se_obj A `SummarizedExperiment` object containing:
#'   - **Assay data**: A gene expression matrix with genes (Entrez IDs) as rows
#'     and samples as columns.
#'     - Expression values must be **positive** (e.g., FPKM or log2(FPKM+1)).
#'     - Values should not be gene-centered or globally scaled.
#'
#' @return A character vector of intrinsic subtype predictions assigned to each
#'   sample using the AIMS method.
#'
#' @references
#' Paquet ER, Hallett MT.
#' *Absolute assignment of breast cancer intrinsic molecular subtype.*
#' Journal of the National Cancer Institute. 2015;107(1):dju357.
#' https://doi.org/10.1093/jnci/dju357
#'
#' @examples
#' ## Example using SummarizedExperiment input
#' data("OSLO2EMIT0obj")
#' res <- BS_AIMS(
#'     se_obj = OSLO2EMIT0obj$data_input$se_SSP
#' )
#'
#' @export

BS_AIMS <- function(se_obj) {
    ## loading datasets

    data_env <- new.env(parent = emptyenv())
    data("AIMSmodel", envir = data_env, package = "BreastSubtypeR")
    AIMSmodel <- data_env[["AIMSmodel"]]
    data("BreastSubtypeRobj", envir = data_env, package = "BreastSubtypeR")
    BreastSubtypeRobj <- data_env[["BreastSubtypeRobj"]]


    # Check that input is a SummarizedExperiment object
    if (!inherits(se_obj, "SummarizedExperiment")) {
        stop("Input must be a SummarizedExperiment object.")
    }

    # Extract gene expression matrix
    gene_expr <- assay(se_obj)

    # Extract AIMS-specific genes
    genes <- as.character(
        BreastSubtypeRobj$genes.signature$EntrezGene.ID[
            which(BreastSubtypeRobj$genes.signature$AIMS == "Yes")
        ]
    )
    gene_expr <- gene_expr[rownames(gene_expr) %in% genes, ]

    arguments <- rlang::dots_list(
        eset = as.matrix(gene_expr),
        EntrezID = rownames(gene_expr)
    )

    call <- rlang::call2(applyAIMS_AIMS, !!!arguments)

    res_AIMS <- eval(call)
    return(res_AIMS)
}


#' Intrinsic Subtyping using SSPBC (BS_sspbc)
#'
#' @name BS_sspbc
#' @description
#' Implements **SSPBC (Single Sample Predictor for Breast Cancer)**, a
#' refinement of the original AIMS methodology trained on the large,
#' population-based SCAN-B RNA-seq cohort. SSPBC provides robust
#' single-sample predictions, independent of cohort composition, and supports
#' multiple model variants for different applications.
#'
#' @param se_obj A `SummarizedExperiment` object containing:
#'   - **Assay data**: A gene expression matrix with genes (Entrez IDs) as rows
#'     and samples as columns.
#'     - Expression values must be **positive** (e.g., FPKM or log2(FPKM+1)).
#'     - Values should not be gene-centered or globally scaled.
#'
#' @param ssp.name Character. Specifies the SSPBC model to use:
#'   - `"ssp.pam50"`: Predicts PAM50-based intrinsic subtypes.
#'   - `"ssp.subtype"`: Predicts Prosigna-like subtypes (four subtypes, excluding Normal-like).
#'
#' @return A character vector of intrinsic subtype predictions for each sample,
#'   as estimated by the SSPBC method.
#'
#' @references
#' Staaf J, Häkkinen J, Hegardt C, Saal LH, Kimbung S, Hedenfalk I, et al.
#' *RNA sequencing-based single sample predictors of molecular subtype and risk of recurrence for clinical assessment of early-stage breast cancer.*
#' NPJ Breast Cancer. 2022;8(1):27.
#' https://doi.org/10.1038/s41523-022-00465-3
#'
#' @examples
#' ## Example using SSPBC with the PAM50 model
#' data("OSLO2EMIT0obj")
#' res <- BS_sspbc(
#'     se_obj = OSLO2EMIT0obj$data_input$se_SSP,
#'     ssp.name = "ssp.pam50"
#' )
#'
#' @export

BS_sspbc <- function(se_obj, ssp.name = "ssp.pam50") {
    data_env <- new.env(parent = emptyenv())
    data("sspbc.models", envir = data_env, package = "BreastSubtypeR")
    sspbc.models <- data_env[["sspbc.models"]]
    data("sspbc.models.fullname", envir = data_env, package = "BreastSubtypeR")
    sspbc.models.fullname <- data_env[["sspbc.models.fullname"]]

    # Check that input is a SummarizedExperiment object
    if (!inherits(se_obj, "SummarizedExperiment")) {
        stop("Input must be a SummarizedExperiment object.")
    }

    # Extract gene expression matrix
    gene_expr <- assay(se_obj)

    arguments <- rlang::dots_list(
        gex = gene_expr,
        id = rownames(gene_expr),
        id.type = "EntrezGene",
        ssp.name = ssp.name,
        .homonyms = "last"
    )

    call <- rlang::call2(applySSP, !!!arguments)
    res_sspbc <- eval(call)
    return(res_sspbc)
}

#' Intrinsic Subtyping with Multiple Approaches (BS_Multi)
#'
#' @name BS_Multi
#' @description
#' Executes multiple intrinsic molecular subtyping methods in parallel.
#' Users can either specify a set of classifiers directly, or enable the
#' **AUTO mode**, which dynamically selects methods based on cohort composition
#' (e.g., ER/HER2 distribution, subtype purity, subgroup size).
#' AUTO reduces misclassification in skewed or subtype-specific cohorts by
#' disabling methods whose assumptions are violated, but does not perform
#' consensus voting—subtypes are still returned per method.
#'
#' @param data_input The output from the [`Mapping()`] function, containing
#'   processed gene expression data prepared for subtyping.
#'
#' @param methods Character vector specifying the subtyping methods to run.
#'   Available options include:
#'   - `"parker.original"`: Original PAM50 (Parker et al., 2009).
#'   - `"genefu.scale"`: PAM50 (scaled version; Gendoo et al., 2016).
#'   - `"genefu.robust"`: PAM50 (robust version; Gendoo et al., 2016).
#'   - `"cIHC"`: Conventional ER-balancing with immunohistochemistry (Ciriello et al., 2015).
#'   - `"cIHC.itr"`: Iterative ER-balancing (Curtis et al., 2012).
#'   - `"PCAPAM50"`: PCA-based PAM50 using ESR1 balancing (Raj-Kumar et al., 2019).
#'   - `"ssBC"`: Subgroup-specific gene-centering (Zhao et al., 2015).
#'   - `"ssBC.v2"`: Updated subgroup-specific centering (Fernandez-Martinez et al., 2020).
#'   - `"AIMS"`: Absolute Intrinsic Molecular Subtyping (Paquet & Hallett, 2015).
#'   - `"sspbc"`: SSPBC, a large-cohort SSP trained on SCAN-B (Staaf et al., 2022).
#'   - `"AUTO"`: Cohort-aware selection of compatible methods (must be the only entry).
#'
#'   **Notes:**
#'   - If `"AUTO"` is selected, it must be the sole value in `methods`.
#'   - Otherwise, at least **two** methods must be specified.
#'
#' @param Subtype Logical. If `TRUE`, returns four subtypes (Luminal A, Luminal B,
#'   HER2-enriched, Basal-like), excluding Normal-like.
#'
#' @param hasClinical Logical. If `TRUE`, incorporates clinical data from
#'   `colData(se_obj)`. Required columns:
#' - "TSIZE": Tumor size (0 = \eqn{\le 2}{<= 2} cm; 1 = \eqn{> 2}{> 2} cm).
#' - "NODE": Lymph node status (0 = negative; \eqn{\ge 1}{>= 1} = positive).
#'
#' @return A list containing per-method subtype assignments for each sample.
#'
#' @references
#' Yang Q, Hartman J, Sifakis EG.
#' *BreastSubtypeR: A Unified R/Bioconductor Package for Intrinsic Molecular Subtyping in Breast Cancer Research.*
#' NAR Genomics and Bioinformatics. 2025. https://doi.org/10.1093/nargab/lqaf131. Selected as Editor’s Choice.
#'
#' Parker JS, Mullins M, Cheung MCU, Leung S, Voduc D, et al.
#' *Supervised risk predictor of breast cancer based on intrinsic subtypes.*
#' J Clin Oncol. 2009;27(8):1160-1167. https://doi.org/10.1200/JCO.2008.18.1370
#'
#' Gendoo DMA, Ratanasirigulchai N, Schröder MS, Paré L, Parker JS, Prat A, et al.
#' *Genefu: An R/Bioconductor package for computation of gene expression-based signatures in breast cancer.*
#' Bioinformatics. 2016;32(7):1097-1099. https://doi.org/10.1093/bioinformatics/btv693
#'
#' Ciriello G, Gatza ML, Beck AH, Wilkerson MD, Rhie SK, Pastore A, et al.
#' *Comprehensive molecular portraits of invasive lobular breast cancer.*
#' Cell. 2015;163(2):506-519. https://doi.org/10.1016/j.cell.2015.09.033
#'
#' Curtis C, Shah SP, Chin SF, Turashvili G, Rueda OM, Dunning MJ, et al.
#' *The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups.*
#' Nature. 2012;486(7403):346-352. https://doi.org/10.1038/nature10983
#'
#' Zhao X, Rodland EA, Tibshirani R, Plevritis S.
#' *Molecular subtyping for clinically defined breast cancer subgroups.*
#' Breast Cancer Res. 2015;17(1):29. https://doi.org/10.1186/s13058-015-0520-4
#'
#' Fernandez-Martinez A, Krop IE, Hillman DW, Polley MY, Parker JS, Huebner L, et al.
#' *Survival, pathologic response, and genomics in CALGB 40601 (Alliance), a neoadjuvant Phase III trial of paclitaxel–trastuzumab with or without lapatinib in HER2-positive breast cancer.*
#' J Clin Oncol. 2020;38(36):4184-4197. https://doi.org/10.1200/JCO.20.01276
#'
#' Paquet ER, Hallett MT.
#' *Absolute assignment of breast cancer intrinsic molecular subtype.*
#' J Natl Cancer Inst. 2015;107(1):dju357. https://doi.org/10.1093/jnci/dju357
#'
#' Staaf J, Häkkinen J, Hegardt C, Saal LH, Kimbung S, Hedenfalk I, et al.
#' *RNA sequencing-based single sample predictors of molecular subtype and risk of recurrence for clinical assessment of early-stage breast cancer.*
#' NPJ Breast Cancer. 2022;8(1):27. https://doi.org/10.1038/s41523-022-00465-3
#'
#' @examples
#' ## Example: run multiple methods
#' data("OSLO2EMIT0obj")
#' methods <- c("parker.original", "genefu.scale", "genefu.robust")
#' res.test <- BS_Multi(
#'     data_input = OSLO2EMIT0obj$data_input,
#'     methods = methods,
#'     Subtype = FALSE,
#'     hasClinical = FALSE
#' )
#'
#' @export

BS_Multi <- function(
        data_input,
        methods = "AUTO",
        Subtype = FALSE,
        hasClinical = FALSE) {
    valid_methods <- c(
        "parker.original", "genefu.scale", "genefu.robust",
        "ssBC", "ssBC.v2", "cIHC", "cIHC.itr", "PCAPAM50",
        "AIMS", "sspbc", "AUTO"
    )

    invalid_methods <- methods[!methods %in% valid_methods]

    if (length(invalid_methods) > 0) {
        stop("Invalid method(s) specified: ", paste(invalid_methods, collapse = ", "))
    }

    if (length(methods) == 1 && methods[1] == "AUTO") {
        message("Running AUTO mode for subtyping.")
    } else if (length(methods) < 2) {
        stop("Select at least two methods or set method to 'AUTO'.")
    }

    ## extract pheno table
    pheno <- colData(data_input$se_NC) %>% data.frame()
    if (!"PatientID" %in% colnames(pheno)) {
        stop("colData(se_NC) must have a 'PatientID' column.")
    }
    rownames(pheno) <- pheno$PatientID

    # Check ER and HER2 columns in pheno
    if (!("ER" %in% colnames(pheno)) && any(methods %in% c("ssBC", "ssBC.v2", "cIHC", "cIHC.itr", "PCAPAM50"))) {
        stop("The 'ER' column is required for selected methods.")
    }
    if (!("HER2" %in% colnames(pheno)) && "ssBC.v2" %in% methods) {
        stop("The 'HER2' column is required for the 'ssBC.v2' method.")
    }

    ## AUTO mode
    # methods = "AUTO"
    cohort.select <- "ERpos"
    samples_ER.icd <- NULL
    samples_ERHER2.icd <- NULL
    if (length(methods) == 1 && methods[1] == "AUTO") {
        AUTO.output <- get_methods(pheno)

        samples_ER.icd <- AUTO.output$samples_ER.icd
        samples_ERHER2.icd <- AUTO.output$samples_ERHER2.icd
        methods <- AUTO.output$methods
        cohort.select <- AUTO.output$cohort.select
    }

    ## run each method
    results <- lapply(methods, function(method) {
        ## try NC-based
        if (method == "parker.original") {
            message(method, " is running!")
            return(
                BS_parker(
                    data_input$se_NC,
                    calibration = "Internal",
                    internal = "medianCtr", # default to medianCtr in BS_parker()
                    Subtype = Subtype,
                    hasClinical = hasClinical
                )
            )
        }
        if (method == "genefu.scale") {
            message(method, " is running!")
            return(
                BS_parker(
                    data_input$se_NC,
                    calibration = "Internal",
                    internal = "meanCtr",
                    Subtype = Subtype,
                    hasClinical = hasClinical
                )
            )
        }
        if (method == "genefu.robust") {
            message(method, " is running!")
            return(
                BS_parker(
                    data_input$se_NC,
                    calibration = "Internal",
                    internal = "qCtr",
                    Subtype = Subtype,
                    hasClinical = hasClinical
                )
            )
        }
        if (method == "cIHC") {
            ## try conventional IHC
            message(method, " is running!")
            return(
                BS_cIHC(
                    data_input$se_NC,
                    Subtype = Subtype,
                    hasClinical = hasClinical
                )
            )
        }

        if (method == "cIHC.itr") {
            ## try iterative IHC
            message(method, " is running!")
            return(
                BS_cIHC.itr(
                    data_input$se_NC,
                    Subtype = Subtype,
                    hasClinical = hasClinical
                )
            )
        }

        if (method == "PCAPAM50") {
            message(method, " is running!")
            result <- tryCatch(
                {
                    # Attempt to run PCAPAM50
                    BS_PCAPAM50(
                        data_input$se_NC,
                        Subtype = Subtype,
                        hasClinical = hasClinical
                    )
                },
                error = function(e) {
                    # Error handling
                    warning("PCAPAM50 failed in this iteration. Error: ")
                    return(NULL) # Return NULL or a dummy tibble with NAs
                }
            )

            # Return result (or NULL if failed)
            return(result)
        }

        if (method == "ssBC") {
            message(method, " is running!")

            if ("TN" %in% colnames(pheno) && cohort.select == "TNBC") {
                res_ssBC <- BS_ssBC(
                    data_input$se_NC,
                    s = "TN",
                    Subtype = Subtype,
                    hasClinical = hasClinical
                )
            } else {
                if (!is.null(samples_ER.icd) & length(samples_ERHER2.icd) < nrow(pheno)) {
                    res_ssBC <- BS_ssBC(
                        data_input$se_NC[, samples_ER.icd],
                        s = "ER",
                        Subtype = Subtype,
                        hasClinical = hasClinical
                    )
                } else {
                    res_ssBC <- BS_ssBC(
                        data_input$se_NC,
                        s = "ER",
                        Subtype = Subtype,
                        hasClinical = hasClinical
                    )
                }
            }

            ## Keep patients' results for AUTO mode
            if (!is.null(samples_ER.icd) &&
                length(samples_ER.icd) < nrow(pheno)) {
                unprocessed_patients <- base::setdiff(pheno$PatientID, samples_ER.icd)

                # Create NA-filled data frame for unprocessed patients with matching structure
                na_df <- data.frame(
                    PatientID = unprocessed_patients,
                    matrix(
                        nrow = length(unprocessed_patients),
                        ncol = ncol(res_ssBC$BS.all) - 1, # Subtract 1 for PatientID column
                        dimnames = list(NULL, colnames(res_ssBC$BS.all)[-1])
                    ),
                    row.names = unprocessed_patients
                )

                res_ssBC$BS.all <- rbind(res_ssBC$BS.all, na_df)
                res_ssBC$BS.all <- res_ssBC$BS.all[pheno$PatientID, ]
            }

            return(res_ssBC)
        }

        if (method == "ssBC.v2") {
            message(method, " is running!")

            if ("TN" %in% colnames(pheno) && cohort.select == "TNBC") {
                res_ssBC.v2 <- BS_ssBC(
                    data_input$se_NC,
                    s = "TN.v2",
                    Subtype = Subtype,
                    hasClinical = hasClinical
                )
            } else {
                if (!is.null(samples_ERHER2.icd) & length(samples_ERHER2.icd) < nrow(pheno)) {
                    res_ssBC.v2 <- BS_ssBC(
                        data_input$se_NC[, samples_ERHER2.icd],
                        s = "ER.v2",
                        Subtype = Subtype,
                        hasClinical = hasClinical
                    )
                } else {
                    res_ssBC.v2 <- BS_ssBC(
                        data_input$se_NC,
                        s = "ER.v2",
                        Subtype = Subtype,
                        hasClinical = hasClinical
                    )
                }
            }

            ## Keep patients' results for AUTO mode
            if (!is.null(samples_ERHER2.icd) &&
                length(samples_ERHER2.icd) < nrow(pheno)) {
                unprocessed_patients <- base::setdiff(pheno$PatientID, samples_ERHER2.icd)


                # Create NA-filled data frame for unprocessed patients with matching structure
                na_df <- data.frame(
                    PatientID = unprocessed_patients,
                    matrix(
                        nrow = length(unprocessed_patients),
                        ncol = ncol(res_ssBC.v2$BS.all) - 1, # Subtract 1 for PatientID column
                        dimnames = list(NULL, colnames(res_ssBC.v2$BS.all)[-1])
                    ),
                    row.names = unprocessed_patients
                )

                res_ssBC.v2$BS.all <- rbind(res_ssBC.v2$BS.all, na_df)
                res_ssBC.v2$BS.all <- res_ssBC.v2$BS.all[pheno$PatientID, ]
            }


            return(res_ssBC.v2)
        }

        if (method == "AIMS") {
            message(method, " is running!")
            data_env <- new.env(parent = emptyenv())
            data("BreastSubtypeRobj", envir = data_env, package = "BreastSubtypeR")
            BreastSubtypeRobj <- data_env[["BreastSubtypeRobj"]]


            res_AIMS <- BS_AIMS(data_input$se_SSP)

            res_AIMS$BS.all <- data.frame(
                PatientID = rownames(res_AIMS$cl),
                BS = res_AIMS$cl[, 1]
            )

            if (Subtype) {
                res_AIMS$BS.all$BS.Subtype <- res_AIMS$BS.all$BS
            }
            return(res_AIMS)
        }

        if (method == "sspbc") {
            message(method, " is running!")

            res_sspbc <- BS_sspbc(
                se_obj = data_input$se_SSP,
                ssp.name = "ssp.pam50"
            )

            BS.all <- data.frame(
                PatientID = rownames(res_sspbc),
                BS = res_sspbc[, 1],
                row.names = rownames(res_sspbc)
            )

            if (Subtype) {
                res_sspbc.Subtype <- BS_sspbc(
                    se_obj = data_input$se_SSP,
                    ssp.name = "ssp.subtype"
                )
                BS.all$BS.Subtype <- res_sspbc.Subtype[, 1]
            }

            return(list(BS.all = BS.all))
        }

        stop("Unknown method: ", method)
    })

    ## Map method names and align outputs by PatientID (robust to order/missing)
    names(results) <- methods

    # Prefer se_NC samples; fall back to se_SSP if needed
    samples_NC <- if (!is.null(data_input$se_NC)) colnames(SummarizedExperiment::assay(data_input$se_NC)) else character(0)
    samples_SSP <- if (!is.null(data_input$se_SSP)) colnames(SummarizedExperiment::assay(data_input$se_SSP)) else character(0)
    samples <- if (length(samples_NC)) samples_NC else samples_SSP

    # Hold per-method calls; start empty and fill by matching PatientID
    res_subtypes <- data.table::data.table(row_id = samples)
    if (Subtype) {
        res_subtypes.Subtype <- data.table::data.table(row_id = samples)
    }

    for (method in methods) {
        x <- results[[method]]
        # default (all NA) if the method failed or returned nothing usable
        vec5 <- rep(NA_character_, length(samples))
        vec4 <- rep(NA_character_, length(samples))

        if (!is.null(x) && is.data.frame(x$BS.all) && "PatientID" %in% names(x$BS.all)) {
            i <- match(samples, x$BS.all$PatientID) # align by PatientID
            # 5-class column (may be absent for some SSPs or error cases)
            if ("BS" %in% names(x$BS.all)) {
                v <- x$BS.all$BS
                # coerce to character to avoid factor levels surprises
                if (is.factor(v)) v <- as.character(v)
                vec5 <- v[i]
            }
            # 4-class column (may not exist for some methods or when Subtype = FALSE)
            if (Subtype && "BS.Subtype" %in% names(x$BS.all)) {
                v <- x$BS.all$BS.Subtype
                if (is.factor(v)) v <- as.character(v)
                vec4 <- v[i]
            }
        }

        data.table::set(res_subtypes, j = method, value = vec5)
        if (Subtype) data.table::set(res_subtypes.Subtype, j = method, value = vec4)
    }

    # Convert to data.frame and drop ONLY row_id once
    res_subtypes <- as.data.frame(res_subtypes, stringsAsFactors = FALSE, check.names = FALSE)
    rownames(res_subtypes) <- res_subtypes$row_id
    res_subtypes$row_id <- NULL

    if (Subtype) {
        res_subtypes.Subtype <- as.data.frame(res_subtypes.Subtype, stringsAsFactors = FALSE, check.names = FALSE)
        rownames(res_subtypes.Subtype) <- res_subtypes.Subtype$row_id
        res_subtypes.Subtype$row_id <- NULL
    }

    ## entropy index (compute on method columns; don't drop any)
    if (ncol(res_subtypes) > 0) {
        res_subtypes$entropy <- apply(res_subtypes, 1, get_entropy)
    } else {
        res_subtypes$entropy <- NA_real_
    }

    if (Subtype) {
        # AIMS is 5-class; for 4-class summary set Normal to NA only if the column exists
        if ("AIMS" %in% colnames(res_subtypes.Subtype)) {
            ix <- which(res_subtypes.Subtype$AIMS == "Normal")
            if (length(ix)) res_subtypes.Subtype$AIMS[ix] <- NA
        }
        if (ncol(res_subtypes.Subtype) > 0) {
            res_subtypes.Subtype$entropy <- apply(res_subtypes.Subtype, 1, get_entropy)
        } else {
            res_subtypes.Subtype$entropy <- NA_real_
        }
    }

    if (Subtype) {
        res <- list(
            res_subtypes = res_subtypes,
            res_subtypes.Subtype = res_subtypes.Subtype,
            results = results
        )
    } else {
        res <- list(res_subtypes = res_subtypes, results = results)
    }
    return(res)
}
