#' @title OSLO2EMIT0obj: Example dataset (OSLO2-EMIT0 cohort subset)
#'
#' @description
#' Example object derived from the OSLO2-EMIT0 cohort (Staaf et al., 2022).
#' Includes a subset of normalized expression data, clinical metadata, feature
#' annotations, and example outputs from \code{Mapping()} and \code{BS_Multi()}.
#'
#' @docType data
#' @usage data("OSLO2EMIT0obj")
#'
#' @format A list with:
#' \describe{
#'   \item{\code{se_obj}}{A \code{SummarizedExperiment} containing a subset of the
#'   log2-transformed, normalized expression matrix (log2(FPKM+0.1)) with \code{colData} clinical
#'   metadata and row-level feature annotations.}
#'   \item{\code{data_input}}{Example output structure produced by \code{Mapping()}.}
#'   \item{\code{res}}{Example results from \code{BS_Multi()} run in \emph{AUTO} mode.}
#' }
#'
#' @references
#' Staaf J, Häkkinen J, Hegardt C, Saal LH, Kimbung S, Hedenfalk I, et al.
#' RNA sequencing-based single sample predictors of molecular subtype and risk of
#' recurrence for clinical assessment of early-stage breast cancer.
#' \emph{NPJ Breast Cancer}. 2022;8(1):27. https://doi.org/10.1038/s41523-022-00465-3
#'
#' @examples
#' library(BreastSubtypeR)
#' data("OSLO2EMIT0obj")
"OSLO2EMIT0obj"


#' @title TCGABRCAobj: Example dataset (TCGA-BRCA subset)
#'
#' @description
#' Example object derived from TCGA-BRCA. Includes a subset of normalized metadata
#' + raw counts (as a SummarizedExperiment), and example outputs from \code{Mapping()}
#' and \code{BS_Multi()} to facilitate runnable examples.
#'
#' @docType data
#' @usage data("TCGABRCAobj")
#'
#' @format A list with:
#' \describe{
#'   \item{\code{se_obj}}{A \code{SummarizedExperiment} containing the integer raw-count matrix
#'   (top 5,000 variable genes), \code{rowData} with \code{probe}, \code{SYMBOL}, \code{ENTREZID},
#'   \code{Length}, and \code{colData} with \code{PatientID}, \code{ER}, \code{PR}, \code{HER2}.}
#'   \item{\code{data_input}}{Example \code{Mapping()} output created from \code{se_obj}.}
#'   \item{\code{res}}{Example \code{BS_Multi()} results (e.g., run in \emph{AUTO} mode).}
#' }
#'
#' @source
#' The Cancer Genome Atlas (TCGA) BRCA via GDC; counts summarized with \code{recount3};
#' clinical data retrieved with \code{TCGAbiolinks}.
#'
#' @references
#' The Cancer Genome Atlas Network.
#' Comprehensive molecular portraits of human breast tumors.
#' \emph{Nature}. 2012;490(7418):61–70. https://doi.org/10.1038/nature11412
#'
#' Colaprico A, Silva TC, Olsen C, Garofano L, Cava C, Garolini D, et al.
#' TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data.
#' \emph{Nucleic Acids Res}. 2016;44(8):e71. https://doi.org/10.1093/nar/gkv1507
#'
#' Collado-Torres L, Nellore A, Kammers K, Ellis SE, Taub MA, Hansen KD, et al.
#' Reproducible RNA-seq analysis using recount2.
#' \emph{Nat Biotechnol}. 2017;35(4):319–321. https://doi.org/10.1038/nbt.3838
#'
#' @examples
#' library(BreastSubtypeR)
#' data("TCGABRCAobj")
#' names(TCGABRCAobj)
#' # str(TCGABRCAobj$se_obj); head(colData(TCGABRCAobj$se_obj))
"TCGABRCAobj"


#' @title BreastSubtypeRobj: Resources for NC-based methods
#'
#' @description
#' List of reference resources required by nearest-centroid (NC) subtyping
#' methods: platform medians, centroids, signatures, subgroup quantiles, and
#' metadata from the UNC232 training cohort.
#'
#' @docType data
#' @usage data("BreastSubtypeRobj")
#'
#' @format A list with:
#' \describe{
#'   \item{\code{medians}}{Matrix/data frame of platform-specific medians
#'   for \strong{11} expression/sequencing platforms, derived as described in
#'   Picornell et al. (2019). Platform columns include:
#'   \code{nCounter},
#'   \code{totalRNA.FFPE.20151111}, \code{RNAseq.Freeze.20120907}, \code{RNAseq.V2}, \code{RNAseq.V1},
#'   \code{GC.4x44Kcustom}, \code{Agilent_244K}, \code{commercial_1x44k_postMeanCollapse_WashU}, \code{commercial_4x44k_postMeanCollapse_WashU_v2},
#'   \code{htp1.5_WU_update}, \code{arrayTrain_postMeanCollapse}.}
#'   \item{\code{centroid}}{PAM50 centroids used by \code{parker.original}.}
#'   \item{\code{genes.sig50}}{Data frame of the 50 PAM50 genes with a proliferation flag.}
#'   \item{\code{ssBC.subgroupQuantile}}{Original and updated subgroup-specific
#'   quantiles used by \code{ssBC} and \code{ssBC.v2}; see Details.}
#'   \item{\code{genes.signature}}{Marker genes used across NC- and SSP-based methods.}
#'   \item{\code{UNC232}}{Summary data for the UNC232 training cohort.}
#'   \item{\code{platform.UNC232}}{Platform annotation for UNC232.}
#' }
#'
#' @details
#' The original `ER` and `TN` methods retain the reference quantiles for
#' subgroup-specific gene centering introduced by Zhao et al. (2015).
#' Fernandez-Martinez et al. (2020) described a HER2/ER-specific application.
#' The five updated columns match the public `SIGMA.txt` reference associated
#' with Fernandez-Martinez et al. (2024), across all 50 PAM50 genes:
#' `ER.v2` uses `ERneg_HER2neg`, `ERpos_HER2neg`, `HER2pos_ERneg` and
#' `HER2pos_ERpos`; `TN.v2` uses `TNBC`.
#'
#' @source
#' [Perou Lab normalization repository](https://github.com/perou-lab/PAM50_ER_HER2_normalization).
#' The updated columns were verified against
#' [SIGMA.txt at commit 9ba84eb](https://github.com/perou-lab/PAM50_ER_HER2_normalization/blob/9ba84ebc131a58a18446e235c9422d9802f56671/SIGMA.txt).
#'
#' @references
#' Parker JS, Mullins M, Cheung MCU, Leung S, Voduc D, et al.
#' Supervised risk predictor of breast cancer based on intrinsic subtypes.
#' \emph{J Clin Oncol}. 2009;27(8):1160–1167. https://doi.org/10.1200/JCO.2008.18.1370
#'
#' Zhao X, Rødland EA, Tibshirani R, Plevritis S.
#' Molecular subtyping for clinically defined breast cancer subgroups.
#' \emph{Breast Cancer Res}. 2015;17:29. https://doi.org/10.1186/s13058-015-0520-4
#'
#' Fernandez-Martinez A, Krop IE, Hillman DW, et al.
#' Survival, Pathologic Response, and Genomics in CALGB 40601 (Alliance),
#' a Neoadjuvant Phase III Trial of Paclitaxel-Trastuzumab With or Without
#' Lapatinib in HER2-Positive Breast Cancer.
#' \emph{J Clin Oncol}. 2020;38(35):4184-4193. https://doi.org/10.1200/JCO.20.01276
#'
#' Fernandez-Martinez A, Rediti M, Tang G, et al.
#' Tumor Intrinsic Subtypes and Gene Expression Signatures in Early-Stage
#' ERBB2/HER2-Positive Breast Cancer: A Pooled Analysis of CALGB 40601,
#' NeoALTTO, and NSABP B-41 Trials.
#' \emph{JAMA Oncol}. 2024;10(5):603-611. https://doi.org/10.1001/jamaoncol.2023.7304
#'
#' Picornell AC, Echavarria I, Alvarez E, López-Tarruella S, Jerez Y, Hoadley K, et al.
#' Breast cancer PAM50 signature: correlation and concordance between RNA-seq and digital multiplexed gene expression technologies in a TNBC series.
#' \emph{BMC Genomics}. 2019;20(1):452. https://doi.org/10.1186/s12864-019-5849-0
#'
#' @examples
#' library(BreastSubtypeR)
#' data("BreastSubtypeRobj")
"BreastSubtypeRobj"


#' @title AIMSmodel: Model object for AIMS
#'
#' @description
#' Model definition for AIMS consisting of 100 pairwise rules and a Naive Bayes
#' classifier (via \pkg{e1071}) as described by Paquet & Hallett (2015).
#'
#' @details
#' The 100 rules are of the form “EntrezID gene A < EntrezID gene B”. A subset
#' of \code{k} rules (typically 20) is used within a Naive Bayes classifier to
#' assign subtypes (Basal-like, HER2-enriched, LumA, LumB, Normal-like) on a
#' per-sample basis.
#'
#' @docType data
#' @usage data("AIMSmodel")
#'
#' @return
#' \item{all.pairs}{Character vector of the 100 AIMS rules (EntrezID comparisons).}
#' \item{k}{Integer; number of optimal rules (commonly 20).}
#' \item{one.vs.all.tsp}{Naive Bayes classifier object used with the rules.}
#' \item{selected.pairs.list}{Rules ranked by discriminative power per subtype.}
#'
#' @references
#' Paquet ER, Hallett MT.
#' Absolute assignment of breast cancer intrinsic molecular subtype.
#' \emph{J Natl Cancer Inst}. 2015;107(1):dju357. https://doi.org/10.1093/jnci/dju357
#'
#' @examples
#' library(BreastSubtypeR)
#' data("AIMSmodel")
"AIMSmodel"


#' @title sspbc.models: Short names for 11 SSPBC predictors
#'
#' @description
#' List of 11 single-sample predictor (SSP) models from Staaf et al. (2022),
#' indexed by short names used by \code{sspbc}.
#'
#' @docType data
#' @usage data("sspbc.models")
#'
#' @details
#' Names correspond to short model identifiers. The contents are identical to
#' \code{sspbc.models.fullname}, which uses full model names.
#'
#' @return
#' \item{sspbc.models}{Named list of 11 SSP models used by \code{sspbc}.}
#'
#' @references
#' Staaf J, Häkkinen J, Hegardt C, Saal LH, Kimbung S, Hedenfalk I, et al.
#' RNA sequencing-based single sample predictors of molecular subtype and risk of
#' recurrence for clinical assessment of early-stage breast cancer.
#' \emph{NPJ Breast Cancer}. 2022;8(1):27. https://doi.org/10.1038/s41523-022-00465-3
#'
#' @examples
#' library(BreastSubtypeR)
#' data("sspbc.models")
"sspbc.models"


#' @title sspbc.models.fullname: Full names for 11 SSPBC predictors
#'
#' @description
#' List of the same 11 SSP models (Staaf et al., 2022) indexed by full model names.
#'
#' @docType data
#' @usage data("sspbc.models.fullname")
#'
#' @details
#' Identical content to \code{sspbc.models} but with full model names as list keys.
#'
#' @return
#' \item{sspbc.models.fullname}{Named list of 11 SSP models used by \code{sspbc}.}
#'
#' @references
#' Staaf J, Häkkinen J, Hegardt C, Saal LH, Kimbung S, Hedenfalk I, et al.
#' RNA sequencing-based single sample predictors of molecular subtype and risk of
#' recurrence for clinical assessment of early-stage breast cancer.
#' \emph{NPJ Breast Cancer}. 2022;8(1):27. https://doi.org/10.1038/s41523-022-00465-3
#'
#' @examples
#' library(BreastSubtypeR)
#' data("sspbc.models.fullname")
"sspbc.models.fullname"


#' @title Gene.ID.ann: Gene annotation table
#'
#' @description
#' Annotation table for GENCODE Human Release 27 genes (\code{Gene.ID}) used by
#' StringTie summarisation. Includes HGNC, EntrezGene, and RefSeq identifiers
#' derived from GENCODE v27 metadata.
#'
#' @docType data
#' @usage data("Gene.ID.ann")
#'
#' @details
#' Used by internal SSP application functions to translate identifiers prior to
#' classification with SSP models.
#'
#' @return
#' \item{Gene.ID.ann}{Data frame of annotations for GENCODE v27 genes.}
#'
#' @references
#' Staaf J, Häkkinen J, Hegardt C, Saal LH, Kimbung S, Hedenfalk I, et al.
#' RNA sequencing-based single sample predictors of molecular subtype and risk of
#' recurrence for clinical assessment of early-stage breast cancer.
#' \emph{NPJ Breast Cancer}. 2022;8(1):27. https://doi.org/10.1038/s41523-022-00465-3
#'
#' @examples
#' library(BreastSubtypeR)
#' data("Gene.ID.ann")
"Gene.ID.ann"
