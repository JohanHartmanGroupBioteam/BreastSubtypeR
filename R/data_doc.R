#' @title OSLO2EMIT0obj: Example dataset (OSLO2-EMIT0 cohort subset)
#'
#' @description
#' Example object derived from the OSLO2-EMIT0 cohort (Staaf et al., 2022).
#' Includes protocol-adjusted expression, historical clinical metadata and
#' feature annotations, and regenerated \code{Mapping()} and \code{BS_Multi()}
#' results. It retains the existing public object name and principal list schema.
#'
#' @docType data
#' @usage data("OSLO2EMIT0obj")
#'
#' @format A list with:
#' \describe{
#'   \item{\code{se_obj}}{A \code{SummarizedExperiment} with 569 unique model
#'   genes and 102 samples. The assay is \code{log2(FPKM + 1)} calculated from
#'   library-protocol-adjusted FPKM. The original 356 rows, all 36,312 values,
#'   feature annotations and sample metadata/order are retained exactly;
#'   213 genes are appended in adjusted-source order. Clinical metadata are in
#'   \code{colData}, feature annotations in \code{rowData}, and reconstruction
#'   provenance in \code{metadata(se_obj)$issue199_reconstruction}.}
#'   \item{\code{data_input}}{Regenerated \code{Mapping()} output containing
#'   \code{se_NC} (50 genes by 102 samples) and \code{se_SSP}
#'   (339 genes by 102 samples).}
#'   \item{\code{res}}{Regenerated \code{BS_Multi()} AUTO results, containing
#'   \code{res_subtypes} and \code{results}; the summary uses five classes
#'   and clinical ROR scoring is disabled.}
#' }
#'
#' @details
#' The 569-gene union covers the bundled PAM50, AIMS and SSPBC
#' \code{ssp.subtype}, \code{ssp.pam50} and \code{ssp.ror} model identifiers.
#' Default \code{Mapping()} still selects its existing NC/SSP registries;
#' this object does not add \code{se_SSP_ROR} or establish complete native
#' SSP-ROR preparation. SSPBC models were developed using unadjusted StringTie
#' FPKM. Prediction equivalence with this protocol-adjusted input has not
#' been established. Current nonraw SSP preprocessing uses \code{2^x},
#' retaining the \code{+1} offset for this input. No unadjusted matrix is included.
#'
#' Complete results were regenerated under BreastSubtypeR 1.5.12 at source
#' commit \code{cb7a7c70e0f94e1690d15ac2ec2b464a40934f36}.
#' Mapping used \code{RawCounts = FALSE}, \code{method = "max"} and
#' \code{impute = TRUE}, with seed 1992026. Prediction used
#' \code{methods = "AUTO"}, \code{Subtype = FALSE},
#' \code{hasClinical = FALSE} and \code{er_subsample = "target"},
#' with seed 20260831 and unchanged method defaults. RNG kinds were
#' Mersenne-Twister, Inversion and Rejection.
#'
#' Under that same implementation, expansion leaves both mapped inputs and
#' every complete method result identical to recomputation from the historical
#' 356-row input. Subtype calls and entropy also match the historical stored
#' object. Some stored historical ROR, Confidence and output-row details differ
#' from current recomputation, including the original proliferation-score
#' standardization correction (#198); these are not expansion effects.
#' The exact original save-event environment and seed remain unknown.
#' Reconstruction scripts and versioned manifests are in
#' \code{data-raw/issue199} in the package source repository.
#'
#' @source
#' \code{OSLO2EMIT0.103.genematrix_noNeg.Rdata}, the verified
#' library-protocol-adjusted FPKM artifact in Mendeley Data, version 4,
#' file identifier \code{216a4cc7-e688-4c72-9513-2ad081e2c801}.
#' \doi{10.17632/yzxtxn4nmd.4}
#' The source contains 19,675 features and 103 sample identifiers. The original
#' 102 sample identifiers/order are retained; \code{OSLO2EMIT0.073} is absent.
#' Its historical omission rationale and the original feature-selection
#' recipe are unresolved. Historical values were copied without recomputation:
#' independently transforming the source differs by one ULP at 51 original
#' cells. The new rows use the verified source transform consistently.
#' MYO10 retains the verified Entrez identifier 4651.
#'
#' @references
#' Vallon-Christersson J (2025). RNA Sequencing-Based Single Sample Predictors
#' of Molecular Subtype and Risk of Recurrence for Clinical Assessment of
#' Early-Stage Breast Cancer. Mendeley Data, V4.
#' \doi{10.17632/yzxtxn4nmd.4}
#'
#' Staaf J, Hakkinen J, Hegardt C, Saal LH, Kimbung S, Hedenfalk I, et al.
#' RNA sequencing-based single sample predictors of molecular subtype and risk of
#' recurrence for clinical assessment of early-stage breast cancer.
#' \emph{npj Breast Cancer}. 2022;8:94. https://doi.org/10.1038/s41523-022-00465-3
#'
#' @examples
#' library(BreastSubtypeR)
#' data("OSLO2EMIT0obj")
#' dim(OSLO2EMIT0obj$se_obj) # 569 genes, 102 samples
#' names(OSLO2EMIT0obj) # se_obj, data_input, res
"OSLO2EMIT0obj"


#' @title TCGABRCAobj: Example dataset (TCGA-BRCA subset)
#'
#' @description
#' Example object containing recount3 TCGA-BRCA base-pair coverage sums,
#' historical clinical metadata and feature annotations, and regenerated
#' \code{Mapping()} and \code{BS_Multi()} outputs. The existing public object
#' name and principal list schema are retained.
#'
#' @docType data
#' @usage data("TCGABRCAobj")
#'
#' @format A list with:
#' \describe{
#'   \item{\code{se_obj}}{A \code{SummarizedExperiment} with 5,350 versioned
#'   Ensembl features and 100 samples. The integer assay retains the recount3
#'   G026 base-pair coverage-sum convention, not transformed read counts, TPM
#'   or FPKM. All original 5,000 rows and 500,000 values are unchanged, followed
#'   by 350 missing model genes in C-locale lexicographic Entrez-ID order.
#'   \code{rowData} contains \code{probe}, \code{SYMBOL}, \code{ENTREZID} and
#'   \code{Length}; \code{colData} contains \code{PatientID}, \code{ER},
#'   \code{PR} and \code{HER2}. Original annotations, sample metadata and
#'   order are preserved. Reconstruction provenance is in
#'   \code{metadata(se_obj)$issue199_reconstruction}.}
#'   \item{\code{data_input}}{Regenerated \code{Mapping()} output containing
#'   \code{se_NC} (50 genes by 100 samples) and \code{se_SSP}
#'   (339 genes by 100 samples).}
#'   \item{\code{res}}{Regenerated \code{BS_Multi()} AUTO results, containing
#'   \code{res_subtypes} and \code{results}; the summary uses five classes
#'   and clinical ROR scoring is disabled.}
#' }
#'
#' @details
#' The historical 5,000-gene selection is retained without variance reranking.
#' Added genes complete the verified 569-gene PAM50/AIMS/SSPBC model union.
#' All 535,000 counts agree with the original G026 source. Gene lengths come
#' from the G026 gene-sum GTF score field, not genomic spans. Versioned
#' Ensembl mappings are retained, including the documented Entrez 9906
#' recovery to \code{ENSG00000215790.6} (SLC35E2; length 3264).
#' Historical missing or repeated Entrez annotations remain unchanged;
#' the 5,350 versioned feature identifiers are unique.
#' No \code{recount3::transform_counts()} conversion was applied.
#'
#' Complete results were regenerated under BreastSubtypeR 1.5.12 at source
#' commit \code{cb7a7c70e0f94e1690d15ac2ec2b464a40934f36}.
#' Mapping used \code{RawCounts = TRUE}, \code{method = "max"} and
#' \code{impute = TRUE}, with seed 1992026. Prediction used
#' \code{methods = "AUTO"}, \code{Subtype = FALSE},
#' \code{hasClinical = FALSE} and \code{er_subsample = "target"},
#' with seed 20260831 and unchanged method defaults. RNG kinds were
#' Mersenne-Twister, Inversion and Rejection.
#'
#' Expansion changes full-input library sizes and upper-quartile factors
#' before Mapping selects its existing registries. Relative to recomputation
#' of the historical input under the same implementation, 64 method calls
#' change across 40 samples: genefu.robust 8, ssBC 4, ssBC.v2 6, cIHC 5,
#' cIHC.itr 6, PCAPAM50 8, AIMS 16 and sspbc 11. This comparison separates
#' dataset effects from earlier implementation and summary-schema changes.
#' The original stored object has additional differences, including
#' 15 ssBC.v2 calls and one PCAPAM50 call on current-code recomputation;
#' exclusive historical build/seed attribution is unresolved.
#'
#' The historical incomplete mapped input (37 NC and 178 SSP genes) remains
#' an immutable test fixture. The current input maps to all 50 PAM50 genes
#' and both proliferation panels (11 and 18 genes). Complete source-gene
#' coverage does not establish native SSP-ROR preparation or equivalence to
#' unadjusted StringTie FPKM training input. No \code{se_SSP_ROR} is provided.
#' Reconstruction scripts and versioned manifests are in
#' \code{data-raw/issue199} in the package source repository.
#'
#' @source
#' TCGA-BRCA recount3 G026 gene sums
#' (\code{tcga.gene_sums.BRCA.G026.gz}) and gene-sum annotation
#' (\code{human.gene_sums.G026.gtf.gz}). The complete count source contains
#' 63,856 features and 1,256 samples; the original 100-sample selection is
#' preserved, with source UUIDs resolved against historical count vectors.
#' Clinical metadata originally retrieved with \code{TCGAbiolinks} are
#' retained without reacquisition or recoding.
#'
#' @references
#' The Cancer Genome Atlas Network.
#' Comprehensive molecular portraits of human breast tumors.
#' \emph{Nature}. 2012;490(7418):61-70. https://doi.org/10.1038/nature11412
#'
#' Colaprico A, Silva TC, Olsen C, Garofano L, Cava C, Garolini D, et al.
#' TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data.
#' \emph{Nucleic Acids Res}. 2016;44(8):e71. https://doi.org/10.1093/nar/gkv1507
#'
#' recount3 project: \url{https://rna.recount.bio/}
#'
#' @examples
#' library(BreastSubtypeR)
#' data("TCGABRCAobj")
#' dim(TCGABRCAobj$se_obj) # 5350 genes, 100 samples
#' names(TCGABRCAobj) # se_obj, data_input, res
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
#' Parker JS, Mullins M, Cheang MCU, Leung S, Voduc D, et al.
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
#' The 100 rules are of the form \code{EntrezID gene A < EntrezID gene B}.
#' The bundled \code{k = 20} identifies the fitted-model setting; its
#' \code{all.pairs} vector and fitted naive Bayes tables contain all 100
#' rule predictors, involving 151 unique genes. Thus \code{k} is not the
#' total rule denominator. The model assigns Basal-like, HER2-enriched,
#' LumA, LumB and Normal-like subtypes independently for each sample.
#' See \code{\link{InputCoverage}} for descriptive input reporting.
#'
#' @docType data
#' @usage data("AIMSmodel")
#'
#' @return
#' \item{all.pairs}{Character vector of the 100 AIMS rules (EntrezID comparisons).}
#' \item{k}{Integer fitted-model setting (20 in the bundled model), not the total rule count.}
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
#' Staaf J, Hakkinen J, Hegardt C, Saal LH, Kimbung S, Hedenfalk I, et al.
#' RNA sequencing-based single sample predictors of molecular subtype and risk of
#' recurrence for clinical assessment of early-stage breast cancer.
#' \emph{npj Breast Cancer}. 2022;8:94. https://doi.org/10.1038/s41523-022-00465-3
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
#' Staaf J, Hakkinen J, Hegardt C, Saal LH, Kimbung S, Hedenfalk I, et al.
#' RNA sequencing-based single sample predictors of molecular subtype and risk of
#' recurrence for clinical assessment of early-stage breast cancer.
#' \emph{npj Breast Cancer}. 2022;8:94. https://doi.org/10.1038/s41523-022-00465-3
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
#' Staaf J, Hakkinen J, Hegardt C, Saal LH, Kimbung S, Hedenfalk I, et al.
#' RNA sequencing-based single sample predictors of molecular subtype and risk of
#' recurrence for clinical assessment of early-stage breast cancer.
#' \emph{npj Breast Cancer}. 2022;8:94. https://doi.org/10.1038/s41523-022-00465-3
#'
#' @examples
#' library(BreastSubtypeR)
#' data("Gene.ID.ann")
"Gene.ID.ann"
