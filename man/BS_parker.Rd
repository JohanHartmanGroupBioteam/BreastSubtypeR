% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BreastSubtypeR_methods.R
\name{BS_parker}
\alias{BS_parker}
\title{Original Parker Intrinsic Subtyping (BS_parker)}
\usage{
BS_parker(
  se_obj,
  calibration = "None",
  internal = NA,
  external = NA,
  medians = NA,
  Subtype = FALSE,
  hasClinical = FALSE
)
}
\arguments{
\item{se_obj}{A \code{SummarizedExperiment} object containing:
\itemize{
\item \strong{Assay data}: A gene expression matrix with genes (Gene Symbol) in rows
and samples in columns. This matrix should be log-transformed.
\item \strong{Column metadata} (\code{colData}): Optional clinical information.
}}

\item{calibration}{The calibration method to use. Options include:
\itemize{
\item \code{"None"}: No calibration is applied.
\item \code{"Internal"}: Calibration is performed using internal strategies.
See the \code{internal} parameter for details.
\item \code{"External"}: Calibration is performed using external medians.
See the \code{external} parameter for details.
}}

\item{internal}{The internal calibration strategy to apply when \code{calibration = "Internal"}. Options are:
\itemize{
\item \code{"medianCtr"} (default): Median-centered calibration.
\item \code{"meanCtr"}:
Mean-centered calibration (aligned with \code{genefu.scale} function).
\item \code{"qCtr"}:
Quantile-based calibration (aligned with \code{genefu.robust} function).
}}

\item{external}{Specify the platform name (i.e., column name) for external
medians calculated from the training cohort. To use user-provided medians,
set this parameter to \code{"Given.mdns"} and provide the values via the
\code{medians} parameter.}

\item{medians}{A matrix or table of user-provided medians. Required if
\code{external = "Given.mdns"}. The first column should list 50 genes, and the
second column should provide the corresponding median values.}

\item{Subtype}{Logical. If \code{TRUE}, the function predicts four subtypes by
excluding the Normal-like subtype.}

\item{hasClinical}{Logical. If \code{TRUE}, the function uses clinical data from
the \code{pheno} table. Required columns include:
\itemize{
\item \code{"TSIZE"}: Tumor size.
\item \code{"NODE"}: Lymph node status.
}}
}
\value{
A list containing the intrinsic subtypes assigned using the
Parker-based method.
}
\description{
This function predicts breast cancer intrinsic subtypes using
the original Parker et al., 2019 method, as well as some variations of the
original approach.
}
\examples{
data("OSLO2EMIT0obj")
res <- BS_parker(
    se_obj = OSLO2EMIT0obj$data_input$se_NC,
    calibration = "Internal",
    internal = "medianCtr",
    Subtype = FALSE,
    hasClinical = FALSE
)

}
\references{
\itemize{
\item Bernard PS, Parker JS, Mullins M, Cheung MCU, Leung S, Voduc D, et al.
\emph{Supervised risk predictor of breast cancer based on intrinsic subtypes}.
Journal of Clinical Oncology. 2009;27(8).
https://doi.org/10.1200/JCO.2008.18.1370
\item Gendoo DMA, Ratanasirigulchai N, Schröder MS, Paré L, Parker JS, Prat A, et
al. \emph{Genefu: An R/Bioconductor package for computation of gene
expression-based signatures in breast cancer}. Bioinformatics. 2016;32(7).
https://doi.org/10.1093/bioinformatics/btv693
}
}
