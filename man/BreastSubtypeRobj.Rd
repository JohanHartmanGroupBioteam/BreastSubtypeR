% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_doc.R
\docType{data}
\name{BreastSubtypeRobj}
\alias{BreastSubtypeRobj}
\title{BreastSubtypeRobj: Data for NC-based Methods}
\format{
A list with the following elements:
\describe{
\item{\code{medians}}{A matrix of medians prepared for nine sequencing platforms.}
\item{\code{centroid}}{The centroids provided by the \code{parker.original} method.}
\item{\code{genes.sig50}}{A data frame of 50 genes used in NC-based methods, including their proliferation information.}
\item{\code{ssBC.subgroupQuantile}}{Subgroup medians prepared by the \code{ssBC} method.}
\item{\code{genes.signature}}{A collection of genes used in both NC-based and SSP-based methods.}
\item{\code{UNC232}}{Data from the UNC232 training cohort.}
\item{\code{platform.UNC232}}{The sequencing platform used for the UNC232 training cohort.}
}
}
\usage{
data("BreastSubtypeRobj")
}
\description{
A list object containing the data required for nearest-centroid (NC)-based
molecular subtyping methods. This includes platform medians, centroids, gene
signatures, and subgroup quantiles, as well as information from the UNC232
training cohort.
}
\examples{
library(BreastSubtypeR)
data("BreastSubtypeRobj")

}
\references{
\itemize{
\item Bernard PS, Parker JS, Mullins M, Cheung MCU, Leung S, Voduc D, et al.
\emph{Supervised risk predictor of breast cancer based on intrinsic subtypes}.
Journal of Clinical Oncology. 2009;27(8).
https://doi.org/10.1200/JCO.2008.18.1370
\item Zhao X, Rodland EA, Tibshirani R, Plevritis S. \emph{Molecular
subtyping for clinically defined breast cancer subgroups.} Breast Cancer
Research. 2015;17(1). https://doi.org/10.1186/s13058-015-0520-4
\item Fernandez-Martinez A, Krop IE, Hillman DW, Polley MY, Parker JS, Huebner L,
et al. \emph{Survival, pathologic response, and genomics in CALGB 40601
(Alliance), a neoadjuvant Phase III trial of paclitaxel-trastuzumab with or
without lapatinib in HER2-positive breast cancer.} Journal of Clinical
Oncology. 2020. https://doi.org/10.1200/JCO.20.01276
\item Picornell AC, Echavarria I, Alvarez E, López-Tarruella S, Jerez Y, Hoadley K,
et al. Breast cancer PAM50 signature: Correlation and concordance between
RNA-Seq and digital multiplexed gene expression technologies in a triple
negative breast cancer series. BMC Genomics. 2019;20(1).
https://doi.org/10.1186/s12864-019-5849-0
}
}
\keyword{datasets}
