% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runShinyBreastSubtypeR.R
\name{iBreastSubtypeR}
\alias{iBreastSubtypeR}
\title{iBreastSubtypeR}
\usage{
iBreastSubtypeR()
}
\value{
A table with subtyping and ROR score
}
\description{
Starts an interactive BreastSubtypeR shiny web app.

BreastSubtypeR integrates intrinsic molecular subtyping methods for breast
cancer, including nearest-centroid (NC-based) and single-sample predictor
(SSP-based) approaches. It employs standardized input and output formats,
providing a unified framework that is highly compatible with other R
packages in the gene expression profiling field.

The \code{iBreastSubtypeR()} function launches an interactive Shiny web
application. This app enables users to configure the arguments of subtyping
functions and execute subtyping on their local computer. For detailed
descriptions of the arguments, including their default and alternative
values, please refer to the manual pages of the respective functions.

Step 1:

The input data can be loaded from the user's workspace or by selecting a
CSV/text file for the expression data, a CSV/text file for clinical
information, and a CSV/text file for feature annotations.

After loading the necessary files, users can click the "Map Now" button
once and wait for the notification. If the Mapping() function runs
successfully, a message will appear stating, "You may now proceed to Step
2."

Step 2:

Users can select the desired subtyping method and adjust the relevant
parameters to conduct their analysis. Once the analysis is complete, a
message will indicate, "Analysis is complete." Two visualizations will be
displayed, and you will have the option to download the results as a text
file. If you wish to continue your analysis, you can directly run another
method without needing to repeat Step 1.
}
\examples{

library(BreastSubtypeR)

# This will open your browser with the BreastSubtypeR shiny web app
\donttest{
iBreastSubtypeR()
}

}
\keyword{BreastSubtypeR}
\keyword{Shiny}
