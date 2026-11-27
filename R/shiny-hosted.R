# Local launchers leave this option unset. The hosted wrapper requires a profile.
.shiny_hosted_profile <- function(profile = getOption("BreastSubtypeR.profile", "local"),
                                  deployment = FALSE) {
    if (nzchar(system.file("HOSTED-ONLY", package = "BreastSubtypeR")) &&
        (!identical(profile, "example-only") ||
         !identical(getOption("BreastSubtypeR.workerAddress"), "127.0.0.1")))
        stop("Start the hosted application through its supervised gateway.", call. = FALSE)
    if (identical(profile, "example-only")) return(TRUE)
    if (!deployment && identical(profile, "local")) return(FALSE)
    if (identical(profile, "research-uploads"))
        stop("Visitor uploads are closed pending verified admission and hosting controls.", call. = FALSE)
    stop("Missing or unsupported BreastSubtypeR operating profile.", call. = FALSE)
}

.shiny_iteration_limit <- function() {
    limit <- getOption("BreastSubtypeR.maxIterations", 100L)
    if (!identical(limit, 100L) && !identical(limit, 100))
        stop("The example pilot requires a fixed ceiling of 100 iterations.", call. = FALSE)
    100L
}

.shiny_allowed_method <- function(value) {
    is.character(value) && !is.object(value) && is.null(dim(value)) &&
        length(value) == 1L && !is.na(value) &&
        value %in% c("PAM50", "cIHC", "cIHC.itr", "PCAPAM50", "ssBC", "AIMS", "sspbc", "AUTO Mode")
}

.shiny_allowed_iterations <- function(value, limit) {
    typeof(value) %in% c("double", "integer") && !is.object(value) && is.null(dim(value)) &&
        length(value) == 1L && is.finite(value) && value >= 10 &&
        value <= limit && value == floor(value)
}
