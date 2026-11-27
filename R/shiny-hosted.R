# Local launchers leave this option unset. The hosted wrapper requires a profile.
.shiny_hosted_profile <- function(profile = getOption("BreastSubtypeR.profile", "local"),
                                  deployment = FALSE) {
    if (identical(profile, "example-only")) return(TRUE)
    if (!deployment && identical(profile, "local")) return(FALSE)
    if (identical(profile, "research-uploads"))
        stop("Visitor uploads are closed pending verified admission and hosting controls.", call. = FALSE)
    stop("Missing or unsupported BreastSubtypeR operating profile.", call. = FALSE)
}

.shiny_iteration_limit <- function() {
    limit <- getOption("BreastSubtypeR.maxIterations", 100L)
    if (!is.numeric(limit) || length(limit) != 1L || !is.finite(limit) ||
        limit != floor(limit) || limit < 10 || limit > 100)
        stop("Invalid hosted iteration limit (expected an integer from 10 to 100).", call. = FALSE)
    as.integer(limit)
}
