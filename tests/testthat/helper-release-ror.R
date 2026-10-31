r143_capture <- function(expr) {
    warnings <- list()
    value <- withCallingHandlers(force(expr), warning = function(w) {
        warnings[[length(warnings) + 1L]] <<- w
        invokeRestart("muffleWarning")
    })
    list(value = value, warnings = warnings)
}
r143_count <- function(x) sum(vapply(x$warnings, inherits, logical(1),
    what = "bsr_experimental_ror_pc_warning"))
r143_input <- function() {
    x <- new.env()
    data("OSLO2EMIT0obj", package = "BreastSubtypeR", envir = x)
    x$OSLO2EMIT0obj$data_input
}
r143_run <- function(method, input, clinical, subtype) {
    args <- list(se_obj = input$se_NC, hasClinical = clinical, Subtype = subtype)
    if (method == "BS_cIHC.itr") args$iteration <- 3L
    if (method == "BS_ssBC") args$s <- "ER"
    do.call(getExportedValue("BreastSubtypeR", method), args)
}
r143_fields <- c("ROR-PC (Subtype + Clinic + Prolif.Subtype)",
    "ROR-PC Group (Subtype + Clinic + Prolif.Subtype)")
