# Retained public-base wrapper: f6dbdbb874785008d8511155b1ba8b72f7c0f22d.
# Positive aligned-input oracle only. The unchanged native calculations and
# manual dispatch are reused in a private lexical environment; no namespace
# binding is changed. Expected values do not pass through the corrected guard.
container_baseline_run <- local({
    oracle <- new.env(parent = asNamespace("BreastSubtypeR"))
    oracle$BS_ssBC <- function(se_obj, s, Subtype = FALSE, hasClinical = FALSE) {
    if (!inherits(se_obj, "SummarizedExperiment")) {
        stop("Input must be a SummarizedExperiment object.")
    }
    gene_expr <- assay(se_obj)
    pheno <- .nc_pheno(se_obj, hasClinical)
    if (!"PatientID" %in% colnames(pheno)) {
        stop("The 'colData' of 'se_obj' must include a 'PatientID' column.")
    }
    rownames(pheno) <- pheno$PatientID
    if (hasClinical) {
        req <- c("TSIZE", "NODE")
        miss <- setdiff(req, colnames(pheno))
        if (length(miss))
            stop("When hasClinical = TRUE, colData(se_obj) must include: TSIZE and NODE.\n                               Missing: ",
                paste(miss, collapse = ", "), ".")
        if (!is.numeric(pheno$NODE))
            stop("colData(se_obj)$NODE must be numeric.")
    }
    required_columns <- switch(s, ER = c("ER"), ER.v2 = c("ER", "HER2"), TN = c("TN"), TN.v2 = c("TN"),
        stop("Invalid value for 's'. Must be one of 'ER', 'ER.v2', 'TN', or 'TN.v2'."))
    missing_columns <- setdiff(required_columns, colnames(pheno))
    if (length(missing_columns) > 0) {
        stop(paste("The following required columns are missing from 'colData':", paste(missing_columns,
            collapse = ", ")))
    }
    arguments <- rlang::dots_list(mat = gene_expr, df.cln = pheno, s = s, Subtype = Subtype,
        hasClinical = hasClinical, receptor_data = colData(se_obj), .homonyms = "last")
    call <- rlang::call2(makeCalls.ssBC, !!!arguments)
    res_ssBC <- eval(call)
    .warn_experimental_ror_pc(res_ssBC)
    return(res_ssBC)
}
    environment(oracle$BS_ssBC) <- oracle
    oracle$BS_Multi <- get("BS_Multi", asNamespace("BreastSubtypeR"))
    environment(oracle$BS_Multi) <- oracle
    function(d, route, four = FALSE) {
        ids <- colnames(d$se_NC); cd <- SummarizedExperiment::colData(d$se_NC)
        stopifnot(sum(names(cd) == "PatientID") == 1L,
            identical(as.character(cd[["PatientID"]]), ids),
            identical(as.character(data.frame(cd)[["PatientID"]]), ids))
        warnings <- list()
        result <- withCallingHandlers(switch(route,
            ER = oracle$BS_ssBC(d$se_NC, "ER", Subtype = four, hasClinical = TRUE),
            ER.v2 = oracle$BS_ssBC(d$se_NC, "ER.v2", Subtype = four, hasClinical = TRUE),
            multi = oracle$BS_Multi(d, methods = c("ssBC", "ssBC.v2"), Subtype = four, hasClinical = TRUE)),
            warning = function(w) {
                warnings[[length(warnings) + 1L]] <<- list(class = class(w), ids = w$patient_ids, message = conditionMessage(w))
                invokeRestart("muffleWarning")
            }, message = function(m) {
                stopifnot(grepl("ssBC.*running", conditionMessage(m)))
                invokeRestart("muffleMessage")
            })
        list(result = result, warnings = warnings)
    }
})
