# Display labels are separate from the legacy data keys. Keep result schemas
# and scoring code unchanged; these helpers only identify and explain outputs.
.experimental_ror_pc_fields <- function() {
    c(
        "ROR-PC (Subtype + Clinic + Prolif.Subtype)" =
            "Experimental 18-gene ROR-PC",
        "ROR-PC Group (Subtype + Clinic + Prolif.Subtype)" =
            "Experimental 18-gene ROR-PC group"
    )
}

.experimental_ror_pc_notice <- function() {
    paste0(
        "Experimental 18-gene ROR-PC: this research calculation differs ",
        "from the published Lien implementation. Its prognostic performance ",
        "and risk-category thresholds have not been independently validated. ",
        "It is not a clinical Prosigna result and must not guide patient care."
    )
}

.has_experimental_ror_pc <- function(result) {
    scores <- if (is.data.frame(result)) result else if (is.list(result)) {
        result[["score.ROR"]]
    } else NULL
    is.data.frame(scores) && nrow(scores) > 0L &&
        all(names(.experimental_ror_pc_fields()) %in% names(scores))
}

# Signal only at the completed public-call boundary, never inside RORgroup().
# BS_Multi collects this specific condition from children and signals once for
# its returned results; no session-global state or broad warning suppression.
.warn_experimental_ror_pc <- function(result) {
    if (.has_experimental_ror_pc(result)) {
        warning(structure(
            list(message = .experimental_ror_pc_notice(), call = NULL),
            class = c("bsr_experimental_ror_pc_warning", "warning", "condition")
        ))
    }
    invisible(result)
}
