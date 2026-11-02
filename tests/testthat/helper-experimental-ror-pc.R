## Handle only the approved experimental condition in existing scientific
## regressions. Other warnings retain their usual testthat/R behavior.
ror221_notice <- paste0(
    "Experimental 18-gene ROR-PC: this research calculation differs from the ",
    "published Lien implementation. Its prognostic performance and ",
    "risk-category thresholds have not been independently validated. ",
    "It is not a clinical Prosigna result and must not guide patient care."
)
ror221_fields <- c("ROR-PC (Subtype + Clinic + Prolif.Subtype)",
    "ROR-PC Group (Subtype + Clinic + Prolif.Subtype)")

expect_experimental_ror_pc <- function(expr, n = 1L) {
    count <- 0L
    result <- withCallingHandlers(force(expr),
        bsr_experimental_ror_pc_warning = function(w) {
            count <<- count + 1L
            expect_identical(conditionMessage(w), ror221_notice)
            expect_null(conditionCall(w))
            invokeRestart("muffleWarning")
        })
    expect_identical(count, as.integer(n))
    result
}
