# Supplemental observation of the unchanged installed test expressions. The
# normal installed suites run separately and supply the primary test results.
library(BreastSubtypeR)
library(testthat)
root <- system.file("tests/testthat", package="BreastSubtypeR")
path <- file.path(root,"test-issue199-posterior-bounds.R")
here <- dirname(sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1L]))
source(file.path(here,"evidence-utils.R"))
env <- new.env(parent=globalenv())
sys.source(file.path(root,"helper-issue199-portability.R"),env)
observations <- list()
env$expect_error <- function(object, regexp=NULL, ...) {
    caller <- parent.frame()
    scalar <- function(name) {
        v <- get0(name,envir=caller,inherits=FALSE)
        if(is.atomic(v) && length(v)==1L) as.character(unname(v)) else ""
    }
    expression <- paste(deparse(substitute(object)),collapse=" ")
    seen <- character()
    value <- testthat::expect_error(withCallingHandlers(object,
        error=function(e) seen <<- c(seen,conditionMessage(e))), regexp, ...)
    observations[[length(observations)+1L]] <<- data.frame(expression,expected=regexp,
        observed=if(length(seen))seen[1L] else "NO ERROR",method=scalar("method"),
        dataset=scalar("dataset"),change=scalar("change"),size=scalar("size"),
        fraction=scalar("fraction"))
    value
}
for (expr in parse(path)) eval(expr,env)
obs <- do.call(rbind,observations)
write_tsv(obs,file.path(Sys.getenv("EVIDENCE_DIR"),"tests","observed-first-safeguards.tsv"))
stopifnot(nrow(obs)>0L,all(mapply(grepl,obs$expected,obs$observed)))
