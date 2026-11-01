args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
candidate <- normalizePath(args[1], winslash="/")
lib <- normalizePath(args[2], winslash="/")
e <- normalizePath(args[3], winslash="/")
cat(R.version.string, "\n")
before <- file.path(e,"before-CITATION")
after <- file.path(candidate,"inst/CITATION")
checks <- character()
ok <- function(value,label) {stopifnot(isTRUE(value));checks <<- c(checks,label);cat("PASS:",label,"\n")}
capture <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(expr,warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
  })
  list(value=value,warnings=warnings)
}
git_value <- function(ref) {
    value <- system2("git", c("-C", shQuote(candidate), "rev-parse", shQuote(ref)), stdout=TRUE)
    stopifnot(is.null(attr(value, "status")))
    value
}
ok(identical(git_value("HEAD"), "99f148f9ede122aa3a120cf2869ced0990bcf740"),
   "Citation checks use the corrected candidate SHA")
ok(identical(git_value("HEAD^{tree}"), "3987d2eac26a283d5e9398e58153c6449251396c"),
   "Citation checks use the corrected candidate tree")
desc0 <- read.dcf(file.path(e,"before-DESCRIPTION"))
desc1 <- read.dcf(file.path(candidate,"DESCRIPTION"))
ok(identical(colnames(desc0),colnames(desc1)),"DESCRIPTION fields unchanged")
ok(identical(desc0[,setdiff(colnames(desc0),"Suggests"),drop=FALSE],
    desc1[,setdiff(colnames(desc1),"Suggests"),drop=FALSE]),"Only Suggests changed")
deps <- function(x)trimws(strsplit(x,",",fixed=TRUE)[[1L]])
ok(identical(deps(desc1[1,"Suggests"]),c("S4Vectors",deps(desc0[1,"Suggests"]))),"Exactly S4Vectors added to Suggests")
p0 <- parse(before,keep.source=FALSE);p1 <- parse(after,keep.source=FALSE)
ok(length(p0)==3L && length(p1)==3L,"CITATION parses")
ok(identical(p0[1:2],p1[1:2]),"Published citation executable expressions unchanged")
ok(identical(p0[[3L]],p1[[3L]][[3L]][[2L]]),"Original Manual expression preserved inside guard")
meta0 <- as.list(desc0[1,]);meta1 <- as.list(desc1[1,])
original <- capture(utils::readCitationFile(before,meta=meta0))
corrected <- capture(utils::readCitationFile(after,meta=meta1))
ok(length(corrected$value)==2L,"Actual DESCRIPTION metadata retains Article and Manual")
ok(identical(original,corrected),"Actual DESCRIPTION citation objects and warnings unchanged")
ok(identical(utils::toBibtex(original$value),utils::toBibtex(corrected$value)),"Actual DESCRIPTION BibTeX unchanged")
if(length(corrected$warnings))cat("Inherited source-metadata warning:",paste(corrected$warnings,collapse="; "),"\n")
none <- capture(utils::readCitationFile(after))
ok(length(none$value)==1L && !length(none$warnings),"No metadata returns one citation with zero warnings")
ok(identical(none$value[[1L]],original$value[[1L]]),"Published article citation unchanged without metadata")
ok(identical(utils::toBibtex(none$value),utils::toBibtex(original$value[1L])),"Published article BibTeX unchanged without metadata")
for(meta in list(NULL,list(),list(Package="BreastSubtypeR"),
    meta1[setdiff(names(meta1),"Version")],meta1[setdiff(names(meta1),"Authors@R")])) {
  x <- capture(utils::readCitationFile(after,meta=meta))
  ok(length(x$value)==1L && !length(x$warnings) &&
      identical(as.list(x$value[[1L]]),as.list(none$value[[1L]])) &&
      identical(utils::toBibtex(x$value),utils::toBibtex(none$value)),
      "Absent/incomplete metadata skips only Manual")
}
legacy <- meta1;legacy$Author <- "Qiao Yang [aut], Emmanouil G. Sifakis [aut]";legacy[["Authors@R"]] <- NULL
x <- capture(utils::readCitationFile(before,meta=legacy))
y <- capture(utils::readCitationFile(after,meta=legacy))
ok(identical(x,y) && identical(utils::toBibtex(x$value),utils::toBibtex(y$value)),"Legacy Author metadata behavior and BibTeX unchanged")
{
 installed_path <- normalizePath(find.package("BreastSubtypeR", lib.loc=lib), winslash="/")
 ok(identical(installed_path, normalizePath(file.path(lib,"BreastSubtypeR"), winslash="/")),
    "Citation uses the R CMD check installation")
 bytes <- function(p)readBin(p,"raw",n=file.info(p)$size)
 ok(identical(bytes(file.path(installed_path,"CITATION")),bytes(after)),
    "Installed CITATION bytes match the corrected source")
 meta <- utils::packageDescription("BreastSubtypeR",lib.loc=lib)
 ok(identical(meta$Version,"1.4.3"),"Installed package is version 1.4.3")
 installed <- capture(utils::citation("BreastSubtypeR",lib.loc=lib))
 expected <- capture(utils::readCitationFile(before,meta=meta))
 direct <- capture(utils::readCitationFile(after,meta=meta))
 ok(length(installed$value)==2L && !length(installed$warnings),"Installed citation() returns Article and Manual with zero warnings")
 ok(identical(installed,expected)&&identical(direct,expected),"Installed citation() and original valid-metadata behavior identical")
 ok(identical(utils::toBibtex(installed$value),utils::toBibtex(expected$value)),"Installed citation BibTeX unchanged")
 ok(identical(installed$value[[1L]],none$value[[1L]]),"Installed published article retained exactly")
 writeLines(as.character(utils::toBibtex(installed$value)),file.path(e,"installed-citation.bib"),useBytes=TRUE)
 saveRDS(meta,file.path(e,"installed-package-metadata.rds"))
}
saveRDS(list(checks=checks,source_original=original,source_corrected=corrected,
    no_metadata=none,R=R.version.string),file.path(e,"citation-results.rds"))
cat("Total:",length(checks),"focused checks passed\n")

writeLines(c("### Canonical citation verification", paste0("PASS: ", length(checks), " checks."),
    paste0("- ", checks), "Original source-metadata warnings (preserved):",
    if(length(original$warnings)) original$warnings else "None."), file.path(e,"summary.txt"))
jsonlite::write_json(list(completed=TRUE, candidate_sha=git_value("HEAD"),
    candidate_tree=git_value("HEAD^{tree}"), R=R.version.string, installed_path=installed_path,
    checks=checks, source_warnings=original$warnings, installed_warnings=installed$warnings),
    file.path(e,"citation-results.json"), pretty=TRUE, auto_unbox=TRUE)
