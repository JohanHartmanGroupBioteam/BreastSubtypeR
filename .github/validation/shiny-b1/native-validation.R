# Run only in a newly installed native macOS ARM64 library; no mocks.
library(BreastSubtypeR)
library(testthat)
stopifnot(identical(unname(Sys.info()["sysname"]),"Darwin"),
          R.version$arch %in% c("aarch64","arm64"))
stopifnot(normalizePath(find.package("BreastSubtypeR")) ==
          normalizePath(file.path(Sys.getenv("R_LIBS"),"BreastSubtypeR")))
target<-system.file("tests/testthat/test-shiny-example-parity.R",package="BreastSubtypeR")
stopifnot(identical(digest::digest(file=target,algo="sha256"),
 "722fc7f3aa394042e834104e88c993523c31ad0832c26d7725f85b8bf787bbae"))
e<-new.env(parent=globalenv())
for(x in parse(target)) if(!(is.call(x)&&identical(x[[1L]],as.name("test_that")))) eval(x,e)
obj<-e$shiny199_object("OSLO2EMIT0obj")
path<-e$shiny199_checked_path("OSLO2EMIT0_GEX_log2.FPKM.txt")
gex<-as.matrix(utils::read.delim(path,row.names=1L,check.names=FALSE))
reference<-as.matrix(SummarizedExperiment::assay(obj$se_obj))
storage.mode(gex)<-storage.mode(reference)<-"double"
profile<-e$shiny199_oslo_profile(gex,reference)
cat(e$shiny199_runtime(),"\nProfile:",profile,"\n")
cat("Observed binary64 fingerprint:",e$shiny199_numeric_fingerprint(gex),"\n")
cat("Binary reference fingerprint:",e$shiny199_numeric_fingerprint(reference),"\n")
cat("Bioconductor:",as.character(BiocManager::version()),"\n")
print(sessionInfo())
res<-test_package("BreastSubtypeR",filter="^shiny-(example-parity|expression-import)$",
 reporter="summary",stop_on_failure=FALSE,stop_on_warning=FALSE)
d<-as.data.frame(res)
write.table(d[c("file","test","nb","passed","failed","error","warning","skipped")],
 "native-blocks.tsv",sep="\t",row.names=FALSE)
print(colSums(d[c("passed","failed","error","warning","skipped")]))
native<-identical(profile,"native-full-table")
stopifnot(nrow(d)==14L,!any(d$failed|d$error|d$warning),
 sum(d$passed)==if(native) 234L else 236L,sum(d$skipped)==as.integer(native))
if(native) stopifnot(identical(d$test[which(d$skipped)],
 "OSLO text is bitwise identical to its frozen binary on the exact profile"))
