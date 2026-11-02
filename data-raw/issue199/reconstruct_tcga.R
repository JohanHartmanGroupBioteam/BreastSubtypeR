# Issue #199 Stage 2: input-only reconstruction, no Mapping or predictions.
# Rscript --vanilla reconstruct_tcga.R --baseline=FILE --archive=DIR --model-dir=DIR --out=NEW_DIR
options(stringsAsFactors=FALSE, digits=17)
a <- commandArgs(TRUE)
stopifnot(all(grepl("^--[^=]+=.+$",a)))
keys <- sub("^--([^=]+)=.*$","\\1",a)
opt <- setNames(sub("^--[^=]+=","",a),keys)
stopifnot(!anyDuplicated(keys),setequal(keys,c("baseline","archive","model-dir","out")))
suppressPackageStartupMessages(library(SummarizedExperiment))
for(p in c("AnnotationDbi","digest","jsonlite")) stopifnot(requireNamespace(p,quietly=TRUE))
sha <- function(p) { con<-file(p,"rb");on.exit(close(con));n<-file.info(p)$size;b<-readBin(con,"raw",n=n);stopifnot(length(b)==n);digest::digest(b,algo="sha256",serialize=FALSE) }
one <- function(p) { e<-new.env();nm<-load(p,e);stopifnot(length(nm)==1L);e[[nm]] }
archive <- normalizePath(opt[["archive"]],winslash="/",mustWork=TRUE)
paths <- c(baseline=opt[["baseline"]],
 counts=file.path(archive,"source-data/tcga/recount3/777c7b814bb5_tcga.gene_sums.BRCA.G026.gz"),
 gtf=file.path(archive,"source-data/tcga/recount3/777c554950c1_human.gene_sums.G026.gtf.gz"),
 samples=file.path(archive,"source-data/tcga/recount3/777c6c8330f_tcga.tcga.BRCA.MD.gz"),
 annotation_db=file.path(archive,"source-data/tcga/annotation/org.Hs.eg.sqlite"),
 added_plan=file.path(archive,"validation/tables/TCGA_added_feature_provenance.tsv"),
 sample_plan=file.path(archive,"validation/tables/TCGA_source_sample_candidates.tsv"),
 expanded_plan=file.path(archive,"provenance/reference-tables/TCGA_expanded_annotations.tsv"),
 archive_manifest=file.path(archive,"SHA256SUMS.txt"))
for(n in c("BreastSubtypeRobj","AIMSmodel","sspbc.models","Gene.ID.ann")) paths[n]<-file.path(opt[["model-dir"]],paste0(n,".RData"))
paths <- vapply(paths,normalizePath,"",winslash="/",mustWork=TRUE)
expected <- c(baseline="b94309c6f915d5f5fcb18368f3876f281b2111acf50eeb25b167f201f499a45d",
 counts="310b604478902e79573e70ff5de261013cc3c9f6ae1f11d262c9835bdfd110c5",
 gtf="2d1bdff8d91e08f1291fd47850d9d61cfa6768cc3a4cb8516775a873b007d0f6",
 samples="d8d5bf6145103283eff4370fe6f9162e8e431d7f2098b00f1e6e335a92aebd65",
 annotation_db="6e4ea433589342d028a65aef45ce6ccd70ac457099f145b145006f312b288d5f",
 BreastSubtypeRobj="dac97591d52816121c6d06ab2c7cbca16385f205616c9fc8fa73ae9793c03096",
 AIMSmodel="52b1190a09535d1df0ff09c06a4e0a180374da81649faaaf00194fd178b07b92",
 sspbc.models="9cf0c03578c4ed7fc5a0b8b2c8a1eb683778df039e6b02d53ddc5d65d500634c",
 Gene.ID.ann="c4291f8495bd639c0775e5cd83c3e5747988e04467587ca094a51ab76c404e0a")
arch_lines<-readLines(paths["archive_manifest"])
for(n in c("added_plan","sample_plan","expanded_plan")) {
 rel<-substring(paths[n],nchar(archive)+2L)
 z<-arch_lines[substring(arch_lines,67)==rel]
 stopifnot(length(z)==1L);expected[n]<-substr(z,1,64)
}
hashes <- vapply(paths,sha,"")
stopifnot(all(hashes[names(expected)]==expected))
out<-opt[["out"]]
stopifnot(!file.exists(out),dir.exists(dirname(out)))
stopifnot(dir.create(out))
out<-normalizePath(out,winslash="/",mustWork=TRUE)
tab<-function(x,n)write.table(x,file.path(out,paste0(n,".tsv")),sep="\t",row.names=FALSE,quote=FALSE,na="NA")
jout<-function(x,n)jsonlite::write_json(x,file.path(out,paste0(n,".json")),pretty=TRUE,auto_unbox=TRUE,na="null",digits=NA)
checks<-data.frame(check=character(),status=character(),detail=character())
check<-function(n,v,d="") { checks[nrow(checks)+1L,]<<-list(n,if(isTRUE(v))"PASS" else "FAIL",d);tab(checks,"construction_checks");if(!isTRUE(v))stop(n,": ",d) }
check("input_hashes",all(hashes[names(expected)]==expected),"Baseline, counts, GTF, sample metadata, annotation DB, models and retained reference tables.")
old<-one(paths["baseline"])$se_obj
rd<-rowData(old);cd<-colData(old);vals<-assay(old)
check("baseline_schema",identical(dim(old),c(5000L,100L))&&identical(assayNames(old),"counts")&&typeof(vals)=="integer")
check("baseline_unique_feature_and_sample_IDs",!anyDuplicated(rownames(old))&&!anyDuplicated(colnames(old))&&!anyNA(rownames(old))&&!anyNA(colnames(old)))
ref<-one(paths["BreastSubtypeRobj"]);aims<-one(paths["AIMSmodel"]);ssp<-one(paths["sspbc.models"]);model_ann<-one(paths["Gene.ID.ann"])
pairs<-lapply(c(list(AIMS=aims),ssp[c("ssp.subtype","ssp.pam50","ssp.ror")]),function(m) {
 p<-strsplit(m$all.pairs,"<",fixed=TRUE);stopifnot(all(lengths(p)==2L))
 data.frame(rule=seq_along(p),gene_a=sub("^e","",vapply(p,"[","",1)),gene_b=sub("^e","",vapply(p,"[","",2)))
})
req<-c(list(PAM50=as.character(ref$genes.sig50$EntrezGene.ID)),lapply(pairs,function(p)unique(c(p$gene_a,p$gene_b))))
union_ids<-sort(unique(unlist(req,use.names=FALSE)),method="radix")
needed<-setdiff(union_ids,as.character(rd$ENTREZID))
check("model_union",length(union_ids)==569L&&length(needed)==350L&&length(req$ssp.ror)==296L&&nrow(pairs$ssp.ror)==361L)
md<-read.delim(gzfile(paths["samples"]),check.names=FALSE,quote="",comment.char="")
candidate_ids<-unique(md$external_id[which(md$tcga_barcode %in% colnames(old))])
hcon<-gzfile(paths["counts"]);header<-readLines(hcon,n=3L);close(hcon)
columns<-strsplit(header[3],"\t",fixed=TRUE)[[1]]
check("full_source_header",identical(header[1],"##annotation=G026")&&length(columns)==1257L&&!anyDuplicated(columns)&&all(candidate_ids %in% columns))
classes<-rep("NULL",length(columns));classes[1]<-"character";classes[match(candidate_ids,columns)]<-"integer"
cat("Reading selected source columns from the complete G026 gzip...\n")
dat<-read.delim(gzfile(paths["counts"]),skip=2,check.names=FALSE,quote="",comment.char="",colClasses=classes)
mat<-as.matrix(dat[,-1,drop=FALSE]);rownames(mat)<-dat[[1]]
check("full_source_features",nrow(mat)==63856L&&!anyDuplicated(rownames(mat))&&typeof(mat)=="integer")
sample_candidates<-do.call(rbind,lapply(seq_len(ncol(old)),function(j) {
 ids<-unique(md$external_id[which(md$tcga_barcode==colnames(old)[j])])
 data.frame(sample_id=colnames(old)[j],external_id=ids,metadata_rows=length(ids),
 all_5000_cells_exact=vapply(ids,function(id)identical(unname(mat[rownames(old),id]),unname(vals[,j])),logical(1)))
}))
check("unique_fingerprint_resolution",all(tapply(sample_candidates$all_5000_cells_exact,sample_candidates$sample_id,sum)==1L))
prior_samples<-read.delim(paths["sample_plan"],check.names=FALSE,colClasses=c("character","character","integer","logical"))
rownames(sample_candidates)<-NULL
check("retained_sample_resolution",identical(sample_candidates,prior_samples),"Table row names reset; all ordered scientific columns must match exactly.")
sel<-sample_candidates[sample_candidates$all_5000_cells_exact,]
sel<-sel[match(colnames(old),sel$sample_id),]
selected_md<-md[match(sel$external_id,md$external_id),]
check("sample_resolution_barcodes",identical(selected_md$tcga_barcode,colnames(old))&&!anyDuplicated(sel$external_id)&&!anyNA(sel$external_id))
mat<-mat[,sel$external_id,drop=FALSE];colnames(mat)<-colnames(old)
check("historical_source_counts",identical(mat[rownames(old),,drop=FALSE],vals),"500,000 raw coverage counts, including integer storage.")
gtf<-read.delim(gzfile(paths["gtf"]),header=FALSE,sep="\t",quote="",comment.char="#",check.names=FALSE)
gtf_ids<-sub('.*gene_id "([^"]+)".*',"\\1",gtf[[9]])
check("G026_feature_identity",length(gtf_ids)==63856L&&!anyDuplicated(gtf_ids)&&setequal(gtf_ids,rownames(mat)))
db<-AnnotationDbi::loadDb(paths["annotation_db"]) # dbFileConnect uses SQLITE_RO.
tab(AnnotationDbi::metadata(db),"annotation_database_metadata")
# Recover historical first mappings using the archived database; never update it.
keys_all<-unique(sub("\\..*$","",rownames(mat)))
emap<-AnnotationDbi::mapIds(db,keys=keys_all,keytype="ENSEMBL",column="ENTREZID",multiVals="first")
smap<-AnnotationDbi::mapIds(db,keys=keys_all,keytype="ENSEMBL",column="SYMBOL",multiVals="first")
keys_source<-sub("\\..*$","",rownames(mat))
ann<-data.frame(probe=rownames(mat),SYMBOL=unname(smap[keys_source]),ENTREZID=unname(emap[keys_source]),
 Length=as.numeric(gtf[[6]][match(rownames(mat),gtf_ids)]),row.names=rownames(mat))
check("historical_annotation_source_identity",identical(ann[rownames(old),,drop=FALSE],as.data.frame(rd)))
selection<-do.call(rbind,lapply(needed,function(id) {
 m<-model_ann[which(as.character(model_ann$EntrezGene)==id),]
 feature<-intersect(m$Gene.ID,rownames(mat));stopifnot(length(feature)==1L)
 d<-ann[feature,,drop=FALSE];hist<-d$ENTREZID
 if(is.na(hist)) { stopifnot(id=="9906",feature=="ENSG00000215790.6");d$ENTREZID<-id;d$SYMBOL<-as.character(m$HGNC[match(feature,m$Gene.ID)]) } else stopifnot(hist==id)
 data.frame(d,historical_ENTREZID=hist,model_Gene_ID=feature,
 annotation_origin=if(is.na(hist))"bundled_Gene.ID.ann_exact_versioned_feature" else "org.Hs.eg.db_3.21.0_first_and_bundled_exact_feature",
 length_origin="recount_G026_GTF_score",source_row=match(feature,rownames(mat)),alternative_source_rows=sum(ann$ENTREZID %in% id))
}))
plan<-read.delim(paths["added_plan"],colClasses="character",check.names=FALSE)
check("added_feature_plan",identical(names(selection),names(plan))&&all(vapply(names(plan),function(n)identical(as.character(selection[[n]]),plan[[n]]),logical(1))),"350 ordered exact versioned feature/Entrez mappings and annotation recovery match retained evidence.")
check("added_unique_genes",nrow(selection)==350L&&!anyDuplicated(selection$probe)&&!anyDuplicated(selection$ENTREZID)&&!anyNA(selection$ENTREZID)&&!any(selection$probe %in% rownames(old)))
check("Entrez_9906_recovery",sum(is.na(selection$historical_ENTREZID))==1L&&selection$probe[is.na(selection$historical_ENTREZID)]=="ENSG00000215790.6"&&selection$ENTREZID[is.na(selection$historical_ENTREZID)]=="9906")
new_rd<-S4Vectors::DataFrame(selection[,names(rd),drop=FALSE],row.names=selection$probe)
extra<-SummarizedExperiment(assays=list(counts=mat[selection$probe,,drop=FALSE]),rowData=new_rd,colData=cd)
candidate<-rbind(old,extra)
prov<-list(stage="issue199_stage2_tcga",source_commit="6b2a4ab16798493f87713a91cee4674b2ae4df8e",package_version="1.5.12",
 input_scale="untransformed recount3 G026 base-pair coverage sums; integer storage",
 source_dimensions=c(63856L,1256L),candidate_dimensions=c(5350L,100L),
 construction="Original 5000 rows copied exactly, then 350 missing model genes in C-locale lexicographic Entrez-ID order; exact bundled versioned features.",
 length_origin="G026 gene-sum GTF score (field 6), not genomic span.",
 sample_selection="Historical barcode + unique exact 5000-count fingerprint among all matching source external IDs.",
 annotation_recovery="Entrez 9906: ENSG00000215790.6, symbol from bundled Gene.ID.ann; archived org.Hs.eg.db mapping absent.",
 contains_predictions=FALSE,transformation="None. No transform_counts, ranking, TPM or FPKM conversion.",
 source_hashes=as.list(hashes))
check("historical_metadata_identity",identical(metadata(candidate),metadata(old)),"Check original metadata before appending the new provenance entry.")
metadata(candidate)$issue199_reconstruction<-prov
check("candidate_schema",identical(class(candidate),class(old))&&identical(dim(candidate),c(5350L,100L))&&identical(assayNames(candidate),assayNames(old)))
check("historical_counts_identity",identical(assay(candidate)[seq_len(5000L),,drop=FALSE],vals))
check("historical_feature_identity",identical(rowData(candidate)[seq_len(5000L),],rd)&&identical(rownames(candidate)[seq_len(5000L)],rownames(old)))
check("sample_and_clinical_identity",identical(colData(candidate),cd)&&identical(colnames(candidate),colnames(old)))
check("all_source_counts_identity",identical(assay(candidate),mat[rownames(candidate),colnames(candidate),drop=FALSE]),"535,000 cells.")
check("all_G026_lengths",identical(as.numeric(rowData(candidate)$Length),as.numeric(gtf[[6]][match(rownames(candidate),gtf_ids)]))&&all(is.finite(rowData(candidate)$Length)&rowData(candidate)$Length>0))
check("count_domain",typeof(assay(candidate))=="integer"&&all(is.finite(assay(candidate)))&&all(assay(candidate)>=0)&&all(assay(candidate)==floor(assay(candidate))))
check("model_coverage",all(union_ids %in% rowData(candidate)$ENTREZID))
expected_ann<-read.delim(paths["expanded_plan"],colClasses=c("character","character","character","numeric"))
check("retained_expanded_annotations",identical(unname(as.matrix(as.data.frame(rowData(candidate)))),unname(as.matrix(expected_ann))))
tab(sample_candidates,"sample_resolution_candidates")
tab(data.frame(position=seq_len(ncol(old)),sample_id=colnames(old),external_id=sel$external_id,rail_id=selected_md$rail_id,source_position=match(sel$external_id,columns)-1L,barcode_candidates=sel$metadata_rows,unique_fingerprint_match=TRUE),"sample_manifest")
feat<-data.frame(position=seq_len(nrow(candidate)),as.data.frame(rowData(candidate)),historical=seq_len(nrow(candidate))<=5000L,source_row=match(rownames(candidate),rownames(mat)),gtf_row=match(rownames(candidate),gtf_ids),length_origin="G026_GTF_score")
tab(feat,"feature_manifest");tab(selection,"added_feature_provenance")
tab(do.call(rbind,lapply(names(req),function(n)data.frame(model=n,ENTREZID=req[[n]]))),"model_requirements")
tab(do.call(rbind,lapply(names(pairs),function(n)cbind(model=n,pairs[[n]]))),"model_rules")
# Explain pre-existing missing/repeated annotations without altering them.
ids<-as.character(rowData(candidate)$ENTREZID)
anomaly<-feat[is.na(ids)|duplicated(ids)|duplicated(ids,fromLast=TRUE),]
anomaly$explanation<-ifelse(is.na(anomaly$ENTREZID),"Preserved historical missing Entrez annotation","Preserved historical many-features-to-one Entrez mapping")
tab(anomaly,"identifier_annotations")
check("no_new_annotation_ambiguity",all(anomaly$historical)&&!anyNA(rowData(candidate)$probe)&&!anyDuplicated(rownames(candidate))&&!anyDuplicated(colnames(candidate)))
file<-file.path(out,"TCGABRCA_stage2_input.rds")
saveRDS(candidate,file,compress="xz",version=3)
check("RDS_roundtrip_identity",identical(readRDS(file),candidate))
check("input_preservation",identical(vapply(paths,sha,""),hashes))
tab(data.frame(role=names(paths),path=unname(paths),bytes=file.info(paths)$size,sha256=unname(hashes)),"input_sha256")
jout(prov,"provenance")
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
ns<-loadedNamespaces()
tab(data.frame(package=ns,version=vapply(ns,function(p)as.character(packageVersion(p)),""),path=vapply(ns,find.package,"")),"package_versions")
jout(list(candidate=list(path=basename(file),sha256=sha(file),bytes=file.info(file)$size,dimensions=dim(candidate),object_bytes=as.numeric(object.size(candidate))),checks=checks,input_hashes=as.list(hashes)),"manifest")
cat("COMPLETE:",nrow(checks),"construction checks PASS\n")
