# Shiny example tables

These six tab-delimited files are synchronized with the expanded packaged
objects staged under issue #199. The app's **Load example data** button loads
OSLO; the TCGA files can be uploaded together after selecting **Raw counts**.

| Files | Expression | Features x samples | Annotation columns | Clinical columns |
|---|---|---|---|---|
| OSLO2EMIT0_* | log2(protocol-adjusted FPKM + 1) | 569 x 102 | probe, SYMBOL, ENTREZID | PatientID, ER, PR, HER2, TSIZE, NODE, IHC, Prosigna.Subtype, Prosigna.ROR, MKI67.corrected, Ki67.estimated, IHC.estimated |
| TCGABRCAsub_* | recount3 gene-level base-pair coverage sums | 5350 x 100 | probe, SYMBOL, ENTREZID, Length | PatientID, ER, PR, HER2 |

Expression is exported with 17 significant decimal digits. All counts and
expression values round-trip exactly to the corresponding packaged se_obj
with the recorded R reader. The first expression column contains the
versioned Ensembl probe ID; its header is blank. All other headers are sample
IDs in the packaged order. Annotation rows follow that same feature order.
TCGA Length is the G026 gene-sum GTF score field, in base pairs. No read-count
conversion, log transform of TCGA counts, or alternative normalization occurs
during export.

All shared clinical fields come from the corresponding candidate. OSLO now
uses the canonical PR+/PR- labels. To preserve the existing 12-column
upload schema, its six additional fields are retained by exact PatientID
match from the original public clinical table at commit
cb7a7c70e0f94e1690d15ac2ec2b464a40934f36; they are not newly inferred candidate
metadata. The original checkout table has SHA-256
00a24cfde1bc7837a5f5cc46d9c2e4ce55a60a624b99caf5ee276005d89b8996.
Its six fields are IHC, Prosigna.Subtype, Prosigna.ROR, MKI67.corrected,
Ki67.estimated and IHC.estimated. They do not enter the app's AUTO calculation.
The OSLO annotation table now contains exactly the 569 candidate features;
the previous 19,675-row background annotation is no longer needed.

The app still calls default Mapping with imputation off and retains its
existing classifier, eligibility, voting, entropy and download behavior.
Complete source-gene coverage does not provide complete native SSP-ROR
preparation or establish compatibility with the unadjusted StringTie FPKM
used for SSPBC training. Issue #131 and normalization research #208 remain
separate.

Recreate the tables with the source recipe
data-raw/issue199/export-shiny-examples.R. It requires the exact validated
candidate RData files and the hash-pinned historical OSLO clinical TSV. It
accepts relocatable paths and emits an optional SHA-256 manifest. Original
large expression sources are not needed for this export step. Publication
remains subject to the data-distribution review described in the dataset
provenance documentation.

