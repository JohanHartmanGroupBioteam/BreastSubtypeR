# Issue #199: versioned example reconstruction

The two public objects retain their names and list schema. The local Stage 3C
integration copies the approved Stage 3A files exactly; it does not regenerate
them to substitute a different artifact. Complete source-gene coverage does not
establish native SSP-ROR preparation or training-input equivalence.

## Reproduction

Reconstruction requires the approved generation implementation
`cb7a7c70e0f94e1690d15ac2ec2b464a40934f36`, R 4.5.2 on Windows and the recorded
[dependency stack](manifests/reconstruction-software-v1.tsv), plus the retained
external sources and original complete historical objects. Current main is not
the approved generation implementation: its fingerprint rejection is intentional.
The separately pinned f878f58c replay and post-#219 instructions below do not
relax this verification.

Use `/generation-checkout` for the pinned cb7a7c70 source installation and
`/integration-checkout` for this integration's reconstruction wrapper, approved
comparison objects and hash-pinned models. Set `--package-source` explicitly
to `/integration-checkout`; see
[Reproduction after merged #219](#reproduction-after-merged-219).
Install the pinned generation checkout into a **new private library**:

~~~text
R CMD INSTALL --install-tests --library=/new/library /generation-checkout
~~~

Set `R_LIBS` to that library and your pinned dependency library. No script
installs dependencies or downloads source data. Run from any working directory:

~~~text
Rscript --vanilla /integration-checkout/data-raw/issue199/rebuild-datasets.R --package-source=/integration-checkout --archive=/retained/2026-10-26_expanded-inputs --historical-dir=/retained/original-complete-objects --out=/new/reproduction
~~~

The archive path is external and relocatable. Preserve its internal layout and
`SHA256SUMS.txt`. The wrapper refuses existing output directories; the copied
Stage 1/2 constructors check source/model inputs against pinned SHA-256 values.
The new outputs are compared component by component with the approved package
objects using R `identical()`. The script never overwrites packaged data.
Serialization hashes can differ for identical R values; the approved artifact
hashes identify the actual integrated bytes.

Compact historical fixtures are bundled in
`tests/testthat/fixtures/historical-cb7a7c70` for regression reconstruction,
not as additional public datasets. The original complete historical objects
remain archived externally and are not included in the package.
Do not run the historical frozen-fixture generators.

## Required external sources

The retained archive supplies:

- `source-data/oslo/OSLO2EMIT0.103.genematrix_noNeg.Rdata`: protocol-adjusted
  linear FPKM from Mendeley Data V4, file
  `216a4cc7-e688-4c72-9513-2ad081e2c801`. The unadjusted matrix is not consumed.
- `source-data/oslo/Gene.ID.ann.Rdata`,
  `provenance/reference-tables/OSLO_exact_discrepancies.tsv` and
  `OSLO_planned_features.tsv` in the same reference-table directory.
- `source-data/tcga/recount3/777c7b814bb5_tcga.gene_sums.BRCA.G026.gz`:
  the complete 63,856-gene by 1,256-sample count source.
- `source-data/tcga/recount3/777c554950c1_human.gene_sums.G026.gtf.gz` and
  `777c6c8330f_tcga.tcga.BRCA.MD.gz` in the same recount3 directory.
- `source-data/tcga/annotation/org.Hs.eg.sqlite`: the pinned 383,188,992-byte
  annotation database, not a substitute from a current online annotation.
- `validation/tables/TCGA_added_feature_provenance.tsv`,
  `TCGA_source_sample_candidates.tsv` in the same validation-table directory,
  and `provenance/reference-tables/TCGA_expanded_annotations.tsv`.

Large original sources and the annotation database are intentionally not
bundled. See [source-inputs.tsv](manifests/source-inputs.tsv) for hashes and
external archive paths. Access to sources and redistribution clearance are
separate from local computational reproducibility; see the
[data attribution and distribution review](../../inst/DATA-PROVENANCE.md).

## Scientific recipe

Input model files originated at
`6b2a4ab16798493f87713a91cee4674b2ae4df8e` and are unchanged at
`cb7a7c70e0f94e1690d15ac2ec2b464a40934f36` (both package version 1.5.12).
Approved-artifact generation used the latter implementation; supported current-main replay is documented below.
Function fingerprints reject incompatible installed implementations.

OSLO preserves all original 356 rows, 36,312 values and 102 samples. The 213
additional rows use `log2(protocol-adjusted FPKM + 1)`, consistent with the
existing assay scale. MYO10 retains verified Entrez ID 4651. The 51 one-ULP
differences from independently recomputing historical rows are preserved rather
than overwritten. Sample `OSLO2EMIT0.073` remains excluded to preserve the
historical cohort; its original exclusion rationale is unresolved.

TCGA preserves all original 5,000 rows, 500,000 integer counts and 100 samples,
then appends 350 verified features in the approved order. Versioned Ensembl
identifiers and G026 GTF score-field lengths are retained; Entrez 9906 is
recovered from pinned model annotation. Ambiguous sample barcode
`TCGA-AC-A2QH-01B-04R-A22O-07` resolves uniquely to UUID
`6aa999cb-5d78-463b-9068-720c09eb9d34` by agreement with all 5,000 historical
counts. Competing source records are not interchangeable. The historical 230
missing annotations and eight rows sharing Entrez IDs are preserved. No
`transform_counts()`, reranking, TPM or FPKM conversion is applied to the
stored recount3 base-pair coverage sums.

Mapping uses `RawCounts = FALSE` for OSLO and `TRUE` for TCGA,
`method = "max"`, `impute = TRUE`, `verbose = FALSE`, and seed 1992026.
Default Mapping behavior is unchanged; no `se_SSP_ROR` is added.
`BS_Multi()` uses `methods = "AUTO"`, `Subtype = FALSE`,
`hasClinical = FALSE`, `er_subsample = "target"`, and seed 20260831.
RNG kinds are Mersenne-Twister/Inversion/Rejection. Method-internal defaults
remain unchanged. These seeds describe the reproducible generation, not an
unproven original save event.

## Versioned evidence

`manifests/` contains approved hashes, source hashes, feature manifests and
function fingerprints. Comparisons distinguish historical saved objects (A),
current-code historical inputs (B), and current-code expanded inputs (C).
OSLO B and C have identical mapped data and full method results; historical A
has correction-related ROR/Confidence differences. TCGA B-to-C changes 64 method
calls across 40 samples. Exact samples and before/after calls are retained as a
test oracle in `tests/testthat/fixtures`. Missing original build/seed evidence
prevents exclusive attribution of every A-to-B difference.

For Shiny table reproduction, follow the usage in
[export-shiny-examples.R](export-shiny-examples.R). It preserves public column
schemas. Six ancillary OSLO clinical fields come from the explicitly
hash-pinned historical Shiny table; the six shared fields come from candidate
`colData`. No additional clinical fields are attached to the public
`SummarizedExperiment`.

PR #228 merged the #199 integration after the scoped review was completed.
`DESCRIPTION` remains at 1.5.12; release and version-target decisions are separate.

The completed public-source and field-level
review has resolved the redistribution gate for these #199 candidates; source
terms and attribution requirements continue to apply. This is not institutional
approval or blanket clearance for future data, uses or releases. The historical
derivation limitation for the three unchanged OSLO Shiny fields is recorded in
the [data attribution and distribution review](../../inst/DATA-PROVENANCE.md).

CI, Rscript launcher verification, PDF tooling and version-target decisions
remain separate technical and release checks.

The complete-generation dependency versions and build versions are recorded in manifests/reconstruction-software-v1.tsv (R 4.5.2 on Windows). This records the mixed dependency stack that reproduced the approved values; it is not a claim of a clean current Bioconductor environment. Source/model/artifact preflight alone can be run with the same command plus --verify-only=true.

Archived public source URLs and the exact org.Hs.eg.db 3.21.0 annotation snapshot identifiers are retained in manifests/source-origins-v1.tsv, joined to the verified hashes. These are historical retrieval origins, not a promise of current download availability; automatic downloads are intentionally absent. Local reconstruction requires access to the retained archive and pinned dependency stack.

## Stage 3C: original inputs and current main

Dataset reconstruction now requires `--historical-dir` containing the original
`OSLO2EMIT0obj.RData` and `TCGABRCAobj.RData` files. Their SHA-256 values remain
43393c2bf06f0ea44b262d5525b75154ab99162fac5de898433556b6805450ec and
b94309c6f915d5f5fcb18368f3876f281b2111acf50eeb25b167f201f499a45d.
Use retained originals or the hash-verified historical Git sources identified
in the source-origin manifest. Compact test fixtures are not substituted for
these archived source files. Manifests resolve `historical/` against this
explicit directory, `archive/` against `--archive`, and `package/` against the
source checkout. No script downloads sources automatically.

The approved artifacts were generated at cb7a7c70. Stage 3C also validates
default Mapping and complete predictions under current main
`f878f58c1c12ca6f1cd5bf4489667935bf25969b`. Versioned function manifests identify
both supported implementations; the new manifest pins all 92 namespace
functions. Each reconstruction records its actual source commit and the
approved generation commit separately.

Merged #131 adds an opt-in Mapping input; the approved stored `data_input`
continues to contain only `se_NC` and `se_SSP`. On f878f58c, the expanded inputs
prepare all 296 SSP-ROR genes and 361 rules with `include_ssp_ror = TRUE`, while
the legacy components remain identical. OSLO preparation is `2^x` from adjusted
log2(FPKM+1), hence adjusted FPKM+1; TCGA preparation is linear FPKM calculated
by current Mapping from the original coverage sums. Coverage is not evidence
of equivalence to native unadjusted StringTie training input.

Use `compact-historical-fixtures.R` to reproduce compact test containers from
the original complete objects and approved public candidates; see that
script's usage and the fixture README. Historical metadata, mapped inputs and
results remain independently frozen. Only identical assay prefixes are shared
with public objects; independently frozen assay and complete-object digests
reject drift. Reconstructed complete objects retain every original value and
all scientific expectations. The original full files remain archived.

## Reproduction after merged #219

The merged #219 implementation at
`9110054b2032d0d560297c2c0166e39cb7776946` changes Mapping donor preparation.
It is not a supported reconstruction fingerprint in `rebuild-datasets.R`.
Compatibility tests against saved objects do not authorize a new generation
implementation or replacement artifact. Keep the existing manifests unchanged.

To reproduce the approved generation, create a separate source checkout at
`cb7a7c70e0f94e1690d15ac2ec2b464a40934f36`, install that checkout into the new
private library described above, and use the retained R 4.5.2 Windows dependency
stack in `manifests/reconstruction-software-v1.tsv`. Run this integration
checkout's reconstruction wrapper with
`--package-source=/path/to/this-integration-checkout`; retain the same external archive,
historical originals and new output-directory arguments. The wrapper and its
manifests, approved comparison objects and hash-pinned unchanged models come
from this integration checkout; only the installed implementation comes from
the explicitly pinned generation checkout. The original checkout's historical
packaged objects must not be substituted for the approved comparison objects.

The existing f878f58c replay remains separately pinned. Do not install current
main and bypass a fingerprint rejection. Function serialization fingerprints
and exact-value replay also depend on the recorded environment; passing Linux
compatibility tests does not establish cross-platform reconstruction identity.
