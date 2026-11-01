# Example data attribution and distribution review (1 November 2026)

This local integration is not distribution clearance. Publication, public push, release and external hosting remain blocked pending author/institutional sign-off on the combined expression and clinical-annotation derivatives. No additional institutional approval is asserted to exist.

The package code is GPL-3. That software license does not replace third-party data terms.

| Component | Evidence and attribution | Remaining release gate |
|---|---|---|
| OSLO adjusted expression and public annotation | Johan Vallon-Christersson, Mendeley Data V4, DOI 10.17632/yzxtxn4nmd.4; CC BY 4.0. Cite Staaf et al., npj Breast Cancer 2022, DOI 10.1038/s41523-022-00465-3. Retain attribution, license link and modification notice (subset, annotation correction, uniform adjusted log2 transformation, expansion). | Confirm the exact clinical fields and inherited ancillary Shiny fields fall within the source license/consent and that institutional redistribution review is complete. |
| TCGA recount3 coverage sums | Johns Hopkins University recount3 resource; registry declares Public Domain. Cite Wilks et al., Genome Biology 22:323 (2021), DOI 10.1186/s13059-021-02533-6; acknowledge the TCGA Research Network. | Confirm the retained clinical metadata origin and that all shipped fields are open-access derivatives; recount3's declaration does not independently clear every linked annotation. |
| Annotation/model mappings | Pinned G026 GTF, archived org.Hs.eg.db SQLite and package model annotation. Preserve respective source attribution; large database and originals remain external. | Verify notices for any redistributed annotation derivatives and institutional requirements. No original-source/database license is inferred from GPL-3. |
| Historical regression fixtures | Unmodified complete objects already used by the package, preserved solely in tests with exact hashes. | Same clinical/attribution review applies; prior packaging is not independent proof of permission. |

OSLO is log2(protocol-adjusted FPKM+1), not unadjusted StringTie FPKM. SSPBC development used unadjusted StringTie FPKM; prediction equivalence is unestablished. TCGA stores original recount3 base-pair coverage sums, not fragment counts or transformed read counts. The historical values and cohorts remain fixed; expansions and regenerated mapped/prediction components are documented in the dataset help and NEWS.

NCI asks for TCGA Research Network acknowledgment and reports that former publication moratoria are lifted. Those publication statements are not a substitute for checking access conditions or clinical-field redistribution.

Sources checked 1 November 2026:
- [Mendeley dataset and CC BY 4.0 designation](https://data.mendeley.com/datasets/yzxtxn4nmd/4).
- [CC BY 4.0 license](https://creativecommons.org/licenses/by/4.0/).
- [recount3 registry: Public Domain and citation](https://registry.opendata.aws/recount/).
- [TCGA acknowledgment and publication guidance](https://www.cancer.gov/ccg/research/genome-sequencing/tcga/using-tcga-data/citing).
- [Bioconductor package-data guidance](https://contributions.bioconductor.org/data.html).
- [Current BiocCheck size guidance](https://bioconductor.org/packages/devel/bioc/vignettes/BiocCheck/inst/doc/BiocCheck.html).

Small examples may be bundled with documented sources and compression. Current BiocCheck guidance specifies a 10 MB software source-tarball limit and 5 MB individual-file limit; older 5 MB-total guidance must not be confused with the current limit. Large data belong in ExperimentHub; git-lfs is not supported by the Bioconductor build system. Measure the final tarball including self-contained tests and Shiny exports before review. ExperimentHub would require cleared rights, metadata, versioned accessions and a separately reviewed loader/API design; it is not introduced by this integration.

A/B/C provenance distinctions remain essential: historical saved A, current-code historical input B, and expanded current-code C. The original save-event build/seed is incomplete. Dataset coverage alone does not establish SSP-ROR preparation completeness or equivalence to its training input; #131 and #208 remain separate.
