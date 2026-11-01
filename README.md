# BreastSubtypeR <a href='https://github.com/yqkiuo/BreastSubtypeR'><img src="inst/ShinyBreastSubtypeR/logo.svg" alt="BreastSubtypeR logo" align="right" height="110"/></a>

<!-- badges: start -->
[![Bioconductor Release](https://bioconductor.org/shields/years-in-bioc/BreastSubtypeR.svg)](https://bioconductor.org/packages/BreastSubtypeR)
[![Bioconductor Devel](https://bioconductor.org/shields/build/devel/bioc/BreastSubtypeR.svg)](https://bioconductor.org/packages/devel/bioc/html/BreastSubtypeR.html)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://opensource.org/licenses/GPL-3.0)
[![Paper: NAR Genomics & Bioinformatics](https://img.shields.io/badge/Paper-NAR%20Genomics%20%26%20Bioinformatics-0a7)](https://doi.org/10.1093/nargab/lqaf131)
<!-- badges: end -->

**BreastSubtypeR** is an assumption-aware, multi-method R/Bioconductor package with a local Shiny app.
It consolidates published intrinsic subtyping methods under one API and lets you run multiple classifiers at once.
**AUTO** inspects cohort diagnostics to select compatible methods and reduce misclassification.

> *Research use only; in clinical practice, intrinsic molecular subtyping is standardized via approved diagnostics (e.g., Prosigna®).*

📄 **Publication:** *NAR Genomics and Bioinformatics* (2025), **Editor’s Choice** → [doi:10.1093/nargab/lqaf131](https://doi.org/10.1093/nargab/lqaf131)

<details>
<summary><strong>How to cite</strong> (plain text &amp; BibTeX)</summary>

**Plain text**

Yang Q., Hartman J., Sifakis E.G. (2025). BreastSubtypeR: a unified R/Bioconductor package for intrinsic molecular subtyping in breast cancer research. *NAR Genomics and Bioinformatics*. <https://doi.org/10.1093/nargab/lqaf131>

**BibTeX**
```bibtex
@article{Yang2025BreastSubtypeR,
  author  = {Yang, Qiao and Hartman, Johan and Sifakis, Emmanouil G.},
  title   = {{BreastSubtypeR}: a unified R/Bioconductor package for intrinsic molecular subtyping in breast cancer research},
  journal = {NAR Genomics and Bioinformatics},
  year    = {2025},
  volume  = {7},
  number  = {4},
  pages   = {lqaf131},
  doi     = {10.1093/nargab/lqaf131},
  url     = {https://doi.org/10.1093/nargab/lqaf131}
}

```
</details>

------------------------------------------------------------------------

## Features

- **Unified interface for published methods:** consolidates PAM50 variants, AIMS, ssBC/sspbc, and others under one consistent API.
- **Run multiple methods at once (`BS_Multi`):** execute several classifiers in a single call and compare results side by side.
- **AUTO (cohort-aware selection):** checks ER/HER2 distribution, subtype purity, and subgroup sizes; disables incompatible classifiers.
- **Method-specific pre-processing:** automatically routes raw RNA-seq counts, precomputed FPKM, or log2-processed microarray/nCounter matrices.
- **Robust mapping:** Entrez ID–based gene mapping with conflict resolution.
- **Local Shiny app (`iBreastSubtypeR`):** point-and-click analysis; data stay on your machine.
- **Reproducibility:** Bioconductor distribution, unit tests, vignettes, and `SummarizedExperiment` compatibility.

------------------------------------------------------------------------

### Methods included (single-method implementations)

| **Method id** | **Short description** | **Group** | **Reference** |
|------------------|-------------------|------------------|------------------|
| `parker.original` | Original PAM50 by Parker et al., 2009 | NC-based | [Parker et al., 2009](https://doi.org/10.1200/JCO.2008.18.1370) |
| `genefu.scale` | PAM50 implementation as in the genefu R package (scaled version) | NC-based | [Gendoo et al., 2016](https://doi.org/10.1093/bioinformatics/btv693) |
| `genefu.robust` | PAM50 implementation as in the genefu R package (robust version) | NC-based | [Gendoo et al., 2016](https://doi.org/10.1093/bioinformatics/btv693) |
| `cIHC` | Conventional ER-balancing using immunohistochemistry (IHC) | NC-based | [Ciriello et al., 2015](https://doi.org/10.1016/j.cell.2015.09.033) |
| `cIHC.itr` | Iterative version of cIHC | NC-based | [Curtis et al., 2012](https://doi.org/10.1038/nature10983) |
| `PCAPAM50` | Selects IHC-defined ER subsets, then uses Principal Component Analysis (PCA) to create ESR1 expression-based ER-balancing | NC-based | [Raj-Kumar et al., 2019](https://doi.org/10.1038/s41598-019-44339-4) |
| `ssBC` | Subgroup-specific gene-centering PAM50 | NC-based | [Zhao et al., 2015](https://doi.org/10.1186/s13058-015-0520-4) |
| `ssBC.v2` | Updated subgroup-specific gene-centering PAM50 with refined quantiles | NC-based | [Fernandez-Martinez et al., 2020](https://doi.org/10.1200/JCO.20.01276) |
| `AIMS` | Absolute Intrinsic Molecular Subtyping (AIMS) method | SSP-based | [Paquet & Hallett, 2015](https://doi.org/10.1093/jnci/dju357) |
| `sspbc` | Single-Sample Predictors for Breast Cancer (AIMS adaptation) | SSP-based | [Staaf et al., 2022](https://doi.org/10.1038/s41523-022-00465-3) |

(See the vignette for implementation details.)

------------------------------------------------------------------------

## Installation

Install the released version from Bioconductor:

``` r
# Requires R >= 4.5.0
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("BreastSubtypeR")

# Devel:
BiocManager::install("BreastSubtypeR", version = "devel")
```

Or install from GitHub:

``` r
if (!require("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("yqkiuo/BreastSubtypeR")
```

## Release 1.4.3: experimental ROR interpretation

This is the Bioconductor 3.23 release branch. Version 1.4.3 adds interpretation
warnings while preserving the 1.4.2 calculations and output schema.

**Experimental 18-gene ROR-PC** and **Experimental 18-gene ROR-PC group** are
display labels for the unchanged legacy columns:

```text
ROR-PC (Subtype + Clinic + Prolif.Subtype)
ROR-PC Group (Subtype + Clinic + Prolif.Subtype)
```

This research calculation differs from the published Lien implementation.
Its prognostic performance and risk-category thresholds have not been
independently validated. It is not a clinical Prosigna result and must not
guide patient care. The groups are not validated patient-risk assignments.
The original 11-gene Parker ROR-PC, stored under
`ROR-PC (Subtype + Clinic + Prolif)`, remains a distinct output.

Relevant public NC calls emit `bsr_experimental_ror_pc_warning` once when
both experimental fields are returned; `BS_Multi()` aggregates this class
once across its returned results. Clinical-off calls, SSP-only outputs and
`Mapping()` do not emit it. With `options(warn = 2)`, the new warning becomes
an error. Expert workflows can handle this class specifically with
`withCallingHandlers()` and `invokeRestart("muffleWarning")`; unrelated
warnings should remain visible.

Shiny displays a persistent notice when corresponding Full metrics results
are available and selected. Full metrics exports append the interpretation
warning to an existing comment line, preserving the TSV schema, filenames
and number of header lines. Calls-only and AUTO exports are unchanged.
The release's Spearman correlations, calibration, both proliferation terms,
coefficients, thresholds, NODE rules and legacy edge cases are retained.
Development fixes #198 and #206 are not included.

## Quick start

> *These examples use datasets shipped with the package. For your own data, provide a `SummarizedExperiment` **with clinical metadata in `colData`** (e.g., `PatientID`, `ER`, `HER2`; for ROR: `TSIZE`, `NODE`).*

``` r
library(BreastSubtypeR)

# Example data
data("BreastSubtypeRobj")
data("OSLO2EMIT0obj")
```

**1) Map & prepare (method-specific pre-processing + mapping)**

``` r
data_input <- Mapping(OSLO2EMIT0obj$se_obj, RawCounts = FALSE, method = "max", impute = TRUE)
```

**2) Multi-method run (user-defined)**

``` r
res <- BS_Multi(data_input = data_input, methods = c("parker.original","PCAPAM50","sspbc"))
head(res$res_subtypes, 5)
```

**3) AUTO mode (cohort-aware selection) + visualize**

``` r
res_auto <- BS_Multi(data_input = data_input, methods = "AUTO")
Vis_Multi(res_auto$res_subtypes)
```

***AUTO logic (clarifications)***

- **ER/HER2-defined cohorts** (ER+/HER2−, ER−/HER2−, ER+/HER2+, ER−/HER2+): 
**ssBC.v2 only** + SSP (AIMS, sspbc).
- **ER-only** (ER+ or ER−) and **TNBC** (size permitting): 
**ssBC and/or ssBC.v2** + SSP.
- ER balance gate (simulation-based): `lower_ratio = 0.39`, `upper_ratio = 0.69`.
- Minimum sizes (defaults): ER+ = 15, ER− = 18, TN = 18. 
Subgroups use half of ER totals (rounded):
    -   ER+: `n_ERposHER2pos_threshold = n_ERposHER2neg_threshold = round(15 / 2)`
    -   ER−: `n_ERnegHER2pos_threshold = n_ERnegHER2neg_threshold = round(18 / 2)`

*Provenance*: ER+/ER− minimums are simulation-based defaults. 
TN currently uses the ER− minimum (18). 
Subgroup gates and TN minimum may be refined in future releases 
as additional simulations become available.

**4) Launch the local Shiny app**

``` r
BreastSubtypeR::iBreastSubtypeR() # interactive GUI (local)
```

***Notes:***

- The app runs locally; no data leave your machine.

- If you see a missing UI dependency:

``` r
install.packages(c("shiny","bslib"))
```

***Example data (for Shiny & scripts)***

You can explore the app without preparing files:

- In **iBreastSubtypeR**, click **“Load example data…”** to pre-fill
  the UI with a small demo cohort. Then click **Preprocess & map** and proceed to Step 2.

- Programmatically, the same files are shipped inside the package:
```r
  exdir   <- system.file("RshinyTest", package = "BreastSubtypeR")
  gex     <- file.path(exdir, "OSLO2EMIT0_GEX_log2.FPKM.txt")
  clin    <- file.path(exdir, "OSLO2EMIT0_clinical.txt")
  anno    <- file.path(exdir, "OSLO2EMIT0_anno.txt")
  stopifnot(all(file.exists(gex, clin, anno)))
```

## Vignette & documentation

A comprehensive usage guide (input types, AUTO details, full method descriptions) is included as a **vignette**.

See function help pages for specifics (e.g., `?BS_Multi`, `?Mapping`, `?iBreastSubtypeR`).

## Contributing & issues

**Canonical source**  
Bioconductor package page: <https://bioconductor.org/packages/BreastSubtypeR>  
Bioconductor DOI: <https://doi.org/10.18129/B9.bioc.BreastSubtypeR>  
Mirrors: <https://github.com/yqkiuo/BreastSubtypeR> (personal), <https://github.com/JohanHartmanGroupBioteam/BreastSubtypeR> (org)

**Support & bugs**  
Bugs/PRs: <https://github.com/yqkiuo/BreastSubtypeR/issues>

## License

GPL-3
