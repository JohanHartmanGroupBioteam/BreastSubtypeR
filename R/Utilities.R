#'
#' Functions for BreastSubtypeR package
#' @keywords internal
#' @import ComplexHeatmap
#' @import RColorBrewer
#' @import circlize
#' @import ggrepel
#' @import magrittr
#' @import impute
#' @importFrom dplyr mutate across everything
#' @importFrom ggplot2 ggplot aes geom_point geom_bar geom_text geom_line geom_vline geom_hline
#' @importFrom ggplot2 geom_boxplot scale_x_continuous scale_color_manual scale_fill_manual
#' @importFrom ggplot2 coord_polar labs theme theme_classic theme_minimal theme_void element_text
#' @importFrom SummarizedExperiment SummarizedExperiment colData rowData assay
#'
#' @noRd
NULL
NULL

#' Impute missing values using k-nearest neighbors
#'
#' @param x Numeric matrix with missing values (NA).
#' @param verbose Logical; if TRUE, prints missing probes/samples.
#' @return Imputed matrix.
#' @noRd
impute_missing <- function(x, verbose = FALSE) {
    if (anyNA(x)) {
        if (verbose) {
            probeid_NA <- rownames(x)[rowSums(is.na(x)) > 0]
            sample_NA <- colnames(x)[colSums(is.na(x)) > 0]
            message("Missing probes: ", paste(probeid_NA, collapse = ", "))
            message("Missing samples: ", paste(sample_NA, collapse = ", "))
        }
        x <- impute::impute.knn(x)$data
    }
    x
}


#' Collapse duplicate genes by a summary statistic
#'
#' @param x Numeric matrix (rows = probes/genes, cols = samples).
#' @param y Dataframe with `probe` and `ENTREZID` columns.
#' @param method Method for resolving duplicates ("max", "mean", etc.).
#' @return Deduplicated matrix (rows = unique Entrez IDs).
#' @noRd
duplicate_genes <- function(x, y, method) {
    probeid <- rownames(x)
    entrezid <- y$ENTREZID
    names(entrezid) <- y$probe

    ## process probeid in input data
    entrezid <- entrezid[probeid]
    ## remove NA
    entrezid <- entrezid[!(is.na(entrezid))]
    x <- x[names(entrezid), ]
    entrezid <- factor(entrezid, levels = unique(entrezid))
    ## names are unique probeid and content are redundant entrezid

    ## This is for probeID or transcriptID
    ## split expression matrix
    split_mat <- split.data.frame(as.data.frame(x), entrezid, drop = FALSE)

    # function to calculate the desired statistic
    calculate_stat <- function(mat, method) {
        switch(method,
            "mean" = apply(mat, 2, mean, na.rm = TRUE),
            "median" = apply(mat, 2, median, na.rm = TRUE),
            "max" = mat[which.max(rowSums(mat)), ],
            "stdev" = {
                temp <- mat
                stdevs <- apply(mat, 1, sd, na.rm = TRUE)
                vals <- temp[match(max(stdevs), stdevs), ]
                return(vals)
            },
            "iqr" = {
                temp <- mat
                iqrs <- apply(mat, 1, function(x) {
                    quantile(x, 0.75, na.rm = TRUE) -
                        quantile(x, 0.25, na.rm = TRUE)
                })
                vals <- temp[match(max(iqrs), iqrs), ]
                return(vals)
            }
        )
    }

    ## keep processed x
    x <- mapply(
        calculate_stat,
        split_mat,
        MoreArgs = list(method = method),
        SIMPLIFY = TRUE,
        USE.NAMES = TRUE
    )
    x <- apply(x, 1, unlist)

    return(x)
}


#' Collapse duplicate genes by a summary statistic
#'
#' @param x Numeric matrix (rows = probes/genes, cols = samples).
#' @param y Dataframe with `probe` and `ENTREZID` columns.
#' @return Deduplicated matrix (rows = unique Entrez IDs).
#' @noRd
prepare_nc_matrix <- function(x, genes.sig50, samplenames, verbose) {
    ## print necessary information
    ## Parker
    missing_ID_parker <- setdiff(genes.sig50$EntrezGene.ID, rownames(x))

    if (length(missing_ID_parker) == 0 & verbose) {
        message("Genes used in NC-based methods are covered.")
    } else if (verbose) {
        message("These genes are missing for NC-based methods:")
        message(paste(missing_ID_parker, collapse = "\n"))
    }

    ## get matrix for NC (symbol as rows, sample as col)
    genes_nc <- genes.sig50$EntrezGene.ID
    x_NC <- x[na.omit(match(genes_nc, rownames(x))), ]
    rownames(x_NC) <- genes.sig50$Symbol[match(rownames(x_NC), genes_nc)]
    x_NC <- data.frame(x_NC)
    colnames(x_NC) <- samplenames

    return(x_NC)
}


#' Prepare matrix for SSP-based methods
#'
#' @param x Numeric matrix (rows = Entrez IDs, cols = samples).
#' @param genes.s Dataframe with gene signatures.
#' @param RawCounts Logical; if TRUE, input is raw counts.
#' @param verbose Logical; print missing genes.
#' @return Processed matrix for SSP.
#' @noRd
prepare_ssp_matrix <- function(x, genes.s, RawCounts, samplenames, verbose) {
    ## AIMS
    genes_ssp <- genes.s$EntrezGene.ID[genes.s$SSP_based == "Yes"]
    missing <- setdiff(genes_ssp, rownames(x))

    if (verbose && length(missing) > 0) {
        message("Missing SSP genes: ", paste(missing, collapse = ", "))
    } else if (verbose) {
        message("Genes used in SSP-based methods are covered.")
    }

    ## get matrix for AIMS (entrezID as colnames)
    x_SSP <- x[intersect(genes_ssp, rownames(x)), , drop = FALSE]
    colnames(x_SSP) <- samplenames

    ## exponential transformation for microarray/RNAseq
    if (!RawCounts) {
        x_SSP <- 2^x_SSP
    }

    return(x_SSP)
}


#' Map Gene IDs and Handle missing data
#'
#' @param method A string specifying the method for resolving duplicate probes in microarray or RNA-seq data. Options include:
#'   - `"iqr"`: Selects the probe with the highest interquartile range (IQR), typically used for short-oligo arrays (e.g., Affymetrix).
#'   - `"mean"`: Chooses the probe with the highest average expression, commonly used for long-oligo arrays (e.g., Agilent, Illumina).
#'   - `"max"`: Retains the probe with the highest expression value, often used for RNA-seq data.
#'   - `"stdev"`: Selects the probe with the highest standard deviation.
#'   - `"median"`: Chooses the probe with the highest median expression value.
#' @noRd

domapping <- function(se_obj,
    RawCounts = FALSE,
    method = "max",
    impute = TRUE,
    verbose = TRUE) {
    ## 1. Input raw counts
    if (RawCounts && !"Length" %in% colnames(rowData(se_obj))) {
        stop(
            "Missing gene length information: for RawCounts=TRUE, ",
            "rowData(se_obj)$Length must be provided (in base pairs) ",
            "to enable FPKM calculation."
        )
    }

    # 2. Load gene signatures
    data_env <- new.env(parent = emptyenv())
    data("BreastSubtypeRobj", envir = data_env, package = "BreastSubtypeR")
    BreastSubtypeRobj <- data_env[["BreastSubtypeRobj"]]
    genes.s <- BreastSubtypeRobj$genes.signature

    # 3. Normalize raw counts (if applicable)
    if (RawCounts) {
        # Pull counts from SummarizedExperiment and build a DGEList
        counts <- SummarizedExperiment::assay(se_obj)
        y_all <- edgeR::DGEList(counts = round(as.matrix(counts)))

        # Light positive-only subset (very gentle): CPM > 0 in ≥30% samples
        prop_pos <- rowMeans(edgeR::cpm(y_all) > 0.1)
        keep_sub <- prop_pos >= 0.30
        if (!any(keep_sub)) stop("Cannot form a positive subset for UQ normalization.")
        y_sub <- y_all[keep_sub, , keep.lib.sizes = TRUE]

        # Estimate UQ factors on the subset (p = 0.75)
        y_sub <- edgeR::calcNormFactors(y_sub, method = "upperquartile", p = 0.75)

        # Apply those factors to the full data
        y_all$samples$norm.factors <- y_sub$samples$norm.factors

        # NC-based input for downstream: log-CPM on the FULL matrix (uses effective lib sizes)
        x <- edgeR::cpm(y_all, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)

        # FPKM from raw counts + gene lengths (in base pairs)
        gl <- as.numeric(SummarizedExperiment::rowData(se_obj)$Length)
        names(gl) <- rownames(se_obj)
        gl <- gl[rownames(y_all)]
        counts.fpkm <- edgeR::rpkm(y_all,
            gene.length = gl,
            normalized.lib.sizes = FALSE, log = FALSE
        )
    } else {
        x <- SummarizedExperiment::assay(se_obj)
    }

    # 4. Check feature data
    ## Extract data from SummarizedExperiment
    y <- rowData(se_obj) %>% data.frame()
    ## check feature data
    if (length(y) == 0) {
        stop("Please provide feature annotation to do probeID mapping ")
    }
    ## change data type
    y$ENTREZID <- as.integer(y$ENTREZID)
    samplenames <- colnames(x)

    # 5. Filter by signature genes and impute
    ## filter by ENTREZID
    y <- y[y$ENTREZID %in% genes.s$EntrezGene.ID, ]
    x <- x[y$probe, ]
    if (impute && anyNA(x)) x <- impute_missing(x, verbose)
    if (RawCounts && impute && anyNA(x)) {
        counts.fpkm <- impute_missing(counts.fpkm, verbose)
    }

    # 6. Deduplicate by EntrezID
    x <- duplicate_genes(x, y, method)
    if (RawCounts) {
        counts.fpkm <- duplicate_genes(counts.fpkm, y, method)
    }

    # 7. Prepare output (NC and SSP matrices)
    x_NC <- prepare_nc_matrix(x, BreastSubtypeRobj$genes.sig50, samplenames, verbose)
    x_SSP <- prepare_ssp_matrix(
        if (RawCounts) counts.fpkm else x,
        genes.s,
        RawCounts,
        samplenames,
        verbose
    )

    result <- list(x_NC = x_NC, x_SSP = x_SSP)

    return(result)
}

#' Function implementing AUTO decision logic
#' @noRd
get_methods <- function(pheno) {
    #### AUTO mode
    cohort.select <- "ERpos"
    samples_ER.icd <- NULL
    samples_ERHER2.icd <- NULL
    methods <- NULL

    if (ncol(pheno) == 0) {
        message("The pheno table has not been detected.")
        message("Running methods: AIMS, & sspbc")
        methods <- c("AIMS", "sspbc")
    } else if (!all(c("ER", "HER2") %in% colnames(pheno))) {
        stop("The 'AUTO' mode requires both 'ER' and 'HER2' columns in the 'pheno' dataframe.")
    } else {
        # Calculate sample sizes
        sample_counts <- with(pheno, table(ER, HER2))
        n_ERpos <- sum(pheno$ER == "ER+")
        n_ERneg <- sum(pheno$ER == "ER-")

        n_ERnegHER2pos <- sum(pheno$ER == "ER-" & pheno$HER2 == "HER2+")
        n_ERnegHER2neg <- sum(pheno$ER == "ER-" & pheno$HER2 == "HER2-")
        n_ERposHER2pos <- sum(pheno$ER == "ER+" & pheno$HER2 == "HER2+")
        n_ERposHER2neg <- sum(pheno$ER == "ER+" & pheno$HER2 == "HER2-")

        # Set thresholds
        n_ERpos_threshold <- 15 # simulations-based cut-off
        n_ERneg_threshold <- 18 # simulations-based cut-off
        n_ERposHER2pos_threshold <- round(n_ERpos_threshold / 2) # simulations ongoing (available on the next release)
        n_ERposHER2neg_threshold <- round(n_ERpos_threshold / 2) # simulations ongoing (available on the next release)
        n_ERnegHER2pos_threshold <- round(n_ERneg_threshold / 2) # simulations ongoing (available on the next release)
        n_ERnegHER2neg_threshold <- round(n_ERneg_threshold / 2) # simulations ongoing (available on the next release)
        upper_ratio <- 0.69 # simulations-based cut-off
        lower_ratio <- 0.39 # simulations-based cut-off

        ## main panel
        if ("TN" %in% colnames(pheno)) {
            message("A TNBC cohort has been detected.")
            cohort.select <- "TNBC"

            if (!("TN" %in% colnames(pheno))) {
                stop("Provide \"TN\" in pheno for: ssBC(TN) & ssBC.v2 (TN)")
            }

            message("Running methods: ssBC (TN), ssBC.v2 (TN), AIMS & sspbc")
            methods <- c("ssBC", "ssBC.v2", "AIMS", "sspbc")
        } else if (n_ERposHER2neg == 0 && n_ERnegHER2neg == 0) {
            message("A HER2+ cohort has been detected.")
            cohort.select <- "HER2pos"

            if (n_ERposHER2pos < n_ERposHER2pos_threshold && n_ERnegHER2pos < n_ERnegHER2pos_threshold) {
                message("A small HER2+ cohort has been detected.")
                message("Running methods: AIMS, & sspbc")
                methods <- c("AIMS", "sspbc")
            } else if (n_ERposHER2pos >= n_ERposHER2pos_threshold && n_ERnegHER2pos < n_ERnegHER2pos_threshold) {
                message("A ER+/HER2+ cohort has been detected.")
                message("Running methods: ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC.v2", "AIMS", "sspbc")
            } else if (n_ERposHER2pos < n_ERposHER2pos_threshold && n_ERnegHER2pos >= n_ERnegHER2pos_threshold) {
                message("A ER-/HER2+ cohort has been detected.")
                message("Running methods: ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC.v2", "AIMS", "sspbc")
            } else {
                message("Running methods: ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC.v2", "AIMS", "sspbc")
            }
        } else if (n_ERpos < n_ERpos_threshold && n_ERneg < n_ERneg_threshold) {
            message("A small number of ER-/ER+ samples has been detected.")
            message("Running methods: AIMS & sspbc")
            methods <- c("AIMS", "sspbc")
        } else if (n_ERpos >= n_ERpos_threshold && n_ERneg < n_ERneg_threshold) {
            if (n_ERposHER2pos >= n_ERposHER2pos_threshold && n_ERposHER2neg >= n_ERposHER2neg_threshold) {
                message("Running methods for ER+ samples:
                ssBC, ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC", "ssBC.v2", "AIMS", "sspbc")
            } else if (n_ERposHER2pos >= n_ERposHER2pos_threshold && n_ERposHER2neg < n_ERposHER2neg_threshold) {
                message("Running methods for ER+/HER2+ samples:
                ssBC, ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC", "ssBC.v2", "AIMS", "sspbc")
            } else if (n_ERposHER2pos < n_ERposHER2pos_threshold && n_ERposHER2neg >= n_ERposHER2neg_threshold) {
                message("Running methods for ER+/HER2- samples:
                        ssBC, ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC", "ssBC.v2", "AIMS", "sspbc")
            }
        } else if (n_ERpos < n_ERpos_threshold && n_ERneg >= n_ERneg_threshold) {
            if (n_ERnegHER2pos >= n_ERnegHER2pos_threshold && n_ERnegHER2neg >= n_ERnegHER2neg_threshold) {
                message("Running methods for ER- samples:
                ssBC, ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC", "ssBC.v2", "AIMS", "sspbc")
            } else if (n_ERnegHER2pos >= n_ERnegHER2pos_threshold && n_ERnegHER2neg < n_ERnegHER2neg_threshold) {
                message("Running methods for ER-/HER2+ samples:
                ssBC, ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC", "ssBC.v2", "AIMS", "sspbc")
            } else if (n_ERnegHER2pos < n_ERnegHER2pos_threshold && n_ERnegHER2neg >= n_ERnegHER2neg_threshold) {
                message("Running methods for ER-/HER2- samples:
                ssBC, ssBC.v2, AIMS, & sspbc")
                methods <- c("ssBC", "ssBC.v2", "AIMS", "sspbc")
            }
        } else if (n_ERpos >= n_ERpos_threshold && n_ERneg >= n_ERneg_threshold) {
            ## for other NC-based methods
            ratio_ER <- n_ERpos / (n_ERneg + n_ERpos)

            if (ratio_ER > lower_ratio && ratio_ER < upper_ratio) {
                message(
                    "Running methods:
                    parker.original, genefu.scale, genefu.robust, ssBC, ssBC.v2, cIHC, cIHC.itr, PCAPAM50, AIMS & sspbc"
                )
                methods <- c(
                    "parker.original",
                    "genefu.scale",
                    "genefu.robust",
                    "ssBC",
                    "ssBC.v2",
                    "cIHC",
                    "cIHC.itr",
                    "PCAPAM50",
                    "AIMS",
                    "sspbc"
                )
            } else {
                message(
                    "The ER+ ratio in the current dataset differs from that observed in the UNC232 training cohort."
                )
                message("Running methods:
                        genefu.robust, ssBC, ssBC.v2, cIHC, cIHC.itr, PCAPAM50, AIMS & sspbc")
                methods <- c(
                    "genefu.robust",
                    "ssBC",
                    "ssBC.v2",
                    "cIHC",
                    "cIHC.itr",
                    "PCAPAM50",
                    "AIMS",
                    "sspbc"
                )
            }
        }

        if (cohort.select != "TNBC") {
            ## subsetting samples for ssBC and ssBC.v2
            # Handle ssBC & ssBC.v2 method for imbalanced subtypes
            ERHER2_counts <- c(
                n_ERpos,
                n_ERneg,
                n_ERnegHER2pos,
                n_ERnegHER2neg,
                n_ERposHER2pos,
                n_ERposHER2neg
            )
            names(ERHER2_counts) <- c(
                "ERpos",
                "ERneg",
                "ERnegHER2pos",
                "ERnegHER2neg",
                "ERposHER2pos",
                "ERposHER2neg"
            )

            er_idx <- c(
                ERpos = ERHER2_counts["ERpos"] > n_ERpos_threshold,
                ERneg = ERHER2_counts["ERneg"] > n_ERneg_threshold
            )
            samples_ER <- names(ERHER2_counts)[seq(1, 2)][er_idx]

            erher2_idx <- c(
                ERposHER2pos = ERHER2_counts["ERposHER2pos"] > n_ERposHER2pos_threshold,
                ERposHER2neg = ERHER2_counts["ERposHER2neg"] > n_ERposHER2neg_threshold,
                ERnegHER2pos = ERHER2_counts["ERnegHER2pos"] > n_ERnegHER2pos_threshold,
                ERnegHER2neg = ERHER2_counts["ERnegHER2neg"] > n_ERnegHER2neg_threshold
            )
            samples_ERHER2 <- names(ERHER2_counts)[seq(3, 6)][erher2_idx]

            if (cohort.select != "HER2pos" && cohort.select != "TNBC") {
                if (length(samples_ER) > 0) {
                    message("ssBC for samples: ", paste(samples_ER, collapse = ", "))
                    samples_ER.icd <- unlist(lapply(samples_ER, function(subtype) {
                        subtype <- str_replace_all(subtype, "pos", "+")
                        subtype <- str_replace_all(subtype, "neg", "-")
                        ER_status <- subtype
                        rownames(pheno)[pheno$ER == ER_status]
                    }))
                }
            }

            if (cohort.select != "TNBC") {
                if (length(samples_ERHER2) > 0) {
                    message(
                        "ssBC.v2 for samples: ",
                        paste(samples_ERHER2, collapse = ", ")
                    )

                    samples_ERHER2.icd <- unlist(
                        lapply(samples_ERHER2, function(subtype) {
                            subtype <- subtype |>
                                str_replace_all("pos", "+") |>
                                str_replace_all("neg", "-")
                            ER_sts <- substr(subtype, 1, 3)
                            HER2_sts <- substr(subtype, 4, 8)
                            rownames(pheno)[pheno$ER == ER_sts & pheno$HER2 == HER2_sts]
                        })
                    )
                }
            }
        }
    }

    if (length(samples_ER.icd) == 0) {
        samples_ER.icd <- NULL
    }
    if (length(samples_ERHER2.icd) == 0) {
        samples_ERHER2.icd <- NULL
    }

    return(list(samples_ER.icd = samples_ER.icd, samples_ERHER2.icd = samples_ERHER2.icd, methods = methods, cohort.select = cohort.select))
}

#' Function for consensus subtype
#' @noRd
get_consensus_subtype <- function(patient_row) {
    patient_row <- unlist(patient_row, use.names = FALSE)
    counts <- table(patient_row)
    max_subtype <- names(counts)[which.max(counts)]
    return(max_subtype)
}

#' Function for entropy calculation
#' @noRd
get_entropy <- function(patient_row) {
    freq <- table(patient_row)
    prob <- freq / sum(freq)
    entropy <- -sum(prob * log2(prob))
    return(entropy)
}


#' Function to get the average correlation and ROR
#' @noRd
get_average_subtype <- function(res_ihc_iterative, consensus_subtypes) {
    ## correlation and ROR to be averaged
    sum_colnames <- c("Basal", "Her2", "LumA", "LumB", "Normal")

    all_patients <- names(res_ihc_iterative[[1]]$predictions)

    sum_cols_list <- mapply(
        function(res_ihc) {
            # res_ihc = as.data.frame( res_ihc_iterative[[1]])

            res_ihc$distances <- as.data.frame(res_ihc$distances)

            ## if FALSE, make the cell as NULL
            keep <- res_ihc$predictions == consensus_subtypes
            res_ihc$distances[!keep, ] <- as.list(rep(NA, ncol(res_ihc$distances)))
            
            ## convert every column to numeric (preserving NA), same as before
            res <- dplyr::mutate(
              res_ihc$distances,
              dplyr::across(
                dplyr::everything(),
                ~ ifelse(!is.na(.x), as.numeric(as.character(.x)), NA_real_)
              )
            )

            return(res)
        },
        res_ihc_iterative,
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE
    )


    ## count_na
    count_weight_save <- Reduce(`+`, lapply(sum_cols_list, function(x) {
        x[!is.na(x)] <- 1
        x[is.na(x)] <- 0
        return(x)
    }))

    ## sum all for each cell
    sum_cols_save <- Reduce(`+`, lapply(sum_cols_list, function(x) {
        ## change NA cell to 0 cell
        x[is.na(x)] <- 0

        return(x)
    }))

    ## get the mean for each cell
    ## only when subtype is supported by consensus_subtypes
    ## for each iteration and each patient
    mean_cols_save <- sum_cols_save / count_weight_save

    ## get the mean testdata
    sum_cols_list.testdata <- mapply(
        function(res_ihc) {
            res_ihc$testData
        },
        res_ihc_iterative,
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE
    )

    ## sum all for each cell
    mean_cols_save.testdata <- Reduce(
        `+`,
        sum_cols_list.testdata
    ) / length(sum_cols_list.testdata)


    ## distances.Subtype for ROR
    sum_cols_list.Subtype <- mapply(
        function(res_ihc) {
            # res_ihc = res_ihc_iterative[[1]]

            res_ihc$distances.Subtype <- data.frame(res_ihc$dist.RORSubtype)

            ## if FALSE, make the cell as NULL
            keep <- res_ihc$predictions == consensus_subtypes
            res_ihc$distances.Subtype[!keep, ] <- as.list(rep(NA, ncol(res_ihc$distances.Subtype)))
            
            res <- dplyr::mutate(
              res_ihc$distances.Subtype,
              dplyr::across(
                dplyr::everything(),
                ~ ifelse(!is.na(.x), as.numeric(as.character(.x)), NA_real_)
              )
            )

            return(res)
        },
        res_ihc_iterative,
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE
    )

    ## count_na
    count_weight_save.Subtype <- Reduce(
        `+`,
        lapply(sum_cols_list.Subtype, function(x) {
            x[!is.na(x)] <- 1
            x[is.na(x)] <- 0
            return(x)
        })
    )


    ## sum all for each cell
    sum_cols_save.Subtype <- Reduce(
        `+`,
        lapply(sum_cols_list.Subtype, function(x) {
            ## change NA cell to 0 cell
            x[is.na(x)] <- 0
            return(x)
        })
    )

    ## get the mean for each cell
    ## only when subtype is supported by consensus_subtypes
    ## for each iteration and each patient
    sum_cols_save.Subtype <- sum_cols_save.Subtype / count_weight_save.Subtype

    res <- list(
        mean_distance = mean_cols_save,
        mean_distance.Subtype = sum_cols_save.Subtype,
        testdata = mean_cols_save.testdata
    )


    return(res)
}

#### Visualization  ####

#' Boxplot of Correlation per Subtype
#'
#' @name Vis_boxplot
#' @description This function generates a boxplot to visualize the correlation
#' distribution between different subtypes of breast cancer, based on the
#' provided correlation table and subtype information.
#'
#' @param out A data frame containing the columns `"PatientID"` and `"Subtype"`.
#'   The `"PatientID"` column should have unique identifiers for each patient,
#'   and the `"Subtype"` column should specify the assigned subtype for each
#'   patient.
#' @param correlations A data frame or matrix containing the correlation values
#'   computed from NC-based methods.
#'
#' @return A `ggplot` object representing the boxplot visualization of the
#'   correlation distributions across the different subtypes.
#'
#' @examples
#' data("OSLO2EMIT0obj")
#' res <- OSLO2EMIT0obj$res
#'
#' # Prepare data: Subtype information and correlation matrix
#' out <- data.frame(
#'     PatientID = res$results$genefu.robust$BS.all$PatientID,
#'     Subtype = res$results$genefu.robust$BS.all$BS
#' )
#'
#' correlations <- res$results$genefu.robust$outList$distances
#'
#' # Generate the boxplot
#' p <- Vis_boxplot(out, correlations)
#' plot(p)
#'
#' @export

Vis_boxplot <- function(out, correlations) {
    df <- data.frame(
        predictions = out$Subtype,
        cor = apply(correlations, 1, max)
    )

    plot <- ggplot(df, aes(x = .data$predictions, y = .data$cor)) +
        geom_boxplot() +
        labs(x = "", y = "Correlation") +
        theme_classic() +
        theme(axis.text = element_text(size = 12))

    return(plot)
}


#' Heatmap Visualization of Gene Expression by Subtype
#'
#' @name Vis_heatmap
#' @description This function generates a heatmap to visualize gene expression
#' patterns across breast cancer subtypes, based on the provided gene expression
#' matrix and subtype information.
#'
#' @param x A gene expression matrix, where genes are rows and samples are
#'   columns. The data should be log2 transformed.
#' @param out A data frame containing two columns: `"PatientID"` and
#'   `"Subtype"`. The `"PatientID"` column should contain unique patient
#'   identifiers, and the `"Subtype"` column should specify the assigned subtype
#'   for each patient.
#'
#' @return A `ggplot` or `heatmap` object (depending on implementation)
#'   representing the heatmap of gene expression across different subtypes.
#'
#' @examples
#' library(SummarizedExperiment)
#' data("OSLO2EMIT0obj")
#' res <- OSLO2EMIT0obj$res
#'
#' # Prepare data: Gene expression matrix and subtype information
#' x <- assay(OSLO2EMIT0obj$data_input$se_NC)
#' out <- data.frame(
#'     PatientID = res$results$genefu.robust$BS.all$PatientID,
#'     Subtype = res$results$genefu.robust$BS.all$BS
#' )
#'
#' # Generate the heatmap
#' p <- Vis_heatmap(x, out)
#' plot(p)
#'
#' @export

Vis_heatmap <- function(x, out) {
    scaled_mat <- t(scale(t(x[, out$PatientID])))

    ## color
    scaled_dat <- c(min(scaled_mat), 0, max(scaled_mat))
    col_fun <- circlize::colorRamp2(scaled_dat, c("green", "black", "red"))

    ## column annotation
    col_anno <- data.frame(
        row.names = out$PatientID,
        Subtype = out$Subtype
    )
    anno_col <- ComplexHeatmap::HeatmapAnnotation(
        df = col_anno,
        show_legend = FALSE,
        col = list(
            Subtype = c(
                "Basal" = "red",
                "Her2" = "hotpink",
                "LumA" = "darkblue",
                "LumB" = "skyblue",
                "Normal" = "green"
            )
        )
    )

    heatmap <- ComplexHeatmap::Heatmap(
        scaled_mat,
        name = "Expr",
        col = col_fun,
        ## annotation
        top_annotation = anno_col,

        ## clustering
        ## as original heatmap plot
        clustering_distance_rows = "spearman",
        clustering_method_rows = "average",
        show_row_dend = FALSE,
        cluster_column_slices = TRUE,
        column_split = col_anno$Subtype,
        clustering_distance_columns = "spearman",
        clustering_method_columns = "average",

        ## general
        show_column_names = FALSE,
        show_heatmap_legend = TRUE,
        row_names_gp = grid::gpar(fontsize = 6)
    )

    return(heatmap)
}




#' PCA Plot Visualization of Gene Expression by Subtype
#'
#' @name Vis_PCA
#' @description This function generates a PCA plot to visualize the principal
#' components of gene expression data, colored by the assigned subtypes.
#' Optionally, it can display a scree plot of eigenvalues to evaluate the
#' explained variance.
#'
#' @param x A gene expression matrix, where genes are rows and samples are
#'   columns. The data should be log2 transformed.
#' @param out A data frame containing two columns: `"PatientID"` and
#'   `"Subtype"`. The `"PatientID"` column should contain unique patient
#'   identifiers, and the `"Subtype"` column should specify the assigned subtype
#'   for each patient.
#' @param Eigen Logical. If `TRUE`, the function will display a scree plot
#'   showing the eigenvalues of the principal components.
#'
#' @return A `ggplot` object representing the PCA plot, colored by subtype. If
#'   `Eigen` is set to `TRUE`, a scree plot of the eigenvalues is also included.
#'
#' @examples
#' library(SummarizedExperiment)
#' data("OSLO2EMIT0obj")
#' res <- OSLO2EMIT0obj$res
#'
#' # Prepare data: Gene expression matrix and subtype information
#' x <- assay(OSLO2EMIT0obj$data_input$se_NC)
#' out <- data.frame(
#'     PatientID = res$results$genefu.robust$BS.all$PatientID,
#'     Subtype = res$results$genefu.robust$BS.all$BS
#' )
#'
#' # Generate the PCA plot
#' p <- Vis_PCA(x = x, out = out)
#' plot(p)
#'
#' # Generate PCA plot with scree plot of eigenvalues
#' p_with_eigen <- Vis_PCA(x = x, out = out, Eigen = TRUE)
#' plot(p_with_eigen)
#'
#' @export

Vis_PCA <- function(x, out, Eigen = FALSE) {
    Subtype.color <- c(
        "Basal" = "red",
        "Her2" = "hotpink",
        "LumA" = "darkblue",
        "LumB" = "skyblue",
        "Normal" = "green"
    )
    pca <- prcomp(t(x), center = TRUE, scale. = TRUE)

    if (Eigen) {
        # Scree plot
        variance <- pca$sdev^2 / sum(pca$sdev^2) * 100
        scree_data <- data.frame(PC = seq_along(variance), Variance = variance)

        screeplot <- ggplot(
            scree_data[seq(1, 10), ],
            aes(x = .data$PC, y = .data$Variance)
        ) +
            geom_bar(
                stat = "identity",
                width = 0.5,
                fill = "steelblue"
            ) +
            geom_line() +
            geom_point(size = 2) +
            scale_x_continuous(breaks = seq(1, 10, 1)) +
            geom_text(
                aes(
                    x = .data$PC,
                    label = paste0(round(.data$Variance, 2), "%")
                ),
                nudge_y = 1
            ) +
            labs(
                x = "Principal Component",
                y = "Percentage of variance Explained"
            ) +
            theme_classic() +
            theme(
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 14)
            )

        return(screeplot)
    } else {
        ## PCA plot
        scores <- as.data.frame(pca$x)
        scores$PatientID <- rownames(scores)
        scores <- dplyr::left_join(scores, out, by = "PatientID")
        rownames(scores) <- scores$PatientID

        pc1_var <- round(100 * summary(pca)$importance[2, 1], 2)
        pc2_var <- round(100 * summary(pca)$importance[2, 2], 2)

        pcaplot <- ggplot(
            data = scores,
            aes(x = .data$PC1, y = .data$PC2, color = .data$Subtype)
        ) +
            geom_point(size = 2) +
            geom_vline(xintercept = 0, linetype = "dashed") +
            geom_hline(yintercept = 0, linetype = "dashed") +
            scale_color_manual(name = "Subtype", values = Subtype.color) +
            labs(
                x = paste0("PC1 (", pc1_var, "% variance)"),
                y = paste0("PC2 (", pc2_var, "% variance)")
            ) +
            theme_minimal() +
            theme(
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 14),
                legend.title = element_text(size = 14),
                legend.text = element_text(size = 12)
            )

        return(pcaplot)
    }
}

#' Pie Chart Visualization of Subtype Distribution
#'
#' @name Vis_pie
#' @description This function generates a pie chart to visualize the
#' distribution of breast cancer subtypes in a cohort, based on the provided
#' `Subtype` data.
#'
#' @param out A data frame containing two columns: `"PatientID"` and
#'   `"Subtype"`. The `"PatientID"` column should contain unique patient
#'   identifiers, and the `"Subtype"` column should specify the assigned subtype
#'   for each patient.
#'
#' @return A `ggplot` object representing a pie chart showing the proportion of
#'   each subtype in the dataset.
#'
#' @examples
#' data("OSLO2EMIT0obj")
#' res <- OSLO2EMIT0obj$res
#'
#' # Prepare data: Subtype information
#' out <- data.frame(
#'     PatientID = res$results$genefu.robust$BS.all$PatientID,
#'     Subtype = res$results$genefu.robust$BS.all$BS
#' )
#'
#' # Generate the pie chart
#' p <- Vis_pie(out = out)
#' plot(p)
#'
#' @export

Vis_pie <- function(out) {
    data <- data.frame(table(out$Subtype))
    data <- data %>%
        dplyr::mutate(
            perc = round(.data$`Freq` / sum(.data$`Freq`) * 100, 2),
            csum = rev(cumsum(rev(.data$perc))),
            pos = .data$perc / 2 + dplyr::lead(.data$csum, 1),
            pos = dplyr::if_else(is.na(.data$pos), .data$perc / 2, .data$pos)
        )


    Subtype.color <- c(
        "Basal" = "red",
        "Her2" = "hotpink",
        "LumA" = "darkblue",
        "LumB" = "skyblue",
        "Normal" = "green"
    )

    plot <- ggplot(data, aes(x = "", y = .data$perc, fill = .data$Var1)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(name = "Subtype", values = Subtype.color) +
        geom_text_repel(
            data = data,
            aes(y = .data$pos, label = paste0(.data$Freq, " (", .data$perc, "%", ")")),
            size = 4.5,
            nudge_x = 0.6,
            show.legend = FALSE,
            segment.color = NA
        ) +
        coord_polar("y", direction = -1) +
        theme_void() +
        theme(
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12)
        )


    return(plot)
}


#' Multi-Method Subtype Heatmap Visualization
#'
#' @name Vis_Multi
#' @description This function generates a heatmap to visualize breast cancer subtypes
#' classified by multiple subtyping methods. It helps users compare how different methods
#' assign subtypes to the same set of samples.
#'
#' @param data Output of the \code{\link{BS_Multi}} function.
#'
#' @return Returns a heatmap visualizing the subtype classifications across multiple methods.
#'
#' @examples
#' data("OSLO2EMIT0obj")
#'
#' # Assuming `OSLO2EMIT0obj$res$res_subtypes` contains multi-method subtype results
#' p <- Vis_Multi(OSLO2EMIT0obj$res$res_subtypes)
#' plot(p)
#'
#' @export

Vis_Multi <- function(data) {
    ## order by subtype
    data <- data[order(data[, ncol(data) - 1]), ]
    ## order by entropy
    data <- data[order(data[, ncol(data)], decreasing = FALSE), ]
    mat <- data[, -ncol(data)]

    Labels <- unique(as.vector(as.matrix(mat)))

    ## preset
    ctg <- data.frame(
        Cat = rep(c("NC-based", "SSP-based"), c(8, 2)),
        row.names = c(
            "parker.original",
            "genefu.scale",
            "genefu.robust",
            "cIHC",
            "cIHC.itr",
            "PCAPAM50",
            "ssBC",
            "ssBC.v2",
            "AIMS",
            "sspbc"
        )
    )
    ctg$Cat <- factor(ctg$Cat, levels = c("NC-based", "SSP-based"))

    ## color
    Subtype.color <- c(
        "Basal" = "red",
        "Her2" = "hotpink",
        "LumA" = "darkblue",
        "LumB" = "skyblue",
        "Normal" = "green"
    )
    col_cat <- setNames(c("#984ea3", "#ff7f00"), c("NC-based", "SSP-based"))

    ## make row annotation
    row_anno <- data.frame(
        Category = ctg[colnames(mat), ],
        row.names = colnames(mat)
    )
    row_anno <- HeatmapAnnotation(
        df = row_anno,
        which = c("row"),
        col = list(Category = col_cat),
        annotation_legend_param = list(
            title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
            gap = unit(2, "points"),
            labels_gp = grid::gpar(fontsize = 10),
            border = "white"
        )
    )
    # ## column annotation
    col_anno <- HeatmapAnnotation(
        which = c("column"),
        Entropy = anno_barplot(data[rownames(data), ncol(data)], bar_width = 1)
    )

    Subtype.color <- Subtype.color[names(Subtype.color) %in% Labels]
    p <- ComplexHeatmap::Heatmap(
        t(as.matrix(mat)),
        name = "Subtypes",
        col = Subtype.color,
        row_names_gp = grid::gpar(fontsize = 10, fontface = "bold"),
        right_annotation = row_anno,
        top_annotation = col_anno,
        show_column_names = FALSE,
        heatmap_legend_param = list(
            title = "Intrinsic Subtype",
            labels = Labels[match(names(Subtype.color), Labels)],
            title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
            gap = unit(2, "points"),
            labels_gp = grid::gpar(fontsize = 10),
            border = "white"
        )
    )

    return(p)
}
