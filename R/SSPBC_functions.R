## adapted from original sspbc codes (Staaf et al., 2022)

#' Function sspbc
#' @import methods
#' @import stats
#' @import grid
#' @importFrom e1071 naiveBayes
#' 
#' @description
#' This function assign classes to breast cancer samples using a selection of provided models. It works on raw gene expression data.
#' Provided models were developed using gene expression data from mRNAseq generated using HiSat/StringTie.
#' For recommended usage see details and examples.
#'
#' @noRd

applySSP <- function(
        tsv,
        ssp = "",
        plot = FALSE,
        txt = FALSE,
        report = FALSE,
        add.is.num = TRUE,
        mylas = 1,
        ssp.name = "",
        gex = NULL,
        id = NULL,
        id.type = "Gene.ID",
        full.out = FALSE,
        output) {
    data_env <- new.env(parent = emptyenv())
    data("Gene.ID.ann", envir = data_env, package = "BreastSubtypeR")
    Gene.ID.ann <- data_env[["Gene.ID.ann"]]
    data("sspbc.models", envir = data_env, package = "BreastSubtypeR")
    sspbc.models <- data_env[["sspbc.models"]]
    data("sspbc.models.fullname", envir = data_env, package = "BreastSubtypeR")
    sspbc.models.fullname <- data_env[["sspbc.models.fullname"]]


    if (report) {
        # report on Gene.ID.ann
        message("Gene.ID.ann nrow:", nrow(Gene.ID.ann), "\n")
    }

    # load gex data as matrix
    # if gex=NULL read data from tsv file
    if (is.null(gex)) {
        # read tsv from single sample gene.tsv file
        mymatrix <- read_StringTie_tsv_FPKM(tsv, id = Gene.ID.ann$Gene.ID)
        # head(mymatrix)
    } else {
        # if data is a matrix make sure it is numeric
        if (is.matrix(gex)) {
            if (is.numeric(gex)) {
                mymatrix <- gex
                rownames(mymatrix) <- id
            } else {
                stop("provided gex matrix is not numeric")
            }
        } else {
            # if data is a vector convert it into matrix
            if (is(gex, "numeric")) {
                mymatrix <- as.matrix(gex)
                rownames(mymatrix) <- id
            }
        }
    }

    # translate gene id
    myid <- translate_id2entrez(
        id = rownames(mymatrix),
        ann = Gene.ID.ann,
        id.type = id.type,
        e = TRUE
    )

    # load ssp model
    if ((ssp.name == "") & (ssp %in% names(sspbc.models.fullname))) {
        aims.gs <- sspbc.models.fullname[[ssp]]
        # loaded ssp object must be named aims.gs
    } else {
        # ssp.name <- "ssp.cc15"
        if (ssp.name %in% names(sspbc.models)) {
            aims.gs <- sspbc.models[[ssp.name]]
        } else {
            stop(
                "Specified ssp.name not among names in sspbc.models: ",
                paste(names(sspbc.models), collapse = ", ")
            )
        }
    }

    # fix issue with $isnumeric
    # add logi (aims.gs$ all.pairs), i.e, one TRUE for each aims.gs$ all.pairs
    if (add.is.num) {
        one_vs_all_tsp <- aims.gs$one.vs.all.tsp[[aims.gs$k]]
        one_vs_all_tsp$isnumeric[names(one_vs_all_tsp$tables)] <- TRUE
    }

    # apply ssp model
    resultslist <- applyAIMS(mymatrix, myid, aims.gs)

    names(resultslist)

    # handle results when a single sample is classified
    if (ncol(mymatrix) == 1) {
        cl <- resultslist$cl[[1]]
        all.probs <- c(resultslist$all.probs[[1]])
        names(all.probs) <- colnames(resultslist$all.probs[[1]])

        if (report) {
            message("assigned class:", cl, "\n")
        }

        # plot only available when a single sample is classified
        if (plot) {
            pdf(
                file = paste(output, "all.probs.pdf", sep = "/"),
                width = 2.53,
                height = 1.19,
                pointsize = 6
            )
            op <- par(
                oma = c(2, 2, 0, 0),
                mar = c(2, 2, 1, 0)
            )
            barplot(all.probs, las = mylas)
            mtext(
                "posterior probability",
                side = 2,
                line = 3,
                cex = 0.9
            )
            # mtext("class", side=1, line=2, cex=0.9)
            par(op)
            dev.off()
        } # end if plot

        # text only available when a single sample is classified
        if (txt) {
            write.table(
                data.frame(all.probs),
                file = paste(output, "all.probs.txt", sep = "/"),
                col.names = FALSE,
                quote = FALSE,
                sep = "\t"
            )
            write.table(
                cl,
                file = paste(output, "class.txt", sep = "/"),
                col.names = FALSE,
                row.names = FALSE,
                quote = FALSE,
                sep = "\t"
            )
        } # end if txt
    }

    # handle results when a matrix with multiple samples are classified
    if (ncol(mymatrix) > 1) {
        cl <- resultslist$cl
        colnames(cl) <- "class"
        all.probs <- resultslist$all.probs[[1]]
        row.names(all.probs) <- row.names(cl)
    }

    if (full.out) {
        return(resultslist)
    } else {
        return(cl)
    }
} # end function




###########################################
#   Johan Vallon-Christersson             #
#   johan.vallon-christersson@med.lu.se   #
###########################################


####
# #' function: read_StringTie_tsv_FPKM
#
# INPUT: gene.tsv file from StingTie. Must contain column holding gene id with
# column name 'Gene ID' or 'Gene.ID' and column 'FPKM' with geneexpression data.
# vector of unique gene id to return data for (genes to be included in returned
# matrix). Must be same type of id found in column 'Gene.ID' in gene.tsv file.
# OUTPUT: returns genematrix with one column (FPKM data sum on gene id) and rows
# for gene id (row names are id)
#
# v1: first implementation the gene.tsv is read using read.delim() and invalid
# characters in column names are translated to ".". That is, column name 'Gene
# ID' will be translated to 'Gene.ID'.


# # examples and manual stuff
#
# # specify input gene.tsv gene.tsv <- "gene.tsv"
#
# # load gene annotation load("Gene.ID.ann.Rdata")
#
# # some gene id to collect some.gene.id <-  Gene.ID.ann$Gene.ID[1:5]
#
# # run function mymatrix <- read_StringTie_tsv_FPKM(tsv=gene.tsv,
# id=some.gene.id, report=TRUE)
#
# head(mymatrix)


#####################################
# function read_StringTie_tsv_FPKM
#####################################

# function
read_StringTie_tsv_FPKM <- function(tsv, id, report = FALSE) {
    # number of columns and rows in genematrix
    my.ncol <- 1
    my.nrow <- length(id)

    # test verify that id are unique
    if (length(id) == length(unique(id))) {
        if (report) {
            message("provided id(s) are unique, ")
        }
    } else {
        stop("mfkr id")
    }

    # create genematrix
    genematrix <- matrix(nrow = my.nrow, ncol = my.ncol)

    # rownames and colname
    row.names(genematrix) <- id
    colnames(genematrix) <- "FPKMsum"

    # head(genematrix)

    # read gene.tsv
    if (report) {
        message("creating matrix, ")
    }

    # read from file
    read.data <- read.delim(file = tsv, as.is = TRUE)
    # head(read.data)

    # aggregate fpkm on gene
    sum_fpkm_gene <- aggregate(FPKM ~ Gene.ID, data = read.data, sum)
    # head(sum_fpkm_gene)

    # find id index in sum_fpkm_gene$Gene.ID
    id.i <- match(id, sum_fpkm_gene$Gene.ID)

    # test verify that all provided id are found
    if (length(which(is.na(id.i) == TRUE)) == 0) {
        # add id.i FPKM to genematrix
        genematrix[, 1] <- sum_fpkm_gene$FPKM[id.i]
    } else {
        stop("provided id not found in tsv")
    }

    if (report) {
        message("done")
    }
    return(genematrix)
} # end function

#####################################
#####################################


###########################################
#    Johan Vallon-Christersson            #
#    johan.vallon-christersson@med.lu.se  #
###########################################


####
# function: translate_id2entrez
#
# INPUT:
# - vector of gene id to translate (id of type used in our StringTie pipeline,
# i.e., 'Gene.ID')
# - for v1.2, other identifiers supported. id is
# expected to be inputted as.character or will be transformed as.caharcter gene
# annotation object uncluding columns 'Gene.ID' and 'EntrezGene' (the gene
# annotation file created for our StringTie pipeline).
# - for v1.2, to use other
# identifier as input these must also be included in the gene annotation object.
# - vector of unique gene id to return data for (genes to be included in
# returned matrix). Must be same type of id found in column 'Gene.ID'
# in gene.tsv file.
# OUTPUT:
# - vector of EntrezGene or vector with EntrezGeneMust with appended
# prefix 'e' (as required by our SSP models).
#
# v1: first implementation Returned vector will have NA for id not found in gene
# annotation object. Note that in gene annotation object all id (Gene.ID)
# without an EntrezGene have character "NA" in column EntrezGene and these will
# be returned
#
# v1.2: added parameter 'id.type' (default Gene.ID) to specify type for inputted
# ID allowed vales: Gene.ID (default), Gene.Name, HGNC, or EntrezGene.
#
# # examples and manual stuff
#
# # load gene annotation
# load("Gene.ID.ann.Rdata")
#
# # some gene id to translate examples
# some.gene.id <-  Gene.ID.ann$Gene.ID[c(1:5)]
# # some.gene.id <-  Gene.ID.ann$Gene.ID[c(1:5, 19625)]
# # some.gene.id <- c(some.gene.id, "gaga", some.gene.id)
# # type <- "HGNC"
# some.gene.id <-  Gene.ID.ann$Gene.Name[c(1:5, 19625)]
# # some.gene.id <- c(some.gene.id, "gaga", some.gene.id)
# # type <- "EntrezGene"
# some.gene.id <-  Gene.ID.ann$EntrezGene[c(1:5, 19625)]
# # some.gene.id <- c(some.gene.id, "gaga", some.gene.id)
#
# # run function
# myid <- translate_id2entrez(
# id=some.gene.id,
# ann=Gene.ID.ann,
# e=TRUE,
# report=TRUE)
#
# head(myid)


#####################################
# function translate_id2entrez
#####################################

# function
translate_id2entrez <- function(
        id,
        ann,
        id.type = "Gene.ID",
        e = FALSE,
        report = FALSE) {
    if (id.type %in% c("Gene.ID", "Gene.Name", "HGNC", "EntrezGene")) {
        # find first index for id in ann
        id.i <- match(as.character(id), ann[, id.type])
    } else {
        stop("specified type not supported")
    }

    if (report) {
        # test verify that all provided id are found
        if (length(which(is.na(id.i) == TRUE)) == 0) {
            message("all id found in ann, ")
        } else {
            message(length(which(is.na(id.i) == TRUE)), "id not found in ann, ")
        }
    }

    # get EntrezGene and assign id as names
    entrez <- ann$EntrezGene[id.i]
    names(entrez) <- id

    entrezNA <- length(which(entrez[is.na(entrez) == FALSE] == "NA"))

    if (report) {
        if (entrezNA > 0) {
            message(entrezNA, "found id have EntrezGene NA, ")
        } else {
            message("all found id have an EntrezGene, ")
        }
    }

    # add prefix
    if (e) {
        if (report) {
            message("adding prefix e, ")
        }

        entrez[is.na(entrez) == FALSE &
            entrez != "NA"] <- paste("e", entrez[is.na(entrez) == FALSE &
            entrez != "NA"], sep = "")
    }

    if (report) {
        message("done")
    }

    return(entrez)
} # end function

#####################################
#####################################

#### SSP_functions.R ####
.smaller <- function(x, y) {
    x < y
}

## Functions necessary to run AIMS
.comp.sel.pairs <- function(dataset, sel.pairs, func = .smaller) {
    to.ret <- matrix(NA, nrow = length(sel.pairs), ncol(dataset$exprs))
    for (spi in seq(1, length(sel.pairs), 1)) {
        ss.cur <- strsplit(sel.pairs[spi], "<")[[1]]
        gene1 <- which(dataset$GeneName == ss.cur[1])
        gene2 <- which(dataset$GeneName == ss.cur[2])
        # stopifnot(length(gene1) == 1 & length(gene2) == 1)
        if (length(gene1) == 1 & length(gene2) == 1) {
            to.ret[spi, ] <- func(
                dataset$exprs[gene1, ],
                dataset$exprs[gene2, ]
            )
        } else {
            message(
                "You are missing the pair or have more than one",
                sel.pairs[spi],
                "in",
                dataset$name
            )
        }
    }
    to.ret <- apply(to.ret, 2, as.numeric)
    rownames(to.ret) <- sel.pairs
    to.ret
}

.one.vs.all.tsp <- function(D, GeneName, one.vs.all.tsp) {
    ## Need to add some QC
    ## First compute
    train.pairs <- .comp.sel.pairs(
        list(
            exprs = D,
            GeneName = GeneName
        ),
        one.vs.all.tsp$all.pairs
    )

    classes <- matrix("", ncol = length(one.vs.all.tsp$k), nrow = ncol(D))
    prob <- matrix(0, ncol = length(one.vs.all.tsp$k), nrow = ncol(D))
    colnames(classes) <- colnames(prob) <- as.character(one.vs.all.tsp$k)
    rownames(classes) <- rownames(prob) <- colnames(D)
    nb.d <- data.frame(t(train.pairs))
    all.probs <- list()
    for (ki in one.vs.all.tsp$k) {
        message(sprintf("Current k = %d", ki))
        prob.train <- predict(
            one.vs.all.tsp$one.vs.all.tsp[[ki]],
            nb.d,
            type = "raw"
        )
        cur.cl <- apply(prob.train, 1, function(prob.cur) {
            colnames(prob.train)[which.max(prob.cur)]
        })
        cur.prob <- apply(prob.train, 1, function(prob.cur) {
            max(prob.cur)
        })
        prob[, as.character(ki)] <- cur.prob
        classes[, as.character(ki)] <- cur.cl
        all.probs[[as.character(ki)]] <- prob.train
    }
    invisible(list(
        cl = classes,
        prob = prob,
        all.probs = all.probs,
        rules.matrix = train.pairs
    ))
}

.get.all.pairs.genes <- function(all.pairs) {
    genes <- c()
    for (cp in strsplit(all.pairs, "<")) {
        genes <- c(genes, cp)
    }
    unique(genes)
}
## Remove the duplicated Entrez by keeping the most higly expressed This is
## closer to the single sample selection D is a raw gene expression matrix rows
## == genes and columns patients
.removeDuplicatedEntrezPerPatients <- function(D, EntrezID, probes) {
    ## Maybe we have nothing to do already
    if (all(!duplicated(EntrezID))) {
        return(list(dataset = D, EntrezID = EntrezID))
    } else {
        uniqEntrez <- sort(unique(EntrezID))
        newD <- matrix(0, nrow = length(uniqEntrez), ncol = ncol(D))
        if (!missing(probes)) {
            sel.probes <- matrix("", nrow = length(uniqEntrez), ncol = ncol(D))
        }
        for (i in seq(1, ncol(D), 1)) {
            curD <- D[, i]
            curEntrez <- EntrezID
            oi <- order(curD, decreasing = TRUE) ## order by raw expression
            curD <- curD[oi]
            curEntrez <- curEntrez[oi]
            cur.sel <- !duplicated(curEntrez) ## remove duplicated
            curD <- curD[cur.sel]
            curEntrez <- curEntrez[cur.sel]
            reorder <- match(uniqEntrez, curEntrez)
            newD[, i] <- curD[reorder]
            if (!missing(probes)) {
                sel.probes[, i] <- probes[oi][cur.sel][reorder]
            }
        }
        colnames(newD) <- colnames(D)
        if (!missing(probes)) {
            colnames(sel.probes) <- colnames(D)
            return(list(
                dataset = newD,
                EntrezID = uniqEntrez,
                probes = sel.probes
            ))
        } else {
            return(list(dataset = newD, EntrezID = uniqEntrez))
        }
    }
}

.apply.nbc <- function(D, EntrezID, sel.nbc) {
    ## Verify the number of rows of D and EntrezIDs have the same length
    if (nrow(D) != length(EntrezID)) {
        stop(
            sprintf(
                "You need the same number of rows and EntrezID.
                Right now nrow(D) = %d and length(EntrezID) = %d",
                nrow(D),
                length(EntrezID)
            )
        )
    }
    ## AIMS needs to be applied on expression values > 0.
    ## Require 95% of the values to be > 0
    if (!all(apply(D, 2, function(x) {
        (sum(x < 0, na.rm = TRUE) / length(x)) < 0.05
    }))) {
        stop(
            "AIMS needs to be applied on expressionn values > 0.
            Did you gene-centered your matrix D? You should not."
        )
    }

    ## Verify if D is a numeric matrix
    if (!all(apply(D, 2, is.numeric))) {
        stop(
            "Verify D is a numeric matrix.
            apply(D,2,as.numeric) could do the job
            or verify the first column doesn't contain probe ids"
        )
    }

    D <- apply(D, 2, as.numeric)
    EntrezID <- as.character(EntrezID)
    sel.ids.nb <- .get.all.pairs.genes(sel.nbc$all.pairs)
    sel.gn <- EntrezID %in% sel.ids.nb
    D <- D[sel.gn, , drop = FALSE]
    Entrez <- EntrezID[sel.gn]
    col.D.test <- .removeDuplicatedEntrezPerPatients(D, Entrez)
    pred.test <- .one.vs.all.tsp(
        D = col.D.test$dataset,
        GeneName = col.D.test$EntrezID,
        sel.nbc
    )

    ## Add more information to the output variable
    pred.test$data.used <- col.D.test$dataset
    pred.test$EntrezID.used <- col.D.test$EntrezID

    invisible(pred.test)
}

applyAIMS <- function(eset, EntrezID, AIMSmodel) {
    D <- NA
    if (!is(eset, "ExpressionSet") & !is(eset, "matrix")) {
        stop(
            "eset argument should be either an
        ExpressionSet (Biobase) or a numerical matrix"
        )
    }

    if (is(eset, "ExpressionSet")) {
        D <- exprs(eset)
    } else if (is(eset, "matrix")) {
        D <- eset
    }

    .apply.nbc(D, EntrezID, AIMSmodel)
}
