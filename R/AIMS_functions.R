## adjusted from AIMS v1.38.0
.smaller_AIMS <- function(x, y) {
    x < y
}

## Functions necessary to run AIMS
## stopifnot(length(gene1) == 1 & length(gene2) == 1)
.comp.sel.pairs_AIMS <- function(dataset, sel.pairs, func = .smaller_AIMS) {
    to_ret <- matrix(NA, nrow = length(sel.pairs), ncol(dataset$exprs))

    for (spi in seq(1, length(sel.pairs), 1)) {
        ss.cur <- strsplit(sel.pairs[spi], "<")[[1]]
        gene1 <- which(dataset$GeneName == ss.cur[1])
        gene2 <- which(dataset$GeneName == ss.cur[2])
        if (length(gene1) == 1 & length(gene2) == 1) {
            to_ret[spi, ] <- func(
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
    to_ret <- apply(to_ret, 2, as.numeric)
    rownames(to_ret) <- sel.pairs
    to_ret
}

.one.vs.all.tsp_AIMS <- function(D, GeneName, one.vs.all.tsp) {
    ## Need to add some QC
    ## First compute
    train.pairs <- .comp.sel.pairs_AIMS(
        list(exprs = D, GeneName = GeneName),
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

        ## need to add this for the new version of e1071 (1.7-1)
        ## there is now a new isnumeric field in the naive bayes object
        ## we need to update our model consequently.
        isnumeric <- list()
        for (ni in names(one.vs.all.tsp$one.vs.all.tsp[[ki]]$tables)) {
            isnumeric[[ni]] <- TRUE
        }
        one.vs.all.tsp$one.vs.all.tsp[[ki]]$isnumeric <- isnumeric

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
        ## Count observed comparisons (both 0 and 1) used by this model.
        model.rules <- match(
            names(one.vs.all.tsp$one.vs.all.tsp[[ki]]$tables), names(nb.d)
        )
        evaluable <- colSums(!is.na(train.pairs[model.rules, , drop = FALSE]))
        no.rules <- evaluable == 0L
        cur.cl[no.rules] <- NA_character_
        cur.prob[no.rules] <- NA_real_
        prob.train[no.rules, ] <- NA_real_
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

.get.all.pairs.genes_AIMS <- function(all.pairs) {
    genes <- c()
    for (cp in strsplit(all.pairs, "<")) {
        genes <- c(genes, cp)
    }
    unique(genes)
}

## Remove the duplicated Entrez by keeping the most higly expressed.
## This is closer to the single sample selection.
## D is a raw gene expression matrix rows == genes and columns patients
.removeDuplicatedEntrezPerPatients_AIMS <- function(D, EntrezID, probes) {
    ## Maybe we have nothing to do already
    if (all(!duplicated(EntrezID))) {
        return(list(dataset = D, EntrezID = EntrezID))
    } else {
        uniqEntrez <- sort(unique(EntrezID))
        newD <- matrix(0, nrow = length(uniqEntrez), ncol = ncol(D))
        if (!missing(probes)) {
            sel.probes <- matrix("",
                nrow = length(uniqEntrez),
                ncol = ncol(D)
            )
        }
        for (i in seq(1, ncol(D), 1)) {
            curD <- D[, i]
            curEntrez <- EntrezID
            oi <- order(curD, decreasing = TRUE)
            ## order by raw expression
            curD <- curD[oi]
            curEntrez <- curEntrez[oi]
            cur.sel <- !duplicated(curEntrez)
            ## remove duplicated
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

.apply.nbc_AIMS <- function(D, EntrezID, sel.nbc) {
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
    if (nrow(D) > 0L && !all(apply(D, 2, function(x) {
        (sum(x < 0, na.rm = TRUE) / length(x)) < 0.05
    }))) {
        stop(
            "AIMS needs to be applied on expressionn values > 0.
            Did you gene - centered your matrix D ? You should not."
        )
    }

    ## Verify if D is a numeric matrix
    if (!all(apply(D, 2, is.numeric))) {
        stop(
            "Verify D is a numeric matrix.
            apply(D,2,as.numeric) could do the job or
            verify the first column doesn't contain probe ids"
        )
    }

    ## Keep zero/one-row inputs as matrices when model genes are absent.
    D <- matrix(as.numeric(D), nrow = nrow(D), ncol = ncol(D),
        dimnames = list(NULL, colnames(D)))
    EntrezID <- as.character(EntrezID)
    sel.ids.nb <- .get.all.pairs.genes_AIMS(sel.nbc$all.pairs)
    sel.gn <- EntrezID %in% sel.ids.nb
    D <- D[sel.gn, , drop = FALSE]
    Entrez <- EntrezID[sel.gn]
    col.D.test <- .removeDuplicatedEntrezPerPatients(D, Entrez)
    pred.test <- .one.vs.all.tsp_AIMS(
        D = col.D.test$dataset,
        GeneName = col.D.test$EntrezID,
        sel.nbc
    )

    ## Add more information to the output variable
    pred.test$data.used <- col.D.test$dataset
    pred.test$EntrezID.used <- col.D.test$EntrezID

    invisible(pred.test)
}

applyAIMS_AIMS <- function(eset, EntrezID) {
    D <- NA
    if (!is(eset, "ExpressionSet") & !is(eset, "matrix")) {
        stop("eset should be either an ExpressionSet (Biobase) or a matrix")
    }

    if (is(eset, "ExpressionSet")) {
        D <- exprs(eset)
    } else if (is(eset, "matrix")) {
        D <- eset
    }

    data("AIMSmodel", package = "BreastSubtypeR", envir = environment())
    .apply.nbc_AIMS(D, EntrezID, AIMSmodel)
}
