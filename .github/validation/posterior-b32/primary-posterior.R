# Independent reconstruction from the pristine numeric model, in e1071 table order.
native_oracle <- function(nb, rules) {
    nd <- data.frame(t(rules))
    j <- match(names(nb$tables), names(nd))
    stopifnot(!anyNA(j), !anyDuplicated(j))
    C <- length(nb$apriori)
    result <- t(vapply(seq_len(nrow(nd)), function(i) {
        if (all(is.na(nd[i, j]))) return(rep(NA_real_, C))
        density <- vapply(seq_along(j), function(v) {
            x <- nd[i, j[v]]
            if (is.na(x)) return(rep(1, C))
            tab <- nb$tables[[v]]
            sd <- tab[, 2L]; sd[sd <= 0] <- .001
            p <- dnorm(x, tab[, 1L], sd); p[p <= 0] <- .001
            p
        }, numeric(C))
        L <- log(nb$apriori) + apply(log(density), 1L, sum)
        vapply(L, function(lp) 1/sum(exp(L-lp)), numeric(1L))
    }, numeric(C)))
    dimnames(result) <- list(NULL, names(nb$apriori))
    result
}
posterior_rows <- list()
primary_posterior <- function(actual, reference, dataset) {
    before <- .Random.seed; kind <- RNGkind()
    approved <- issue199_ror_approved(dataset)
    stopifnot(identical(reference, approved$res, num.eq=FALSE, single.NA=FALSE))
    issue199_posterior_arithmetic()
    models <- issue199_posterior_models()
    policy <- issue199_portability_policy()
    for (method in c("AIMS","sspbc")) {
        model <- models[[method]]; key <- as.character(model$k)
        a <- actual$results[[method]]; b <- reference$results[[method]]
        g <- issue199_posterior_method(a, b, model, method, dataset)
        p <- a$all.probs[[key]]; q <- b$all.probs[[key]]
        nb <- model$one.vs.all.tsp[[model$k]]
        wanted <- native_oracle(nb, a$rules.matrix)
        candidate <- issue199_posterior_oracle(nb, a$rules.matrix)
        stopifnot(identical(wanted, p, num.eq=FALSE, single.NA=FALSE),
            identical(candidate, wanted, num.eq=FALSE, single.NA=FALSE),
            identical(a$rules.matrix, b$rules.matrix, num.eq=FALSE),
            identical(names(nb$tables), make.names(rownames(a$rules.matrix))),
            identical(a$BS.all, b$BS.all), identical(a$cl, b$cl),
            identical(dim(p), c(if(dataset=="OSLO2EMIT0obj") 102L else 100L, 5L)))
        path <- paste0("/res/results/",method,"/all.probs/",key)
        rule <- policy[policy$dataset==dataset & policy$path==path, ]
        stopifnot(nrow(rule)==1L)
        delta <- abs(p-q); scale <- pmax(abs(p),abs(q)); normal <- q >= .Machine$double.xmin
        ceiling <- if(method=="AIMS") 2^-35 else 2^-33
        # Record every possible x from both possible rule values, including the floor.
        xs <- unlist(lapply(nb$tables,function(tab) {
            sd <- tab[,2]; sd[sd<=0] <- .001
            c(abs((0-tab[,1])/sd),abs((1-tab[,1])/sd))
        }),use.names=FALSE)
        issue199_posterior_x(xs)
        row <- data.frame(dataset,method,key,path,nrows=nrow(p),ncols=ncol(p),cells=length(p),
            oracle_bitwise=length(p),changed=sum(p!=q),max_abs=max(delta),
            max_relative=max(delta[normal]/scale[normal]),old_rtol=rule$rtol,
            old_exceeded=sum(delta>rule$atol+rule$rtol*scale),
            ceiling,new_exceeded=sum(delta[normal]>ceiling*scale[normal]),
            x_entries=length(xs),max_unfloored_x=max(xs[xs<=32]),floored=sum(xs>=64),
            max_sum=g$max_sum,max_span=g$max_span,min_gap=g$min_gap,
            overflow_margin=g$overflow_margin,nonzero_margin=g$nonzero_margin,
            subnormal_margin=g$subnormal_margin,normal_margin=g$normal_margin,
            zeros=g$zero,subnormal=g$subnormal,actual_hash=numeric_fingerprint(p),
            oracle_hash=numeric_fingerprint(wanted),reference_hash=numeric_fingerprint(q))
        stopifnot(row$new_exceeded==0L,identical(before,.Random.seed),identical(kind,RNGkind()))
        posterior_rows[[paste(dataset,method)]] <<- row
        write_tsv(do.call(rbind,posterior_rows),file.path(out,"b32-routes.tsv"))
    }
}
