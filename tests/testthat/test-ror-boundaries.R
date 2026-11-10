## #172: author-approved inclusive medium intervals for the four original
## research ROR models. Fixed IEEE-754 inputs/targets; no tolerance, rounding,
## reverse scaling, production classifier oracle, or runtime fixture search.
ror172_scores <- c("ROR-S (Subtype Only)", "ROR-P (Subtype + Prolif)",
    "ROR-C (Subtype + Clinic)", "ROR-PC (Subtype + Clinic + Prolif)")
ror172_groups <- c("ROR-S Group (Subtype Only)", "ROR-P Group (Subtype + Prolif)",
    "ROR-C Group (Subtype + Clinic)", "ROR-PC Group (Subtype + Clinic + Prolif)")
ror172_add <- c("ROR-PC (Subtype + Clinic + Prolif.Subtype)",
    "ROR-PC Group (Subtype + Clinic + Prolif.Subtype)")

ror172_cases <- function() {
    data.frame(model = rep(1:4, each = 6),
        boundary = rep(rep(c("lower", "upper"), each = 3), 4),
        side = rep(c("below", "equal", "above"), 8),
        cutoff = rep(c(-.15, .1, -.25, .1, -.1, .2, -.2, .2), each = 3),
        raw = c(
            -0x1.3333333333334p-3, -.15, -0x1.3333333333332p-3,
            0x1.9999999999999p-4, .1, 0x1.999999999999bp-4,
            -0x1.0000000000001p-2, -.25, -0x1.fffffffffffffp-3,
            0x1.9999999999999p-4, .1, 0x1.999999999999bp-4,
            -0x1.999999999999bp-4, -.1, -0x1.9999999999999p-4,
            0x1.9999999999999p-3, .2, 0x1.999999999999bp-3,
            -0x1.999999999999bp-3, -.2, -0x1.9999999999999p-3,
            0x1.9999999999999p-3, .2, 0x1.999999999999bp-3),
        expected = rep(c("low", "med", "med", "med", "med", "high"), 4))
}

ror172_input <- function() {
    ids <- paste0("boundary ", seq_len(24), "/sample")
    genes <- c("ESR1", "ERBB2", "CCNB1", "MKI67", "ANLN", "FOXA1", "PGR")
    td <- matrix(c(-1, 1, 0, 0, 0, 0, 0), nrow = 7, ncol = 24,
        dimnames = list(genes, ids))
    ## These columns have mean 0 and SD 1 exactly: the represented original
    ## panel has mean -1, the separate 18-gene panel has mean 0. All other
    ## columns also have positive SD, with original-panel mean exactly 0.
    td[, c(7:9, 19:21)] <- c(1, 1, -1, -1, 1, -1, 0)
    d <- matrix(0, nrow = 24, ncol = 5,
        dimnames = list(ids, c("Basal", "Her2", "LumA", "LumB", "Normal")))
    d[, "Her2"] <- c(0, 0, 0, 0, 0, -0x1.0624dd2f1a9fcp-10,
        0, 0, 0, 0x1.f5c28f5c28f5cp-2, .5, 0x1.f5c28f5c28f5cp-2,
        rep(0, 6), rep(1, 6))
    d[, "LumA"] <- c(
        0x1.b3ec599dafffdp-2, 0x1.b3ec599dafffbp-2, 0x1.b3ec599dafffap-2,
        -0x1.229d9113caaa7p-2, -0x1.229d9113caaa8p-2, -0x1.22fa51545abebp-2,
        -0x1.dc46c0cff23f8p-1, -0x1.dc46c0cff23fdp-1, -0x1.dc46c0cff24p-1,
        -0x1.638faea5bcc1p-1, -0x1.5fd609d6beb3dp-1, -0x1.638faea5bcc12p-1,
        0x1.889462bbd8d5p-2, 0x1.889462bbd8d4fp-2, 0x1.889462bbd8d4ep-2,
        -0x1.889462bbd8d4ep-1, -0x1.889462bbd8d4fp-1, -0x1.889462bbd8d5p-1,
        -0x1.2f14093fbf214p-1, -0x1.2f14093fbf216p-1, -0x1.2f14093fbf219p-1,
        -0x1.dac96b3c5dbfbp-1, -0x1.dac96b3c5dbfep-1, -0x1.dac96b3c5dc01p-1)
    d[c(10, 12), "LumB"] <- 0x1.44390eddbd6d2p-52
    d[22:24, "LumB"] <- 0x1.999999999999ap-1
    d4 <- d[, 1:4]
    d4[,] <- 0
    d4[, "LumA"] <- seq(-.9, .9, length.out = 24)
    list(out = list(testData = td, distances = d, dist.RORSubtype = d4,
        centroids = matrix(seq_len(35) / 10, nrow = 7,
            dimnames = list(genes, colnames(d))),
        predictions = setNames(rep("LumA", 24), ids),
        predictions.FourSubtype = setNames(rep("LumA", 24), ids)),
        ## Deliberately reverse clinical order to exercise sample matching.
        clinical = data.frame(PatientID = ids, TSIZE = 0,
            NODE = rep(c(0, 1, 3, 4, NA_real_, 0), 4))[24:1, ])
}

## Independent numerical oracle, transcribed from the distributed research
## coefficients. Use mean/sd, not the production standardize/panel helpers.
ror172_raw <- function(input) {
    d <- input$out$distances
    p <- apply(input$out$testData, 2, function(x) {
        mean(((x - mean(x)) / sd(x))[c("CCNB1", "MKI67")])
    })
    t <- input$clinical$TSIZE[match(rownames(d), input$clinical$PatientID)]
    cbind(
        .04210193 * d[, 1] + .12466938 * d[, 2] +
            -.35235561 * d[, 3] + .14213283 * d[, 4],
        -.0009299747 * d[, 1] + .0692289192 * d[, 2] +
            -.0951505484 * d[, 3] + .0493487685 * d[, 4] + .3385116381 * p,
        .0442770 * d[, 1] + .1170297 * d[, 2] +
            -.2608388 * d[, 3] + .1055908 * d[, 4] + .1813751 * t,
        -.009383416 * d[, 1] + .073725503 * d[, 2] +
            -.090436516 * d[, 3] + .053013865 * d[, 4] +
            .131605960 * t + .327259375 * p)
}

test_that("all 24 fixtures reach exact raw cutoffs and immediate IEEE neighbours", {
    cases <- ror172_cases()
    input <- ror172_input()
    raw <- ror172_raw(input)
    expect_identical(unname(raw[cbind(1:24, cases$model)]), cases$raw)
    expect_true(all(apply(input$out$testData, 2, sd) > 0))
    expect_true(all(abs(input$out$distances) <= 1))
    expect_identical(cases$raw[cases$side == "equal"],
        c(-.15, .1, -.25, .1, -.1, .2, -.2, .2))
    ## ULP gaps are fixed independently of the fixture arithmetic. At -0.25
    ## they differ on either side because the exponent changes there.
    gaps <- matrix(c(2^-55, 2^-55, 2^-56, 2^-56,
        2^-54, 2^-55, 2^-56, 2^-56, 2^-56, 2^-56,
        2^-55, 2^-55, 2^-55, 2^-55, 2^-55, 2^-55), nrow = 2)
    for (i in 1:8) {
        triple <- cases$raw[(3 * i - 2):(3 * i)]
        expect_identical(diff(triple), gaps[, i])
    }
    ## Distinct raw inputs can collapse to the same exported double. This
    ## guards against testing/classifying from reverse-scaled output.
    expect_identical(100 * (cases$raw[4] + .35) / .85,
        100 * (cases$raw[5] + .35) / .85)
})

test_that("original ROR equality is med with clinical data on or off", {
    input <- ror172_input()
    snapshot <- serialize(input, NULL)
    cases <- ror172_cases()
    raw <- ror172_raw(input)
    for (clinical in c(FALSE, TRUE)) for (subtype in c(FALSE, TRUE)) {
        ans <- BreastSubtypeR:::RORgroup(input$out, input$clinical,
            hasClinical = clinical, Subtype = subtype)
        models <- if (clinical) 1:4 else 1:2
        for (m in models) {
            at <- which(cases$model == m)
            expect_identical(ans[[ror172_groups[m]]][at], cases$expected[at])
            ## All reported numerical scores, including nontarget models,
            ## retain the exact formula and scaling, without tolerance.
            expect_identical(ans[[ror172_scores[m]]],
                unname(100 * (raw[, m] + .35) / .85))
            expect_type(ans[[ror172_groups[m]]], "character")
            expect_type(ans[[ror172_scores[m]]], "double")
        }
        fields <- as.vector(rbind(ror172_scores[models], ror172_groups[models]))
        expect_identical(names(ans), c("patientID", colnames(input$out$distances),
            "Call", if (subtype) "Call.Subtype", "Confidence", fields,
            if (clinical) ror172_add, "ER", "HER2"))
        expect_identical(rownames(ans), names(input$out$predictions))
        expect_identical(ans$patientID, names(input$out$predictions))
        expect_identical(ans$Call, unname(input$out$predictions))
        if (subtype) expect_identical(ans$Call.Subtype,
            unname(input$out$predictions.FourSubtype))
        expect_identical(as.matrix(ans[, colnames(input$out$distances)]),
            input$out$distances)
        expect_identical(ans$ER, unname(input$out$testData["ESR1", ]))
        expect_identical(ans$HER2, unname(input$out$testData["ERBB2", ]))
        if (clinical) {
            no_clinical <- BreastSubtypeR:::RORgroup(input$out, input$clinical,
                hasClinical = FALSE, Subtype = subtype)
            expect_identical(ans[, names(no_clinical)], no_clinical)
            ## The additional 18-gene panel is MKI67 + ANLN here, whose raw
            ## mean is zero. Its score and NODE policy remain independent.
            additional <- 54.7690 * (-.3172 * input$out$dist.RORSubtype[, 3] +
                .1981 * 0 + .1133 * 0 + .8826)
            expect_identical(ans[[ror172_add[1]]], unname(additional))
            node <- input$clinical$NODE[match(ans$patientID, input$clinical$PatientID)]
            expected <- ifelse(is.na(node), NA_character_, ifelse(node >= 4, "high",
                ifelse(node == 0,
                    ifelse(additional <= 40, "low", ifelse(additional <= 60, "med", "high")),
                    ifelse(additional <= 15, "low", ifelse(additional <= 40, "med", "high")))))
            expect_identical(ans[[ror172_add[2]]], unname(expected))
        }
        expect_identical(serialize(input, NULL), snapshot)
    }
})

test_that("original ROR boundary changes preserve missing and nonfinite results", {
    input <- ror172_input()
    at <- c(2, 8, 14, 20, 4) # four lower equalities and a NODE >= 4 control
    for (clinical in c(FALSE, TRUE)) {
        baseline <- BreastSubtypeR:::RORgroup(input$out, input$clinical,
            hasClinical = clinical)
        for (kind in c("correlation", "proliferation", "size")) {
            for (value in c(NA_real_, NaN, -Inf, Inf)) {
                changed <- input
                if (kind == "correlation") changed$out$distances[at, "Basal"] <- value
                if (kind == "proliferation") changed$out$testData[, at] <- value
                if (kind == "size") changed$clinical$TSIZE[
                    match(names(input$out$predictions)[at], changed$clinical$PatientID)] <- value
                score <- function() BreastSubtypeR:::RORgroup(changed$out, changed$clinical,
                    hasClinical = clinical)
                if (clinical && kind == "size" && is.infinite(value)) {
                    expect_warning(ans <- score(), class = "bsr_invalid_tsize_warning")
                } else {
                    ans <- score()
                }
                affected <- switch(kind, correlation = 1:4, proliferation = c(2, 4), size = 3:4)
                models <- if (clinical) 1:4 else 1:2
                for (m in intersect(affected, models)) {
                    expect_identical(ans[[ror172_scores[m]]][at], rep(NA_real_, length(at)))
                    expect_identical(ans[[ror172_groups[m]]][at], rep(NA_character_, length(at)))
                }
                unaffected <- as.vector(rbind(ror172_scores[setdiff(models, affected)],
                    ror172_groups[setdiff(models, affected)]))
                expect_identical(ans[, unaffected, drop = FALSE],
                    baseline[, unaffected, drop = FALSE])
                expect_identical(ans[-at, ], baseline[-at, ])
                expect_identical(names(ans), names(baseline))
                expect_identical(vapply(ans, typeof, ""), vapply(baseline, typeof, ""))
                if (clinical && kind == "correlation") {
                    expect_identical(ans[, ror172_add], baseline[, ror172_add])
                }
                if (clinical && kind == "proliferation") {
                    expect_identical(ans[[ror172_add[1]]][4], NA_real_)
                    expect_identical(ans[[ror172_add[2]]][4], "high")
                }
            }
        }
        ## Zero within-sample variance still means undefined original
        ## proliferation, never a fallback medium score/group.
        constant <- input
        constant$out$testData[, at] <- 1
        ans <- BreastSubtypeR:::RORgroup(constant$out, constant$clinical,
            hasClinical = clinical)
        for (m in if (clinical) c(2, 4) else 2) {
            expect_identical(ans[[ror172_scores[m]]][at], rep(NA_real_, length(at)))
            expect_identical(ans[[ror172_groups[m]]][at], rep(NA_character_, length(at)))
        }
    }
})
