# Freeze both alignment calls and standardization together with sspPredict.
# Use a private environment with base and stats bindings, independent of any
# production classifiers, batching helpers, or mocked package bindings.
ssp_scalar_reference <- local({
    ref <- new.env(parent = asNamespace("stats"))
    sys.source(testthat::test_path("fixtures", "ssp-scalar-17fb9d0.R"), ref)
    ref$sspPredict
})

ssp_strict_equal <- function(x, y) {
    identical(x, y, num.eq = FALSE, single.NA = FALSE, attrib.as.set = FALSE)
}

ssp_rng_state <- function() {
    list(seed = get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
         kind = RNGkind())
}

ssp_restore_rng <- function(state) {
    if (!identical(RNGkind(), state$kind))
        do.call(RNGkind, as.list(state$kind))
    if (!is.null(state$seed)) {
        assign(".Random.seed", state$seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
    }
}

# The same invocation expression is used for both functions so an error raised
# at the outermost call also has comparable call text. Inner condition calls
# are preserved as language objects.
ssp_capture <- function(fn, serialized_inputs, std, distm, Subtype, state, info) {
    inputs <- unserialize(serialized_inputs)
    trace <- list()
    record <- function(type, condition) {
        trace[[length(trace) + 1L]] <<- list(
            type = type, classes = class(condition),
            message = conditionMessage(condition), call = conditionCall(condition))
    }
    ssp_restore_rng(state)
    value <- tryCatch(withCallingHandlers(
        fn(inputs$x, inputs$y, std = std, distm = distm, Subtype = Subtype),
        warning = function(w) {
            record("warning", w)
            # With warn=2, allow R's normal warning-to-error conversion.
            if (getOption("warn") < 2L) invokeRestart("muffleWarning")
        },
        message = function(m) {
            record("message", m)
            invokeRestart("muffleMessage")
        }),
        error = function(e) {
            record("error", e)
            NULL
        })
    after <- ssp_rng_state()
    testthat::expect_true(ssp_strict_equal(after, state),
        info = paste(info, "RNG state after call"))
    testthat::expect_true(ssp_strict_equal(inputs, unserialize(serialized_inputs)),
        info = paste(info, "input values and attributes after call"))
    failed <- length(trace) > 0L && identical(trace[[length(trace)]]$type, "error")
    list(value = value, trace = trace, failed = failed)
}

expect_ssp_scalar_equivalence <- function(x, y, std = FALSE, distm = "spearman",
                                          Subtype = TRUE, info = NULL) {
    caller <- ssp_rng_state()
    on.exit(ssp_restore_rng(caller), add = TRUE)
    original <- serialize(list(x = x, y = y), NULL)
    set.seed(713L)
    seeded <- ssp_rng_state()
    absent <- list(seed = NULL, kind = seeded$kind)
    states <- list(seeded = seeded, absent = absent)
    for (name in names(states)) {
        state <- states[[name]]
        label <- paste(info, "seed", name)
        reference <- ssp_capture(ssp_scalar_reference, original, std, distm, Subtype,
                                 state, paste(label, "reference"))
        actual <- ssp_capture(BreastSubtypeR:::sspPredict, original, std, distm, Subtype,
                              state, paste(label, "candidate"))
        testthat::expect_true(ssp_strict_equal(actual$trace, reference$trace),
            info = paste(label, "ordered condition classes, messages and calls"))
        testthat::expect_true(ssp_strict_equal(actual$value, reference$value),
            info = paste(label, "complete return value"))
        if (name == "seeded") seeded_result <- actual
    }
    invisible(seeded_result)
}
