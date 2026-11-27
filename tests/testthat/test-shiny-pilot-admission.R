test_that("hosted pilot admits only scalar method and iteration arguments and recovers", {
    withr::local_options(BreastSubtypeR.profile="example-only", BreastSubtypeR.maxIterations=100L,
                         BreastSubtypeR.workerAddress="127.0.0.1")
    limit <- BreastSubtypeR:::.shiny_iteration_limit
    expect_identical(limit(), 100L)
    withr::with_options(list(BreastSubtypeR.maxIterations=50L), expect_error(limit(), "fixed ceiling"))
    expect_equal(formals(BreastSubtypeR::BS_cIHC.itr)$iteration, 100)
    app <- new.env(parent=asNamespace("shiny"))
    sys.source(system.file("ShinyBreastSubtypeR/server.R", package="BreastSubtypeR", mustWork=TRUE), app)
    received <- list()
    fit <- BreastSubtypeR::BS_cIHC.itr
    testthat::local_mocked_bindings(BS_cIHC.itr=function(se_obj, iteration, ...) {
        received[[length(received)+1L]] <<- iteration
        fit(se_obj=se_obj, iteration=iteration, ...)
    }, .package="BreastSubtypeR")
    shiny::testServer(app$server, {
        i <- 0L
        run <- function(...) {i <<- i+1L; session$setInputs(...,run=i)}
        session$setInputs(map=1)
        expect_identical(reactive_files$error,"Select a packaged example.")
        session$setInputs(example_choice="oslo",confirm_example=1,k_subtypes="5",ratio=54/64)
        run(BSmethod="cIHC.itr",iteration=100)
        expect_match(reactive_files$error,"Preprocess")
        session$setInputs(map=2)
        for (value in list(101,1e12,"100",NULL,NA_real_,NaN,Inf,-Inf,10.5,9,c(10,100),matrix(100),c(10,1e6),"50",1e6)) {
            run(BSmethod="cIHC.itr", iteration=value)
            expect_length(received,0L)
            expect_match(reactive_files$error,"whole number")
            expect_null(current_prediction())
        }
        for (value in list(c("cIHC.itr","PAM50"),c("cIHC.itr","cIHC.itr"),c("PAM50","cIHC.itr"),NULL,NA_character_,"unknown",42L,factor("PAM50"))) {
            run(BSmethod=value,iteration=100)
            expect_length(received,0L)
            expect_match(reactive_files$error,"one supported")
            expect_null(current_prediction())
        }
        run(BSmethod="cIHC.itr",iteration=100L)
        expect_identical(received,list(100L))
        expect_identical(current_prediction()$method,"cIHC.itr")
        session$setInputs(BSmethod="PAM50")
        do.call(session$setInputs, setNames(list("External","nCounter"),
            c(calibration_ns()("calibration"),calibration_ns()("external"))))
        run(BSmethod="PAM50")
        expect_identical(current_prediction()$method,"PAM50")
        file <- tempfile()
        write_download(file)
        expect_equal(nrow(read.delim(file,comment.char="#")),102L)
        session$setInputs(load_example_step1=1,example_choice="tcga")
        expect_identical(current_prediction()$method,"PAM50")
        session$setInputs(confirm_example=2)
        expect_null(current_prediction())
        expect_error(write_download(file),"No current prediction",fixed=TRUE)
        expect_identical(safe_detail(simpleError("SYNTHETIC_PRIVATE_VALUE")),
                         "Review the settings or use the local BreastSubtypeR application.")
    })
})
