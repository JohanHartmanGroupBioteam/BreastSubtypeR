# State transitions and download callbacks of the installed Shiny app.
shiny_app_env <- function() {
    skip_if_not_installed("shiny")
    skip_if_not_installed("bslib")
    env <- new.env(parent = asNamespace("shiny"))
    server_path <- system.file("ShinyBreastSubtypeR/server.R",
                               package = "BreastSubtypeR")
    eval(parse(file = server_path, encoding = "UTF-8", keep.source = FALSE),
         envir = env)
    env
}
shiny_defaults <- function(session, method = "PAM50") {
    session$setInputs(BSmethod = method, k_subtypes = "5",
                     is_raw_counts = "norm", hasClinical = FALSE,
                     export_kind = "calls", calibration = "Internal",
                     internal = "medianCtr", external = "Given.mdns",
                     s = "ER.v2", ratio = 54 / 64, iteration = 10)
    session$setInputs(example_choice = "oslo", confirm_example = 1)
    session$setInputs(map = 1)
}
shiny_upload <- function(name) {
    path <- system.file("RshinyTest", name, package = "BreastSubtypeR")
    data.frame(name = name, size = file.info(path)$size,
               type = "text/plain", datapath = path)
}
