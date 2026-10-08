# CI-only helper for .github/workflows/BiocCheck.yaml.
#
# Not part of the package: .Rbuildignore excludes .github, and the workflow
# verifies that the built tarball contains no .github files. Every mode runs
# in its own R process, because BiocCheck keeps its results in a
# reference-class singleton that separate runs must not share.
#
# The process exit status is the gate for every check mode:
#   0  the check completed and reported no ERROR and no WARNING
#   1  the check completed and reported at least one ERROR or WARNING
#   2  the check did not complete, its result was missing or malformed, or
#      the identity of the checked tarball/installation was not established
# NOTEs never change the exit status; they are recorded and summarised.
#
# Usage: Rscript bioccheck-ci.R <mode> <args...>
#   environment SOURCE_DIR OUT_DIR
#   gitclone    SOURCE_DIR OUT_DIR
#   rcmdcheck   RCHECK_DIR OUT_DIR
#   bioccheck   TARBALL RCHECK_DIR BUILD_OUTPUT OUT_DIR
#   summary     EVIDENCE_DIR
#   selftest    OUT_DIR
#   fixture     KIND OUT_DIR    (synthetic BiocCheck results, for selftest)

EXIT_FINDINGS <- 1L
EXIT_INCOMPLETE <- 2L
EXPECTED_BIOC <- "3.24"

abort <- function(...) stop(paste0(...), call. = FALSE)

write_result <- function(x, path) {
    jsonlite::write_json(
        x, path, auto_unbox = TRUE, pretty = TRUE, null = "null",
        na = "null", digits = NA
    )
}

read_result <- function(path) {
    if (!file.exists(path))
        return(NULL)
    jsonlite::read_json(path, simplifyVector = FALSE)
}

log_command <- function(...) {
    evidence <- Sys.getenv("EVIDENCE_DIR")
    if (nzchar(evidence) && dir.exists(evidence))
        cat(paste0(...), "\n", sep = "",
            file = file.path(evidence, "commands.txt"), append = TRUE)
}

require_file <- function(path, what) {
    if (!file.exists(path))
        abort("missing mandatory evidence (", what, "): ", path)
    normalizePath(path, winslash = "/")
}

same_path <- function(a, b)
    identical(normalizePath(a, winslash = "/"), normalizePath(b, winslash = "/"))

gate_status <- function(counts) {
    if (counts[["error"]] + counts[["warning"]] > 0L) EXIT_FINDINGS else 0L
}

print_counts <- function(stage, counts) {
    message(sprintf(
        "%s: %d ERROR(s), %d WARNING(s), %d NOTE(s)", stage,
        counts[["error"]], counts[["warning"]], counts[["note"]]
    ))
}

# BiocCheck result object -------------------------------------------------

# Reads conditions through the documented fields (error, warning, note, log)
# and getNum() of the 'BiocCheck' reference class returned by BiocCheck() and
# BiocCheckGitClone(). Anything else is treated as a malformed result.
bioccheck_findings <- function(res) {
    if (!methods::is(res, "BiocCheck"))
        abort("expected a 'BiocCheck' result object, got: ",
            paste(class(res), collapse = "/"))
    conditions <- c("error", "warning", "note")
    counts <- res$getNum(conditions)
    if (!is.integer(counts) || length(counts) != 3L || anyNA(counts) ||
        !identical(names(counts), conditions))
        abort("malformed counts from BiocCheck getNum()")
    log <- res$log
    findings <- list()
    for (cond in conditions) {
        items <- res[[cond]]
        if (length(items) != counts[[cond]])
            abort("BiocCheck '", cond, "' list does not match getNum()")
        for (i in seq_along(items)) {
            el <- items[[i]]
            txt <- unlist(el, use.names = FALSE)
            if (!length(txt) || !is.character(txt))
                abort("malformed BiocCheck '", cond, "' entry")
            in_check <- vapply(log, function(entries) {
                any(vapply(entries, identical, logical(1L), el))
            }, logical(1L))
            findings[[length(findings) + 1L]] <- list(
                condition = toupper(cond),
                check = if (any(in_check))
                    names(log)[which(in_check)[1L]] else NA_character_,
                debug = names(items)[i] %||% NA_character_,
                message = txt[1L],
                details = I(txt[-1L])
            )
        }
    }
    list(counts = as.list(counts), findings = findings)
}

# environment -------------------------------------------------------------

run_environment <- function(source_dir, out_dir) {
    out <- file.path(out_dir, "environment.json")
    write_result(list(stage = "environment", completed = FALSE), out)
    pkg <- read.dcf(require_file(
        file.path(source_dir, "DESCRIPTION"), "DESCRIPTION"
    ))[, "Package"]

    r_version <- getRversion()
    bioc_version <- packageVersion("BiocVersion")
    bioc_reported <- as.character(BiocManager::version())

    fields <- c("RemoteType", "RemotePkgRef", "RemoteRepos")
    ip <- installed.packages(fields = fields)
    pkgs <- data.frame(
        Package = ip[, "Package"], Version = ip[, "Version"],
        LibPath = ip[, "LibPath"], RemotePkgRef = ip[, "RemotePkgRef"],
        RemoteRepos = ip[, "RemoteRepos"], row.names = NULL
    )
    utils::write.table(
        pkgs, file.path(out_dir, "installed-packages.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE, na = ""
    )
    writeLines(
        utils::capture.output(utils::sessionInfo()),
        file.path(out_dir, "sessionInfo.txt")
    )

    # Bioconductor provenance as recorded by pak at installation time.
    is_bioc <- grepl("bioconductor", pkgs$RemoteRepos, ignore.case = TRUE)
    repo_version <- ifelse(is_bioc, sub(
        ".*/([0-9]+\\.[0-9]+)/(bioc|data/annotation|data/experiment|workflows|books)/?$",
        "\\1", pkgs$RemoteRepos
    ), NA_character_)
    bioc_pkgs <- pkgs[is_bioc, c("Package", "Version", "RemoteRepos")]
    bioc_pkgs$repo_version <- repo_version[is_bioc]

    key <- c(
        "BiocCheck", "BiocVersion", "BiocManager", "BiocBaseUtils",
        "BiocFileCache", "biocViews", "httr2", "gert", "jsonlite", "rcmdcheck",
        "knitr", "rmarkdown", "BiocStyle", "testthat", "callr", "cli"
    )
    key_pkgs <- lapply(key, function(p) {
        row <- pkgs[pkgs$Package == p, , drop = FALSE]
        list(
            package = p,
            version = I(row$Version),
            library = I(row$LibPath),
            repository = I(ifelse(is.na(row$RemoteRepos), "", row$RemoteRepos))
        )
    })
    pandoc <- tryCatch(
        as.character(rmarkdown::pandoc_version()),
        error = function(e) paste("unavailable:", conditionMessage(e))
    )

    problems <- character()
    if (r_version < "4.6.0" || r_version >= "4.7.0")
        problems <- c(problems, paste("expected R 4.6.x, found R", r_version))
    if (bioc_version < EXPECTED_BIOC || bioc_version >= "3.25")
        problems <- c(problems, paste(
            "expected BiocVersion 3.24.x, found", bioc_version
        ))
    if (!identical(bioc_reported, EXPECTED_BIOC))
        problems <- c(problems, paste(
            "BiocManager::version() reports", bioc_reported
        ))
    if (!identical(Sys.getenv("R_BIOC_VERSION"), EXPECTED_BIOC))
        problems <- c(problems, "R_BIOC_VERSION is not 3.24")
    for (p in c("BiocCheck", "BiocVersion", "gert", "jsonlite", "rcmdcheck"))
        if (!p %in% pkgs$Package)
            problems <- c(problems, paste(p, "is not installed"))
    mixed <- bioc_pkgs$Package[
        is.na(bioc_pkgs$repo_version) | bioc_pkgs$repo_version != EXPECTED_BIOC
    ]
    if (length(mixed))
        problems <- c(problems, paste(
            "Bioconductor packages not installed from a 3.24 repository:",
            paste(mixed, collapse = ", ")
        ))
    for (p in c("BiocCheck", "BiocVersion"))
        if (!p %in% bioc_pkgs$Package[bioc_pkgs$repo_version %in% EXPECTED_BIOC])
            problems <- c(problems, paste(
                p, "has no recorded Bioconductor 3.24 repository"
            ))
    preinstalled <- find.package(pkg, lib.loc = .libPaths(), quiet = TRUE)
    if (length(preinstalled))
        problems <- c(problems, paste(
            pkg, "is already installed (would mask the checked tarball):",
            paste(preinstalled, collapse = ", ")
        ))

    rec <- list(
        stage = "environment",
        completed = TRUE,
        problems = I(problems),
        R = R.version.string,
        R_home = R.home(),
        platform = R.version$platform,
        BiocVersion = as.character(bioc_version),
        BiocManager_version = bioc_reported,
        BiocCheck = as.character(packageVersion("BiocCheck")),
        pandoc = pandoc,
        lib_paths = I(.libPaths()),
        R_LIBS_USER = Sys.getenv("R_LIBS_USER"),
        R_LIBS_SITE = Sys.getenv("R_LIBS_SITE"),
        R_BIOC_VERSION = Sys.getenv("R_BIOC_VERSION"),
        repos = as.list(getOption("repos")),
        bioc_repositories = as.list(BiocManager::repositories()),
        bioconductor_repo_versions = as.list(c(table(ifelse(
            is.na(bioc_pkgs$repo_version), "unknown", bioc_pkgs$repo_version
        )))),
        bioconductor_packages = nrow(bioc_pkgs),
        installed_packages = nrow(pkgs),
        key_packages = key_pkgs,
        git_head = system2("git", c("-C", shQuote(source_dir), "rev-parse", "HEAD"),
            stdout = TRUE),
        github = list(
            event = Sys.getenv("GITHUB_EVENT_NAME"),
            ref = Sys.getenv("GITHUB_REF"),
            sha = Sys.getenv("GITHUB_SHA"),
            pr_head_sha = Sys.getenv("PR_HEAD_SHA"),
            run_id = Sys.getenv("GITHUB_RUN_ID"),
            run_attempt = Sys.getenv("GITHUB_RUN_ATTEMPT")
        )
    )
    write_result(rec, out)
    message(R.version.string, "; BiocVersion ", bioc_version,
        "; BiocManager::version() ", bioc_reported,
        "; BiocCheck ", rec$BiocCheck)
    message("Library paths: ", paste(.libPaths(), collapse = " | "))
    message("Bioconductor packages by repository version: ",
        paste(names(rec$bioconductor_repo_versions),
            unlist(rec$bioconductor_repo_versions), sep = " x ",
            collapse = ", "))
    if (length(problems))
        abort("environment verification failed:\n  ",
            paste(problems, collapse = "\n  "))
    0L
}

# BiocCheckGitClone -------------------------------------------------------

run_gitclone <- function(source_dir, out_dir) {
    out <- file.path(out_dir, "gitclone-result.json")
    write_result(list(stage = "BiocCheckGitClone", completed = FALSE), out)
    source_dir <- normalizePath(source_dir, winslash = "/", mustWork = TRUE)
    if (!requireNamespace("gert", quietly = TRUE))
        abort("gert is required so that BiocCheckGitClone() reads git-tracked files")
    head <- gert::git_info(source_dir)$commit
    log_command("BiocCheck::BiocCheckGitClone(", deparse(source_dir), ")")
    res <- BiocCheck::BiocCheckGitClone(source_dir)
    found <- bioccheck_findings(res)
    write_result(c(list(
        stage = "BiocCheckGitClone",
        completed = TRUE,
        source_dir = source_dir,
        git_head = head,
        bioccheck_version = as.character(packageVersion("BiocCheck"))
    ), found), out)
    print_counts("BiocCheckGitClone", found$counts)
    gate_status(found$counts)
}

# R CMD check -------------------------------------------------------------

status_count <- function(status_line, word) {
    m <- regmatches(status_line,
        regexpr(paste0("[0-9]+ ", word, "s?"), status_line))
    if (length(m)) as.integer(sub(" .*", "", m)) else 0L
}

run_rcmdcheck <- function(rcheck_dir, out_dir) {
    out <- file.path(out_dir, "rcmdcheck-result.json")
    write_result(list(stage = "R CMD check", completed = FALSE), out)
    run_outcome <- Sys.getenv("RCMDCHECK_RUN_OUTCOME", "unknown")
    rcheck_dir <- normalizePath(rcheck_dir, winslash = "/", mustWork = TRUE)
    check_log <- require_file(file.path(rcheck_dir, "00check.log"), "00check.log")
    require_file(file.path(rcheck_dir, "00install.out"), "00install.out")

    lines <- readLines(check_log, warn = FALSE)
    status_line <- grep("^Status: ", lines, value = TRUE)
    if (length(status_line) != 1L)
        abort("00check.log has no single 'Status:' line; R CMD check did not complete")
    parsed <- rcmdcheck::parse_check(check_log)
    counts <- list(
        error = length(parsed$errors),
        warning = length(parsed$warnings),
        note = length(parsed$notes)
    )
    from_status <- list(
        error = status_count(status_line, "ERROR"),
        warning = status_count(status_line, "WARNING"),
        note = status_count(status_line, "NOTE")
    )
    if (!identical(counts, from_status))
        abort("parsed findings do not match '", status_line, "'")

    test_out <- list.files(file.path(rcheck_dir, "tests"),
        pattern = "\\.Rout(\\.fail)?$", full.names = TRUE)
    test_lines <- unlist(lapply(test_out, readLines, warn = FALSE))
    test_summary <- grep("\\[ FAIL [0-9]+ \\| WARN [0-9]+ \\| SKIP [0-9]+ \\| PASS [0-9]+ \\]",
        test_lines, value = TRUE)

    findings <- c(
        lapply(parsed$errors, function(x) list(condition = "ERROR", message = x)),
        lapply(parsed$warnings, function(x) list(condition = "WARNING", message = x)),
        lapply(parsed$notes, function(x) list(condition = "NOTE", message = x))
    )
    write_result(list(
        stage = "R CMD check",
        completed = TRUE,
        run_outcome = run_outcome,
        rcheck_dir = rcheck_dir,
        status_line = status_line,
        counts = counts,
        findings = findings,
        tests = I(basename(test_out)),
        test_summary = I(utils::tail(test_summary, 1L)),
        vignette_lines = I(grep("vignette", lines, value = TRUE, ignore.case = TRUE))
    ), out)
    print_counts("R CMD check", counts)
    status <- gate_status(counts)
    if (status == 0L && !identical(run_outcome, "success"))
        abort("R CMD check exited with outcome '", run_outcome,
            "' although 00check.log reports no ERROR or WARNING")
    status
}

# BiocCheck ---------------------------------------------------------------

run_bioccheck <- function(tarball, rcheck_dir, build_output, out_dir) {
    out <- file.path(out_dir, "bioccheck-result.json")
    write_result(list(stage = "BiocCheck", completed = FALSE), out)
    tarball <- require_file(tarball, "source tarball")
    rcheck_dir <- normalizePath(rcheck_dir, winslash = "/", mustWork = TRUE)
    build_output <- require_file(build_output, "R CMD build output")
    install_out <- require_file(file.path(rcheck_dir, "00install.out"), "00install.out")

    # Identity: one tarball, the installation R CMD check made from it, and
    # no other installation of the package anywhere on the library path.
    untar_dir <- tempfile("tarball-description-")
    files <- utils::untar(tarball, list = TRUE)
    desc_entry <- grep("^[^/]+/DESCRIPTION$", files, value = TRUE)
    if (length(desc_entry) != 1L)
        abort("tarball does not contain exactly one top-level DESCRIPTION")
    utils::untar(tarball, files = desc_entry, exdir = untar_dir)
    desc <- read.dcf(file.path(untar_dir, desc_entry))
    pkg <- desc[, "Package"]
    version <- desc[, "Version"]
    if (!identical(basename(tarball), sprintf("%s_%s.tar.gz", pkg, version)))
        abort("tarball name does not match its DESCRIPTION: ", basename(tarball))
    if (!identical(basename(rcheck_dir), paste0(pkg, ".Rcheck")))
        abort("check directory does not belong to ", pkg, ": ", rcheck_dir)
    installed_desc <- read.dcf(require_file(
        file.path(rcheck_dir, pkg, "DESCRIPTION"), "installed DESCRIPTION"
    ))
    if (!identical(unname(installed_desc[, "Version"]), unname(version)) ||
        !identical(unname(installed_desc[, "Packaged"]), unname(desc[, "Packaged"])))
        abort("the installation in ", rcheck_dir, " was not made from ", tarball)
    elsewhere <- find.package(pkg, lib.loc = .libPaths(), quiet = TRUE)
    if (length(elsewhere))
        abort(pkg, " is also installed in ", paste(elsewhere, collapse = ", "))
    if (grepl(":", install_out, fixed = TRUE))
        abort("BiocCheck splits 'install' on ':'; path is unusable: ", install_out)
    .libPaths(c(rcheck_dir, .libPaths()))
    if (!same_path(.libPaths()[1L], rcheck_dir) ||
        !nzchar(system.file(package = pkg, lib.loc = rcheck_dir)))
        abort("could not put the R CMD check installation first on .libPaths()")
    tarball_sha256 <- unname(tools::sha256sum(tarball))

    # BiocCheck 1.49.x: install = "check:<file>" skips installation and uses
    # .libPaths()[1]; libloc is the library for the man-page checks. Its own
    # copy of <file> into the .BiocCheck folder is lost to a later on.exit(),
    # so 00install.out is preserved here instead.
    file.copy(install_out, file.path(out_dir, "00install.out"), overwrite = TRUE)
    args <- list(
        package = tarball,
        debug = TRUE,
        install = paste0("check:", install_out),
        libloc = rcheck_dir,
        `build-output-file` = build_output
    )
    log_command(".libPaths(", deparse(.libPaths()), ")")
    log_command(paste(deparse(as.call(c(
        quote(BiocCheck::BiocCheck), args
    ))), collapse = " "))
    res <- do.call(BiocCheck::BiocCheck, args)
    found <- bioccheck_findings(res)

    meta <- res$metadata
    if (!identical(meta$Package, unname(pkg)) ||
        !identical(meta$PackageVersion, unname(version)) ||
        !isTRUE(meta$isTarBall) || !same_path(meta$installDir, rcheck_dir))
        abort("BiocCheck metadata does not match the checked tarball/installation")
    bioccheck_log <- require_file(
        file.path(meta$BiocCheckDir, "00BiocCheck.log"), "00BiocCheck.log"
    )
    file.copy(bioccheck_log, file.path(out_dir, "00BiocCheck.log"), overwrite = TRUE)

    write_result(c(list(
        stage = "BiocCheck",
        completed = TRUE,
        mode = "existing package (new-package not set); debug = TRUE",
        package = unname(pkg),
        version = unname(version),
        tarball = tarball,
        tarball_sha256 = tarball_sha256,
        rcheck_dir = rcheck_dir,
        install_dir = meta$installDir,
        install_log = install_out,
        build_output = build_output,
        bioccheck_dir = meta$BiocCheckDir,
        bioccheck_version = meta$BiocCheckVersion,
        bioc_version = meta$BiocVersion,
        bioc_devel_password_set = nzchar(Sys.getenv("BIOC_DEVEL_PASSWORD"))
    ), found), out)
    print_counts("BiocCheck", found$counts)
    gate_status(found$counts)
}

# Job summary -------------------------------------------------------------

strip_ansi <- function(x) gsub("\033\\[[0-9;]*[A-Za-z]", "", x)

one_line <- function(x) {
    x <- gsub("[[:space:]]+", " ", strip_ansi(paste(x, collapse = " ")))
    gsub("|", "\\|", trimws(x), fixed = TRUE)
}

annotate <- function(level, title, msg) {
    esc <- function(s) gsub("\n", "%0A", gsub("\r", "%0D", gsub("%", "%25", s)))
    cat(sprintf("::%s title=%s::%s\n", level, esc(title), esc(msg)))
}

stage_status <- function(outcome, result) {
    if (outcome %in% c("", "skipped"))
        return("NOT EXECUTED")
    if (identical(outcome, "cancelled"))
        return("CANCELLED")
    if (is.null(result))
        return(if (identical(outcome, "success")) "INCOMPLETE (no result)" else "FAIL (no result recorded)")
    if (!isTRUE(result$completed))
        return("INCOMPLETE (did not finish)")
    if (identical(outcome, "success")) "PASS" else "FAIL"
}

counts_text <- function(result) {
    if (is.null(result$counts))
        return("")
    sprintf("%s ERROR, %s WARNING, %s NOTE",
        result$counts$error, result$counts$warning, result$counts$note)
}

# Account checks run inside BiocCheck(); outcomes reported only as console
# messages (not counted) are read from the step transcript.
account_checks <- function(result, console) {
    text <- one_line(console)
    has <- function(p) grepl(p, text, fixed = TRUE)
    if (is.null(result) || !isTRUE(result$completed)) {
        na <- "NOT EXECUTED (BiocCheck did not complete)"
        return(list(devel = na, support = na, tag = na))
    }
    in_check <- function(check)
        Filter(function(f) identical(f$check, check), result$findings)
    devel <- in_check("Checking for bioc-devel mailing list subscription...")
    devel <- if (isTRUE(result$bioc_devel_password_set))
        "UNEXPECTED: admin credentials were present"
    else if (length(devel) && any(grepl("requires admin credentials", vapply(devel, `[[`, "", "message"))))
        "NOT VERIFIED - requires Bioconductor admin credentials, which this workflow does not use (BiocCheck records a NOTE)"
    else
        "NOT EXECUTED (no result found)"
    reg <- in_check("Checking for support site registration...")
    reg_msgs <- vapply(reg, function(f) paste(f$condition, f$message), "")
    support <- if (any(grepl("Unable to retrieve email info", reg_msgs)))
        "NOT VERIFIED - Support Site lookup failed (network/service, not a package defect); BiocCheck WARNING fails this job"
    else if (any(grepl("Register your email", reg_msgs)))
        "FAIL - maintainer email is not registered on the Support Site"
    else if (has("Maintainer is registered at support site"))
        "PASS - maintainer is registered on the Support Site"
    else
        "NOT EXECUTED (no result found)"
    tag <- if (any(grepl("Watched Tags", reg_msgs)))
        "FAIL - package tag is not in the maintainer's Watched Tags"
    else if (has("is already in your 'Watched Tags'"))
        "PASS - package tag is in the maintainer's Watched Tags"
    else if (has("Unable to retrieve 'Watched Tags' profile"))
        "NOT VERIFIED - Support Site lookup failed; BiocCheck reports this only as a message, so it does not change the counts"
    else
        "NOT EXECUTED (registration lookup did not succeed)"
    list(devel = devel, support = support, tag = tag)
}

findings_table <- function(findings, with_check = TRUE) {
    if (!length(findings))
        return("None.\n")
    hdr <- if (with_check) "| Level | Check | Finding | Function |\n|---|---|---|---|"
        else "| Level | Finding |\n|---|---|"
    rows <- vapply(findings, function(f) {
        msg <- one_line(c(f$message, unlist(f$details)))
        if (with_check)
            sprintf("| %s | %s | %s | `%s` |", f$condition,
                one_line(f$check %||% ""), msg, one_line(f$debug %||% ""))
        else
            sprintf("| %s | %s |", f$condition, msg)
    }, "")
    paste(c(hdr, rows, ""), collapse = "\n")
}

run_summary <- function(evidence) {
    outcome <- function(name) Sys.getenv(paste0("OUTCOME_", name), "")
    env <- read_result(file.path(evidence, "environment", "environment.json"))
    selftest <- read_result(file.path(evidence, "selftest", "selftest.json"))
    gitclone <- read_result(file.path(evidence, "gitclone", "gitclone-result.json"))
    rcheck <- read_result(file.path(evidence, "check", "rcmdcheck-result.json"))
    bioc <- read_result(file.path(evidence, "bioccheck", "bioccheck-result.json"))
    console_file <- file.path(evidence, "bioccheck", "bioccheck-console.log")
    console <- if (file.exists(console_file)) readLines(console_file, warn = FALSE) else ""
    sha_file <- file.path(evidence, "build", "tarball.sha256")
    build_log <- file.path(evidence, "build", "build-output.txt")

    problems <- character()
    claim <- function(name, ok, what)
        if (identical(outcome(name), "success") && !ok)
            problems <<- c(problems, what)
    claim("ENVIRONMENT", isTRUE(env$completed) && !length(env$problems),
        "environment step succeeded without a complete environment record")
    claim("SELFTEST", isTRUE(selftest$completed) && isTRUE(selftest$all_passed),
        "self-test step succeeded without a complete passing record")
    for (st in list(list("GITCLONE", gitclone, "BiocCheckGitClone"),
        list("RCMDCHECK", rcheck, "R CMD check"),
        list("BIOCCHECK", bioc, "BiocCheck"))) {
        r <- st[[2L]]
        claim(st[[1L]], isTRUE(r$completed) &&
            identical(as.integer(r$counts$error) + as.integer(r$counts$warning), 0L),
            paste(st[[3L]], "step succeeded without a complete, clean result"))
    }
    claim("BUILD", file.exists(sha_file) && file.exists(build_log),
        "build step succeeded without tarball hash and build log")
    if (identical(outcome("RCMDCHECK"), "success"))
        for (f in c("00check.log", "00install.out"))
            claim("RCMDCHECK", file.exists(file.path(rcheck$rcheck_dir %||% "", f)),
                paste("R CMD check succeeded but", f, "is missing"))
    if (identical(outcome("BIOCCHECK"), "success"))
        for (f in c("00BiocCheck.log", "00install.out", "bioccheck-console.log"))
            claim("BIOCCHECK", file.exists(file.path(evidence, "bioccheck", f)),
                paste("BiocCheck succeeded but", f, "is missing"))

    status <- c(
        environment = stage_status(outcome("ENVIRONMENT"), env),
        selftest = stage_status(outcome("SELFTEST"), selftest),
        gitclone = stage_status(outcome("GITCLONE"), gitclone),
        build = if (outcome("BUILD") %in% c("", "skipped")) "NOT EXECUTED"
            else if (identical(outcome("BUILD"), "success") && file.exists(sha_file)) "PASS"
            else "FAIL",
        rcmdcheck = stage_status(outcome("RCMDCHECK"), rcheck),
        bioccheck = stage_status(outcome("BIOCCHECK"), bioc),
        pristine = if (identical(outcome("PRISTINE"), "success")) "PASS"
            else if (outcome("PRISTINE") %in% c("", "skipped")) "NOT EXECUTED" else "FAIL"
    )
    if (identical(outcome("ENVIRONMENT"), "failure") && length(env$problems))
        status[["environment"]] <- paste("FAIL -",
            one_line(paste(unlist(env$problems), collapse = "; ")))
    accounts <- account_checks(bioc, console)

    sha <- if (file.exists(sha_file))
        sub("^\\*", "", strsplit(readLines(sha_file, 1L), "[[:space:]]+")[[1L]])
    else c("", "")
    vignettes_built <- file.exists(build_log) &&
        any(grepl("^\\* creating vignettes \\.\\.\\.", readLines(build_log, warn = FALSE)))
    gh <- env$github %||% list()

    md <- c(
        "## BiocCheck (R 4.6, Bioc 3.24)",
        "",
        "Advisory job: not a required status check. ERRORs and WARNINGs still fail it; NOTEs are listed and do not.",
        "",
        "| Stage | Status | Result |",
        "|---|---|---|",
        sprintf("| Environment setup | %s | %s |", status[["environment"]],
            one_line(c(env$R, "; BiocVersion", env$BiocVersion, "; BiocCheck", env$BiocCheck))),
        sprintf("| Failure-gate self-test (synthetic) | %s | %s |", status[["selftest"]],
            if (is.null(selftest)) "" else sprintf("%s/%s controls behaved as expected",
                selftest$passed, selftest$total)),
        sprintf("| BiocCheckGitClone (pristine checkout) | %s | %s |", status[["gitclone"]], counts_text(gitclone)),
        sprintf("| R CMD build | %s | %s |", status[["build"]],
            if (nzchar(sha[2L])) sprintf("`%s`, vignettes built: %s", sha[2L], if (vignettes_built) "yes" else "NO") else ""),
        sprintf("| R CMD check --no-manual | %s | %s |", status[["rcmdcheck"]],
            one_line(c(counts_text(rcheck), unlist(rcheck$test_summary)))),
        sprintf("| BiocCheck (existing package) | %s | %s |", status[["bioccheck"]], counts_text(bioc)),
        sprintf("| Account: Bioc-devel subscription | %s | |", accounts$devel),
        sprintf("| Account: Support Site registration | %s | |", accounts$support),
        sprintf("| Account: Support Site watched tag | %s | |", accounts$tag),
        sprintf("| Source checkout unmodified | %s | |", status[["pristine"]]),
        "",
        "### Identity",
        "",
        sprintf("- Event: `%s` on `%s`", gh$event %||% "", gh$ref %||% ""),
        sprintf("- Tested commit (checkout HEAD): `%s`", env$git_head %||% "unknown"),
        sprintf("- GITHUB_SHA: `%s`; PR head: `%s`", gh$sha %||% "", if (nzchar(gh$pr_head_sha %||% "")) gh$pr_head_sha else "n/a"),
        sprintf("- Tarball: `%s`", sha[2L]),
        sprintf("- Tarball SHA-256: `%s`", sha[1L]),
        sprintf("- Package version: `%s`", bioc$version %||% ""),
        sprintf("- Installation used by BiocCheck: `%s` (log `%s`)", bioc$install_dir %||% "", bioc$install_log %||% ""),
        sprintf("- %s; BiocVersion %s; BiocManager::version() %s; BiocCheck %s; pandoc %s",
            env$R %||% "", env$BiocVersion %||% "", env$BiocManager_version %||% "",
            env$BiocCheck %||% "", env$pandoc %||% ""),
        sprintf("- Library paths: `%s`", paste(unlist(env$lib_paths), collapse = "`, `")),
        sprintf("- Bioconductor packages by repository version: %s",
            paste(names(env$bioconductor_repo_versions), unlist(env$bioconductor_repo_versions),
                sep = " x ", collapse = ", ")),
        "",
        "### BiocCheck findings",
        "",
        findings_table(bioc$findings),
        "### BiocCheckGitClone findings",
        "",
        findings_table(gitclone$findings),
        "### R CMD check findings",
        "",
        findings_table(rcheck$findings, with_check = FALSE),
        "### Not verified by this job",
        "",
        "- Bioc-devel mailing-list subscription: needs Bioconductor admin credentials. A passing job does not verify it.",
        "- Anything marked NOT VERIFIED or NOT EXECUTED above.",
        "",
        sprintf("Evidence: artifact `BiocCheck-evidence-attempt-%s` (commands, logs, results, tarball and hash).",
            gh$run_attempt %||% Sys.getenv("GITHUB_RUN_ATTEMPT"))
    )
    if (length(problems))
        md <- c(md, "", "### Evidence problems", "", paste("-", problems))

    writeLines(md, file.path(evidence, "summary.md"))
    step_summary <- Sys.getenv("GITHUB_STEP_SUMMARY")
    if (nzchar(step_summary))
        cat(md, sep = "\n", file = step_summary, append = TRUE)
    cat(md, sep = "\n")

    for (st in list(list("BiocCheck", bioc), list("BiocCheckGitClone", gitclone),
        list("R CMD check", rcheck)))
        for (f in st[[2L]]$findings)
            annotate(if (identical(f$condition, "NOTE")) "notice" else "error",
                paste(st[[1L]], f$condition), one_line(c(f$message, unlist(f$details))))
    for (a in unlist(accounts))
        if (startsWith(a, "NOT VERIFIED"))
            annotate("warning", "BiocCheck account check not verified", a)
    if (length(problems)) {
        for (p in problems) annotate("error", "Missing or inconsistent evidence", p)
        return(EXIT_INCOMPLETE)
    }
    0L
}

# Synthetic controls ------------------------------------------------------

# Builds a BiocCheck result with the package's own condition handlers, so the
# synthetic results go through exactly the same reading and gate as real ones.
synthetic_result <- function(kind) {
    ns <- asNamespace("BiocCheck")
    bc <- BiocCheck::.BiocCheck
    bc$zero()
    bc$log <- list()
    raise <- function(condition, msg) get(paste0("handle", condition), envir = ns)(msg)
    get("handleCheck", envir = ns)("Synthetic control check...")
    switch(kind,
        clean = NULL,
        notes = {
            raise("Note", "synthetic note 1")
            raise("Note", "synthetic note 2")
        },
        error = {
            raise("Error", "synthetic error")
            raise("Note", "synthetic note")
        },
        warning = raise("Warning", "synthetic warning"),
        malformed = return(list(error = list(), warning = list(), note = list())),
        fatal = abort("synthetic fatal execution error"),
        abort("unknown fixture kind: ", kind)
    )
    bc
}

run_fixture <- function(kind, out_dir) {
    out <- file.path(out_dir, "fixture-result.json")
    write_result(list(stage = paste("fixture", kind), completed = FALSE), out)
    found <- bioccheck_findings(synthetic_result(kind))
    write_result(c(list(stage = paste("fixture", kind), completed = TRUE), found), out)
    print_counts(paste("fixture", kind), found$counts)
    gate_status(found$counts)
}

run_selftest <- function(out_dir) {
    self <- normalizePath(sub("^--file=", "",
        grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]))
    rscript <- file.path(R.home("bin"), "Rscript")
    out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
    out <- file.path(out_dir, "selftest.json")
    write_result(list(stage = "selftest", completed = FALSE), out)
    # Children inherit the environment; keep their evidence and summaries out
    # of the real ones.
    Sys.setenv(GITHUB_STEP_SUMMARY = "")

    # Run a mode in a child process as the workflow does: on Linux through
    # 'bash -eo pipefail' with '| tee', otherwise directly. Child output goes
    # only to its log, so synthetic findings never become job annotations.
    child <- function(name, ...) {
        dir <- file.path(out_dir, name)
        dir.create(dir, showWarnings = FALSE, recursive = TRUE)
        Sys.setenv(EVIDENCE_DIR = dir)
        args <- c(...)
        args[args == "{dir}"] <- dir
        log <- file.path(dir, "console.log")
        cmd <- paste(shQuote(c(rscript, self, args)), collapse = " ")
        if (identical(.Platform$OS.type, "unix")) {
            system2("bash", c("--noprofile", "--norc", "-eo", "pipefail", "-c",
                shQuote(paste(cmd, "2>&1 | tee", shQuote(log), "> /dev/null"))))
        } else {
            system2(rscript, shQuote(c(self, args)), stdout = log, stderr = log)
        }
    }
    check_header <- c(
        "* using log directory '/tmp/fake.Rcheck'",
        "* using R version 4.6.0 (2026-04-24)",
        "* using platform: x86_64-pc-linux-gnu",
        "* using option '--no-manual'",
        "* checking for file 'fake/DESCRIPTION' ... OK",
        "* this is package 'fake' version '0.0.1'"
    )
    fake_rcheck <- function(name, status_line, body = character(), install = TRUE) {
        dir <- file.path(out_dir, name, "fake.Rcheck")
        dir.create(dir, showWarnings = FALSE, recursive = TRUE)
        if (!is.null(status_line))
            writeLines(c(check_header, body, "* DONE", "", status_line),
                file.path(dir, "00check.log"))
        if (install)
            writeLines("* installing *source* package 'fake' ...", file.path(dir, "00install.out"))
        dir
    }
    result_of <- function(name, file = "fixture-result.json")
        read_result(file.path(out_dir, name, file))

    controls <- list()
    expect <- function(name, expected, observed, ok_result = TRUE) {
        controls[[length(controls) + 1L]] <<- list(
            control = name, expected_exit = expected, observed_exit = observed,
            result_ok = ok_result,
            passed = identical(as.integer(observed), as.integer(expected)) && isTRUE(ok_result)
        )
    }
    counts_are <- function(r, e, w, n)
        isTRUE(r$completed) && identical(c(r$counts$error, r$counts$warning, r$counts$note),
            as.integer(c(e, w, n)))

    for (fx in list(
        list("clean", 0L, c(0, 0, 0)), list("notes", 0L, c(0, 0, 2)),
        list("error", 1L, c(1, 0, 1)), list("warning", 1L, c(0, 1, 0)))) {
        name <- paste0("bioccheck-", fx[[1L]])
        st <- child(name, "fixture", fx[[1L]], "{dir}")
        expect(name, fx[[2L]], st, do.call(counts_are, c(list(result_of(name)), as.list(fx[[3L]]))))
    }
    for (fx in c("malformed", "fatal")) {
        name <- paste0("bioccheck-", fx)
        st <- child(name, "fixture", fx, "{dir}")
        r <- result_of(name)
        expect(name, EXIT_INCOMPLETE, st, !is.null(r) && identical(r$completed, FALSE))
    }

    Sys.setenv(RCMDCHECK_RUN_OUTCOME = "success")
    for (fx in list(
        list("rcmdcheck-ok", "Status: OK", character(), 0L),
        list("rcmdcheck-note", "Status: 1 NOTE",
            c("* checking R code for possible problems ... NOTE", "  synthetic note"), 0L),
        list("rcmdcheck-warning", "Status: 1 WARNING",
            c("* checking Rd files ... WARNING", "  synthetic warning"), 1L),
        list("rcmdcheck-no-status", NULL, character(), 2L))) {
        if (is.null(fx[[2L]])) {
            dir <- fake_rcheck(fx[[1L]], NULL)
            writeLines(c(check_header, "* checking examples ..."),
                file.path(dir, "00check.log"))
        } else {
            dir <- fake_rcheck(fx[[1L]], fx[[2L]], fx[[3L]])
        }
        st <- child(fx[[1L]], "rcmdcheck", dir, "{dir}")
        expect(fx[[1L]], fx[[4L]], st)
    }
    dir <- fake_rcheck("rcmdcheck-missing-log", NULL)
    expect("rcmdcheck-missing-log", 2L, child("rcmdcheck-missing-log", "rcmdcheck", dir, "{dir}"))
    Sys.setenv(RCMDCHECK_RUN_OUTCOME = "failure")
    dir <- fake_rcheck("rcmdcheck-error", "Status: 1 ERROR",
        c("* checking examples ... ERROR", "  synthetic error"))
    expect("rcmdcheck-error", 1L, child("rcmdcheck-error", "rcmdcheck", dir, "{dir}"))
    dir <- fake_rcheck("rcmdcheck-exit-mismatch", "Status: OK")
    expect("rcmdcheck-exit-mismatch", 2L, child("rcmdcheck-exit-mismatch", "rcmdcheck", dir, "{dir}"))
    Sys.unsetenv("RCMDCHECK_RUN_OUTCOME")

    # Summaries: a failing check still yields a summary and a clean exit of
    # the summary step; a step that claims success without its result does not.
    sum_dir <- file.path(out_dir, "summary-after-failure")
    dir.create(file.path(sum_dir, "bioccheck"), recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(out_dir, "bioccheck-error", "fixture-result.json"),
        file.path(sum_dir, "bioccheck", "bioccheck-result.json"))
    Sys.setenv(OUTCOME_BIOCCHECK = "failure")
    st <- child("summary-after-failure", "summary", sum_dir)
    md <- file.path(sum_dir, "summary.md")
    expect("summary-after-failure", 0L, st, file.exists(md) &&
        any(grepl("BiocCheck (existing package) | FAIL", readLines(md), fixed = TRUE)))
    sum_dir <- file.path(out_dir, "summary-missing-result")
    dir.create(sum_dir, showWarnings = FALSE)
    Sys.setenv(OUTCOME_BIOCCHECK = "success")
    st <- child("summary-missing-result", "summary", sum_dir)
    expect("summary-missing-result", 2L, st, file.exists(file.path(sum_dir, "summary.md")))
    Sys.unsetenv("OUTCOME_BIOCCHECK")

    passed <- vapply(controls, `[[`, logical(1L), "passed")
    for (ct in controls)
        message(sprintf("%-28s expected exit %d, observed %d, result ok: %s -> %s",
            ct$control, ct$expected_exit, ct$observed_exit, ct$result_ok,
            if (ct$passed) "PASS" else "FAIL"))
    write_result(list(
        stage = "selftest", completed = TRUE, platform = R.version$platform,
        total = length(controls), passed = sum(passed), all_passed = all(passed),
        controls = controls
    ), out)
    if (!all(passed))
        abort("failure-gate self-test failed")
    0L
}

# Entry point -------------------------------------------------------------

main <- function(args) {
    if (!length(args))
        abort("no mode given")
    nargs <- c(environment = 2L, gitclone = 2L, rcmdcheck = 2L, bioccheck = 4L,
        summary = 1L, selftest = 1L, fixture = 2L)
    mode <- args[1L]
    if (!mode %in% names(nargs) || length(args) - 1L != nargs[[mode]])
        abort("usage: bioccheck-ci.R <mode> <args>; got: ", paste(args, collapse = " "))
    a <- args[-1L]
    switch(mode,
        environment = run_environment(a[1L], a[2L]),
        gitclone = run_gitclone(a[1L], a[2L]),
        rcmdcheck = run_rcmdcheck(a[1L], a[2L]),
        bioccheck = run_bioccheck(a[1L], a[2L], a[3L], a[4L]),
        summary = run_summary(a[1L]),
        selftest = run_selftest(a[1L]),
        fixture = run_fixture(a[1L], a[2L])
    )
}

status <- tryCatch(
    main(commandArgs(trailingOnly = TRUE)),
    error = function(e) {
        message("ERROR (incomplete): ", conditionMessage(e))
        EXIT_INCOMPLETE
    }
)
quit(save = "no", status = status)
