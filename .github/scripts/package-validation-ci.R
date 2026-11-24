# CI-only helper for .github/workflows/package-validation.yaml.
#
# Derived from .github/scripts/bioccheck-ci.R, which stays unchanged for the
# existing BiocCheck and release-validation workflows until the migration is
# complete. Differences from that helper: the expected R and Bioconductor
# versions come from the CI lane (.github/ci/lanes.json, resolved by
# resolve-ci-lane.sh and passed as environment variables) instead of being
# fixed in the code; the lane, the source tree and the R build are recorded
# and checked; and the summary describes the consolidated Package validation
# job, in which BiocCheck findings gate the job like R CMD check findings.
#
# Not part of the package: .Rbuildignore excludes .github, and the workflow
# verifies that the built tarball contains no .github files. Every mode runs
# in its own R process, because BiocCheck keeps its results in a
# reference-class singleton that separate runs must not share.
#
# The process exit status is the gate for every check mode:
#   0  the check completed and reported no ERROR and no WARNING
#   1  the check completed and reported at least one ERROR or WARNING
#   2  the check did not complete, its result was missing or malformed, the
#      CI lane was missing or invalid, or the identity of the checked
#      source/tarball/installation was not established
# NOTEs never change the exit status; they are recorded and summarised.
#
# Lane variables, required by the lane, environment and bioccheck modes and
# shown by the summary: CI_LANE, CI_TARGET_BRANCH, CI_R_VERSION ("X.Y" or
# "devel"), CI_BIOC_VERSION ("X.Y"), CI_CACHE_NAMESPACE, CI_LANES_SHA256.
#
# Usage: Rscript package-validation-ci.R <mode> <args...>
#   lane
#   environment SOURCE_DIR OUT_DIR
#   gitclone    SOURCE_DIR OUT_DIR
#   rcmdcheck   RCHECK_DIR OUT_DIR
#   bioccheck   TARBALL RCHECK_DIR BUILD_OUTPUT OUT_DIR
#   summary     EVIDENCE_DIR
#   selftest    OUT_DIR
#   fixture     KIND OUT_DIR    (synthetic BiocCheck results, for selftest)

EXIT_FINDINGS <- 1L
EXIT_INCOMPLETE <- 2L
# GitHub CI only: account requirements must be verified separately.
# https://bioconductor.org/packages/devel/bioc/vignettes/BiocCheck/inst/doc/BiocCheck.html
BIOCCHECK_CI_OPTIONS <- list(`no-check-bioc-help` = TRUE)

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

# CI lane -----------------------------------------------------------------

lane_from_env <- function() {
    get <- function(v) Sys.getenv(v, "")
    list(
        id = get("CI_LANE"),
        target_branch = get("CI_TARGET_BRANCH"),
        r_version = get("CI_R_VERSION"),
        bioc_version = get("CI_BIOC_VERSION"),
        cache_namespace = get("CI_CACHE_NAMESPACE"),
        lanes_sha256 = get("CI_LANES_SHA256")
    )
}

lane_problems <- function(lane) {
    p <- character()
    ok <- function(x, re) length(x) == 1L && !is.na(x) && grepl(re, x)
    if (!ok(lane$id, "^[a-z0-9][a-z0-9.-]*$"))
        p <- c(p, "CI_LANE is missing or malformed")
    if (!ok(lane$target_branch, "^[A-Za-z0-9._/-]+$"))
        p <- c(p, "CI_TARGET_BRANCH is missing or malformed")
    if (!ok(lane$r_version, "^(devel|[0-9]+\\.[0-9]+)$"))
        p <- c(p, "CI_R_VERSION must be X.Y or devel")
    if (!ok(lane$bioc_version, "^[0-9]+\\.[0-9]+$"))
        p <- c(p, "CI_BIOC_VERSION must be X.Y")
    if (!ok(lane$cache_namespace, "^[A-Za-z0-9][A-Za-z0-9._-]*$"))
        p <- c(p, "CI_CACHE_NAMESPACE is missing or malformed")
    else if (ok(lane$bioc_version, "^[0-9]+\\.[0-9]+$") &&
        !grepl(paste0("bioc-", lane$bioc_version), lane$cache_namespace, fixed = TRUE))
        p <- c(p, "CI_CACHE_NAMESPACE does not contain bioc-<CI_BIOC_VERSION>")
    if (!ok(lane$lanes_sha256, "^[0-9a-f]{64}$"))
        p <- c(p, "CI_LANES_SHA256 is missing or malformed")
    p
}

require_lane <- function() {
    lane <- lane_from_env()
    problems <- lane_problems(lane)
    if (length(problems))
        abort("invalid CI lane: ", paste(problems, collapse = "; "))
    lane
}

major_minor <- function(x) {
    x <- as.character(x)
    if (length(x) != 1L || is.na(x) || !grepl("^[0-9]+\\.[0-9]+", x))
        return(NA_character_)
    sub("^([0-9]+\\.[0-9]+).*$", "\\1", x)
}

# NULL if the running R matches the lane, otherwise a problem description.
# "X.Y" accepts any X.Y.z build that is not R-devel; "devel" requires R-devel.
r_version_problem <- function(expected, version, status) {
    version <- as.character(version)
    devel <- identical(status, "Under development (unstable)")
    if (identical(expected, "devel")) {
        if (!devel)
            return(sprintf("expected R-devel, found R %s%s", version,
                if (nzchar(status)) sprintf(" (%s)", status) else ""))
        return(NULL)
    }
    parts <- suppressWarnings(as.integer(strsplit(expected, ".", fixed = TRUE)[[1L]]))
    if (length(parts) != 2L || anyNA(parts))
        return(sprintf("invalid expected R version '%s'", expected))
    v <- package_version(version)
    lo <- package_version(sprintf("%d.%d.0", parts[1L], parts[2L]))
    hi <- package_version(sprintf("%d.%d.0", parts[1L], parts[2L] + 1L))
    if (devel || v < lo || v >= hi)
        return(sprintf("expected R %s.x, found R %s%s", expected, version,
            if (nzchar(status)) sprintf(" (%s)", status) else ""))
    NULL
}

# NULL if an installed BiocVersion belongs to the lane, otherwise a problem.
bioc_version_problem <- function(expected, installed) {
    found <- major_minor(installed)
    if (is.na(found) || !identical(found, expected))
        return(sprintf("expected BiocVersion %s.x, found %s", expected,
            paste(as.character(installed), collapse = " ")))
    NULL
}

run_lane <- function() {
    lane <- require_lane()
    message("CI lane ", lane$id, " for ", lane$target_branch, ": R ", lane$r_version,
        ", Bioconductor ", lane$bioc_version, ", cache ", lane$cache_namespace)
    0L
}

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
    lane <- require_lane()
    EXPECTED_BIOC <- lane$bioc_version
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

    git_out <- function(...) {
        out <- tryCatch(system2("git", shQuote(c("-C", source_dir, ...)),
            stdout = TRUE, stderr = FALSE), error = function(e) character())
        if (length(out) == 1L) out else NA_character_
    }
    git_head <- git_out("rev-parse", "HEAD")
    git_tree <- git_out("rev-parse", "HEAD^{tree}")

    problems <- character()
    problems <- c(problems, r_version_problem(
        lane$r_version, r_version, R.version$status
    ))
    problems <- c(problems, bioc_version_problem(EXPECTED_BIOC, bioc_version))
    if (!identical(bioc_reported, EXPECTED_BIOC))
        problems <- c(problems, paste(
            "BiocManager::version() reports", bioc_reported,
            "but the lane requires", EXPECTED_BIOC
        ))
    if (!identical(Sys.getenv("R_BIOC_VERSION"), EXPECTED_BIOC))
        problems <- c(problems, paste("R_BIOC_VERSION is not", EXPECTED_BIOC))
    if (is.na(git_head) || is.na(git_tree))
        problems <- c(problems, "the source revision could not be read with git")
    else if (nzchar(Sys.getenv("GITHUB_SHA")) &&
        !identical(git_head, Sys.getenv("GITHUB_SHA")))
        problems <- c(problems, paste(
            "checkout HEAD", git_head, "is not GITHUB_SHA", Sys.getenv("GITHUB_SHA")
        ))
    for (p in c("BiocCheck", "BiocVersion", "gert", "jsonlite", "rcmdcheck"))
        if (!p %in% pkgs$Package)
            problems <- c(problems, paste(p, "is not installed"))
    mixed <- bioc_pkgs$Package[
        is.na(bioc_pkgs$repo_version) | bioc_pkgs$repo_version != EXPECTED_BIOC
    ]
    if (length(mixed))
        problems <- c(problems, paste(
            "Bioconductor packages not installed from a", EXPECTED_BIOC, "repository:",
            paste(mixed, collapse = ", ")
        ))
    for (p in c("BiocCheck", "BiocVersion"))
        if (!p %in% bioc_pkgs$Package[bioc_pkgs$repo_version %in% EXPECTED_BIOC])
            problems <- c(problems, paste(
                p, "has no recorded Bioconductor", EXPECTED_BIOC, "repository"
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
        lane = lane,
        R = R.version.string,
        R_status = R.version$status,
        R_svn_rev = R.version[["svn rev"]],
        R_nickname = R.version$nickname,
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
        git_head = git_head,
        git_tree = git_tree,
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
    message("CI lane ", lane$id, " (", lane$target_branch, "): R ",
        lane$r_version, ", Bioconductor ", lane$bioc_version)
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
    write_result(list(stage = "BiocCheck", completed = FALSE,
        options = BIOCCHECK_CI_OPTIONS), out)
    lane <- require_lane()
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
    args <- c(list(
        package = tarball,
        debug = TRUE,
        install = paste0("check:", install_out),
        libloc = rcheck_dir,
        `build-output-file` = build_output
    ), BIOCCHECK_CI_OPTIONS)
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
    if (!identical(major_minor(meta$BiocVersion), lane$bioc_version))
        abort("BiocCheck reports Bioconductor ", paste(meta$BiocVersion, collapse = " "),
            " but the lane requires ", lane$bioc_version)
    bioccheck_log <- require_file(
        file.path(meta$BiocCheckDir, "00BiocCheck.log"), "00BiocCheck.log"
    )
    file.copy(bioccheck_log, file.path(out_dir, "00BiocCheck.log"), overwrite = TRUE)

    write_result(c(list(
        stage = "BiocCheck",
        completed = TRUE,
        mode = "existing package (new-package not set); debug = TRUE",
        options = BIOCCHECK_CI_OPTIONS,
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
    # The recorded counts decide FAIL even if the step outcome says success;
    # a clean result never turns a failed step into PASS.
    if (!is.null(result$counts) &&
        as.integer(result$counts$error) + as.integer(result$counts$warning) > 0L)
        return("FAIL")
    if (identical(outcome, "success")) "PASS" else "FAIL (step exited non-zero)"
}

counts_text <- function(result) {
    if (is.null(result$counts))
        return("")
    sprintf("%s ERROR, %s WARNING, %s NOTE",
        result$counts$error, result$counts$warning, result$counts$note)
}

# The CI option excludes all three account checks, even if BiocCheck fails.
# For older evidence without the option, read outcomes from the transcript.
account_checks <- function(result, console) {
    if (isTRUE(result$options[["no-check-bioc-help"]])) {
        skipped <- "NOT CHECKED - disabled in GitHub CI by no-check-bioc-help = TRUE; verify separately"
        return(list(devel = skipped, support = skipped, tag = skipped))
    }
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
    claim("BIOCCHECK", identical(bioc$options, BIOCCHECK_CI_OPTIONS),
        "BiocCheck step succeeded without the recorded CI account-check exclusion")
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
    lane <- lane_from_env()
    lane_ok <- !length(lane_problems(lane))
    lane_text <- if (lane_ok)
        sprintf("lane `%s` for `%s`: R %s, Bioconductor %s", lane$id,
            lane$target_branch, lane$r_version, lane$bioc_version)
    else
        "INVALID CI lane (see the environment stage)"

    md <- c(
        "## Package validation",
        "",
        sprintf("Tested on %s; running %s.", lane_text, env$R %||% "an unrecorded R"),
        "",
        "ERRORs and WARNINGs from BiocCheckGitClone, R CMD check or BiocCheck, an incomplete stage, missing evidence or a modified checkout fail this job; NOTEs are listed and do not. Whether this check is required for merging is set by the repository ruleset, not by this workflow.",
        "",
        if (isTRUE(bioc$completed))
            "Package BiocCheck validation was performed; its result is shown below."
        else
            "Package BiocCheck validation did not complete; no package validation success is claimed.",
        "",
        "| Stage | Status | Result |",
        "|---|---|---|",
        sprintf("| Environment setup | %s | %s |", status[["environment"]],
            one_line(c(if (lane_ok) paste0("lane ", lane$id, ";"), env$R,
                "; BiocVersion", env$BiocVersion, "; BiocCheck", env$BiocCheck))),
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
        sprintf("| Account: BreastSubtypeR watched-tag registration | %s | |", accounts$tag),
        sprintf("| Source checkout unmodified | %s | |", status[["pristine"]]),
        "",
        "### Identity",
        "",
        sprintf("- Event: `%s` on `%s`", gh$event %||% "", gh$ref %||% ""),
        sprintf("- CI lane: `%s`; target branch `%s`; R `%s`; Bioconductor `%s`; dependency cache `%s`; lanes.json SHA-256 `%s`",
            lane$id, lane$target_branch, lane$r_version, lane$bioc_version,
            lane$cache_namespace, lane$lanes_sha256),
        sprintf("- Tested commit (checkout HEAD): `%s`; tree `%s`", env$git_head %||% "unknown",
            env$git_tree %||% "unknown"),
        sprintf("- GITHUB_SHA: `%s`; PR head: `%s`", gh$sha %||% "", if (nzchar(gh$pr_head_sha %||% "")) gh$pr_head_sha else "n/a"),
        sprintf("- Tarball: `%s`", sha[2L]),
        sprintf("- Tarball SHA-256: `%s`", sha[1L]),
        sprintf("- Package version: `%s`", bioc$version %||% ""),
        sprintf("- Installation used by BiocCheck: `%s` (log `%s`)", bioc$install_dir %||% "", bioc$install_log %||% ""),
        sprintf("- %s (status `%s`, svn rev `%s`); BiocVersion %s; BiocManager::version() %s; BiocCheck %s; pandoc %s",
            env$R %||% "", env$R_status %||% "", env$R_svn_rev %||% "",
            env$BiocVersion %||% "", env$BiocManager_version %||% "",
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
        "GitHub CI uses `no-check-bioc-help = TRUE`. These account requirements must be verified separately; a passing package check does not verify them:",
        "",
        "- Bioconductor (bioc-devel) mailing-list subscription: NOT checked.",
        "- Support Site registration: NOT checked.",
        "- BreastSubtypeR watched-tag registration: NOT checked.",
        "",
        "The maintainer must separately verify subscription and Support Site registration with the DESCRIPTION maintainer email, and that BreastSubtypeR is a watched tag.",
        "Anything marked NOT VERIFIED or NOT EXECUTED above also remains unverified.",
        "",
        sprintf("Evidence: artifact `package-validation-evidence-%s-attempt-%s` (lane, commands, logs, results, tarball and hash).",
            lane$id, gh$run_attempt %||% Sys.getenv("GITHUB_RUN_ATTEMPT"))
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
        clean = , `clean-then-fail` = NULL,
        notes = {
            raise("Note", "synthetic note 1")
            raise("Note", "synthetic note 2")
        },
        error = {
            raise("Error", "synthetic error")
            raise("Note", "synthetic note")
        },
        warning = raise("Warning", "synthetic warning"),
        `support-warning` = raise("Warning", paste(
            "Unable to retrieve email info from the Support Site:",
            "Connection timed out after 15001 milliseconds")),
        malformed = return(list(error = list(), warning = list(), note = list())),
        fatal = abort("synthetic fatal execution error"),
        abort("unknown fixture kind: ", kind)
    )
    bc
}

run_fixture <- function(kind, out_dir) {
    out <- file.path(out_dir, "fixture-result.json")
    write_result(list(stage = paste("fixture", kind), completed = FALSE,
        options = BIOCCHECK_CI_OPTIONS), out)
    found <- bioccheck_findings(synthetic_result(kind))
    write_result(c(list(stage = paste("fixture", kind), completed = TRUE,
        options = BIOCCHECK_CI_OPTIONS), found), out)
    print_counts(paste("fixture", kind), found$counts)
    # A producer that fails after writing a clean result: the exit status,
    # not the file, must decide.
    if (identical(kind, "clean-then-fail"))
        return(3L)
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
        list("error", 1L, c(1, 0, 1)), list("warning", 1L, c(0, 1, 0)),
        list("support-warning", 1L, c(0, 1, 0)))) {
        name <- paste0("bioccheck-", fx[[1L]])
        st <- child(name, "fixture", fx[[1L]], "{dir}")
        expect(name, fx[[2L]], st, do.call(counts_are, c(list(result_of(name)), as.list(fx[[3L]]))))
    }
    st <- child("bioccheck-clean-then-fail", "fixture", "clean-then-fail", "{dir}")
    expect("bioccheck-clean-then-fail", 3L, st,
        counts_are(result_of("bioccheck-clean-then-fail"), 0, 0, 0))
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
        any(grepl("BiocCheck (existing package) | FAIL", readLines(md), fixed = TRUE)) &&
        sum(grepl("| NOT CHECKED - disabled", readLines(md), fixed = TRUE)) == 3L)
    # A clean result file next to a failed step is displayed as FAIL.
    sum_dir <- file.path(out_dir, "summary-clean-result-failed-step")
    dir.create(file.path(sum_dir, "bioccheck"), recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(out_dir, "bioccheck-clean-then-fail", "fixture-result.json"),
        file.path(sum_dir, "bioccheck", "bioccheck-result.json"))
    Sys.setenv(OUTCOME_BIOCCHECK = "failure")
    st <- child("summary-clean-result-failed-step", "summary", sum_dir)
    md <- file.path(sum_dir, "summary.md")
    expect("summary-clean-result-failed-step", 0L, st, file.exists(md) &&
        any(grepl("BiocCheck (existing package) | FAIL", readLines(md), fixed = TRUE)))
    # A WARNING in the result is displayed as FAIL even if the step outcome
    # claims success, and the inconsistency fails the summary step.
    sum_dir <- file.path(out_dir, "summary-masked-warning")
    dir.create(file.path(sum_dir, "bioccheck"), recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(out_dir, "bioccheck-warning", "fixture-result.json"),
        file.path(sum_dir, "bioccheck", "bioccheck-result.json"))
    Sys.setenv(OUTCOME_BIOCCHECK = "success")
    st <- child("summary-masked-warning", "summary", sum_dir)
    md <- file.path(sum_dir, "summary.md")
    expect("summary-masked-warning", 2L, st, file.exists(md) &&
        any(grepl("BiocCheck (existing package) | FAIL", readLines(md), fixed = TRUE)))
    sum_dir <- file.path(out_dir, "summary-missing-result")
    dir.create(sum_dir, showWarnings = FALSE)
    Sys.setenv(OUTCOME_BIOCCHECK = "success")
    st <- child("summary-missing-result", "summary", sum_dir)
    md <- file.path(sum_dir, "summary.md")
    expect("summary-missing-result", 2L, st, file.exists(md) &&
        any(grepl("Package BiocCheck validation did not complete", readLines(md), fixed = TRUE)))

    # A complete clean package result must disclose all account exclusions,
    # even if a stale transcript contains messages claiming account success.
    sum_dir <- file.path(out_dir, "summary-account-exclusions")
    bioc_dir <- file.path(sum_dir, "bioccheck")
    dir.create(bioc_dir, recursive = TRUE, showWarnings = FALSE)
    r <- result_of("bioccheck-clean")
    write_result(r, file.path(bioc_dir, "bioccheck-result.json"))
    for (f in c("00BiocCheck.log", "00install.out", "bioccheck-console.log"))
        writeLines(c("synthetic evidence", "Maintainer is registered at support site",
            "BreastSubtypeR is already in your 'Watched Tags'"), file.path(bioc_dir, f))
    st <- child("summary-account-exclusions", "summary", sum_dir)
    md <- readLines(file.path(sum_dir, "summary.md"))
    expect("summary-account-exclusions", 0L, st,
        any(grepl("Package BiocCheck validation was performed", md, fixed = TRUE)) &&
        any(grepl("BiocCheck (existing package) | PASS", md, fixed = TRUE)) &&
        sum(grepl("| NOT CHECKED - disabled", md, fixed = TRUE)) == 3L &&
        sum(grepl(": NOT checked.", md, fixed = TRUE)) == 3L &&
        any(grepl("must be verified separately", md, fixed = TRUE)) &&
        !any(grepl("^\\| Account:.*\\| PASS", md)))

    # Missing policy evidence cannot silently become a passing CI result.
    sum_dir <- file.path(out_dir, "summary-missing-ci-options")
    missing_dir <- file.path(sum_dir, "bioccheck")
    dir.create(missing_dir, recursive = TRUE, showWarnings = FALSE)
    for (f in c("00BiocCheck.log", "00install.out", "bioccheck-console.log"))
        file.copy(file.path(bioc_dir, f), file.path(missing_dir, f))
    r$options <- NULL
    write_result(r, file.path(missing_dir, "bioccheck-result.json"))
    expect("summary-missing-ci-options", 2L,
        child("summary-missing-ci-options", "summary", sum_dir))

    # An interrupted package check is never described as performed/passing.
    sum_dir <- file.path(out_dir, "summary-incomplete-account-exclusions")
    dir.create(file.path(sum_dir, "bioccheck"), recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(out_dir, "bioccheck-fatal", "fixture-result.json"),
        file.path(sum_dir, "bioccheck", "bioccheck-result.json"))
    Sys.setenv(OUTCOME_BIOCCHECK = "failure")
    st <- child("summary-incomplete-account-exclusions", "summary", sum_dir)
    md <- readLines(file.path(sum_dir, "summary.md"))
    expect("summary-incomplete-account-exclusions", 0L, st,
        any(grepl("Package BiocCheck validation did not complete", md, fixed = TRUE)) &&
        sum(grepl("| NOT CHECKED - disabled", md, fixed = TRUE)) == 3L &&
        !any(grepl("| PASS", md, fixed = TRUE)))
    Sys.unsetenv("OUTCOME_BIOCCHECK")

    # Lane controls: the configured lane passes; a missing or malformed lane
    # variable, or a cache namespace shared across Bioconductor versions, is
    # an incomplete result (exit 2), never a pass.
    lane_vars <- c("CI_LANE", "CI_TARGET_BRANCH", "CI_R_VERSION",
        "CI_BIOC_VERSION", "CI_CACHE_NAMESPACE", "CI_LANES_SHA256")
    saved <- Sys.getenv(lane_vars, unset = NA)
    restore <- function() for (v in lane_vars)
        if (is.na(saved[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, as.list(saved[v]))
    expect("lane-configured", 0L, child("lane-configured", "lane"))
    for (bad in list(
        list("lane-missing-bioc", "CI_BIOC_VERSION", ""),
        list("lane-malformed-bioc", "CI_BIOC_VERSION", "3.24.0"),
        list("lane-malformed-r", "CI_R_VERSION", "4.6.2"),
        list("lane-missing-id", "CI_LANE", ""),
        list("lane-foreign-cache", "CI_CACHE_NAMESPACE", "package-validation-bioc-9.99-r-4.6-v1"),
        list("lane-missing-config-hash", "CI_LANES_SHA256", ""))) {
        do.call(Sys.setenv, setNames(list(bad[[3L]]), bad[[2L]]))
        expect(bad[[1L]], EXIT_INCOMPLETE, child(bad[[1L]], "lane"))
        restore()
    }
    # Version matching used by the environment gate (in-process controls).
    version_case <- function(name, problem, expect_ok)
        expect(name, if (expect_ok) 0L else 1L, if (is.null(problem)) 0L else 1L)
    dev <- "Under development (unstable)"
    version_case("r-4.6-accepts-4.6.2", r_version_problem("4.6", "4.6.2", ""), TRUE)
    version_case("r-4.6-rejects-4.7.0", r_version_problem("4.6", "4.7.0", ""), FALSE)
    version_case("r-4.6-rejects-4.5.3", r_version_problem("4.6", "4.5.3", ""), FALSE)
    version_case("r-4.6-rejects-r-devel", r_version_problem("4.6", "4.6.0", dev), FALSE)
    version_case("r-devel-accepts-r-devel", r_version_problem("devel", "4.7.0", dev), TRUE)
    version_case("r-devel-rejects-release", r_version_problem("devel", "4.7.0", ""), FALSE)
    version_case("bioc-3.24-accepts-3.24.1", bioc_version_problem("3.24", "3.24.1"), TRUE)
    version_case("bioc-3.25-rejects-3.24.0", bioc_version_problem("3.25", "3.24.0"), FALSE)
    version_case("bioc-3.24-rejects-3.2", bioc_version_problem("3.24", "3.2"), FALSE)

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
    nargs <- c(lane = 0L, environment = 2L, gitclone = 2L, rcmdcheck = 2L,
        bioccheck = 4L, summary = 1L, selftest = 1L, fixture = 2L)
    mode <- args[1L]
    if (!mode %in% names(nargs) || length(args) - 1L != nargs[[mode]])
        abort("usage: package-validation-ci.R <mode> <args>; got: ", paste(args, collapse = " "))
    a <- args[-1L]
    switch(mode,
        lane = run_lane(),
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
