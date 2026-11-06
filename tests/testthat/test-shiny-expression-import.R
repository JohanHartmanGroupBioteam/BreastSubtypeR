read_expression_fixture <- function(records, patients, probes, sep = "\t", raw = FALSE) {
    path <- withr::local_tempfile()
    writeLines(records, path, useBytes = TRUE)
    BreastSubtypeR:::.read_shiny_expression(path, sep, patients, probes, raw)
}

test_that("gene-ID layouts preserve exact measurements, identifiers and orientation", {
    patients <- c("009 B", "S-01")
    probes <- c("p-2", "p.1")
    expected <- matrix(c(1.2345678901234567, -20, .125, 0), 2L,
                       dimnames = list(probes, patients))
    for (sep in c("\t", ",")) {
        for (corner in c("PROBE", "", "feature", NA_character_)) {
            header <- if (is.na(corner)) patients else c(corner, patients)
            records <- c(paste(header, collapse = sep),
                paste("p-2", "1.2345678901234567", " +.125 ", sep = sep),
                paste("p.1", "-2E+1", "0", sep = sep))
            expect_identical(read_expression_fixture(records, rev(patients), rev(probes), sep),
                             expected)
        }
    }
    expect_identical(read_expression_fixture(c("009 B\tPrObE\tS-01",
        "1.2345678901234567\tp-2\t.125", "-20\tp.1\t0"), patients, probes), expected)
    expect_identical(read_expression_fixture(c("S1\tS2", "2\t101", "4\t102", "1\t103"),
        c("S1", "S2"), as.character(1:3)),
        matrix(c(2, 4, 1, 101, 102, 103), 3L, dimnames = list(as.character(1:3), c("S1", "S2"))))
    expect_error(read_expression_fixture(c("feature\tS1", "p1\t23", "p2\t12"),
        "S1", c("1", "p1", "p2")), "Ambiguous.*explicit probe")
    # An unmatched sample can share numeric values with annotation IDs. Its
    # repeated measurements rule it out as a possible unique probe column.
    expect_identical(read_expression_fixture(c("unused\tS1", "1\t10", "1\t20"), "S1", c("1", "2")),
        matrix(c(1, 1, 10, 20), 2L, dimnames = list(c("1", "2"), c("unused", "S1"))))
    expect_identical(read_expression_fixture(c("feature\tS1", "p1\t23", "p2\t12"),
        "S1", "p2"), matrix(c(23, 12), 2L, dimnames = list(c("p1", "p2"), "S1")))
})

test_that("singletons, quoting and numeric-looking IDs keep their dimensions and spelling", {
    expect_identical(read_expression_fixture(c('"probe","007","S,""A"""',
        '"0004","-0.03125","0x10"'), c("007", 'S,"A"'), "0004", ","),
        matrix(c(-.03125, 16), 1L, dimnames = list("0004", c("007", 'S,"A"'))))
    expect_identical(read_expression_fixture(c("007", "2147483648"), "007", "1", raw = TRUE),
        matrix(2147483648, 1L, dimnames = list("1", "007")))
    expect_identical(read_expression_fixture(c("probe\tA", "001\t5", "01\t6"), "A",
        c("001", "01"), raw = TRUE), matrix(c(5, 6), 2L, dimnames = list(c("001", "01"), "A")))
    expect_identical(read_expression_fixture(c("probe,S", '"p', '1",1.5'), "S", "p\n1", ","),
        matrix(1.5, 1L, dimnames = list("p\n1", "S")))
    expect_identical(read_expression_fixture(c('probe\tS"1', 'p"1\t12'), 'S"1', 'p"1'),
        matrix(12, 1L, dimnames = list('p"1', 'S"1')))
    p <- withr::local_tempfile()
    writeBin(charToRaw(enc2utf8("\ufeffprobe\tS\u03b2\np\u00e9\t2\n")), p)
    expect_identical(BreastSubtypeR:::.read_shiny_expression(p, "\t", "S\u03b2", "p\u00e9", FALSE),
                     matrix(2, 1L, dimnames = list("p\u00e9", "S\u03b2")))
})

test_that("missingness is mode-specific and explicit infinities retain their signs", {
    rows <- c("probe\tS1\tS2", "g1\tNA\t", "g2\tNaN\t", "g3\tInf\t-Infinity", "g4\t\t+inf")
    patients <- c("S1", "S2")
    probes <- paste0("g", 1:4)
    expected <- matrix(c(NA_real_, NaN, Inf, NA_real_, NA_real_, NA_real_, -Inf, Inf),
                       4L, dimnames = list(probes, patients))
    expect_identical(read_expression_fixture(rows, patients, probes), expected)
    for (missing in c("NA", "NaN", "")) {
        expect_error(read_expression_fixture(c("probe\tS1", paste0("g1\t", missing)),
            "S1", "g1", raw = TRUE), "line 2, probe g1, sample S1.*Raw counts")
    }
    expect_error(read_expression_fixture(c("probe,S1,S2", "g1,12,"), patients, "g1", ",", TRUE),
                 "line 2, probe g1, sample S2 = <blank>.*Raw counts")
    for (raw in c(FALSE, TRUE)) {
        expect_identical(read_expression_fixture(c("probe\tS", "g1\tInf\f", "g2\t-Inf"),
            "S", c("g1", "g2"), raw = raw),
            matrix(c(Inf, -Inf), 2L, dimnames = list(c("g1", "g2"), "S")))
        expect_error(read_expression_fixture(c("probe\tS", "g1\t-1e999"),
            "S", "g1", raw = raw), 'line 2, probe g1, sample S = "-1e999" \\(numeric overflow\\)')
    }
})

test_that("malformed measurements identify the original cell and stop the import", {
    for (token in c("n/a", "#N/A", "-", "1,25", "TRUE", "1+2i", "3mg")) {
        for (sep in c("\t", ",")) {
            field <- if (sep == "," && token == "1,25") '"1,25"' else token
            records <- c(paste("probe", "S-01", "009 B", sep = sep),
                paste("p1", "10.5", "25", sep = sep), paste("p2", "8.75", ".5", sep = sep),
                paste("p3", field, "-2e1", sep = sep))
            expect_error(read_expression_fixture(records, c("S-01", "009 B"), paste0("p", 1:3), sep),
                         paste0('line 4, probe p3, sample S-01 = "', token, '"'), fixed = TRUE)
        }
    }
    msg <- tryCatch(read_expression_fixture(c("probe\tA\tB", "p1\t-\tn/a", "",
        "p2\t#N/A\tbad"), c("A", "B"), c("p1", "p2")), error = conditionMessage)
    expect_match(msg, "4 invalid measurements", fixed = TRUE)
    expect_match(msg, "line 4, probe p2, sample A", fixed = TRUE)
    expect_false(grepl("probe p2, sample B", msg, fixed = TRUE))
})

test_that("record widths are checked before header association", {
    patients <- c("S1", "S2")
    rows <- c("probe,S1,S2", "p1,10.5,12", "p2,9,7,20")
    expect_error(read_expression_fixture(rows, patients, c("p1", "p2"), ","),
                 "line 3 has 4 fields, expected 3")
    # A later over-wide record must not spill into the following record.
    rows <- c("probe,S1,S2", paste0("p", 1:8, ",3,4"), "p9,1,25,4", "p10,3,4")
    expect_error(read_expression_fixture(rows, patients, paste0("p", 1:10), ","),
                 "line 10 has 4 fields, expected 3")
    expect_error(read_expression_fixture(c("S1,S2", "2,101", "4,1,02", "1,103"),
        patients, as.character(1:3), ","), "line 3 has 3 fields, expected 2")
    expect_error(read_expression_fixture(c("probe,S1", 'p1,"1.2"3'), "S1", "p1", ","),
                 "line 2 has an invalid quoted field")
    expect_error(read_expression_fixture(c("probe,S1", 'p1,"1.2'), "S1", "p1", ","),
                 "line 2 has an invalid quoted field")
})

test_that("a missing final newline is valid and a NUL cannot truncate a measurement", {
    path <- withr::local_tempfile()
    text <- charToRaw("probe,S\np1,12")
    writeBin(text, path)
    read <- function() BreastSubtypeR:::.read_shiny_expression(path, ",", "S", "p1", FALSE)
    expect_identical(read(), matrix(12, 1L, dimnames = list("p1", "S")))
    writeBin(c(text, as.raw(0L), charToRaw("invalid\n")), path)
    expect_error(read(), "NUL byte", fixed = TRUE)
})

test_that("short records stop at the physical line and explicit missing fields retain their position", {
    patients <- c("S1", "S2", "S3")
    probes <- c("p1", "p2")
    header <- "probe,S1,S2,S3"
    for (row in c("p2,4,6", "p2,4,5")) {
        expect_error(read_expression_fixture(c(header, "p1,1,2,3", "", row), patients, probes, ","),
                     "line 4 has 3 fields, expected 4")
    }
    expect_error(read_expression_fixture(c(header, "p1,1,2", "p2,4,5"), patients, probes, ","),
                 "line 2 has 3 fields, expected 4")
    expect_identical(read_expression_fixture(c(header, "p1,1,,3", "p2,4,5,"), patients, probes, ","),
        matrix(c(1, 4, NA, 5, 3, NA), 2L, dimnames = list(probes, patients)))
})

test_that("only uniformly empty terminal fields can be export padding", {
    patients <- c("S1", "S2")
    probes <- c("p1", "p2")
    for (sep in c(",", "\t")) {
        lines <- c(paste(c("", patients), collapse = sep),
                   paste("p1", '"7.25"', "", sep = sep), paste("p2", 3, 4, sep = sep))
        expected <- matrix(c(7.25, 3, NA, 4), 2L, dimnames = list(probes, patients))
        for (suffix in c(sep, paste0(sep, sep))) {
            expect_identical(read_expression_fixture(paste0(lines, suffix), patients, probes, sep), expected)
        }
        expect_error(read_expression_fixture(c(lines[1], paste0(lines[-1], sep)), patients, probes, sep),
                     "line 2 has 4 fields, expected 3")
        expect_error(read_expression_fixture(c(paste0(lines[1:2], sep), lines[3]), patients, probes, sep),
                     "line 3 has 3 fields, expected 4")
        expect_error(read_expression_fixture(c(paste0(lines[1:2], sep), paste0(lines[3], sep, 99)),
                     patients, probes, sep), "nonempty header")
    }
    # Data-only delimiters must not manufacture a corner column even with
    # numeric IDs that could corroborate that interpretation.
    expect_error(read_expression_fixture(c("S1,S2", "10,2,", "20,4,"), patients, c("10","20"), ","),
                 "explicit probe header")
})

test_that("text metadata cannot hide an unknown sample column", {
    for (name in c("SYMBOL", "Description", "unknown")) {
        header <- paste("probe", name, "S1", sep = ",")
        expect_error(read_expression_fixture(c(header, "p1,text,8"), "S1", "p1", ","),
                     paste0("sample ", name, ".*annotation table"))
        expect_identical(read_expression_fixture(c(header, "p1,3,8"), c(name, "S1"), "p1", ","),
            matrix(c(3, 8), 1L, dimnames = list("p1", c(name, "S1"))))
    }
})

test_that("an omitted corner requires uniformly wider data records", {
    samples <- c("S1", "S2")
    probes <- c("10", "g2", "g3")
    for (sep in c(",", "\t")) {
        record <- function(...) paste(..., sep = sep)
        header <- record("S1", "S2")
        wide <- record(30, 1, 25)
        short <- c(record(10, 100), record(20, 200))
        expect_error(read_expression_fixture(c(header, short, wide), samples, probes, sep),
                     "line 4 has 3 fields, expected 2.*explicit probe header")
        expect_error(read_expression_fixture(c(header, wide, short), samples, probes, sep),
                     "line 2 has 3 fields, expected 2.*explicit probe header")
        # Even uniformly split decimal commas can resemble an R export.
        # A single coincident numeric annotation ID cannot establish it.
        split <- c(record(10, 1, 25), record(20, 2, 25), record(30, 3, 25))
        expect_error(read_expression_fixture(c(header, split), samples, probes, sep),
                     "explicit probe header")

        # Exercise a real R row-name export, including numeric row IDs and
        # an explicitly empty last measurement in an otherwise full record.
        path <- withr::local_tempfile()
        frame <- data.frame(S1 = c(10, 20, 30), S2 = c(100, NA, 1.25),
                            row.names = c("10", "20", "30"))
        utils::write.table(frame, path, sep = sep, col.names = TRUE, row.names = TRUE, na = "")
        expected <- matrix(c(10, 20, 30, 100, NA, 1.25), 3L,
                           dimnames = list(c("10", "20", "30"), samples))
        expect_identical(read_expression_fixture(readLines(path), samples, rownames(frame), sep),
                         expected)
    }
})

test_that("numeric annotation overlap cannot consume an unmatched sample", {
    header <- "unused,S1,S2"
    rows <- c("10,100,7", "20,200,8", "30,300,9")
    samples <- c("S1", "S2")
    for (probes in list(c("10", "g2", "g3"), c("10", "20", "30"))) {
        expect_error(read_expression_fixture(c(header, rows), samples, probes, ","),
                     "Ambiguous first expression column.*explicit probe header")
    }
    expect_error(read_expression_fixture(c(header, rows[1:2], "n/a,300,9"),
        samples, c("10", "g2", "g3"), ","), "Ambiguous.*explicit probe")
    expect_identical(read_expression_fixture(c(header, rows), samples, as.character(1:3), ","),
        matrix(c(10, 20, 30, 100, 200, 300, 7, 8, 9), 3L,
               dimnames = list(as.character(1:3), c("unused", samples))))
    expect_identical(read_expression_fixture(c("probe,S1,S2", rows), samples, "10", ","),
        matrix(c(100, 200, 300, 7, 8, 9), 3L, dimnames = list(c("10", "20", "30"), samples)))
})
