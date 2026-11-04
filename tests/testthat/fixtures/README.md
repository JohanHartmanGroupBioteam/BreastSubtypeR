# Regression fixture provenance

The public examples were expanded for issue #199. Tests of historical behavior
must not silently acquire a new input or expected result when public data changes.

The compact RDS files in `historical-cb7a7c70/` reproduce the complete
pre-#199 public objects present at main commit
`cb7a7c70e0f94e1690d15ac2ec2b464a40934f36`. Their original complete RData
containers remain immutable in the retained Stage 3A/3B archives; their hashes
and sizes are preserved in `original-files.tsv` and `manifest.tsv`.
The snapshot name does not identify their original generation environment.
Original save-event software and RNG provenance remains incomplete.
These are test-only fixtures, not additional public datasets.

Stage 3C removes only duplicate historical assay storage. Each compact fixture
independently retains the original feature annotations, sample metadata,
SummarizedExperiment metadata, mapped inputs and frozen results. The loader
restores the first 356 OSLO or 5,000 TCGA assay rows from the installed public
object, whose historical values were preserved exactly. Row/sample order,
storage type and every assay attribute are checked against frozen original
content before reconstruction. The complete reconstructed historical object
must then match a separate frozen digest. Altering public values or either
fixture cannot silently redefine the historical baseline.

The runtime uses base-R MD5 checks over version-2 XDR serialization, normalizing
only the two writer/minimum-reader version integers in the serialization
header (bytes 7-14). No value, type, attribute or ordering is normalized.
The corresponding SHA-256 hashes, original container/payload hashes and compact
container hashes are in `manifest.tsv`. These digests were computed directly
from the immutable originals, not from a new classifier run. Independent
validation also compares the full reconstructed objects with the original
archives using R `identical()`. The checks use no private archive path or
network access when the package tests run.

- OSLO retains 356 rows, 102 samples, the historical mapped inputs and stored
  results. AUTO calls/entropy, the #198 independent reference checks, PCA-PAM50
  ordering, no-call behavior and visualization baselines use this fixed input.
- TCGA retains 5,000 rows, 100 samples, 37 mapped NC and 178 mapped SSP rows.
  The incomplete proliferation tests still assert 7/11 and 8/18 represented
  genes. Their expectations and independent scoring formulas are unchanged.
- `aims_1.5.6_multisample.rds` and `cihc_itr_1.5.5.rds` retain their
  original bytes. `make-cihc-itr-1.5.5.R` retains its source unchanged
  (Git checkout line endings may differ). Do not regenerate these fixtures
  with current code. The existing version-specific correction exclusions remain.

`issue199-regression-manifest.tsv` records SHA-256 and byte counts of the five
listed files' raw committed Git blob contents; its `file` names resolve under
`tests/testthat/fixtures/`. Issue #235 corrects only the hash/size fields for
`issue199-TCGA-B-C-changes.tsv` and `make-cihc-itr-1.5.5.R`. Their previous values
identify the strict CRLF representation of the pinned LF blobs, not different
scientific values. All three binary entries and all fixture payloads remain
unchanged. The full-repository
[historical correction ledger](../../../data-raw/issue199/manifests/line-ending-corrections-v1.tsv)
preserves old/new hashes, sizes and the comparison commit. The original input
hashes and serialized-object fingerprints in the other manifests are unchanged.
See the source-only [verification instructions](../../../data-raw/issue199/README.md#file-byte-provenance-and-issue-235).
This Git-blob check does not attest to arbitrary checkout or tarball bytes;
ordinary installed-package tests have no new Git or Python requirement.

All 24 test files that previously loaded a public example now explicitly load
the historical snapshot through `load_historical_example()`. Those loader migrations changed no numerical expected value, method setting,
reference formula or tolerance. The 25th direct dependency in the Stage 3A inventory is the protected
old fixture generator, which was not edited or run. Generic Mapping, routing,
identifier and plotting regressions therefore retain exactly their original
input meaning as well as their expected behavior.

`test-expanded-examples.R` separately loads the new public objects and checks
their dimensions, historical value/metadata preservation, current Mapping
outputs, reproducible complete AUTO results, sample order, no-calls and entropy.
It also checks complete TCGA proliferation coverage against the independent
#198 arithmetic rather than changing the old incomplete-panel expectations.

`issue199-cb7a7c70-calls.rds` contains only the B/C summary tables from the
accepted Stage 3A outputs, preserving calls and entropy without duplicating
complete replacement objects. B means current code on the historical input;
C means the same code on the expanded input. `issue199-TCGA-B-C-changes.tsv`
is the accepted exact per-sample/per-method difference oracle: 64 calls in 40
samples (8 genefu.robust, 4 ssBC, 6 ssBC.v2, 5 cIHC, 6 cIHC.itr, 8 PCAPAM50,
16 AIMS, 11 sspbc). The tests reproduce B and C and compare with this independent
retained evidence. OSLO B/C mapped inputs and complete results are identical;
historical A calls and entropy remain unchanged. Historical ROR-P and some
other details predate package corrections and are not asserted to equal C.

Stage 3A used BreastSubtypeR 1.5.12 at the commit above, R 4.5.2,
Mersenne-Twister/Inversion/Rejection, seed 1992026 before Mapping(method="max",
impute=TRUE) and seed 20260831 before BS_Multi(methods="AUTO", Subtype=FALSE,
hasClinical=FALSE, er_subsample="target"). RawCounts is FALSE for OSLO and TRUE
for TCGA. Existing method defaults are preserved. The new complete-result
replay is strict; the separate entropy formula comparison allows 1e-12 for
floating-point evaluation. Historical tests keep their original tolerances.

These files are included under tests/ in the source package and by
`R CMD INSTALL --install-tests`, so installed-package tests do not require
any private archive path. They are not installed in ordinary installations.
The completed #199 public-source and field-level review includes these fixtures
and resolves their scoped redistribution gate; see the dated
[data provenance review](../../../inst/DATA-PROVENANCE.md). Source attribution
and reuse conditions still apply. This is not institutional approval or blanket
clearance for future data, uses or releases; historical derivation limitations
remain documented in that review.

To reproduce the compact fixtures, run the source-only builder from a full
checkout with the original complete objects and a new output directory:

~~~text
Rscript data-raw/issue199/compact-historical-fixtures.R --historical-dir=/retained/original-complete-objects --package-dir=/checkout --out=/new/fixture-build
~~~

The builder verifies original and approved public file SHA-256 values,
constructs no predictions, and verifies full historical object/component
identity. It requires SummarizedExperiment and digest for generation; the
installed test loader requires no new package dependency. Original versioned
AIMS/cIHC fixtures and the B/C call oracle remain byte-identical.
The Stage 3B scripts `migrate_regression.ps1` and
`freeze_regression_oracles.R` document the original extraction and are retained
separately in the Stage 3B archive's `scripts/` directory; they are not tracked
in the current Git tree. Do not rerun
the guarded historical scientific fixture generator.

The Stage 3C packaging assertions check the compact physical containers and
the frozen historical content instead of requiring duplicate gzip containers
inside the package. Scientific expectations and tolerances are unchanged.
New negative tests require failures for modified assay values, row/sample
order, storage type, annotations, metadata, mapped inputs or frozen entropy.
