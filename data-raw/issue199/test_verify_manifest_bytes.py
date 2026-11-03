#!/usr/bin/env python3
"""Maintainer-only verifier tests using synthetic bytes in temporary Git repos.

No package datasets or classifiers are loaded. No index, ref or configuration
in the caller's repository is changed. LF/CRLF checkouts are simulated locally,
not evidence of execution on multiple operating systems.
"""

import csv
import hashlib
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.dont_write_bytecode = True
import verify_manifest_bytes as verifier


SCRIPT = Path(verifier.__file__).resolve()
TEXT_PATH = 'inst/RshinyTest/example.txt'
BINARY_PATH = 'tests/testthat/fixtures/example.rds'
SOURCE_PATH = 'data/OSLO2EMIT0obj.RData'
# A BOM and a non-ASCII UTF-8 character make byte preservation observable.
TEXT = b'\xef\xbb\xbfkey\tvalue\nalpha\tcaf\xc3\xa9\n'
HISTORICAL = b'\xef\xbb\xbfkey\tvalue\r\nalpha\tcaf\xc3\xa9\r\n'


def digest(value):
    return hashlib.sha256(value).hexdigest()


class VerifyBytesTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='bsr235-verifier-')
        self.root = Path(self.temp.name).resolve()

        def cleanup():
            # Bound recursive cleanup to this test's newly allocated temp dir.
            assert self.root.parent == Path(tempfile.gettempdir()).resolve()
            assert self.root.name.startswith('bsr235-verifier-')
            self.temp.cleanup()

        self.addCleanup(cleanup)
        self.repo = self.root / 'repo'
        self.repo.mkdir()
        self.git('init', '--quiet', '--object-format=sha1')
        self.files = {TEXT_PATH: TEXT, BINARY_PATH: b'\x00\xff\r\n\x80binary',
                      SOURCE_PATH: b'\x1f\x8b\x00\x80synthetic-container'}
        self.ref = self.make_commit(self.files)
        self.shiny = self.root / 'shiny.tsv'
        self.regression = self.root / 'regression.tsv'
        self.ledger = self.root / 'ledger.tsv'
        self.shiny_rows = [dict(file='example.txt', object='OSLO2EMIT0obj',
            source_sha256=digest(self.files[SOURCE_PATH]), kind='GEX', rows='1',
            columns='1', units='synthetic', sha256=digest(TEXT), bytes=str(len(TEXT)))]
        self.regression_rows = [dict(file='example.rds', bytes=str(len(self.files[BINARY_PATH])),
                                    sha256=digest(self.files[BINARY_PATH]))]
        self.ledger_rows = [dict(manifest=verifier.SHINY, path=TEXT_PATH,
            old_sha256=digest(HISTORICAL), old_bytes=str(len(HISTORICAL)),
            canonical_sha256=digest(TEXT), canonical_bytes=str(len(TEXT)),
            old_representation='historical-CRLF', canonical_representation='git-blob-LF',
            pinned_commit=self.ref, explanation='Synthetic byte-preservation test')]
        self.write_inputs()

    def git(self, *args, data=None, repo=None):
        target = Path(repo if repo is not None else self.repo).resolve()
        self.assertIn(self.root, target.parents)  # every write stays in this test's temp allocation
        env = dict(verifier.git_environment(), GIT_AUTHOR_NAME='Verifier test', GIT_AUTHOR_EMAIL='test@example.invalid',
                   GIT_COMMITTER_NAME='Verifier test', GIT_COMMITTER_EMAIL='test@example.invalid',
                   GIT_AUTHOR_DATE='2000-01-01T00:00:00Z', GIT_COMMITTER_DATE='2000-01-01T00:00:00Z')
        return subprocess.run(['git', '-C', str(target), *args], input=data, env=env,
                              check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout

    def make_commit(self, files, repo=None, message=b'Synthetic verifier fixture\n'):
        # Plumbing creates only synthetic objects in the isolated test repo;
        # no staging, checkout or branch/history update in the real repository.
        tree = {}
        for path, raw in files.items():
            node = tree
            parts = path.split('/')
            for part in parts[:-1]:
                node = node.setdefault(part, {})
            node[parts[-1]] = self.git('hash-object', '-w', '--stdin', data=raw, repo=repo).strip()

        def make_tree(node):
            entries = []
            for name, value in sorted(node.items()):
                if isinstance(value, dict):
                    mode, kind, oid = b'040000', b'tree', make_tree(value)
                else:
                    mode, kind, oid = b'100644', b'blob', value
                entries.append(mode + b' ' + kind + b' ' + oid + b'\t' + name.encode() + b'\n')
            return self.git('mktree', data=b''.join(entries), repo=repo).strip()

        return self.git('commit-tree', make_tree(tree).decode(), data=message, repo=repo).decode().strip()

    def write_table(self, path, rows):
        with path.open('w', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)

    def write_inputs(self):
        self.write_table(self.shiny, self.shiny_rows)
        self.write_table(self.regression, self.regression_rows)
        self.write_table(self.ledger, self.ledger_rows)

    def run_verifier(self, historical=False, extra=(), env=None):
        args = [sys.executable, '-B', str(SCRIPT), '--repo', str(self.repo), '--ref', self.ref,
                '--shiny-manifest', str(self.shiny), '--regression-manifest', str(self.regression)]
        if historical:
            args += ['--historical-crlf', str(self.ledger)]
        return subprocess.run(args + list(extra), stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)

    def assert_failure(self, text, **kwargs):
        result = self.run_verifier(**kwargs)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(text.encode(), result.stderr)

    def test_raw_binary_and_historical_bytes(self):
        result = self.run_verifier(historical=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(b'CANONICAL PASS: 2 file hash/size pairs; 1 source-container hashes.', result.stdout)
        self.assertIn(b'HISTORICAL-CRLF PASS: 1 pinned representations', result.stdout)
        self.assertIn(self.ref.encode(), result.stdout)

    def test_simulated_checkout_line_endings_do_not_change_git_result(self):
        output = None
        working = self.repo / TEXT_PATH
        working.parent.mkdir(parents=True)
        for raw in (TEXT, HISTORICAL, b'non-newline working-tree mutation'):
            working.write_bytes(raw)
            result = self.run_verifier(historical=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            if output is not None:
                self.assertEqual(result.stdout, output)
            output = result.stdout
            self.assertEqual(working.read_bytes(), raw)  # verifier is read-only

    def test_committed_non_newline_mutation_fails(self):
        changed = dict(self.files, **{TEXT_PATH: TEXT.replace(b'alpha', b'Alpha')})
        self.assertEqual(len(changed[TEXT_PATH]), len(TEXT))
        self.ref = self.make_commit(changed)
        self.assert_failure('SHA-256 mismatch')

    def test_wrong_byte_count_fails_independently_of_hash(self):
        self.shiny_rows[0]['bytes'] = str(len(TEXT) + 1)
        self.write_inputs()
        self.assert_failure('byte-count mismatch')

    def test_missing_target_and_source_blobs_fail(self):
        for path in (TEXT_PATH, BINARY_PATH, SOURCE_PATH):
            with self.subTest(path=path):
                self.ref = self.make_commit({p: raw for p, raw in self.files.items() if p != path})
                self.assert_failure('missing/unreadable committed blob')

    def test_missing_manifest_and_ref_fail(self):
        self.shiny = self.root / 'missing.tsv'
        self.assert_failure('missing.tsv')
        self.shiny = self.root / 'shiny.tsv'
        self.ref = 'not-a-ref'
        self.assert_failure('FAIL:')

    def test_ambiguous_branch_tag_ref_fails(self):
        self.git('update-ref', 'refs/heads/ambiguous', self.ref)
        other = self.make_commit(dict(self.files, **{TEXT_PATH: b'other content\n'}))
        self.git('update-ref', 'refs/tags/ambiguous', other)
        self.ref = 'ambiguous'
        self.assert_failure('ambiguous')
        self.ref = 'refs/heads/ambiguous'
        self.assertEqual(self.run_verifier().returncode, 0)

    def make_caller(self):
        caller = self.root / 'caller'
        caller.mkdir()
        self.git('init', '--quiet', '--object-format=sha1', repo=caller)
        ref = self.make_commit(dict(self.files, **{TEXT_PATH: b'caller-only content\n'}), repo=caller)
        self.git('update-ref', 'refs/heads/selected', ref, repo=caller)
        self.git('read-tree', ref, repo=caller)  # an index whose bytes must survive
        self.git('config', '--local', 'review.caller', 'must-remain-unchanged', repo=caller)
        (caller / '.git' / 'shallow').write_bytes(b'')
        return caller

    def caller_overrides(self, caller):
        gitdir = caller / '.git'
        return dict(GIT_DIR=str(gitdir), GIT_WORK_TREE=str(caller),
                    GIT_COMMON_DIR=str(gitdir), GIT_INDEX_FILE=str(gitdir / 'index'),
                    GIT_OBJECT_DIRECTORY=str(gitdir / 'objects'),
                    GIT_ALTERNATE_OBJECT_DIRECTORIES=str(gitdir / 'objects'),
                    GIT_NAMESPACE='caller-namespace', GIT_SHALLOW_FILE=str(gitdir / 'shallow'),
                    GIT_CONFIG_COUNT='1', GIT_CONFIG_KEY_0='core.worktree',
                    GIT_CONFIG_VALUE_0=str(caller))

    def snapshot(self, directory):
        # Includes refs, loose objects, index, config and all other caller files.
        return {p.relative_to(directory).as_posix(): p.read_bytes()
                for p in directory.rglob('*') if p.is_file()}

    def test_verifier_ignores_repository_overrides(self):
        caller = self.make_caller()
        self.git('update-ref', 'refs/heads/selected', self.ref)
        self.ref = 'refs/heads/selected'
        before = self.snapshot(caller)
        overrides = self.caller_overrides(caller)
        cases = [{key: value} for key, value in overrides.items() if not key.startswith('GIT_CONFIG_')]
        cases += [{key: value for key, value in overrides.items() if key.startswith('GIT_CONFIG_')}, overrides]
        for override in cases:
            with self.subTest(overrides=sorted(override)):
                result = self.run_verifier(historical=True, env=dict(os.environ, **override))
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn(self.ledger_rows[0]['pinned_commit'].encode(), result.stdout)
                self.assertEqual(self.snapshot(caller), before)

    def test_verifier_cannot_read_caller_object_store(self):
        caller = self.make_caller()
        # Identical target blobs in a commit absent from the requested repository.
        self.ref = self.make_commit(self.files, repo=caller, message=b'Caller-only commit\n')
        before = self.snapshot(caller)
        for key in ('GIT_OBJECT_DIRECTORY', 'GIT_ALTERNATE_OBJECT_DIRECTORIES'):
            with self.subTest(variable=key):
                self.assert_failure('FAIL:', env=dict(os.environ, **{key: str(caller / '.git' / 'objects')}))
                self.assertEqual(self.snapshot(caller), before)

    def test_test_helpers_do_not_modify_caller_repository(self):
        caller = self.make_caller()
        before = self.snapshot(caller)
        # Execute the original ambiguity test in a child: its setUp (including
        # git init), object creation and branch/tag writes inherit these values.
        for override in ({'GIT_DIR': str(caller / '.git')}, self.caller_overrides(caller)):
            with self.subTest(overrides=sorted(override)):
                result = subprocess.run([sys.executable, '-B', str(Path(__file__).resolve()),
                    'VerifyBytesTests.test_ambiguous_branch_tag_ref_fails'],
                    env=dict(os.environ, **override), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                self.assertEqual(self.snapshot(caller), before)
                self.assertEqual(result.returncode, 0, result.stderr)

    def test_linked_worktree_with_inherited_overrides(self):
        caller = self.make_caller()
        before = self.snapshot(caller)
        linked = self.root / 'linked'
        self.git('worktree', 'add', '--quiet', '--detach', str(linked), self.ref)
        self.assertTrue((linked / '.git').is_file())
        self.repo = linked
        result = self.run_verifier(historical=True, env=dict(os.environ, **self.caller_overrides(caller)))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.snapshot(caller), before)

    def test_ambiguity_with_disabled_warnings_and_identical_blobs(self):
        first = self.ref
        second = self.make_commit(self.files, message=b'Different commit, identical target blobs\n')
        self.assertNotEqual(first, second)
        self.assertEqual(self.git('rev-parse', first + '^{tree}'), self.git('rev-parse', second + '^{tree}'))
        self.git('update-ref', 'refs/heads/collision', first)
        self.git('update-ref', 'refs/tags/collision', second)
        self.git('config', '--local', 'core.warnAmbiguousRefs', 'false')
        config_before = (self.repo / '.git' / 'config').read_bytes()
        self.ref = 'collision'
        self.assert_failure('ambiguous')
        for ref, expected in (('refs/heads/collision', first), ('refs/tags/collision', second),
                              (first, first), (second, second)):
            with self.subTest(ref=ref):
                self.ref = ref
                result = self.run_verifier()
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn(('Resolved commit: ' + expected).encode(), result.stdout)
        self.assertEqual((self.repo / '.git' / 'config').read_bytes(), config_before)

    def test_duplicates_and_ambiguous_paths_fail(self):
        self.shiny_rows *= 2
        self.write_inputs()
        self.assert_failure('duplicate/case-ambiguous entry')
        self.shiny_rows = [dict(self.shiny_rows[0]), dict(self.shiny_rows[0], file='EXAMPLE.TXT')]
        self.write_inputs()
        self.assert_failure('duplicate/case-ambiguous entry')
        self.shiny_rows = [dict(self.shiny_rows[0], file='../example.txt')]
        self.write_inputs()
        self.assert_failure('ambiguous manifest/path')

    def test_malformed_columns_and_counts_fail(self):
        self.shiny.write_bytes(b'file\tfile\nexample.txt\texample.txt\n')
        self.assert_failure('invalid/duplicate columns')
        self.shiny_rows[0]['bytes'] = '1.0'
        self.write_inputs()
        self.assert_failure('invalid byte-count')

    def test_bad_or_ambiguous_source_hash_fails(self):
        self.shiny_rows[0]['source_sha256'] = '0' * 64
        self.write_inputs()
        self.assert_failure('SOURCE-CONTAINER')
        self.shiny_rows[0]['source_sha256'] = digest(self.files[SOURCE_PATH])
        second = dict(self.shiny_rows[0], file='second.txt', source_sha256='0' * 64)
        self.shiny_rows.append(second)
        self.ref = self.make_commit(dict(self.files, **{'inst/RshinyTest/second.txt': TEXT}))
        self.write_inputs()
        self.assert_failure('ambiguous source_sha256')

    def test_canonical_mode_never_accepts_historical_crlf_digest(self):
        self.shiny_rows[0].update(sha256=digest(HISTORICAL), bytes=str(len(HISTORICAL)))
        self.write_inputs()
        self.assert_failure('SHA-256 mismatch', historical=True)

    def test_historical_hash_and_size_fail_separately(self):
        self.ledger_rows[0]['old_sha256'] = '0' * 64
        self.write_inputs()
        self.assert_failure('SHA-256 mismatch', historical=True)
        self.ledger_rows[0].update(old_sha256=digest(HISTORICAL), old_bytes=str(len(HISTORICAL) + 1))
        self.write_inputs()
        self.assert_failure('byte-count mismatch', historical=True)

    def test_historical_mode_refuses_mixed_line_endings(self):
        mixed = TEXT.replace(b'\n', b'\r\n', 1)
        self.ref = self.make_commit(dict(self.files, **{TEXT_PATH: mixed}))
        self.shiny_rows[0].update(sha256=digest(mixed), bytes=str(len(mixed)))
        self.ledger_rows[0].update(canonical_sha256=digest(mixed), canonical_bytes=str(len(mixed)), pinned_commit=self.ref)
        self.write_inputs()
        self.assert_failure('strict LF input required', historical=True)

    def test_historical_mode_preserves_absent_final_newline(self):
        raw, crlf = b'a\nb', b'a\r\nb'
        self.ref = self.make_commit(dict(self.files, **{TEXT_PATH: raw}))
        self.shiny_rows[0].update(sha256=digest(raw), bytes=str(len(raw)))
        self.ledger_rows[0].update(canonical_sha256=digest(raw), canonical_bytes=str(len(raw)),
            old_sha256=digest(crlf), old_bytes=str(len(crlf)), pinned_commit=self.ref)
        self.write_inputs()
        result = self.run_verifier(historical=True)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_duplicate_ledger_and_options_fail(self):
        self.ledger_rows *= 2
        self.write_inputs()
        self.assert_failure('duplicate/case-ambiguous ledger entry', historical=True)
        self.assert_failure('duplicate option', extra=('--ref', self.ref))


if __name__ == '__main__':
    unittest.main(verbosity=2)
