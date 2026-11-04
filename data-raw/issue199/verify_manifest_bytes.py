#!/usr/bin/env python3
"""Read-only #199 Git-blob verifier; Python 3.9+ standard library and local Git.

Manifest arguments are explicit filesystem inputs supplying expectations.
Targets are read ONLY from the resolved commit, never from the working tree.
No classifiers, export writers, checkout filters, downloads or repairs are run.
"""

import argparse
import csv
import hashlib
import io
import os
from pathlib import Path
import re
import subprocess
import sys


SHINY = 'data-raw/issue199/manifests/shiny-tables.tsv'
REGRESSION = 'tests/testthat/fixtures/issue199-regression-manifest.tsv'
ROOTS = {SHINY: 'inst/RshinyTest/', REGRESSION: 'tests/testthat/fixtures/'}
SHINY_FIELDS = ['file', 'object', 'source_sha256', 'kind', 'rows', 'columns',
                'units', 'sha256', 'bytes']
LEDGER_FIELDS = ['manifest', 'path', 'old_sha256', 'old_bytes', 'canonical_sha256',
                 'canonical_bytes', 'old_representation', 'canonical_representation',
                 'pinned_commit', 'explanation']
SOURCES = {name: 'data/' + name + '.RData' for name in ('OSLO2EMIT0obj', 'TCGABRCAobj')}


class VerificationError(Exception):
    pass


def git_environment():
    # -C alone does not override GIT_DIR, common/object directories, indexes,
    # namespaces or environment-injected config. Clear all inherited Git
    # overrides, including future ones, before reinstating our own controls.
    # Git still discovers .git files normally, including linked worktrees.
    env = {key: value for key, value in os.environ.items()
           if not key.upper().startswith('GIT_')}
    env.update(GIT_NO_REPLACE_OBJECTS='1', GIT_NO_LAZY_FETCH='1',
               GIT_ALLOW_PROTOCOL='', GIT_TERMINAL_PROMPT='0', LC_ALL='C')
    return env


class GitBytes:
    def __init__(self, repo):
        self.repo = str(Path(repo).resolve())
        self.env = git_environment()

    def git(self, *args):
        # Disallow transports even if an older Git ignores GIT_NO_LAZY_FETCH.
        result = subprocess.run(['git', '-C', self.repo, '-c', 'protocol.allow=never',
                                 '-c', 'core.warnAmbiguousRefs=true', *args],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=self.env)
        # rev-parse can exit zero while warning that a branch/tag is ambiguous.
        # Fail closed on Git diagnostics instead of silently choosing that ref.
        if result.returncode or result.stderr:
            raise VerificationError(result.stderr.decode('utf-8', errors='replace').strip())
        return result.stdout  # bytes: never text=True, decode/re-encode or --filters

    def commit(self, ref):
        return self.git('rev-parse', '--verify', '--end-of-options', ref + '^{commit}').decode('ascii').strip()

    def blob(self, commit, path):
        try:
            return self.git('cat-file', 'blob', commit + ':' + path)
        except VerificationError as exc:
            raise VerificationError(f'missing/unreadable committed blob {commit}:{path}: {exc}') from exc


def read_table(path, fields):
    # Only the expectation TSV is decoded. Scientific target bytes are not.
    rows = list(csv.reader(io.StringIO(Path(path).read_bytes().decode('utf-8'), newline=''),
                           delimiter='\t', strict=True))
    if not rows or rows[0] != fields or len(rows) < 2:
        raise VerificationError(f'{path}: invalid/duplicate columns or empty manifest')
    if any(len(row) != len(fields) for row in rows[1:]):
        raise VerificationError(f'{path}: malformed row')
    return [dict(zip(fields, row)) for row in rows[1:]]


def target_path(manifest, name):
    if manifest not in ROOTS or not re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.-]*', name):
        raise VerificationError(f'ambiguous manifest/path: {manifest!r}, {name!r}')
    return ROOTS[manifest] + name


def check_bytes(label, raw, expected_hash, expected_size=None):
    if not re.fullmatch(r'[0-9a-f]{64}', expected_hash):
        raise VerificationError(f'{label}: invalid SHA-256 expectation')
    actual = hashlib.sha256(raw).hexdigest()
    problems = []
    if actual != expected_hash:
        problems.append(f'SHA-256 mismatch: expected {expected_hash}, got {actual}')
    if expected_size is not None:
        if not re.fullmatch(r'0|[1-9][0-9]*', expected_size):
            raise VerificationError(f'{label}: invalid byte-count expectation')
        if len(raw) != int(expected_size):
            problems.append(f'byte-count mismatch: expected {expected_size}, got {len(raw)}')
    if problems:
        raise VerificationError(f'{label}: ' + '; '.join(problems))
    print(f'PASS {label} sha256={actual} bytes={len(raw)}')


def verify_canonical(git, commit, manifests):
    seen, sources = set(), {}
    for manifest, filename, fields in manifests:
        print(f'Expectation input: {Path(filename).resolve()} (schema: {manifest})')
        for row in read_table(filename, fields):
            path = target_path(manifest, row['file'])
            if path.casefold() in seen:
                raise VerificationError(f'duplicate/case-ambiguous entry: {path}')
            seen.add(path.casefold())
            check_bytes(f'CANONICAL {commit}:{path}', git.blob(commit, path), row['sha256'], row['bytes'])
            if manifest == SHINY:
                if row['object'] not in SOURCES:
                    raise VerificationError(f'{path}: unknown source object {row["object"]!r}')
                source = SOURCES[row['object']]
                if source in sources and sources[source] != row['source_sha256']:
                    raise VerificationError(f'ambiguous source_sha256 for {source}')
                sources[source] = row['source_sha256']
    for source, digest in sorted(sources.items()):
        check_bytes(f'SOURCE-CONTAINER {commit}:{source} (SHA-256 only; size observed)',
                    git.blob(commit, source), digest)
    print(f'CANONICAL PASS: {len(seen)} file hash/size pairs; {len(sources)} source-container hashes.')


def verify_historical(git, ledger):
    rows = read_table(ledger, LEDGER_FIELDS)
    seen = set()
    print(f'HISTORICAL-CRLF: in-memory byte substitution only; ledger={Path(ledger).resolve()}')
    for row in rows:
        path, pin = row['path'], row['pinned_commit']
        if target_path(row['manifest'], path.rsplit('/', 1)[-1]) != path:
            raise VerificationError(f'ambiguous ledger path: {path}')
        key = (pin, path.casefold())
        if key in seen:
            raise VerificationError(f'duplicate/case-ambiguous ledger entry: {path}')
        seen.add(key)
        if not re.fullmatch(r'[0-9a-f]{40}|[0-9a-f]{64}', pin) or git.commit(pin) != pin:
            raise VerificationError(f'ledger requires an exact pinned commit: {pin}')
        if (row['old_representation'], row['canonical_representation']) != ('historical-CRLF', 'git-blob-LF'):
            raise VerificationError(f'{path}: unsupported historical representation')
        raw = git.blob(pin, path)
        check_bytes(f'HISTORICAL-BASE {pin}:{path}', raw, row['canonical_sha256'], row['canonical_bytes'])
        if b'\r' in raw or b'\n' not in raw or b'\0' in raw:
            raise VerificationError(f'{path}: strict LF input required; refusing normalization')
        historical = raw.replace(b'\n', b'\r\n')
        check_bytes(f'HISTORICAL-CRLF {pin}:{path}', historical, row['old_sha256'], row['old_bytes'])
        print(f'  inserted_CR_bytes={raw.count(bytes([10]))}; all other bytes and final-newline presence preserved')
    print(f'HISTORICAL-CRLF PASS: {len(rows)} pinned representations; no files rewritten.')


class Once(argparse.Action):
    def __call__(self, parser, namespace, value, option_string=None):
        if getattr(namespace, self.dest) is not None:
            raise argparse.ArgumentError(self, 'duplicate option')
        setattr(namespace, self.dest, value)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ('repo', 'ref', 'shiny-manifest', 'regression-manifest'):
        parser.add_argument('--' + option, required=True, action=Once)
    parser.add_argument('--historical-crlf', metavar='LEDGER', action=Once,
                        help='separately verify the historical ledger at its pinned commits')
    args = parser.parse_args(argv)
    try:
        git = GitBytes(args.repo)
        commit = git.commit(args.ref)
        print(f'Repository: {git.repo}\nResolved commit: {commit}\nRepresentation: raw committed Git blob bytes')
        verify_canonical(git, commit, [(SHINY, args.shiny_manifest, SHINY_FIELDS),
                         (REGRESSION, args.regression_manifest, ['file', 'bytes', 'sha256'])])
        if args.historical_crlf:
            verify_historical(git, args.historical_crlf)
        print('Scope: committed blobs, not working-tree files, generated exports or package tarballs.')
        return 0
    except (VerificationError, OSError, UnicodeError, csv.Error) as exc:
        print(f'FAIL: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
