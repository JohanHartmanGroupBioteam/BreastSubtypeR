#!/usr/bin/env bash
# CI harness only. Never build the harness checkout.
set -euo pipefail
test "$GITHUB_REPOSITORY" = "JohanHartmanGroupBioteam/BreastSubtypeR"
test "$GITHUB_REF" = "refs/heads/ci/validate-1.4.3-bioc3.23"
test "$GITHUB_EVENT_NAME" = "workflow_dispatch"
test "$(git -C "$GITHUB_WORKSPACE" rev-parse HEAD)" = "$GITHUB_SHA"
sha=99f148f9ede122aa3a120cf2869ced0990bcf740
tree=3987d2eac26a283d5e9398e58153c6449251396c
previous=8c319a50220685261c1ef6591fc1077b40168197
baseline=773135c550731713726a73c4aa7d3c90b685bb34
repo=https://github.com/JohanHartmanGroupBioteam/BreastSubtypeR.git
candidate="$RUNNER_TEMP/release143-candidate"
evidence="$RUNNER_TEMP/$1"
deps="$evidence/dependencies"
mkdir -p "$evidence"/{source,environment,selftest,gitclone,build,check,bioccheck,manual,focused,dependencies,citation}
exec > >(tee "$evidence/source/checkout.log") 2>&1
test ! -e "$candidate"
test "$(git ls-remote --heads "$repo" refs/heads/RELEASE_3_23 | cut -f1)" = "$baseline"
test "$(git ls-remote --heads "$repo" refs/heads/backport/221-experimental-ror-pc-1.4.3 | cut -f1)" = "$sha"
test "$(git ls-remote --heads https://git.bioconductor.org/packages/BreastSubtypeR refs/heads/RELEASE_3_23 | cut -f1)" = "$baseline"
git init "$candidate"
git -C "$candidate" fetch --no-tags "$repo" "$sha"
git -C "$candidate" checkout --detach "$sha"
case "$(realpath "$candidate")/" in
  "$(realpath "$GITHUB_WORKSPACE")/"*) echo "Candidate is inside the harness"; exit 2 ;;
esac
test "$(git -C "$candidate" rev-parse HEAD)" = "$sha"
test "$(git -C "$candidate" rev-parse HEAD^{tree})" = "$tree"
git -C "$candidate" merge-base --is-ancestor "$baseline" HEAD
test "$(git -C "$candidate" rev-parse HEAD^)" = "$previous"
test "$(git -C "$candidate" rev-parse "$previous^")" = "$baseline"
test "$(git -C "$candidate" rev-list --count "$baseline..HEAD")" -eq 2
test "$(sed -n 's/^Package:[[:space:]]*//p' "$candidate/DESCRIPTION")" = "BreastSubtypeR"
test "$(sed -n 's/^Version:[[:space:]]*//p' "$candidate/DESCRIPTION")" = "1.4.3"
git -C "$candidate" diff --name-only "$baseline" HEAD | LC_ALL=C sort > "$evidence/source/changed-files.txt"
diff -u "$GITHUB_WORKSPACE/.github/scripts/release143-files.txt" "$evidence/source/changed-files.txt"
test "$(wc -l < "$evidence/source/changed-files.txt")" -eq 19
git -C "$candidate" status --porcelain --untracked-files=all --ignored > "$evidence/source/pristine-before.txt"
test ! -s "$evidence/source/pristine-before.txt"
git -C "$candidate" diff --binary "$baseline" HEAD > "$evidence/source/candidate.diff"
git -C "$candidate" show "$previous:inst/CITATION" > "$evidence/citation/before-CITATION"
git -C "$candidate" show "$previous:DESCRIPTION" > "$evidence/citation/before-DESCRIPTION"
cp "$candidate/DESCRIPTION" "$deps/DESCRIPTION"
cmp "$candidate/DESCRIPTION" "$deps/DESCRIPTION"
cp "$candidate/DESCRIPTION" "$evidence/source/DESCRIPTION"
git -C "$candidate" ls-tree -r HEAD > "$evidence/source/git-tree.txt"
(cd "$candidate" && git ls-files -z | xargs -0 sha256sum) > "$evidence/source/files.sha256"
(cd "$candidate" && sha256sum --check "$GITHUB_WORKSPACE/.github/scripts/release143-source.sha256")
{
  echo "repository=$repo"
  echo "candidate_commit=$sha"
  echo "candidate_tree=$tree"
  echo "release_baseline=$baseline"
  echo "citation_reference_commit=$previous"
  echo "version=1.4.3"
  echo "candidate_directory=$candidate"
  echo "harness_commit=$GITHUB_SHA"
  echo "harness_tree=$(git -C "$GITHUB_WORKSPACE" rev-parse HEAD^{tree})"
  echo "harness_ref=$GITHUB_REF"
  echo "run=$GITHUB_SERVER_URL/$GITHUB_REPOSITORY/actions/runs/$GITHUB_RUN_ID"
} | tee "$evidence/source.txt"
{
  echo "CANDIDATE_DIR=$candidate"
  echo "CANDIDATE_SHA=$sha"
  echo "CANDIDATE_TREE=$tree"
  echo "DEPENDENCY_DIR=$deps"
  echo "EVIDENCE_DIR=$evidence"
  echo "BIOCCHECK_CI=$GITHUB_WORKSPACE/.github/scripts/bioccheck-ci.R"
} >> "$GITHUB_ENV"
echo "git fetch --no-tags $repo $sha; git checkout --detach $sha; verify SHA/tree/parent/version/19 files/pristine" >> "$evidence/commands.txt"
