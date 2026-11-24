#!/usr/bin/env bash
# CI-only: resolve the lane for a target branch from .github/ci/lanes.json.
#
# Usage: resolve-ci-lane.sh LANES_JSON TARGET_BRANCH
#
# On success prints exactly these key=value lines on stdout and exits 0:
#   lane, target_branch, r_version, bioc_version, cache_namespace, lanes_sha256
# Exits 2, with ::error annotations on stderr, if the configuration is missing,
# malformed or inconsistent anywhere (not only for the requested lane), or if
# the target branch is not mapped to a lane. There is no default lane.
set -euo pipefail

fail() {
  echo "::error title=CI lane::$*" >&2
  exit 2
}

[ "$#" -eq 2 ] || fail "usage: resolve-ci-lane.sh LANES_JSON TARGET_BRANCH"
config="$1"
target="$2"
command -v jq > /dev/null || fail "jq is required"
[ -f "$config" ] || fail "lane configuration not found: $config"
[[ "$target" =~ ^[A-Za-z0-9._/-]+$ ]] || fail "invalid target branch name: '$target'"
jq -e 'type == "object"' "$config" > /dev/null 2>&1 ||
  fail "lane configuration is not a JSON object: $config"

problems="$(jq -r '
  def okstr(re): if type == "string" then test(re) else false end;
  def fields: ["bioc_version", "cache_namespace", "r_version"];
  [
    (if (keys | sort) != ["branches", "description", "lanes", "schema"]
     then "top-level keys must be exactly: schema, description, branches, lanes" else empty end),
    (if .schema != 1 then "schema must be 1" else empty end),
    (if (.branches | type) != "object" or (.branches | length) == 0
     then "branches must be a non-empty object" else empty end),
    (if (.lanes | type) != "object" or (.lanes | length) == 0
     then "lanes must be a non-empty object" else empty end),
    (. as $c
     | if ($c.branches | type) == "object" then
         $c.branches | to_entries[] | .key as $b | .value as $l
         | if ($l | okstr("^[a-z0-9][a-z0-9.-]*$") | not)
           then "branch \($b): lane id must match ^[a-z0-9][a-z0-9.-]*$"
           elif (($c.lanes | type) != "object") or ($c.lanes | has($l) | not)
           then "branch \($b): maps to undefined lane \($l)"
           else empty end
       else empty end),
    (if (.lanes | type) == "object" then
       .lanes | to_entries[] | .key as $id | .value as $v
       | if ($id | okstr("^[a-z0-9][a-z0-9.-]*$") | not)
         then "lane \($id): id must match ^[a-z0-9][a-z0-9.-]*$"
         elif ($v | type) != "object"
         then "lane \($id): must be an object"
         elif ($v | keys | sort) != fields
         then "lane \($id): keys must be exactly \(fields | join(", "))"
         else
           (if ($v.r_version | okstr("^(devel|[0-9]+\\.[0-9]+)$") | not)
            then "lane \($id): r_version must be X.Y or devel" else empty end),
           (if ($v.bioc_version | okstr("^[0-9]+\\.[0-9]+$") | not)
            then "lane \($id): bioc_version must be X.Y" else empty end),
           (if ($v.cache_namespace | okstr("^[A-Za-z0-9][A-Za-z0-9._-]*$") | not)
            then "lane \($id): cache_namespace must match ^[A-Za-z0-9][A-Za-z0-9._-]*$"
            elif ($v.bioc_version | okstr("^[0-9]+\\.[0-9]+$"))
                 and ($v.cache_namespace | contains("bioc-" + $v.bioc_version) | not)
            then "lane \($id): cache_namespace must contain bioc-\($v.bioc_version)"
            else empty end)
         end
     else empty end),
    (if (.lanes | type) == "object" then
       [.lanes[] | objects | .cache_namespace | strings] | group_by(.)
       | map(select(length > 1) | .[0])[]
       | "cache_namespace \(.) is used by more than one lane"
     else empty end)
  ] | .[]' "$config")" || fail "could not evaluate lane configuration: $config"

if [ -n "$problems" ]; then
  while IFS= read -r p; do
    echo "::error title=CI lane::$p" >&2
  done <<< "$problems"
  exit 2
fi

lane="$(jq -r --arg b "$target" '.branches[$b] // empty' "$config")"
[ -n "$lane" ] || fail "target branch '$target' is not mapped to a lane in $config"
field() { jq -r --arg l "$lane" --arg f "$1" '.lanes[$l][$f]' "$config"; }
sha="$(sha256sum "$config")"
printf '%s\n' \
  "lane=$lane" \
  "target_branch=$target" \
  "r_version=$(field r_version)" \
  "bioc_version=$(field bioc_version)" \
  "cache_namespace=$(field cache_namespace)" \
  "lanes_sha256=${sha%% *}"
