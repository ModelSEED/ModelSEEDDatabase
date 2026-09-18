#!/bin/bash
# Post the compiled biochemistry JSON to the SOLR cores.
#
# Usage:
#   post_biochemistry.sh              # bare cores: "compounds", "reactions"
#                                     # (new nested schema, new-format JSON)
#   post_biochemistry.sh staging      # cores: "compounds_staging", "reactions_staging"
#                                     # (new nested schema, new-format JSON)
#   post_biochemistry.sh prod         # BARE cores: "compounds", "reactions"
#                                     # (legacy flat schema, legacy-format JSON) —
#                                     # matches the URL the production UI hits
#                                     # today (/solr/compounds/..., /solr/reactions/...
#                                     # with no suffix).
#
# The single-container multi-env layout supports staging + production
# sharing one Solr instance; see entrypoint.sh's SOLR_ENVIRONMENTS
# handling. Env=="prod" keeps the bare (unsuffixed) core names, because that
# is the URL the production UI hits, but since the 2026-09-16 cutover it posts
# the nested (non-legacy) payload like every other env; every
# other env uses the new nested layout with an env-suffixed core name.
#
# Expects the compiled JSONs (both flavours) under
# ${BIOCHEMISTRY_JSON_DIR} (default /data/compilation). Generate them
# by running BOTH compile scripts in Solr/compilation/ against
# the current Biochemistry/*.json.
#
# Idempotent: /update replaces documents matching the same unique key.
# Safe to re-run after regenerating the JSON to refresh the index.

set -euo pipefail

TARGET_ENV="${1:-}"

# env=="prod" → bare core names (matching the production UI's URL) AND
# legacy-format JSON payload. Any other env name → suffixed core name AND
# new-format payload. Unset → bare core name + new-format payload
# (standalone / dev / test default).
if [ "$TARGET_ENV" = "prod" ]; then
    suffix=""
    # CUTOVER 2026-09-16: was "_legacy". Prod now takes the nested payload,
    # matching the nested configset entrypoint.sh creates its cores from. Both
    # must change together: a flat payload into a nested core, or the reverse,
    # will not load correctly. Roll back by restoring "_legacy" here AND the
    # prod branch in entrypoint.sh's configset_for_env.
    json_suffix=""
elif [ -n "$TARGET_ENV" ]; then
    suffix="_${TARGET_ENV}"
    json_suffix=""
else
    suffix=""
    json_suffix=""
fi

SOLR_HOST="${SOLR_HOST:-localhost}"
SOLR_PORT="${SOLR_PORT:-8983}"
SOLR_URL="${SOLR_URL:-http://${SOLR_HOST}:${SOLR_PORT}/solr}"
DATA_DIR="${BIOCHEMISTRY_JSON_DIR:-/data/compilation}"

log() { echo "[post-biochemistry] $*" >&2; }

# Refuse to post a nested payload into a core that was created from a flat
# configset. The configset is bound when a core is created, so changing
# configset_for_env and restarting is NOT enough: entrypoint.sh leaves existing
# cores alone, and Solr would accept the parent documents while silently
# dropping every child. That failure is invisible until the UI returns empty
# nested queries against a core whose doc count looks plausible.
#
# _nest_path_ exists only in the nested configsets, so its presence is a
# reliable proxy for which schema a live core is running.
assert_schema_matches() {
    local core="$1"
    local want_nested="$2"     # "yes" when posting a non-legacy payload
    local fields
    fields=$(curl -fsS "${SOLR_URL}/${core}/schema/fields?wt=json" 2>/dev/null) || {
        log "WARNING: could not read schema for ${core}; skipping the check"
        return 0
    }
    local has_nested=no
    case "$fields" in *'"_nest_path_"'*) has_nested=yes;; esac
    if [ "$want_nested" = "$has_nested" ]; then
        return 0
    fi
    log "ERROR: ${core} is running the $([ "$has_nested" = yes ] && echo nested || echo flat) schema,"
    log "       but the payload for this env is $([ "$want_nested" = yes ] && echo nested || echo flat)."
    log ""
    log "       A core keeps the configset it was created from. To move it:"
    log "         curl \"\${SOLR_URL}/admin/cores?action=UNLOAD&core=${core}&deleteInstanceDir=true\""
    log "         then restart the container so entrypoint.sh recreates it,"
    log "         then re-run this script."
    log ""
    log "       Refusing to post rather than load data the core cannot represent."
    exit 1
}

want_nested=yes
[ "$json_suffix" = "_legacy" ] && want_nested=no

post_core() {
    local core="$1"
    local file="$2"
    if [ ! -f "$file" ]; then
        log "ERROR: ${file} not found — run the appropriate Compile_Biochemistry_for_SOLR*.py first"
        exit 1
    fi
    local bytes
    bytes=$(stat -c%s "$file" 2>/dev/null || wc -c < "$file")
    log "posting ${bytes} bytes to ${core} from ${file} ..."
    curl -fsS -X POST \
        -H 'Content-Type: application/json' \
        --data-binary "@${file}" \
        "${SOLR_URL}/${core}/update?commit=true" \
        > /dev/null
    log "posted ${core}."
}

assert_schema_matches "compounds${suffix}" "$want_nested"
assert_schema_matches "reactions${suffix}" "$want_nested"

post_core "compounds${suffix}"  "${DATA_DIR}/solr_compounds${json_suffix}.json"
post_core "reactions${suffix}"  "${DATA_DIR}/solr_reactions${json_suffix}.json"
# Structures core has one schema shared across all envs; the same JSON
# payload is posted to whichever env-suffixed core is being populated.
post_core "structures${suffix}" "${DATA_DIR}/solr_structures.json"

log "done. Core doc counts:"
for core in "compounds${suffix}" "reactions${suffix}" "structures${suffix}"; do
    n=$(curl -fsS "${SOLR_URL}/${core}/select?q=*:*&rows=0&wt=json" \
        | grep -oE '"numFound":[0-9]+' | head -1 | cut -d: -f2)
    log "  ${core}: ${n} docs (parent + child)"
done
