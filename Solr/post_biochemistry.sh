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
# handling. Env=="prod" routes to the *_legacy configsets AND the
# *_legacy JSON payload AND the bare (unsuffixed) core names; every
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
    json_suffix="_legacy"
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
