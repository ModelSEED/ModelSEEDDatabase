#!/bin/bash
# Post the compiled biochemistry JSON to the SOLR cores.
#
# Usage:
#   post_biochemistry.sh              # bare cores: "compounds", "reactions"
#                                     # (new nested schema, new-format JSON)
#   post_biochemistry.sh staging      # cores: "compounds_staging", "reactions_staging"
#                                     # (new nested schema, new-format JSON)
#   post_biochemistry.sh prod         # cores: "compounds_prod", "reactions_prod"
#                                     # (legacy flat schema, legacy-format JSON)
#
# The single-container multi-env layout supports staging + production
# sharing one Solr instance; see entrypoint.sh's SOLR_ENVIRONMENTS
# handling. Env=="prod" routes to the *_legacy configsets AND the
# *_legacy JSON payload; every other env uses the new nested layout.
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
suffix=""
[ -n "$TARGET_ENV" ] && suffix="_${TARGET_ENV}"

# Legacy-schema envs use the master-era flat JSON payload; everything
# else uses the new nested payload.
if [ "$TARGET_ENV" = "prod" ]; then
    json_suffix="_legacy"
else
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

post_core "compounds${suffix}" "${DATA_DIR}/solr_compounds${json_suffix}.json"
post_core "reactions${suffix}" "${DATA_DIR}/solr_reactions${json_suffix}.json"

log "done. Core doc counts:"
for core in "compounds${suffix}" "reactions${suffix}"; do
    n=$(curl -fsS "${SOLR_URL}/${core}/select?q=*:*&rows=0&wt=json" \
        | grep -oE '"numFound":[0-9]+' | head -1 | cut -d: -f2)
    log "  ${core}: ${n} docs (parent + child)"
done
