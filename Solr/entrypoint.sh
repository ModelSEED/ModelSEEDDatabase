#!/bin/bash
# ModelSEED Biochemistry SOLR — container entrypoint.
#
# Starts Solr in the foreground, waits for it to accept requests, creates
# the compounds + reactions cores from the baked-in configsets if they
# don't already exist, and optionally posts biochemistry JSON.
#
# Posting is idempotent (Solr replaces documents matching the same unique
# key) but slow with 227 MB of payload, so if POST_ON_START=true we still
# skip the post when both cores already contain documents — this makes
# container restarts against a populated /var/solr volume fast.
# Force a re-post via `docker exec ... /scripts/post_biochemistry.sh`.

set -euo pipefail

SOLR_HOST="${SOLR_HOST:-localhost}"
SOLR_PORT="${SOLR_PORT:-8983}"
SOLR_URL="http://${SOLR_HOST}:${SOLR_PORT}/solr"

log() { echo "[modelseed-solr] $*" >&2; }

wait_for_solr() {
    log "waiting for Solr to accept requests at ${SOLR_URL} ..."
    local i=0
    until curl -fsS "${SOLR_URL}/admin/info/system" > /dev/null 2>&1; do
        i=$((i + 1))
        if [ "$i" -gt 60 ]; then
            log "ERROR: Solr didn't become ready after 60s; aborting"
            exit 1
        fi
        sleep 1
    done
    log "Solr is ready."
}

core_exists() {
    local core="$1"
    curl -fsS "${SOLR_URL}/admin/cores?action=STATUS&core=${core}" 2>/dev/null \
        | grep -q "\"name\":\"${core}\""
}

core_doc_count() {
    local core="$1"
    curl -fsS "${SOLR_URL}/${core}/select?q=*:*&rows=0&wt=json" 2>/dev/null \
        | grep -oE '"numFound":[0-9]+' | head -1 | cut -d: -f2 || echo 0
}

create_core() {
    local core="$1"
    if core_exists "$core"; then
        log "core '${core}' already exists — leaving as-is"
        return 0
    fi
    log "creating core '${core}' from configset '${core}' ..."
    solr create -c "$core" -d "$core"
    log "core '${core}' created."
}

# 0. Initialize /var/solr from the image if the persistent volume is
#    empty. The stock image's docker-entrypoint.sh normally does this
#    before invoking solr-foreground; we call the init step directly
#    since our entrypoint takes over the container.
init-var-solr

# 1. Boot Solr in the foreground under the shell as a background job so we
#    can create cores against it. `solr-foreground` is the standard
#    command from the official Solr image.
log "starting Solr ..."
solr-foreground &
SOLR_PID=$!
trap 'log "shutting down (pid ${SOLR_PID}) ..."; kill -TERM "${SOLR_PID}" 2>/dev/null || true; wait "${SOLR_PID}" 2>/dev/null || true; exit 0' TERM INT

# 2. Wait for it to accept requests.
wait_for_solr

# 3. Create cores (idempotent). Environment layout:
#      SOLR_ENVIRONMENTS unset / empty  → bare cores "compounds", "reactions"
#                                          (new nested schema)
#      SOLR_ENVIRONMENTS="staging"      → "compounds_staging", "reactions_staging"
#                                          (new nested schema)
#      SOLR_ENVIRONMENTS="prod"         → "compounds_prod", "reactions_prod"
#                                          (legacy flat schema — master-era layout)
#      SOLR_ENVIRONMENTS="staging prod" → all four cores in one Solr instance,
#                                          each using the configset appropriate
#                                          for its consumer UI
#
#    Configset selection is by env name: env=="prod" uses the *_legacy
#    configsets (matching what the production UI currently expects);
#    everything else (staging, dev, test, ...) uses the new nested
#    schema. Both configset flavours are baked into the image.
if [ -z "${SOLR_ENVIRONMENTS:-}" ]; then
    ENVS=("")
else
    # shellcheck disable=SC2206  # intentional word-splitting on the env list
    ENVS=(${SOLR_ENVIRONMENTS})
fi

# Pick the compounds/reactions configset name for a given env name.
configset_for_env() {
    local env="$1"
    local base="$2"        # "compounds" or "reactions"
    if [ "$env" = "prod" ]; then
        echo "${base}_legacy"
    else
        echo "$base"
    fi
}

for env in "${ENVS[@]}"; do
    suffix=""
    [ -n "$env" ] && suffix="_${env}"
    cfg_cpd=$(configset_for_env "$env" compounds)
    cfg_rxn=$(configset_for_env "$env" reactions)
    if core_exists "compounds${suffix}"; then
        log "core 'compounds${suffix}' already exists — leaving as-is"
    else
        log "creating core 'compounds${suffix}' from configset '${cfg_cpd}' ..."
        solr create -c "compounds${suffix}" -d "$cfg_cpd"
        log "core 'compounds${suffix}' created."
    fi
    if core_exists "reactions${suffix}"; then
        log "core 'reactions${suffix}' already exists — leaving as-is"
    else
        log "creating core 'reactions${suffix}' from configset '${cfg_rxn}' ..."
        solr create -c "reactions${suffix}" -d "$cfg_rxn"
        log "core 'reactions${suffix}' created."
    fi
done

# 4. Optionally post biochemistry data. Skip the post if the target cores
#    are already populated (fast restart on an existing /var/solr volume).
#    When SOLR_ENVIRONMENTS has multiple envs, POST_ON_START only posts
#    into the FIRST one — additional envs are populated explicitly via
#    `docker exec ... /scripts/post_biochemistry.sh <env>`.
if [ "${POST_ON_START:-false}" = "true" ]; then
    first_env="${ENVS[0]}"
    suffix=""
    [ -n "$first_env" ] && suffix="_${first_env}"
    n_cpd=$(core_doc_count "compounds${suffix}")
    n_rxn=$(core_doc_count "reactions${suffix}")
    if [ "${n_cpd:-0}" -gt 0 ] && [ "${n_rxn:-0}" -gt 0 ]; then
        log "cores${suffix:+ (${suffix#_})} already populated (compounds=${n_cpd}, reactions=${n_rxn}); skipping post"
    else
        log "POST_ON_START=true and cores${suffix:+ (${suffix#_})} are empty — posting biochemistry JSON ..."
        /scripts/post_biochemistry.sh "$first_env"
    fi
fi

log "startup complete; core status:"
curl -fsS "${SOLR_URL}/admin/cores?action=STATUS&indexInfo=false" | head -c 500
echo

# 5. Hand control back to Solr.
wait "${SOLR_PID}"
