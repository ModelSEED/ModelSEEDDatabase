#!/usr/bin/env bash
# Deploy the ModelSEED biochemistry SOLR.
#
# Modeled on ModelSEED-UI's deploy_container.sh: an interactive wrapper
# that chooses (1) which environment(s) to deploy — staging, production,
# or both — and (2) whether to run against the Docker container or an
# already-installed local Solr.
#
# Usage:
#   ./deploy_container.sh                          # interactive
#   ./deploy_container.sh staging                  # deploy staging only
#   ./deploy_container.sh production               # deploy production only
#   ./deploy_container.sh both                     # deploy both
#   ./deploy_container.sh <env> --yes              # skip the confirm prompt
#   EXEC_METHOD=local  ./deploy_container.sh both  # use local Solr, not Docker
#   AUTO_YES=true      ./deploy_container.sh       # env auto-detected from branch
#
# Environment vars honored:
#   TARGET_ENV     staging | production | both        (same as first arg)
#   EXEC_METHOD    docker (default) | local
#   AUTO_YES       skip interactive confirmation when true
#   SOLR_URL       for EXEC_METHOD=local; default http://localhost:8983/solr
#
# "production" is the user-facing name; the entrypoint's env token is
# "prod" (bare-name cores + legacy configset). "staging" maps 1:1 to
# the "staging" env token (env-suffixed cores + new nested configset).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

# --- Arg parsing --------------------------------------------------------

TARGET_ENV="${TARGET_ENV:-${1:-}}"
AUTO_YES="${AUTO_YES:-false}"
EXEC_METHOD="${EXEC_METHOD:-}"

# Accept --yes as first or second arg (matches the UI script's UX).
if [[ "${TARGET_ENV}" == "--yes" ]]; then
    TARGET_ENV=""
    AUTO_YES="true"
fi
if [[ "${2:-}" == "--yes" ]]; then
    AUTO_YES="true"
fi

# --- Provenance ---------------------------------------------------------

GIT_BRANCH=$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "unknown")
GIT_COMMIT=$(git rev-parse HEAD 2>/dev/null | cut -c 1-8 || echo "unknown")
DEPLOY_DATE=$(date +"%B %-d, %Y")

# Branch-based default: main/master/production → production; otherwise staging.
if [ -z "${TARGET_ENV}" ]; then
    case "${GIT_BRANCH}" in
        main|master|production) TARGET_ENV="production" ;;
        *)                      TARGET_ENV="staging"    ;;
    esac
fi

# --- Interactive selection ---------------------------------------------

if [ -t 0 ] && [ "${AUTO_YES}" != "true" ]; then
    echo "Select deployment target:"
    echo "  1) Staging     (new nested schema — cores 'compounds_staging', 'reactions_staging')"
    echo "  2) Production  (legacy flat schema — bare cores 'compounds', 'reactions')"
    echo "  3) Both"
    read -r -p "Enter choice [1/2/3, default based on branch: ${TARGET_ENV}]: " env_choice
    case "${env_choice}" in
        1)  TARGET_ENV="staging" ;;
        2)  TARGET_ENV="production" ;;
        3)  TARGET_ENV="both" ;;
        "") ;;   # accept branch-based default
        *)  echo "Invalid choice. Exiting."; exit 1 ;;
    esac

    echo ""
    echo "Select execution method:"
    echo "  1) Docker (containerized Solr — default)"
    echo "  2) Local  (against an already-running Solr at \${SOLR_URL})"
    read -r -p "Enter choice [1 or 2, default 1]: " exec_choice
    case "${exec_choice}" in
        2) EXEC_METHOD="local"  ;;
        *) EXEC_METHOD="docker" ;;
    esac
fi

EXEC_METHOD="${EXEC_METHOD:-docker}"

# --- Validate + normalize ----------------------------------------------

case "${TARGET_ENV}" in
    staging|production|both) ;;
    *)  echo "Invalid target: '${TARGET_ENV}'. Use staging, production, or both."
        exit 1 ;;
esac
case "${EXEC_METHOD}" in
    docker|local) ;;
    *)  echo "Invalid EXEC_METHOD: '${EXEC_METHOD}'. Use docker or local."
        exit 1 ;;
esac

# User-facing "production" → entrypoint token "prod".
declare -a ENV_TOKENS
case "${TARGET_ENV}" in
    staging)    ENV_TOKENS=("staging") ;;
    production) ENV_TOKENS=("prod") ;;
    both)       ENV_TOKENS=("staging" "prod") ;;
esac

# --- Summary + confirm -------------------------------------------------

echo ""
echo "========================================"
echo "Ready to deploy ModelSEED SOLR:"
echo "  Target:      ${TARGET_ENV}   (env tokens: ${ENV_TOKENS[*]})"
echo "  Branch:      ${GIT_BRANCH}"
echo "  Commit:      ${GIT_COMMIT}"
echo "  Date:        ${DEPLOY_DATE}"
echo "  Execution:   ${EXEC_METHOD}"
if [[ "${EXEC_METHOD}" == "local" ]]; then
    echo "  Solr URL:    ${SOLR_URL:-http://localhost:8983/solr}"
fi
echo "========================================"
echo ""

if [ -t 0 ] && [ "${AUTO_YES}" != "true" ]; then
    read -r -p "Proceed? [y/N]: " confirm
    if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
        echo "Aborted. No changes were made."
        exit 0
    fi
fi

# --- Compile step (always) ---------------------------------------------

# Produce both JSON flavours so Docker or Local can pick whichever is
# needed. Cheap (~15s each) and always safe — the JSONs get overwritten.
echo ""
echo ">>> Compiling biochemistry JSON (both flavours)..."
python3 compilation/Compile_Biochemistry_for_SOLR.py
python3 compilation/Compile_Biochemistry_for_SOLR_legacy.py

# --- Docker path -------------------------------------------------------

if [[ "${EXEC_METHOD}" == "docker" ]]; then
    export SOLR_ENVIRONMENTS="${ENV_TOKENS[*]}"

    echo ""
    echo ">>> docker compose up -d --build   (SOLR_ENVIRONMENTS='${SOLR_ENVIRONMENTS}')"
    docker compose up -d --build

    # POST_ON_START auto-posts the FIRST env only. For "both", post the
    # second explicitly. Idempotent — /update replaces existing docs.
    if [ "${#ENV_TOKENS[@]}" -gt 1 ]; then
        echo ""
        echo ">>> waiting for Solr to accept requests..."
        i=0
        until curl -fsS 'http://localhost:8983/solr/admin/info/system' > /dev/null 2>&1; do
            i=$((i + 1))
            [ "$i" -gt 60 ] && { echo "ERROR: Solr didn't become ready after 60s"; exit 1; }
            sleep 1
        done
        for env in "${ENV_TOKENS[@]:1}"; do
            echo ""
            echo ">>> posting biochemistry to '${env}' env..."
            docker compose exec -T solr /scripts/post_biochemistry.sh "${env}"
        done
    fi

    echo ""
    echo ">>> deploy complete. Core status:"
    curl -fsS 'http://localhost:8983/solr/admin/cores?action=STATUS&indexInfo=false' \
        | python3 -c "import json,sys; [print(f'    {n}') for n in sorted(json.load(sys.stdin)['status'])]" || true

# --- Local path --------------------------------------------------------

else
    SOLR_URL="${SOLR_URL:-http://localhost:8983/solr}"

    echo ""
    echo ">>> checking local Solr at ${SOLR_URL}..."
    if ! curl -fsS "${SOLR_URL}/admin/info/system" > /dev/null 2>&1; then
        echo "ERROR: cannot reach Solr at ${SOLR_URL}."
        echo "       Start your local Solr, or set SOLR_URL to point at the right host."
        exit 1
    fi

    # For each env, post to the appropriate core names. This assumes the
    # cores already exist on the local Solr — the docker path takes care
    # of creating cores from configsets, but on bare-metal Solr the
    # sysadmin owns core creation (see Solr/README.md's "Bare-metal Solr 9"
    # section for the one-time `solr create -c NAME -d NAME` incantation).
    for env in "${ENV_TOKENS[@]}"; do
        echo ""
        echo ">>> posting biochemistry to '${env}' env on local Solr..."
        SOLR_URL="${SOLR_URL}" BIOCHEMISTRY_JSON_DIR="${SCRIPT_DIR}/compilation" \
            "${SCRIPT_DIR}/post_biochemistry.sh" "${env}"
    done

    echo ""
    echo ">>> deploy complete."
fi
