# ModelSEED Biochemistry SOLR — Solr 9.x, classic schema + nested docs.
#
# Multi-stage build:
#   Stage 1 (compile) — small Python image runs Compile_Biochemistry_for_SOLR.py
#     against the repo's Biochemistry/ tree and produces the nested-child JSON
#     payload for both cores.
#   Stage 2 (runtime) — Solr 9.8 image + baked-in configsets + baked-in
#     compiled JSON + entrypoint that creates cores and (by default) posts
#     the biochemistry on first start.
#
# TESTING NOTE: for the initial poplar bring-up the biochemistry JSON is
# baked into the image, so `docker compose up` is fully self-contained
# with no external mounts required. Once IT provides a persistent mount
# location for the biochemistry data, the compile stage can be swapped
# for a runtime volume mount + post_biochemistry.sh invocation.
#
# Build context must be the repository root (see docker-compose.yml):
#   docker build -f Solr/Dockerfile -t modelseed-solr .

# -----------------------------------------------------------------------------
# Stage 1 — compile biochemistry JSON into the nested-child Solr payload
# -----------------------------------------------------------------------------
# Non-slim python: rdkit's rdMolDraw2D transitively needs a handful of
# X11 / expat shared libs even for headless SVG output. The slim image
# is missing them (libXext, libXrender, libexpat, ...). Since this is a
# multi-stage build only the compiled JSON copies into the runtime
# image, so using the full python image here has zero effect on the
# final image size.
FROM python:3.12 AS compile

WORKDIR /build

# RDKit for 2D-structure SVG rendering in the structures compile.
RUN pip install --no-cache-dir rdkit

# Mirror the repo layout so the compile scripts' relative paths
# (../../Biochemistry/*.json, ../../Biochemistry/Structures/) resolve
# correctly. Only what the compile scripts actually read is copied —
# the rest of Biochemistry/ (TSV, provenance dumps, thermodynamics
# tables) totals ~3.5 GB and isn't touched.
COPY Biochemistry/compound_*.json Biochemistry/reaction_*.json Biochemistry/
COPY Biochemistry/Structures/Unique_ModelSEED_Structures.txt Biochemistry/Structures/
COPY Solr/compilation/Compile_Biochemistry_for_SOLR.py \
     Solr/compilation/Compile_Biochemistry_for_SOLR_legacy.py \
     Solr/compilation/Compile_Structures_for_SOLR.py \
     Solr/compilation/

# Three compile passes → five JSON files:
#   solr_compounds.json / solr_reactions.json                   (new nested schema, staging)
#   solr_compounds_legacy.json / solr_reactions_legacy.json     (master-era flat, production)
#   solr_structures.json                                        (per-compound smiles/inchi/inchikey + SVG)
# All get baked into the runtime image so the container can post
# any/all without re-running the compile. The structures pass is
# the slowest (~5 min for 30K RDKit renders); build cache reuses it
# when only the biochemistry data (not the SMILES) changes.
RUN cd Solr/compilation && python3 Compile_Biochemistry_for_SOLR.py \
                        && python3 Compile_Biochemistry_for_SOLR_legacy.py \
                        && python3 Compile_Structures_for_SOLR.py

# -----------------------------------------------------------------------------
# Stage 2 — Solr runtime
# -----------------------------------------------------------------------------
FROM solr:9.8

USER root

# Baked-in configsets. Persistent /var/solr volume will pick these up
# when the entrypoint creates cores at startup. Four configsets so a
# single Solr instance can serve staging (new nested schema) and
# production (legacy flat schema) side by side.
COPY --chown=solr:solr Solr/cores/compounds \
     /opt/solr/server/solr/configsets/compounds/conf
COPY --chown=solr:solr Solr/cores/reactions \
     /opt/solr/server/solr/configsets/reactions/conf
COPY --chown=solr:solr Solr/cores/compounds_legacy \
     /opt/solr/server/solr/configsets/compounds_legacy/conf
COPY --chown=solr:solr Solr/cores/reactions_legacy \
     /opt/solr/server/solr/configsets/reactions_legacy/conf
COPY --chown=solr:solr Solr/cores/structures \
     /opt/solr/server/solr/configsets/structures/conf

# Baked-in compiled biochemistry JSON (from the compile stage above).
# Two flavours: `_legacy` files are shaped for the master-era schema
# (deployed today in production); the unsuffixed files are shaped for
# the new nested schema (staging UI).
COPY --from=compile --chown=solr:solr /build/Solr/compilation/solr_compounds.json \
     /data/compilation/solr_compounds.json
COPY --from=compile --chown=solr:solr /build/Solr/compilation/solr_reactions.json \
     /data/compilation/solr_reactions.json
COPY --from=compile --chown=solr:solr /build/Solr/compilation/solr_compounds_legacy.json \
     /data/compilation/solr_compounds_legacy.json
COPY --from=compile --chown=solr:solr /build/Solr/compilation/solr_reactions_legacy.json \
     /data/compilation/solr_reactions_legacy.json
COPY --from=compile --chown=solr:solr /build/Solr/compilation/solr_structures.json \
     /data/compilation/solr_structures.json

# Entrypoint + helper scripts
COPY --chown=solr:solr Solr/entrypoint.sh        /scripts/entrypoint.sh
COPY --chown=solr:solr Solr/post_biochemistry.sh /scripts/post_biochemistry.sh
RUN chmod +x /scripts/entrypoint.sh /scripts/post_biochemistry.sh

USER solr

# Post on first start by default (the biochemistry is baked into the image
# for the poplar test bring-up; flip to false once IT provides a mount and
# post is triggered explicitly).
ENV POST_ON_START=true

# Where the entrypoint looks for the compiled JSON. Baked location by
# default; override to a mount path once IT provides one.
ENV BIOCHEMISTRY_JSON_DIR=/data/compilation

EXPOSE 8983

HEALTHCHECK --interval=30s --timeout=5s --start-period=60s \
  CMD curl -fsS 'http://localhost:8983/solr/admin/cores?action=STATUS' > /dev/null || exit 1

ENTRYPOINT ["/scripts/entrypoint.sh"]
