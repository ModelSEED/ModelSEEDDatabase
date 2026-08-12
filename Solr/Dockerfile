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
FROM python:3.12-slim AS compile

WORKDIR /build

# Mirror the repo layout so the compile script's relative paths
# (../../Biochemistry/*.json) resolve correctly. Only the compound + reaction
# JSON files are needed — the rest of Biochemistry/ (TSV, provenance dumps,
# structure files, thermodynamics tables) totals ~3.5 GB and isn't read by
# the compile step.
COPY Biochemistry/compound_*.json Biochemistry/reaction_*.json Biochemistry/
COPY Solr/compilation/Compile_Biochemistry_for_SOLR.py Solr/compilation/

RUN cd Solr/compilation && python3 Compile_Biochemistry_for_SOLR.py

# -----------------------------------------------------------------------------
# Stage 2 — Solr runtime
# -----------------------------------------------------------------------------
FROM solr:9.8

USER root

# Baked-in configsets. Persistent /var/solr volume will pick these up
# when the entrypoint creates cores at startup.
COPY --chown=solr:solr Solr/cores/compounds \
     /opt/solr/server/solr/configsets/compounds/conf
COPY --chown=solr:solr Solr/cores/reactions \
     /opt/solr/server/solr/configsets/reactions/conf

# Baked-in compiled biochemistry JSON (from the compile stage above).
COPY --from=compile --chown=solr:solr /build/Solr/compilation/solr_compounds.json \
     /data/compilation/solr_compounds.json
COPY --from=compile --chown=solr:solr /build/Solr/compilation/solr_reactions.json \
     /data/compilation/solr_reactions.json

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
