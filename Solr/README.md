# ModelSEED Biochemistry in SOLR

Container-based Solr 9 deployment of the ModelSEED biochemistry database
with a classic-schema layout, nested per-source thermodynamics, and
expanded fields so the UI can render everything a compound or reaction
page needs from a single query.

## Layout

```
Solr/
├── Dockerfile               # Solr 9 image + baked-in configsets + entrypoint
├── docker-compose.yml       # dev/local convenience — mounts + POST_ON_START
├── entrypoint.sh            # start Solr → wait ready → create cores → optional post
├── post_biochemistry.sh     # curl-post the compiled JSON to both cores
├── cores/
│   ├── compounds/           # schema + solrconfig for the compounds core (new nested)
│   │   ├── schema.xml
│   │   ├── schema_types.xml
│   │   └── solrconfig.xml
│   ├── reactions/           # schema + solrconfig for the reactions core (new nested)
│   │   ├── schema.xml
│   │   ├── schema_types.xml
│   │   └── solrconfig.xml
│   ├── compounds_legacy/    # master-era flat schema for production compounds
│   ├── reactions_legacy/    # master-era flat schema for production reactions
│   └── structures/          # smiles / inchi / inchikey / pre-rendered SVG per cpd_id
│       ├── schema.xml
│       ├── schema_types.xml
│       └── solrconfig.xml
├── compilation/
│   ├── Compile_Biochemistry_for_SOLR.py         # new nested compound+reaction JSON
│   ├── Compile_Biochemistry_for_SOLR_legacy.py  # legacy flat compound+reaction JSON
│   ├── Compile_Structures_for_SOLR.py           # structures JSON (RDKit-rendered SVG)
│   ├── solr_compounds.json  / solr_reactions.json          # generated (gitignored)
│   ├── solr_compounds_legacy.json / solr_reactions_legacy.json  # generated (gitignored)
│   └── solr_structures.json                                # generated (gitignored)
├── examples/                # curl / python query samples (unchanged from pre-9)
└── data/                    # persistent Solr index (gitignored)
```

## Prerequisites

- Docker + `docker compose` (v2 syntax).
- Nothing else — Python 3, biochemistry data, and the compile step are
  all inside the image's build.

## Deploy script (wrapper)

`deploy_container.sh` is an interactive wrapper modeled on
ModelSEED-UI's `deploy_container.sh`. It handles both the compile and
the deploy in one command, and prompts (or takes args) for:

- **Which env(s) to deploy** — staging, production, or both
- **How to run** — Docker (containerized Solr) or Local (against an
  already-installed Solr you point at via `SOLR_URL`)

Examples:

```bash
./deploy_container.sh                        # interactive prompts
./deploy_container.sh staging                # deploy staging via Docker
./deploy_container.sh production --yes       # deploy production, skip confirm
./deploy_container.sh both --yes             # both envs, no prompts
EXEC_METHOD=local ./deploy_container.sh both # against a local Solr install
```

Auto-detects the branch — main/master/production → production default,
anything else → staging default (matches the UI script's convention).
Runs the two compile scripts unconditionally, then the appropriate
Docker or Local post steps. See below for the underlying bring-up
sequence the script wraps.

## First-time bring-up (testing mode)

```bash
# From the Solr/ directory:

# 1. Build + start the container. The build:
#    - runs Compile_Biochemistry_for_SOLR.py against ../Biochemistry/
#    - bakes the resulting nested JSON into the image
#    - creates the compounds + reactions cores on first start
#    - posts the biochemistry (~30s)
docker compose up -d --build

# 2. Verify.
curl 'http://localhost:8983/solr/reactions/select?q=id:rxn00001+AND+doc_type:reaction&fl=*,[child limit=100]&wt=json&indent=true'
```

The container listens on `http://localhost:8983`. Admin UI is at
`http://localhost:8983/solr/`.

Restarts of the same container are fast: the entrypoint checks whether
the persistent `./data` volume already has documents and skips the post
when it does. To force a re-post, run
`docker compose exec solr /scripts/post_biochemistry.sh`.

## Refreshing the biochemistry (testing mode)

Whenever `../Biochemistry/*.json` changes, rebuild the image (the compile
step re-runs as part of the build):

```bash
docker compose down
rm -rf data           # clear the persistent volume so post runs again
docker compose up -d --build
```

## Persisting the index on a host path (`SOLR_DATA_DIR`)

By default the Solr index lives in the named volume
`modelseed_solr_data`, so Docker owns the directory and its permissions.
To put the index somewhere durable and inspectable instead — an
IT-provided mount, a project volume — set `SOLR_DATA_DIR` to an absolute
host path:

```bash
# Solr/.env  (gitignored, per-deployment)
SOLR_DATA_DIR=/vol/model-biochem/solr-data
```

Solr runs as **UID 8983** inside the container and will not be able to
write a directory created by your host user. Chown it once before the
first bring-up. If you don't have root on the host but are in the
`docker` group, a throwaway container does the job:

```bash
mkdir -p /vol/model-biochem/solr-data
docker run --rm -v /vol/model-biochem/solr-data:/d alpine \
    chown -R 8983:8983 /d
```

Leaving `SOLR_DATA_DIR` unset keeps the named-volume behaviour, so local
development is unaffected.

## Production mode — external biochemistry mount

The baked-in data pattern is a testing convenience while IT provides a
persistent shared mount for the biochemistry data. To swap to the
production pattern:

1. Remove the two `COPY --from=compile` lines from `Dockerfile` (drops
   the ~227 MB of JSON out of the image).
2. In `docker-compose.yml`, add read-only mounts for the source data and
   a writable mount for the compiled JSON:
   ```yaml
   volumes:
     - /shared/modelseed/Biochemistry:/data/biochemistry:ro
     - /shared/modelseed/compiled:/data/compilation
   ```
3. Run the compile from outside the container against the shared source:
   `python3 compilation/Compile_Biochemistry_for_SOLR.py`.
4. `docker compose exec solr /scripts/post_biochemistry.sh` to load.

The compile script writes to `Solr/compilation/solr_*.json` by default;
edit its `OUT_COMPOUNDS` / `OUT_REACTIONS` to point at the shared mount
instead if you're running the compile from a different working directory.

## Schema highlights

Both cores use nested child documents so the JSON layout of the source
records is preserved in the search index.

**Compounds** (`cpdXXXXX` parent doc, `doc_type=compound`) carry:

- Flat parent fields — `id`, `name`, `abbreviation`, `formula`, `mass`,
  `charge`, `smiles`, `inchikey`, `aliases[]`, `notes[]`, `is_core`,
  `is_obsolete`, `is_cofactor`, `class`, `linked_compound`, ...
- **Denormalized filters** — `has_structure`, `atom_count_C`,
  `atom_count_N`, `atom_count_O`, `atom_count_P`, `atom_count_S`,
  `n_sources_thermodynamics`, `is_electron_carrier`
- **Child docs (`doc_type=thermodynamics`)** — one per source
  (`Group contribution`, `eQuilibrator`, …) with `source_name`, `energy`,
  `error`
- **Child docs (`doc_type=pkas`)** — one per source (`Marvin`, `MolGpKa`,
  …) with `source_name`, `pka_value[]`, `pkb_value[]`
- Legacy flat scalars `deltag`, `deltagerr`, `pka`, `pkb` still populated
  for backwards compatibility

**Reactions** (`rxnXXXXX` parent doc, `doc_type=reaction`) carry:

- Flat parent fields — `id`, `name`, `abbreviation`, `code`, `equation`,
  `definition`, `status`, `aliases[]`, `pathways[]`, `ec_numbers[]`,
  `notes[]`, `is_transport`, `is_obsolete`, `reversibility`, ...
- **Atom mapping** (flattened from the JSON's nested `atom_mapping` dict):
  `atom_mapping_data[]` (multi-valued `cpdA:E#N=cpdB:E#M` strings;
  endpoints may be single atoms or `(x;y;…)` set notation),
  `atom_mapping_confidence` (`clean` | `salvaged`), and
  `atom_mapping_has_symmetry_groups` (`true` when any row uses set
  notation)
- **Denormalized filters** — `has_atom_mapping`,
  `n_sources_thermodynamics`, `sources_agree_direction`, `n_reactants`,
  `n_products`
- **Child docs (`doc_type=thermodynamics`)** — one per source with
  `source_name`, `energy`, `error`, `operator` (`>`, `<`, `=`, `?`)
- **Child docs (`doc_type=stoichiometry`)** — one per participant with
  `compound`, `coefficient`, `compartment`, `participant_charge`,
  `participant_formula`, `participant_name`, `is_reactant`
- Legacy flat scalars `deltag`, `deltagerr`, `reversibility`, flat
  `compound_ids[]` still populated for backwards compatibility

## Querying

### Detail-view: one record with all children

```
GET /solr/reactions/select?q=id:rxn00001+AND+doc_type:reaction&fl=*,[child limit=100]&wt=json
GET /solr/compounds/select?q=id:cpd00002+AND+doc_type:compound&fl=*,[child limit=100]&wt=json
```

Returns the parent plus each nested child under the field name of its
subdocument path (`thermodynamics`, `stoichiometry`, `pkas`), so the UI
can render a full compound or reaction detail page — flat fields,
per-source thermodynamics table, atom-mapping list, stoichiometry table
— from a single request. The `doc_type:...` clause in `q` excludes the
matching child docs (whose IDs start with the parent's, e.g.
`rxn00001::thermo::GC`); the `[child limit=100]` transformer defaults
to 0 children if you omit `limit=`, so it's worth being explicit.

### Filter parents by child conditions

Block-join query — reactions where eQuilibrator's operator is `>` and GC's
operator is `=` (disagreement candidates):

```
q={!parent which="doc_type:reaction"}
   doc_type:thermodynamics AND source_name:eQuilibrator AND operator:">"
```

### Fast flat filters (no join)

Denormalized fields index directly, so common browse-page filters skip
the join entirely:

```
# only atom-mapped reactions, only if all sources agree on direction
fq=has_atom_mapping:true&fq=sources_agree_direction:true

# small central-metabolism compounds only
fq=atom_count_C:[2 TO 6]&fq=has_structure:true

# non-obsolete compounds with at least 2 thermodynamics sources
fq=is_obsolete:false&fq=n_sources_thermodynamics:[2 TO *]
```

## Multi-environment deployment (staging + production)

The single-container layout supports **staging and production data
side by side in the same Solr instance, each using its own schema**.
This is the ModelSEED-UI's expected deployment model (the UI's
`NEXT_PUBLIC_SOLR_*_STAGING` and `_PRODUCTION` env vars each point at
their own core name on the same host) and — critically — lets the
production UI keep working during a staging release without any
UI-side coordination on when to switch fields.

### Three cores per env (five configsets total)

The image bakes in **five configsets** and **three JSON payloads**:

| configset | schema shape | consumed by | JSON payload |
|---|---|---|---|
| `compounds`, `reactions` | new nested (per-source thermodynamics as child docs, `atom_mapping_data` + `atom_mapping_confidence` + `atom_mapping_has_symmetry_groups` fields, denormalized flat filters like `has_atom_mapping`, `sources_agree_direction`, `atom_count_C`, ...) | the staging UI, exercising the 2026-update fields | `solr_compounds.json` / `solr_reactions.json` |
| `compounds_legacy`, `reactions_legacy` | master-era flat (single `deltag`/`deltagerr` scalars, `stoichiometry` as a single semicolon-joined string, no `atom_mapping` / no `thermodynamics` dict) | the production UI (unchanged from what's currently deployed) | `solr_compounds_legacy.json` / `solr_reactions_legacy.json` |
| `structures` | one doc per `compound_id` (all ~45.7K compounds), fields `smiles` / `inchi` / `inchikey` / `svg` are all optional (nulls = absent, not `has_*` booleans) | staging + production UI (same schema everywhere — nothing to be backward-compatible with) | `solr_structures.json` (single payload, shared across envs) |

Three Python scripts under `Solr/compilation/` produce the three
payloads from `Biochemistry/*.json` (+ `Biochemistry/Structures/
Unique_ModelSEED_Structures.txt` for the structures core). All three
run automatically during `docker build` so the runtime image has all
five JSONs ready to post.

The structures compile is the slowest (~5 min for ~30K RDKit
`MolDraw2DSVG` renders at 300 x 300, capped at 100 KB per SVG). It
skips wildcard SMILES (`*` atoms render as literal `*`, misleading)
and lets the UI regenerate SVG on the fly for the skipped compounds.

### Env-name → configset and core-name mapping

The entrypoint has two per-env helpers (`configset_for_env` and
`suffix_for_env`) that together decide, for each env, which configset
its cores use AND what the cores are named:

- `env == "prod"` → **bare core names** (`compounds`, `reactions`,
  matching the URL the production UI hits today) + `*_legacy`
  configsets (matching the schema the production UI expects)
- `env == "staging"` (or any other named env) → env-suffixed core
  names (`compounds_staging`, `reactions_staging`) + new nested
  configsets
- `env` unset (standalone dev) → bare core names + new nested
  configsets

Set `SOLR_ENVIRONMENTS` to a space-separated list of env names — the
entrypoint creates one **triple** of cores per env (compounds +
reactions + structures):

| `SOLR_ENVIRONMENTS` | cores created | configset used per core |
|---|---|---|
| unset / empty (dev) | `compounds`, `reactions`, `structures` | new nested (`structures` for the structures core) |
| `staging` | `compounds_staging`, `reactions_staging`, `structures_staging` | new nested (`structures` for the structures core) |
| `prod` | `compounds`, `reactions`, `structures` (bare) | legacy flat + `structures` |
| `staging prod` | all six above cores | new nested for `_staging`, legacy flat for the bare `compounds`/`reactions`, shared `structures` for both |

The structures core uses the SAME `structures` configset in every env
— there's no legacy variant since the production UI doesn't have this
core today. Only the core name gets the env suffix; the schema and
payload are identical everywhere.

The bare-name-for-prod choice means production URLs stay stable:
`/solr/compounds/select?q=...` continues to hit the production data
regardless of when we cut over the staging schema underneath.

Populate a specific env (auto-picks the matching JSON payload):

```bash
docker compose exec modelseed-solr /scripts/post_biochemistry.sh staging
docker compose exec modelseed-solr /scripts/post_biochemistry.sh prod
```

`POST_ON_START=true` auto-posts into the FIRST env in the list;
subsequent envs need an explicit post. The idempotent-post check on
startup only looks at the first-env cores, so bringing the container
back up after populating both envs is fast.

### Cutover choreography

The staging + prod split means UI cutover happens on the SOLR side,
not the UI side. The production UI's query URL (`/solr/compounds`,
`/solr/reactions`) doesn't need to change — what changes is which
configset + payload sits behind those bare core names.

1. Staging UI evolves against the `_staging` cores using the new
   schema. Production UI keeps querying `/solr/compounds` and
   `/solr/reactions` (bare) which are populated from the legacy
   payload under the legacy schema — no change visible to production
   users.
2. When the new UI is ready to become production, the current
   bare-name cores get replaced with the new-schema configset +
   payload. Two options for the swap:
   - **In-place**: delete the bare cores, recreate them under the
     new configset, re-post the new-format JSON. Brief interruption
     during the recreate + post window.
   - **Blue/green**: leave bare cores untouched, promote the
     `_staging` cores to bare names via configset swap +
     `CORE RENAME`. Zero-interruption; needs a small entrypoint
     script change to skip re-creation.
3. The old legacy JSON payload stays in the image as an
   instant-rollback lane until the new UI has been in production for
   a while.

Both configset flavours + both JSON payloads stay in the image
indefinitely; retiring the legacy schema is a separate cleanup that
happens after the new UI has been in production for a while.

## Alternate deployments

### Just the container, no compose

```bash
docker build -t modelseed-solr Solr/
docker run -d --name modelseed-solr \
    -p 8983:8983 \
    -v $(pwd)/Solr/data:/var/solr \
    -v $(pwd)/Solr/compilation:/data/compilation \
    modelseed-solr
docker exec modelseed-solr /scripts/post_biochemistry.sh
```

### Bare-metal Solr 9 (no container)

The `Solr/cores/{compounds,reactions}/` directories are drop-in configsets
for a standard Solr 9 install. Copy the two directories into
`$SOLR_HOME/configsets/`, then create the cores:

```bash
solr create -c compounds -d compounds
solr create -c reactions -d reactions
```

## Provenance

- Solr version: **9.8** (base image `solr:9.8`). Bump the tag in
  `Dockerfile` when upgrading; also update `luceneMatchVersion` in
  each core's `solrconfig.xml` to match.
- Migration from Solr 7.4.0: the schemas were reworked to move all
  `Trie*Field` types to `Point*Field` (required by Solr 9), and both
  cores gained nested-document plumbing (`_root_`, `_nest_path_`,
  `_nest_parent_`, `doc_type`) plus the denormalized flat fields
  documented above.
