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
│   ├── compounds/           # schema + solrconfig for the compounds core
│   │   ├── schema.xml
│   │   ├── schema_types.xml
│   │   └── solrconfig.xml
│   └── reactions/           # schema + solrconfig for the reactions core
│       ├── schema.xml
│       ├── schema_types.xml
│       └── solrconfig.xml
├── compilation/
│   ├── Compile_Biochemistry_for_SOLR.py    # produces the nested JSON payload
│   ├── solr_compounds.json                  # generated (gitignored)
│   └── solr_reactions.json                  # generated (gitignored)
├── examples/                # curl / python query samples (unchanged from pre-9)
└── data/                    # persistent Solr index (gitignored)
```

## Prerequisites

- Docker + `docker compose` (v2 syntax).
- Nothing else — Python 3, biochemistry data, and the compile step are
  all inside the image's build.

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

### Two schemas ship in the same container

The image bakes in **four configsets** and **two flavours of the
biochemistry JSON payload**:

| configset | schema shape | consumed by | JSON payload |
|---|---|---|---|
| `compounds`, `reactions` | new nested (per-source thermodynamics as child docs, `atom_mapping_data` + `atom_mapping_confidence` + `atom_mapping_has_symmetry_groups` fields, denormalized flat filters like `has_atom_mapping`, `sources_agree_direction`, `atom_count_C`, ...) | the staging UI, exercising the 2026-update fields | `solr_compounds.json` / `solr_reactions.json` |
| `compounds_legacy`, `reactions_legacy` | master-era flat (single `deltag`/`deltagerr` scalars, `stoichiometry` as a single semicolon-joined string, no `atom_mapping` / no `thermodynamics` dict) | the production UI (unchanged from what's currently deployed) | `solr_compounds_legacy.json` / `solr_reactions_legacy.json` |

Two Python scripts under `Solr/compilation/` produce the two payloads
from the same source `Biochemistry/*.json`. Both scripts run
automatically during `docker build` so the runtime image has all four
JSONs ready to post.

### Env-name → configset mapping

The entrypoint's `configset_for_env` helper picks the appropriate
configset per env name:

- `env == "prod"` → uses `compounds_legacy` / `reactions_legacy`
- `env == "staging"` (or anything else) → uses `compounds` / `reactions`

Set `SOLR_ENVIRONMENTS` to a space-separated list of env names — the
entrypoint creates one pair of cores per env:

| `SOLR_ENVIRONMENTS` | cores created | schema used per core |
|---|---|---|
| unset / empty (default) | `compounds`, `reactions` | new nested |
| `staging` | `compounds_staging`, `reactions_staging` | new nested |
| `prod` | `compounds_prod`, `reactions_prod` | legacy flat |
| `staging prod` | all four | new nested for `_staging` cores, legacy flat for `_prod` cores |

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

The staging + prod split means UI cutover happens on the UI side, not
the SOLR side:

1. Staging UI evolves against the `_staging` cores using the new
   schema. Production UI keeps querying the `_prod` cores under the
   legacy schema — no change visible to users.
2. When the new UI is ready to become production, the ModelSEED-UI
   deployment swaps its `NEXT_PUBLIC_SOLR_*_PRODUCTION` env vars to
   point at what were the `_staging` cores (or the UI is redeployed
   as the new production version pointing at a fresh set of cores).
3. The old `_prod` cores stay populated as an instant-rollback lane
   until the new UI is proven stable.

Both configsets stay in the image indefinitely; retiring the legacy
schema is a separate cleanup that happens after the new UI has been
in production for a while.

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
