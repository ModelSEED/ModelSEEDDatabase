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
- **Atom mapping** — `atom_mapping[]` (multi-valued `cpdA:E#N=cpdB:E#M`
  strings) and `atom_mapping_confidence` (`clean` | `salvaged`)
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
side by side in the same Solr instance**, which is the ModelSEED-UI's
expected deployment model (the UI's `NEXT_PUBLIC_SOLR_*_STAGING` and
`_PRODUCTION` env vars each point at their own core name on the same
host).

Set `SOLR_ENVIRONMENTS` to a space-separated list of environment names.
For each env, the entrypoint creates a suffixed pair of cores from
the same baked-in configsets:

| `SOLR_ENVIRONMENTS` | cores created |
|---|---|
| unset / empty (default) | `compounds`, `reactions` |
| `staging` | `compounds_staging`, `reactions_staging` |
| `staging prod` | all four: `compounds_staging`, `reactions_staging`, `compounds_prod`, `reactions_prod` |

Populate a specific env:

```bash
docker compose exec solr /scripts/post_biochemistry.sh staging
docker compose exec solr /scripts/post_biochemistry.sh prod
```

`POST_ON_START=true` auto-posts into the FIRST env in the list;
subsequent envs need an explicit post. The idempotent-post check on
startup only looks at the first-env cores, so bringing the container
back up after populating both envs is fast.

Cutover from staging → prod is a re-post to the prod cores; no
rebuild, no downtime, no core deletion. Roll back by re-posting the
previous compiled JSON.

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
