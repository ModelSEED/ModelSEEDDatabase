# ModelSEED Atom Mappings

Reaction-level atom-atom mappings produced by the **UniversalRDT/ModelSEED**
pipeline (Sebastian Huhn), which wraps the [Reaction Decoder Tool
(RDT)](https://github.com/asad/ReactionDecoder) around every ModelSEED reaction.

Source repository: <https://github.com/sebahu/UniversalRDT/tree/main/ModelSEED>

The raw RDT output (`all_mapping.txt`) is filtered locally by
`Scripts/Structures/Rebuild_AtomMappings_from_raw.py` to produce
`all_mapping_no_problem.txt`; see the *Local filter* section below for the
divergence from upstream's `unite_and_filter_mappings.sh`.

## Files

| File | Rows | Reactions | Description |
|---|---:|---:|---|
| `all_mapping_no_problem.txt` | 1,329,095 | 32,877 | Clean set — filtered from `all_mapping.txt` by the local row-level filter. **This is what is ingested into `reaction_*.json`'s `atom_mapping` field.** |
| `all_mapping.txt` | 1,456,896 | 33,931 | Raw superset from RDT — kept for reproducibility and re-filtering. |
| `rxns_no_problems.txt` | 32,877 | 32,877 | Just the IDs of the clean reactions (equivalent to unique keys in `all_mapping_no_problem.txt`). |
| `rxns_confidence.tsv` | 32,878 | 32,877 | Per-reaction confidence tag (`clean` or `salvaged`) written by the local filter. `clean` = every raw RDT row was already a canonical single-pair same-element row; `salvaged` = at least one raw row was a run-on chain / dangling orphan / cross-element pair / malformed and this reaction's kept pairs are a strict subset of the raw output. Ingested into `reaction_*.json` as `atom_mapping_confidence`. |
| `rxns_with_cpds_without_structure.txt` | 18,621 | 18,621 | IDs of reactions RDT could not attempt because one or more of their compounds lacks a SMILES structure in `Unique_ModelSEED_Structures.txt`. |
| `compounds_without_structure.txt` | 12,318 | — | The compound IDs (cpdXXXXX) referenced by the reactions above but missing a structure. |
| `all_rxns_with_joker.txt` | 11,993 | 11,993 | IDs of reactions whose SMILES contains a wildcard (`*`) atom — fundamentally unmappable by RDT without picking a concrete placeholder atom. Ingested here for reference; excluded from clean set by construction. |

## Row format

Every row in `all_mapping*.txt` is one atom-pair:

```
rxnXXXXX cpdAAAAA:E#N=cpdBBBBB:E#M
```

- `rxnXXXXX` — reaction ID
- `cpdAAAAA:E#N` — atom `N` of element `E` in reactant `cpdAAAAA`, indexed
  per element per compound in InChI canonical order (so `C#1` is the first
  carbon, `C#2` the second, etc.)
- `cpdBBBBB:E#M` — the corresponding atom on the product side
- Hydrogen atoms are excluded from mappings
- Elements are matched (C to C, N to N, etc.); rows with element mismatches
  are filtered out of the clean set

Example, for `rxn00001` (H2O + PPi → 2 Pi):

```
rxn00001 cpd00001:O#1=cpd00009:O#2
rxn00001 cpd00012:O#1=cpd00009:O#1
rxn00001 cpd00012:P#1=cpd00009:P#1
...
```

Note: because the format uses compound identity (not per-molecule instance),
when a reaction produces multiple copies of the same compound (like 2 Pi
here), the reader must consult the reaction's stoichiometry to distribute
target atoms across molecule instances.

## Integration into reaction_*.json

Each reaction present in `all_mapping_no_problem.txt` gets two new fields:

```json
"atom_mapping": [
  "cpd00001:O#1=cpd00009:O#2",
  "cpd00012:O#1=cpd00009:O#1",
  ...
],
"atom_mapping_confidence": "clean"
```

`atom_mapping_confidence` is either `"clean"` or `"salvaged"`:

- **`clean`** — every raw RDT row for this reaction was already a
  canonical single-pair same-element row. Nothing was salvaged;
  nothing had to be dropped.
- **`salvaged`** — at least one raw row was a run-on chain, dangling
  orphan, cross-element pair, or malformed. The kept pairs are a strict
  subset of the raw RDT output for that reaction. Sebastian's
  observation: *"as soon as there is one problem, there are likely
  more (hidden) problems"* — so a `salvaged` mapping may carry
  correct-looking rows that are subtly wrong elsewhere. Reachability /
  neighborhood-level use cases are largely fine with `salvaged`;
  mechanism-level tracing (¹³C flux, exact atom fate) should filter
  to `clean` only.

Reactions not in the clean set do not carry either field. The
population script lives at `Scripts/Structures/Populate_Atom_Mappings.py`.

## Coverage

| Set | Count | % of 56,012 total ModelSEED reactions |
|---|---:|---:|
| Clean mapping present in JSON | 32,877 | 59% |
| Raw rows unrecoverable even at row level | 1,054 | 2% |
| Unmapable — compound(s) lack SMILES | 18,621 | 33% |
| Wildcard SMILES (`*` atom) — permanent | 1,725 | 3% |
| Not attempted by pipeline | 1,735 | 3% |

### Priority scope (v7.0 ModelSEEDTemplates ∪ PlantSEED_v3 Roles)

Of the 9,125 reactions used by the v7.0 templates and PlantSEED_v3 role
assignments (the union), 7,378 (80.9%) currently carry a clean atom mapping.
Breakdown of the 1,747 gap:

| Bucket | Count | Notes |
|---|---:|---|
| Compound(s) lack SMILES | 1,611 | blocked on structure curation |
| RDT ran but no salvageable rows | 57 | all rows malformed or element-mismatched — filter cannot rescue |
| Wildcard SMILES (`*` atom) | 42 | permanent unless placeholder atoms substituted |
| Not attempted by pipeline | 37 | RDT timed out even with extended per-reaction budget |

For the Athaliana_TAIR10 reconstruction in plantseed-v3 specifically (782
unique base reaction IDs across 1,218 modelreactions): 728 mapped (93.1%),
46 blocked on SMILES, 3 RDT-flagged, 1 wildcard, 4 pipeline-timeout.

### PlantSEED biomass reachability (carbon-atom trace from CO₂)

Undirected carbon-atom graph over PlantSEED reactions, seeded at
`cpd00011` (CO₂) + `cpd00242` (bicarbonate): **74 of 76** carbon-containing
biomass components reachable. The two remaining (Glucotropaeolin, Sinalbin)
are a PlantSEED curation gap — no producing reaction exists for the benzenic
glucosinolate branch — not an atom-mapping gap.

## Local filter

`Scripts/Structures/Rebuild_AtomMappings_from_raw.py` reads
`all_mapping.txt` and rewrites `all_mapping_no_problem.txt` under a
row-level filter that recovers pairs the upstream shell filter
throws away wholesale. Upstream's `unite_and_filter_mappings.sh` rejects
the *entire reaction* if any single row is malformed; the local filter
drops only the bad rows and keeps the reaction's valid pairs.

Recoveries the local filter makes over the shell filter:

- **Run-on chains** (`A=B=C=D`) — split into adjacent pairs, keep the
  same-element ones. Frequent in decarboxylations that collapse a
  carboxyl to CO₂: RDT emits `cpd00516:C#6=cpd00516:O#1=cpd00516:O#2=cpd00011:O#1`
  instead of two clean pair rows.
- **Two-letter elements** (`Cl`, `Fe`, `Mg`, `Zn`, `Hg`, `Br`, `Se`, …) —
  the shell filter's element slot is a single char; the local filter
  accepts 1–2 chars.
- **Dangling orphans** (`cpd00011:C#1` with no partner) — dropped at the
  row level rather than killing the reaction. The dangling row's
  intended partner is not inferred.

Element-pair whitelist is enforced by construction (only same-element
pairs are emitted); no separate whitelist file is needed.

The local rewrite raises clean coverage from ~24k to ~33k reactions and
closes six of eight biomass reachability gaps in PlantSEED
(Biotin, Leucine, Lysine, Phosphopantetheine, Thiamin diphosphate,
UDP-Xylose).

## Provenance and regeneration

- Tool: Reaction Decoder Tool v4.0.0 (Rahman et al.)
- Wrapper: <https://github.com/sebahu/UniversalRDT/tree/main/ModelSEED>
- Runtime: ~3 days on an 8-core machine
- Per-reaction timeout: 30s
- Regeneration path:
    1. Rerun `prepareRDT.sh` + `run_rdt.sh` + `unite_and_filter_mappings.sh`
       in the source repo (Sebastian's local copy on this server is at
       `/scratch/seaver/Claude_Projects/MSD_Structures/UniversalRDT/UniversalRDT/`
       with the RDT jar prebuilt at
       `../vendor/ReactionDecoder/target/rdt-4.0.0-jar-with-dependencies.jar`).
       This produces the raw `all_mapping.txt` and Sebastian's own filtered
       clean set.
    2. Copy the raw `all_mapping.txt` (plus any refreshed metadata files) into
       this directory. Ignore Sebastian's `all_mapping_no_problem.txt` and
       `rxns_no_problems.txt` — the local filter regenerates them.
    3. Rerun `Scripts/Structures/Rebuild_AtomMappings_from_raw.py` to
       produce `all_mapping_no_problem.txt` + `rxns_no_problems.txt` from
       the raw file.
    4. Rerun `Scripts/Structures/Populate_Atom_Mappings.py` to update the
       `atom_mapping` field on every `reaction_*.json`.

## Future work

- Fill the ~18K unmapable-because-no-SMILES gap by feeding placeholder
  structures where possible or documenting why (abstract compounds,
  polymers, etc.).
- Contribute the row-level filter back to
  `unite_and_filter_mappings.sh` upstream so downstream consumers benefit
  without needing to run the local rebuild step.
- Add a second mapper (rxnmapper — IBM RXN transformer) as a
  cross-verification source.
