# ModelSEED Atom Mappings

Reaction-level atom-atom mappings produced by the **UniversalRDT/ModelSEED**
pipeline (Sebastian Huhn), which wraps the [Reaction Decoder Tool
(RDT)](https://github.com/asad/ReactionDecoder) around every ModelSEED reaction.

Source repository: <https://github.com/sebahu/UniversalRDT/tree/main/ModelSEED>

## Files

| File | Rows | Reactions | Description |
|---|---:|---:|---|
| `all_mapping_no_problem.txt` | 1,067,343 | 23,913 | Clean set — RDT ran, no parse errors, no element mismatches. **This is what is currently ingested into `reaction_*.json`'s `atom_mapping` field.** |
| `all_mapping.txt` | 1,434,396 | 33,479 | Raw superset — includes reactions with parse errors or element mismatches (e.g., C mapped to N). Kept for review; not currently ingested. |
| `rxns_no_problems.txt` | 23,913 | 23,913 | Just the IDs of the clean reactions (equivalent to unique keys in `all_mapping_no_problem.txt`). |
| `rxns_with_cpds_without_structure.txt` | 18,621 | 18,621 | IDs of reactions RDT could not attempt because one or more of their compounds lacks a SMILES structure in `Unique_ModelSEED_Structures.txt`. |
| `compounds_without_structure.txt` | 12,318 | — | The compound IDs (cpdXXXXX) referenced by the reactions above but missing a structure. |

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

Each reaction present in `all_mapping_no_problem.txt` gets a new field:

```json
"atom_mapping": [
  "cpd00001:O#1=cpd00009:O#2",
  "cpd00012:O#1=cpd00009:O#1",
  ...
]
```

Reactions not in the clean set do not carry the field. The population script
lives at `Scripts/Structures/Populate_Atom_Mappings.py`.

## Coverage

| Set | Count | % of 56,012 total ModelSEED reactions |
|---|---:|---:|
| Clean mapping present in JSON | 23,913 | 43% |
| RDT ran but flagged as problematic (raw only) | 9,566 | 17% |
| Unmapable — compound(s) lack SMILES | 18,621 | 33% |
| Not attempted (in DB, not in Sebastian's runs) | ~4,000 | 7% |

### Priority scope (v7.0 ModelSEEDTemplates + PlantSEED_v3 Roles)

Of the 9,125 reactions used by the v7.0 templates and PlantSEED_v3 role
assignments (the union), 5,614 (61.5%) currently carry a clean atom mapping.
Breakdown of the 3,511 gap:

| Bucket | Count | Notes |
|---|---:|---|
| Compound(s) lack SMILES | 1,611 | blocked on structure curation |
| RDT ran but flagged | 1,369 | multi-target chains, duplicate mappings, element mismatches (O→S, C→N) — inherent RDT graph-alignment limits |
| Not attempted by pipeline | 531 | valid, mass-balanced reactions the UniversalRDT/ModelSEED wrapper silently skipped; 455 are non-obsolete + non-transport + status=OK |

For the Athaliana_TAIR10 reconstruction in plantseed-v3 specifically (782
unique base reaction IDs across 1,218 modelreactions): 571 mapped (73.0%),
46 blocked on SMILES, 126 RDT flagged, 39 pipeline-skipped.

## Provenance and regeneration

- Tool: Reaction Decoder Tool v4.0.0 (Rahman et al.)
- Wrapper: <https://github.com/sebahu/UniversalRDT/tree/main/ModelSEED>
- Runtime: ~3 days on an 8-core machine
- Per-reaction timeout: 30s
- Regeneration path: rerun the three scripts (`prepareRDT.sh` →
  `run_rdt.sh` → `unite_and_filter_mappings.sh`) in the source repo, then
  copy the four `.txt` outputs into this directory and rerun
  `Scripts/Structures/Populate_Atom_Mappings.py`.

## Future work

- Ingest raw+flagged set (add per-entry status) after user review of the
  ~11K problematic reactions.
- Fill the ~18K unmapable-because-no-SMILES gap by feeding placeholder
  structures where possible or documenting why (abstract compounds,
  polymers, etc.).
- Add a second mapper (rxnmapper — IBM RXN transformer) as a
  cross-verification source.
