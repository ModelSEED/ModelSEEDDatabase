# Modular Updates

This directory holds **declarative update manifests** — small YAML
files that describe a single curated change to the database. Each
manifest is a self-contained record of *what* changed, *why*, *who*
made the change, and *which downstream stages must be re-run*. Manifests
are the recommended way to make non-trivial edits to the structure
data, curation overrides, ignore lists, alias tables, and per-source
pKa/protonation files.

The motivation is to replace ad-hoc multi-file edits with one reviewable
artifact. A manifest:

- gives reviewers a single PR diff that explains the curation rationale
  alongside the data change;
- carries enough metadata for `Apply_Manifest.py` to make the right
  edits in the right files;
- lets `Refresh_Pipeline.py` figure out which downstream stages need
  re-running (avoiding "did you remember to run List, then Update,
  then Reprint?");
- preserves audit trail in `applied/<date>-<topic>.yaml` once
  applied, so curation history is queryable without `git log`
  archaeology.

## Workflow

```
1. Curator writes <manifest>.yaml in this directory                            
2. Curator runs:                                                                  
    python Scripts/Updates/Apply_Manifest.py <manifest>.yaml [--dry-run]      
3. Apply_Manifest reads the manifest, dispatches to the right handler,        
   edits the right per-source/Curation files, and prints the list of         
   pipeline stages that need to re-run                                       
4. Curator runs:                                                              
    python Scripts/Updates/Refresh_Pipeline.py --since=<manifest>.yaml       
   to cascade the edit through Print -> List -> Update_Compound_*           
   -> Reprint -> provenance -> FAISS-validator                              
5. On success, Apply_Manifest moves the manifest into applied/ with a       
   YYYY-MM-DD prefix; the curator commits the lot as one PR.                
```

Use `--dry-run` to preview which files would change without writing
anything.

## Manifest types (currently handled)

Each manifest must declare its `type:` field, which dispatches to a
handler in `Scripts/Updates/handlers/`. The minimum-viable handlers
covered by the current `Apply_Manifest.py`:

| `type:` | What it changes | Downstream cascade |
|---|---|---|
| `structure_update` | a row in `Biochemistry/Structures/<source>/{inchi,smiles,inchikey}.tsv` or in `protonations/*.tsv` | Print → List → Update_Cpd_Structures_Formulas_Charge → Reprint → Build_Provenance → FAISS-Tier1 |
| `protonation_replace` | replaces or adds a `protonations/<tool>_<ver>_ph<n>.tsv` file | List → Update_Cpd_Structures_Formulas_Charge → Reprint → Build_Provenance → FAISS-Tier1 |
| `override_add` | appends to `Biochemistry/Curation/overrides/acps_formula_charge.tsv` | Update_Cpd_Structures_Formulas_Charge → Reprint → Build_Provenance |
| `ignore_add` | appends to `Biochemistry/Curation/ignores/<file>.txt` | List → Update_Cpd_Structures_Formulas_Charge → Reprint → Build_Provenance |
| `alias_add` / `alias_remove` | edits `Biochemistry/Aliases/Unique_ModelSEED_Compound_Aliases.txt` | Update_Cpd_Aliases → Reprint → Build_Provenance |
| `pka_replace` | replaces `Biochemistry/Structures/<source>/pkas/<file>.tsv` | Update_Cpd_pKas → Reprint → Build_Provenance |

The cascade stages live in `Scripts/Updates/cascade.py` as a small
DAG; adding a new manifest type is one handler module plus one
cascade entry.

## Manifest schema

Every manifest has the same envelope:

```yaml
type:    <one of the types above>
title:   <short human-readable title>
author:  <github handle or email>
date:    YYYY-MM-DD
reason:  |
  Free text explaining why this change is needed. Cite data sources,
  tickets, papers — whatever a future reader will need to understand
  the rationale.
target:  <type-specific target descriptor>
change:  <type-specific change descriptor>
expected_effects:                # optional self-validation
  affected_compounds: <int>      # rough count; Apply_Manifest warns
                                 # if reality deviates by >10%
  fields_changed:                # which compound TSV columns
    - formula
    - charge
```

The `target:` and `change:` fields are type-specific. See
`examples/` for one of each:

- `examples/structure_update.yaml` — fix a single compound's KEGG InChI
- `examples/protonation_replace.yaml` — drop in a new Marvin 24.0 protonation bundle
- `examples/override_add.yaml` — add an ACP entry
- `examples/ignore_add.yaml` — flag a structure as ignored
- `examples/alias_add.yaml` — add a cross-reference

## Why this matters

Without manifests, an edit like "change cpd05123's protonated SMILES,
because Marvin 24 corrected a tautomerism error" requires the curator
to:

1. Find the right row in `MetaCyc/protonations/marvin_23.4_ph7.tsv`
2. Edit the SMILES, formula, charge columns by hand
3. Remember to re-run `Print_Structure_Formula_Charge.py`
4. Remember to re-run `List_ModelSEED_Structures.py`
5. Remember to re-run `Update_Compound_Structures_Formulas_Charge.py`
6. Remember to re-run `Reprint_Biochemistry.py`
7. Remember to re-run `Build_Compound_Field_Provenance.py`
8. Remember to re-run `Validate_FAISS_Outputs.py`
9. Hope the PR reviewer can reconstruct the rationale from the diff

With manifests, the curator writes a 15-line YAML and runs two
commands. The PR diff includes the manifest itself, so the reviewer
sees the rationale and the data change side-by-side. The cascade
guarantees nothing is forgotten.

## Adding a new manifest type

1. Pick a name (`type: my_new_thing`).
2. Add a handler in `Scripts/Updates/handlers/my_new_thing.py`
   exposing `apply(manifest, dry_run=False) -> list[str]` that
   returns the stage names that need to cascade.
3. Add the type → cascade mapping in `Scripts/Updates/cascade.py`.
4. Add an example manifest in `examples/my_new_thing.yaml`.
5. Update the table above.

That's the entire integration surface — no other code changes.
