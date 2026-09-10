# Marvin 26.1 pKa regeneration

Generated 2026-09-10 by `Scripts/Structures/Run_Marvin_pKas.py`.

Regenerates the per-source Marvin pKa bundles against the current
`<source>/inchi.tsv` and `<source>/smiles.tsv`, replacing the retired `pKaMol.java` path (which needed
MarvinBeans on the classpath) with the `cxcalc` CLI from Marvin Desktop Suite.

## Invocation

```
cxcalc --ignore-error -i x pka --na 20 --nb 20 <structure-list>
```

Calculator 26.1.2 (ChemAxon), one process per input representation per source,
batch mode. `--na/--nb 20` because the widest existing 23.4 row carries 19
tokens. Licensed under the Protonation Plugin Group.

## Two engines, one Marvin

InChI structures go through the **cxcalc CLI**. The SMILES-only remainder --
polymers carrying `*` attachment points, plus organometallics -- goes through the
**Java `chemaxon.calculations.PkaPlugin`** via JPype, because cxcalc refuses those:
*"pka: Calculation result is not defined for query molecules"*. A JRE is enough;
nothing is compiled. This is the route `pKaMol.java` took before the API was
renamed in 26.1 and its committed `.class` stopped running.

They are one engine, verified rather than assumed. Over 300 ChEBI compounds the
plugin reproduces cxcalc to a median |delta| of **0.0000**, max **0.000**, 100%
within 0.01, with 1 site-count difference in 567 sets and identical atom ordering
in 562. Rows that came from the plugin are exactly the ids absent from
`inchi.tsv`.

The plugin does not invent sites on the wildcards: across the Rhea polymers, 0 of
154 predicted sites sat on an atom bonded to a `*`, and site counts match 23.4 in
40 of 40 comparable sets.

## Coverage

| source | InChI (cxcalc) | SMILES-only (plugin) | recovered | total with pKa | rows | 23.4 |
|---|---:|---:|---:|---:|---:|---:|
| ChEBI | 9473 | 2072 | 2038 | 11177 | 21010 | 11186 |
| KEGG | 15319 | 957 | 931 | 15530 | 28079 | 15967 |
| MetaCyc | 19298 | 5775 | 5644 | 24249 | 45162 | 24247 |
| Rhea | 207 | 30 | 30 | 235 | 444 | 235 |
| **total** | **44297** | **8834** | **8643** | **51191** | **94695** | **51635** |

**8643 compounds recovered** through the plugin path. The shortfall against
23.4 drops from 9,121 ids to **485**: Rhea now matches 23.4 exactly,
MetaCyc exceeds it, and 448 of the remaining 485 are KEGG's
SRU exclusions -- absent from both structure files and listed in
`KEGG/KEGG_SRU_041020.txt`.

Zero InChI structures failed to parse. The plugin declined 191 SMILES across all
sources, chiefly dative-bond notation its parser rejects.

## 23.4 vs 26.1

| source | 23.4 ids | 26.1 ids | 23.4-only | site-count differs | median abs delta | within 0.5 |
|---|---:|---:|---:|---:|---:|---:|
| ChEBI | 11186 | 11177 | 14 | 15.1% | 0.000 | 75.2% |
| KEGG | 15967 | 15530 | 448 | 12.6% | 0.000 | 76.5% |
| MetaCyc | 24247 | 24249 | 23 | 14.1% | 0.000 | 75.9% |
| Rhea | 235 | 235 | 0 | 9.5% | 0.000 | 83.6% |

Pairing: per shared `(external_id, kind)`, both value sets sorted and paired by
rank up to the shorter of the two -- **no length filter**. Sort direction is
irrelevant; ascending and descending zips give the same multiset of differences.
Across all 321,345 paired values: median **0.000**, 75.9% within 0.5.

### The recovered population on its own

| source | comparable sets | site-count match | median abs delta | within 0.5 |
|---|---:|---:|---:|---:|
| ChEBI | 3829 | 3485 (91%) | 0.000 | 86.0% |
| KEGG | 1728 | 1596 (92%) | 0.000 | 84.6% |
| MetaCyc | 10529 | 9438 (90%) | 0.000 | 80.6% |
| Rhea | 59 | 58 (98%) | 0.000 | 92.1% |

Across 62,559 paired values on recovered compounds: median
**0.000**, 82.3% within 0.5 -- in line with the InChI
population, so the plugin path is not a lower-quality source, only a
differently-reachable one.

## Encoding

Values are the three-field `<fragment>:<atom>:<pKa>` form these per-source files
keep, semicolon-joined, in cxcalc's own significance order (`apKa1` is the most
significant acidic value) -- deliberately not re-sorted by atom index. These are
microscopic values, per `pka_encoding.py`: per-site predictions on one protonation
state, not a ladder, so significance is the only ordering that carries
information.

Atom indices are in Marvin's atom space, not ModelSEED's. Quoting
`pka_encoding.py`: "Marvin reorders atoms on import, so its indices never
described our structures." They are kept because this file *is* the provenance
record; nothing in the energy path reads the atom slot.

Re-running the pipeline reproduces identical pKa **values** and identical site
counts, but atom indices may permute among symmetry-equivalent sites -- observed
for KEGG `C04477` and MetaCyc `CPD-8254`, where the same eight values were
reassigned among equivalent phosphate oxygens. That is empirical confirmation of
the line above: the atom slot is provenance, not data.

## Not done here

- Compound records are unchanged. `Update_Compound_pKas.py` is the step that
  rewrites them, and it will touch nearly every compound carrying a pKa: 23.4 was
  atom-index ordered and 26.1 is significance ordered, so the stored *strings*
  differ even where the values agree.
- Protonations (`majorms`) were not regenerated; `protonations/` still holds only
  the 23.4 bundle. The plugin route documented here should apply there too, since
  the CLI will refuse the same polymer structures.
- ChEBI keeps its `CHEBI_` prefix and Rhea its `POLYMER_` prefix. `sources.yaml`
  asks future regenerations to drop them; kept here so the id-format migration
  stays a separate reviewable change.
