# Marvin 26.1 pKa regeneration

Generated 2026-09-10 by `Scripts/Structures/Run_Marvin_pKas.py`.

Regenerates the per-source Marvin pKa bundles against the current
`<source>/inchi.tsv`, replacing the retired `pKaMol.java` path (which needed
MarvinBeans on the classpath) with the `cxcalc` CLI from Marvin Desktop Suite.

## Invocation

```
cxcalc --ignore-error -i x pka --na 20 --nb 20 <structure-list>
```

Calculator 26.1.2 (ChemAxon), one process per input representation per source,
batch mode. `--na/--nb 20` because the widest existing 23.4 row carries 19
tokens. Licensed under the Protonation Plugin Group.

## Coverage

| source | structures in | with pKa | rows |
|---|---:|---:|---:|
| ChEBI | 9473 | 9139 | 17178 |
| KEGG | 15319 | 14599 | 26351 |
| MetaCyc | 19298 | 18605 | 34596 |
| Rhea | 207 | 205 | 385 |
| **total** | **44297** | **42548** | **78510** |

Zero InChI structures failed to parse in any source.

## 23.4 vs 26.1

| source | 23.4 ids | 26.1 ids | shared | new | 23.4-only | median abs delta | within 0.5 |
|---|---:|---:|---:|---:|---:|---:|---:|
| ChEBI | 11186 | 9139 | 9134 | 5 | 2052 | 0.000 | 72.2% |
| KEGG | 15967 | 14599 | 14588 | 11 | 1379 | 0.000 | 75.9% |
| MetaCyc | 24247 | 18605 | 18587 | 18 | 5660 | 0.000 | 74.3% |
| Rhea | 235 | 205 | 205 | 0 | 30 | 0.000 | 82.4% |

Pairing: for each shared `(external_id, kind)`, both value sets are sorted and
paired by rank, up to the shorter of the two. **No length filter** -- sets where
the two releases report a different number of sites are included, paired as far
as they go. Sort direction is irrelevant: ascending and descending zips give the
same multiset of differences.

Across all 258,786 rank-paired values: median **0.000**, 74.4% within 0.5.

### Why the site count itself moved

| source | shared sets | same count | different count | % differing |
|---|---:|---:|---:|---:|
| ChEBI | 17120 | 14302 | 2818 | 16.5% |
| KEGG | 26224 | 22834 | 3390 | 12.9% |
| MetaCyc | 34418 | 29187 | 5231 | 15.2% |
| Rhea | 384 | 343 | 41 | 10.7% |
| **total** | **78146** | **66666** | **11480** | **14.7%** |

14.7% of shared sets have a different number of predicted
sites in 26.1 than in 23.4 -- the release changed which atoms it considers
ionizable, not only the values it assigns them. Restricting the comparison to
equal-count sets would flatter the `within 0.5` column by conditioning on the
sets that changed least; the table above does not do this.

## What the 23.4-only column actually is

Not compounds pruned since the 2024-01 snapshot. Three buckets:

| source | 23.4-only | polymer/organometallic (SMILES-only) | in neither structure file | other |
|---|---:|---:|---:|---:|
| ChEBI | 2052 | 2039 | 0 | 13 |
| KEGG | 1379 | 931 | 403 | 45 |
| MetaCyc | 5660 | 5637 | 0 | 23 |
| Rhea | 30 | 30 | 0 | 0 |
| **total** | **9121** | **8637** | **403** | **81** |

**8637 are compounds with a SMILES and no InChI**, and cxcalc cannot
process any of them -- the reason they lack an InChI is the same reason cxcalc
refuses them. Across all four sources the 8834 SMILES-only compounds break
down as 8,728 carrying `*` attachment points (structural repeating units), which
Marvin reads as QUERY molecules and cxcalc declines with *"pka: Calculation
result is not defined for query molecules"*; 105 organometallics (Mg-porphyrins,
chlorophylls) that parse but return no pKa and an empty `atoms` column; and 1
dative-bond SMILES the parser rejects at the `<` character.

Their 23.4 values came from `pKaMol.java` -- the Java `pKaPlugin` API, which
accepted `*`-bearing structures that the cxcalc CLI will not. That API was
renamed in 26.1 (`chemaxon.calculations.PkaPlugin`) and the committed
`pKaMol.class` no longer runs against it. **Until it is rewritten,
`Compounds.loadPerSourcePkas`' accumulation across bundles is the only mechanism
preserving pKas for these compounds** -- it globs every TSV in `pkas/` and lets
the last-sorted file win per `(ext_id, kind)`, so 26.1 overrides 23.4 on shared
ids while 23.4-only ids keep their 23.4 values. That accumulation is load-bearing
here, not incidental.

The 403 in neither structure file are dominated by KEGG's 403, all of
which appear in `KEGG/KEGG_SRU_041020.txt` -- structural repeating units filtered
out of the current structure files. The remaining 81 are simply not in
the current `inchi.tsv` or `smiles.tsv`.

## InChI and SMILES are not interchangeable

Measured on 600 ChEBI compounds that have both representations, running cxcalc on
the InChI and on the SMILES:

| SMILES charge | sets | count mismatch | within 0.05 | within 0.5 |
|---|---:|---:|---:|---:|
| neutral | 836 | 42 (5.0%) | 80.7% | 87.3% |
| charged | 225 | 61 (27.1%) | 53.1% | 66.0% |

A charged SMILES is an already-deprotonated species, so Marvin is answering a
different question than it is for the neutral InChI parent. InChI therefore wins
wherever it exists and SMILES is used only to fill the gap, so no InChI-derived
value is displaced. `--structures {inchi,smiles,both}` exposes the choice.

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
  the 23.4 bundle. Note that bundle carries SMILE rows for the polymer structures
  too, so a 26.1 protonation run should be expected to hit the same query-molecule
  wall documented above.
- ChEBI keeps its `CHEBI_` prefix and Rhea its `POLYMER_` prefix. `sources.yaml`
  asks future regenerations to drop them; kept here so the id-format migration
  stays a separate reviewable change.
