# Marvin 26.1 pKa regeneration

Generated 2026-09-07 by `Scripts/Structures/Run_Marvin_pKas.py`.

Regenerates the per-source Marvin pKa bundles against the current
`<source>/inchi.tsv`, replacing the retired `pKaMol.java` path (which needed
MarvinBeans on the classpath) with the `cxcalc` CLI from Marvin Desktop Suite.

## Invocation

```
cxcalc --ignore-error -i x pka --na 20 --nb 20 <inchi-list>
```

Calculator 26.1.2 (ChemAxon), one process per source, batch mode.
`--na/--nb 20` because the widest existing 23.4 row carries 19 tokens.
Licensed under the Protonation Plugin Group.

## Coverage

| source | structures in | with pKa | rows |
|---|---:|---:|---:|
| ChEBI | 9473 | 9139 | 17178 |
| KEGG | 15319 | 14599 | 26351 |
| MetaCyc | 19298 | 18605 | 34596 |
| Rhea | 207 | 205 | 385 |
| **total** | **44297** | **42548** | **78510** |

Zero structures failed to parse in any source.

## 23.4 vs 26.1

| source | 23.4 ids | 26.1 ids | shared | new | 23.4-only | median abs delta | within 0.5 | max delta |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| ChEBI | 11186 | 9139 | 9134 | 5 | 2052 | 0.000 | 83.3% | 23.87 |
| KEGG | 15967 | 14599 | 14588 | 11 | 1379 | 0.000 | 82.7% | 35.16 |
| MetaCyc | 24247 | 18605 | 18587 | 18 | 5660 | 0.000 | 83.1% | 32.76 |
| Rhea | 235 | 205 | 205 | 0 | 30 | 0.000 | 89.3% | 15.53 |

Values compared as sorted sets per (id, kind) where both versions report the
same count. The median is 0.000 in every source: the two Marvin releases agree
on the great majority of sites. The tail beyond 0.5 is genuine model change
across three releases, not parse error.

The `23.4-only` column is compounds pruned from `inchi.tsv` since the
2024-01 snapshot. They are not lost: `Compounds.loadPerSourcePkas` accumulates
across every TSV in `pkas/` (`out.setdefault((db, ext_id), {})[kind] = value`)
and lets the last-sorted file win per key, so 26.1 overrides 23.4 on shared ids
while 23.4-only ids keep their 23.4 values. Verified: 51,635 ids in 23.4, 34 new
in 26.1, and the loader returns exactly 51,669.

## The id-drift hazard

Under `--ignore-error`, cxcalc *silently drops* molecules it cannot parse **and
renumbers the surviving rows**, so its default `id` column is a position in the
output, not in the input. A single unparseable structure would shift every
subsequent id onto the wrong compound.

Passing `-i` makes cxcalc emit `idError[N]`, where `N` is the true 1-based input
line. The script anchors on that and asserts the indices are strictly increasing
and within range before writing anything.

## Encoding

Values are the three-field `<fragment>:<atom>:<pKa>` form these per-source files
keep, semicolon-joined, in cxcalc's own significance order (`apKa1` is the most
significant acidic value) — deliberately not re-sorted by atom index. These are
microscopic values, per `pka_encoding.py`: per-site predictions on one protonation
state, not a ladder. Significance is the only ordering that carries information.

Atom indices are in Marvin's atom space, not ModelSEED's. Quoting
`pka_encoding.py`: "Marvin reorders atoms on import, so its indices never
described our structures." They are retained because this file is the provenance
record; nothing in the energy path reads the atom slot.

## Not done here

- Compound records are unchanged. `Update_Compound_pKas.py` is the step that
  rewrites them, and it will touch nearly every compound carrying a pKa: 23.4 was
  atom-index ordered and 26.1 is significance ordered, so the stored *strings*
  differ even where the values agree.
- Protonations (`majorms`) were not regenerated; `protonations/` still holds only
  the 23.4 bundle.
- ChEBI keeps its `CHEBI_` prefix and Rhea its `POLYMER_` prefix. `sources.yaml`
  asks future regenerations to drop them; kept here so the id-format migration
  stays a separate reviewable change rather than a side effect of this upgrade.
