# Marvin 26.1 protonation regeneration

Generated 2026-09-18 by `Scripts/Structures/Run_Marvin_Protonations.py`.

Regenerates the per-source pH 7 protonation bundles against the current
`<source>/inchi.tsv` and `<source>/smiles.tsv`, replacing the retired
`ChargeMol.java` path (which needed MarvinBeans on the classpath) with
`chemaxon.calculations.MajorMicrospeciesPlugin` driven through JPype.

This is the companion to `Marvin_26.1_pKa_Run.md`, and the bundle that report
listed under "Not done here".

## The tautomer step is missing, and that is the headline

`ChargeMol.java` ran five steps per fragment:

```
 i)   convert explicit hydrogens to implicit
 ii)  aromatize
 iii) find the dominant tautomer at pH 7      <-- NOT REPRODUCED
 iv)  find the major microspecies at pH 7
 v)   convert implicit hydrogens back to explicit
```

Step (iii) needs `TautomerizationPlugin`, licensed under the **Isomers Plugin
Group**. The licence these regenerations run on covers the **Protonation Plugin
Group** only — pKa and major microspecies. The refusal was verified, not
assumed:

```
$ cxcalc "CC(=O)CC(=O)C" majortautomer -H 7
chemaxon.license.api.LicenseException: No valid license has been found.
Product name: Isomers Plugin Group
```

So this run performs (i), (ii), (iv), (v) and skips (iii).

**The effect is not small.** Across the 53,021 compounds both bundles cover, the
net charge at pH 7 is unchanged for 79.7% and different for 20.3%, skewed
toward deprotonation:

| Δ net charge (26.1 − 23.4) | compounds | share |
|---|---:|---:|
| −2 or lower | 2,173 | 4.10% |
| −1 | 6,234 | 11.76% |
| **0** | **42,252** | **79.69%** |
| +1 | 1,887 | 3.56% |
| +2 or higher | 381 | 0.72% |

That 20.3% is an **upper bound on the tautomer effect, not a measurement of
it**. It also contains the genuine 23.4 → 26.1 engine change, which the pKa
regeneration already showed to be substantial: 26.1 disagrees with 23.4 about
the *number* of ionizable sites in roughly one shared set in seven. A pKa that
moves across 7 flips the protonation state at pH 7 on its own, with no tautomer
involved. The two causes cannot be separated without an Isomers licence, and
this report does not claim otherwise.

Restoring the step, if that licence is obtained, is a few lines: insert a
`TautomerizationPlugin` call between the aromatize and the microspecies call in
`protonate()` and regenerate. Nothing else changes.

## Marvin's InChI and SMILES importers lose stereochemistry

This is the finding most likely to matter to anyone else driving Marvin through
JPype, and it is not documented anywhere obvious.

Round-tripping the source InChI with **no protonation at all** — import it,
export it, compare:

| import path | reproduces the source InChI |
|---|---:|
| `MolImporter.importMol(<InChI>)` | 61% |
| `MolImporter.importMol(<SMILES>)` | 39% |
| `MolImporter.importMol(<RDKit molblock>)` | **99.6%** |
| RDKit alone, InChI → mol → InChI (control) | 99.7% |

The loss is on **import**, not export: both Marvin writers score the same ~60%
behind a lossy import, and the RDKit control proves the InChI carries
everything needed. What is lost is the `/t` and `/m` layers — stereo parity.

So every structure is bridged through RDKit into an MDL molblock before Marvin
sees it. That is what `ChargeMol.java` was doing when it read `args[0]` as a
*mol file* rather than taking a string; the significance of that detail is easy
to miss. 2D coordinates are not required — the molblock's parity flags suffice.

Export uses `inchi:AuxNone,SAbs`. The `SAbs` flag is load-bearing too: without
it the same bridged molecule reproduces the source InChI 30.9% of the time,
against 99.6% with it.

**Measured against the benchmark that matters** — compounds where protonation
changed nothing, so the output should equal the input:

| source | no-op compounds | 26.1 reproduces the source | 23.4 on the same set |
|---|---:|---:|---:|
| ChEBI | 7,128 | 7,055 (99.0%) | 8,618/8,918 (96.6%) |
| KEGG | 8,346 | 8,263 (99.0%) | 8,833/9,031 (97.8%) |
| MetaCyc | 13,738 | 13,551 (98.6%) | 16,748/17,448 (96.0%) |
| Rhea | 135 | 133 (98.5%) | 194/197 (98.5%) |

An unbridged run scores 55–61% here and would have shipped a stereochemistry
regression across tens of thousands of rows while every coverage count still
looked perfect.

The InChIKey is **derived from the InChI string actually written** rather than
exported separately, because Marvin's `inchikey` export does not honour `SAbs`
and so disagrees with its own `inchi` export. For Rhea `POLYMER_10033` the two
describe the same molecule but hash differently: Marvin's key is
`...-ARPYZQPTNA-N`, the key of the InChI it just wrote is `...-VFUOTHLCSA-N`,
and the latter is what 23.4 shipped. Deriving makes the columns consistent by
construction — 0 disagreements across all 44,294 InChI rows.

## Coverage

| source | compounds | SMILE | InChI | InChIKey | rows | 23.4 rows | 23.4 compounds |
|---|---:|---:|---:|---:|---:|---:|---:|
| ChEBI | 11,544 | 11,544 | 9,472 | 9,472 | 30,488 | 30,488 | 11,544 |
| KEGG | 16,275 | 16,275 | 15,318 | 15,318 | 46,911 | 46,911 | 16,275 |
| MetaCyc | 25,071 | 25,071 | 19,297 | 19,297 | 63,665 | 63,665 | 25,071 |
| Rhea | 237 | 237 | 207 | 207 | 651 | 651 | 237 |
| **total** | **53,127** | **53,127** | **44,294** | **44,294** | **141,715** | **141,715** | **53,127** |

Row counts match 23.4 exactly in every source. ChEBI, KEGG and Rhea cover an
identical set of compounds; MetaCyc differs by one in each direction — see
Failures.

The layout reproduces 23.4's, which is one row per *representation*, not per
compound: a SMILE row for every compound in `smiles.tsv`, and InChI plus
InChIKey rows only for compounds that also appear in `inchi.tsv`. The asymmetry
is chemistry, not omission — the compounds with no InChI are polymers carrying
`*` attachment points and organometallics, and InChI cannot represent a query
molecule.

## Agreement with 23.4

| source | shared ids | same net charge | identical InChI | identical InChIKey |
|---|---:|---:|---:|---:|
| ChEBI | 11,544 | 9,026 (78.2%) | 6,795/9,472 (71.7%) | 6,795/9,472 (71.7%) |
| KEGG | 16,275 | 13,636 (83.8%) | 12,195/15,318 (79.6%) | 12,195/15,318 (79.6%) |
| MetaCyc | 25,070 | 19,435 (77.5%) | 13,534/19,297 (70.1%) | 13,534/19,297 (70.1%) |
| Rhea | 237 | 155 (65.4%) | 127/207 (61.4%) | 127/207 (61.4%) |
| **total** | **53,126** | **42,252 (79.5%)** | **32,651/44,294 (73.7%)** | **32,651/44,294 (73.7%)** |

Net charge is the honest headline rather than the formula string: 23.4's
formula and charge columns were re-derived by `Print_Structure_Formula_Charge.py`
with RDKit/OpenBabel while these are Marvin's, so the two can spell the same
molecule differently — R-group counts on the polymers especially — without
disagreeing about chemistry.

InChI and InChIKey agreement track each other exactly, which is the internal
consistency check passing: the key always hashes the string beside it.

## One engine, and why not the CLI

Every compound goes through the Java plugin, including the ones `cxcalc majorms`
could have handled, so that all three representations come from the same
in-memory molecule.

The two agree on chemistry — across 300 ChEBI compounds formula and charge are
identical in 296. **All four exceptions are organometallics** (Mg, Fe and Ni
porphyrins and chlorophylls) where `cxcalc majorms` silently drops the metal
fragment: it returns `C35H34N4O5` for a compound whose protonated form is
`C35H34MgN4O5`. The per-fragment route keeps it. That is the same reason
`ChargeMol.java` fragmented and re-fused rather than protonating the molecule
whole, and it is why the CLI is not used for the bulk.

## Query molecules

Two structure classes need care, and both are recovered rather than dropped:

- `aromatize()` on a molecule containing query atoms leaves bonds no SMILES
  writer can express. Un-aromatized, ChEBI `26432` protonates and writes
  cleanly, and its `[NH+]` matches what 23.4 reported.
- MetaCyc spells its protein-bound cofactors with `~` (any-bond), as in
  `*~O=C/C=C(C)/...` for the retinal-binding proteins. Marvin's SMILES writer
  throws on those rather than degrading them; going out through a molblock and
  letting RDKit write the SMILES keeps the query bond. On `11C-Retinal-RALBPs`
  that route returns the 23.4 string character for character. Without it, 156
  MetaCyc compounds would have been silently lost.

`Run_Marvin_Protonations.protonate_best()` walks a ladder — each representation
aromatized then un-aromatized, InChI before SMILES — with Marvin's own writer
across every rung *before* any molblock rescue, so a degraded early result can
never win over a clean later one.

## Charge on wildcard-adjacent atoms is inherited, not introduced

On the Rhea glycans, the peptide-backbone nitrogen hidden behind `*` is
protonated to `[NH2+]`, because Marvin cannot see the amide bond the wildcard
stands for. In the real polymer that nitrogen is an amide and would be neutral
at pH 7.

This is not new. 23.4 wrote `*[NH2+][C@@H](...` for `POLYMER_12621` and this run
writes the same. Correcting it would be a curation decision about what `*`
means, not a regeneration, so it is recorded here and left alone.

## Failures

Four compounds across 53,131 structures yielded no row, and only one is a
failure:

- **H⁺** — ChEBI `15378`, KEGG `C00080`, MetaCyc `PROTON`, all
  `InChI=1S/p+1`. No heavy atom, so no protonated form. 23.4 carries no row for
  these either; agreeing with the old bundle is the correct outcome.
- **MetaCyc `HypC-Dimer-Fe-CO2`** (`*N~[Fe+2](~S*)~C(=O)=O`) — a genuine
  failure. `PkaPlugin.run()` throws `ArrayIndexOutOfBoundsException: Index -1
  out of bounds for length 8`. 23.4 has one row for it; this is a loss of one
  row out of 141,715. The id is recorded rather than only the count.

MetaCyc also *gains* `A-DNA-WITH-OPPOSING-AP-SITE`, which 23.4 did not cover,
so its compound total is unchanged at 25,071.

Zero structures were unparseable, zero produced a non-standard InChI, and no
compound that carries an InChI in the source failed to get an InChI row back.

## Reproducibility

Re-running a source rewrites the bundle **byte for byte**. This is a stronger
guarantee than the pKa regeneration could give — there, re-runs reproduced
identical values and site counts but atom indices could permute among
symmetry-equivalent sites. Nothing in this bundle carries an atom index, so
there is nothing left to permute.

## Not done here

- Compound records are unchanged. `Update_Compound_Structures_Formulas_Charge.py`
  is the step that rewrites them.
- `Print_Structure_Formula_Charge.py` has not been re-run. It re-derives the
  formula and charge columns of this file in place with RDKit/OpenBabel; the
  values committed here are Marvin's own and agree with 23.4's on spot checks.
  Pyruvate (`C00022`) is representative: identical InChI
  (`InChI=1S/C3H4O3/...../p-1`), identical InChIKey
  (`LCTONWCANYUPML-UHFFFAOYSA-M`), identical formula and charge (`C3H3O3`, −1).
  Its SMILE differs as a *string* — `CC(=O)C([O-])=O` here against
  `CC(=O)C(=O)[O-]` in 23.4 — because 23.4's SMILES column was written by
  RDKit/OpenBabel from Marvin's molfile while these come from Marvin's own
  writer. Same molecule, different canonical atom order; the InChIKey agreeing
  is the proof. No attempt was made to re-canonicalise 53,127 SMILES strings to
  match the old writer's output.
- The tautomer step, pending an Isomers Plugin Group licence.
- ChEBI keeps bare ids in this bundle, matching the 23.4 protonation file —
  note that the *pKa* bundles use a `CHEBI_` prefix. The id-format migration
  stays a separate, reviewable change.
