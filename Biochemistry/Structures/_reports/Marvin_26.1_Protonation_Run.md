# Marvin 26.1 protonation regeneration

Generated 2026-09-22 by `Scripts/Structures/Run_Marvin_Protonations.py`.

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
| −2 or lower | 2,260 | 4.26% |
| −1 | 6,234 | 11.76% |
| **0** | **42,251** | **79.69%** |
| +1 | 1,888 | 3.56% |
| +2 or higher | 388 | 0.73% |

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
| ChEBI | 7,135 | 7,058 (98.9%) | 8,618/8,918 (96.6%) |
| KEGG | 8,351 | 8,263 (98.9%) | 8,833/9,031 (97.8%) |
| MetaCyc | 13,743 | 13,554 (98.6%) | 16,748/17,448 (96.0%) |
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
| MetaCyc | 25,070 | 19,434 (77.5%) | 13,534/19,297 (70.1%) | 13,534/19,297 (70.1%) |
| Rhea | 237 | 155 (65.4%) | 127/207 (61.4%) | 127/207 (61.4%) |
| **total** | **53,126** | **42,251 (79.5%)** | **32,651/44,294 (73.7%)** | **32,651/44,294 (73.7%)** |

Net charge is the honest headline rather than the formula string, because a
formula moves whenever the hydrogen count does and so mixes the protonation
change with everything else. Both columns are now derived the same way in both
bundles — see The formula column below.

InChI and InChIKey agreement track each other exactly, which is the internal
consistency check passing: the key always hashes the string beside it.

## The formula column, and the R groups it nearly lost

The first cut of this bundle wrote Marvin's own `getFormula()` into the formula
column and deferred `Print_Structure_Formula_Charge.py` to a follow-up step.
That was wrong and it shipped. **Marvin omits wildcard atoms from a formula;
this repository renders them as R**, in one line of that script:

```python
formula = re.sub(r'\*', 'R', formula)
```

The result was that SMILE rows whose formula contains R fell from 23.4's 8,704
to **zero**. `Stearoyl-ACPs` went from `C32H60N3O9PR2S` to `C32H60N3O9PS`, and
downstream `Update_Compound_Structures_Formulas_Charge.py` propagated that into
6,052 compound records — turning `cpd00049` "carboxylic acid" from `CHO2R` into
`CHO2`, which is a generic compound quietly ceasing to be generic.

**No structure was affected by THIS defect.** 8,727 SMILE structures carry a
`*` in both bundles, identically; the R groups were lost from the column alone.
That is exactly why it survived the original validation suite — coverage,
compound sets, InChI agreement, InChIKey consistency and stereochemical
fidelity all inspect *structures*, and on this axis the structures were right.
A column-only defect was invisible to all of it. The table below is the check
that was missing.

(Structures were not uniformly safe, though. A separate defect — InChI-first
input selection meeting InChI's metal disconnection — shattered 461 metal
compounds in the same first cut. See "InChI disconnects metals" below. The two
defects are unrelated in mechanism and were found by different checks.)

| source | SMILE rows with R (26.1) | 23.4 | identical formula | identical charge |
|---|---:|---:|---:|---:|
| ChEBI | 2,071 | 2,071 | 16,129/21,016 (76.7%) | 16,138/21,016 (76.8%) |
| KEGG | 957 | 957 | 26,345/31,593 (83.4%) | 26,371/31,593 (83.5%) |
| MetaCyc | 5,669 | 5,646 | 33,723/44,367 (76.0%) | 33,788/44,367 (76.2%) |
| Rhea | 30 | 30 | 282/444 (63.5%) | 282/444 (63.5%) |
| **total** | **8,727** | **8,704** | **76,479/97,420 (78.5%)** | **76,579/97,420 (78.6%)** |

Count R with `R(?![a-z])`, not a substring test for `"R"`: the elements Ru, Rb,
Rh, Re and Rn match the naive test and inflate every figure in this table by 8.
That is not hypothetical — the first version of this section reported 8,735 and
8,712 for exactly that reason, and it makes the regression look milder than it
was, since the eight survivors it appeared to leave were all ruthenium and
rubidium compounds rather than R groups.

Formula and charge now come from `parse_structure` — this repository's own
function, imported rather than reimplemented, computed per row from that row's
structure string. The ~21% that differ from 23.4 are the protonation-state
changes described above; a formula moves whenever the hydrogen count does.

The 26.1 total exceeds 23.4's because the convention is now enforced as an
invariant — **a SMILE structure carrying `*` gets an R in its formula** — which
23.4 violates 23 times and this bundle violates 0 times. Two parsers had to be
reconciled to get there: RDKit renders a dummy atom as `*`, so the substitution
above catches it, but OpenBabel and Marvin both omit dummy atoms entirely, so on
those paths there is nothing for it to rewrite. Structures like `ISOCITHASE-P`'s
`*OP(=O)(=O)=O` have deliberately invalid valences that only OpenBabel will
read, and were losing their R that way.

Counting wildcards to enforce this has its own trap, recorded in the script: a
molblock `R` atom — which is how RDKit writes *every* dummy atom, and therefore
how every structure arrives through the import bridge — reads back from Marvin
as symbol `R#`, not `R`. A wildcard set without `R#` counts zero on a molecule
that plainly has them.

Note when checking this yourself that an InChI string may contain a literal `*`
that is **not** a wildcard: InChI uses `n*` for repeated components, so
`InChI=1S/Mn.2H2O/h;2*1H2` is manganese with two waters. A naive grep reports
255 false positives here and 217 in 23.4. The invariant is meaningful on SMILE
rows only.

## InChI disconnects metals, and that shattered 461 compounds

The first cut of this bundle took InChI-first unconditionally, inheriting the
rule from `Run_Marvin_pKas.py` where it is correct. On metal compounds it is
not, because **InChI disconnects metal–ligand bonds by design**. Triphenyltin
chloride is stored in `smiles.tsv` as the intact molecule, but its InChI is

```
InChI=1S/3C6H5.ClH.Sn/c3*1-2-4-6-5-3-1;;/h3*1-5H;1H;/q;;;;+1/p-1
```

— three phenyl radicals, HCl and a tin atom, five separate components. Marvin
was handed an already-shattered molecule, protonated each piece, and the bundle
shipped `[Cl-].[SnH3+].[c]1ccccc1.[c]1ccccc1.[c]1ccccc1` where 23.4 shipped the
intact structure.

**487 compounds** have an InChI more fragmented than their SMILES — cobalamins,
Ni/Fe/Mg porphyrins, molybdenum cofactors, organotins, iron–sulfur clusters —
and **461 shipped fragmented**.

The rule is now InChI-first *unless the InChI is the more fragmented of the
two*. Comparing fragment counts rather than screening for metals keeps it
general: whatever the reason an InChI has taken a molecule apart, the
representation that keeps it together is the better input. 493 compounds take
the SMILES route on that test.

Two subtleties made this harder than it sounds, both recorded in the script:

- **The InChI is demoted, not merely reordered.** Sharing one ladder lets a
  writable-but-wrong rung beat a correct one — on KEGG `C18384` the SMILES
  yields the right dative-bonded magnesium propionate, which Marvin's SMILES
  writer refuses, while the disconnected InChI yields a writable three-fragment
  answer and would win.
- **The InChI ROW is written from the InChI-derived molecule.** Each column
  should carry what its own representation can express, which is exactly what
  23.4 did: for `CPD-18407` it shipped a connected 8-iron cluster in SMILE and
  the disconnected `InChI=1S/C.8Fe.6HS.3S/...` in InChI. This is also load
  bearing for stability — asking Marvin to write an InChI for a *connected*
  metal cluster aborts the JVM outright with `free(): double free detected in
  tcache 2` from `InChINativeGenerateInChICall`, a native fault no Python or
  Java handler can catch.

| fragmented vs 23.4 | count |
|---|---:|
| first cut (InChI-first everywhere) | 461 |
| after the fragment-count rule | 16 |
| after demotion + per-column sourcing | **8** |

The remaining 8 are not an input-selection problem. Marvin 26.1's microspecies
plugin breaks metal coordination bonds during protonation itself, from either
input: ferrocene (`CPD-21742`, `CPD-21743`) splits into iron plus two
cyclopentadienyls, and `C12862` sheds both ammines. Confirmed by protonating
the connected SMILES directly. That is an engine behaviour change, recorded
rather than worked around.

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
- `Print_Structure_Formula_Charge.py` has not been re-run over these files as a
  whole — but it no longer needs to be for the formula and charge columns, which
  are computed with its own `parse_structure` as the bundle is written. Running
  it should be a no-op on those two columns; it was treating it as a deferred
  follow-up that produced the R-group regression above.
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
