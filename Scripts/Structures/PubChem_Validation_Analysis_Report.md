# PubChem Structure Validation — Analysis Report

**Date:** 2026-06-01
**Script:** `Scripts/Structures/Validate_PubChem_Structures.py`
**Outputs:** `Biochemistry/Structures/{All,Unique}_ModelSEED_Structures_updated.txt`
**Result:** 68 compound structures corrected; 40 new compounds resolved; 3 removed.

---

## 1. Purpose
Validate every structured ModelSEED compound against PubChem, classify the
differences, auto-apply only the safe corrections, and route ambiguous cases to
review files. Inputs are `All_ModelSEED_Structures.txt`, the compound name and
ChEBI/KEGG alias tables, and ChemAxon pKa data; outputs are the corrected
`*_updated` structure files plus a set of review TSVs in
`structure_review_output/`.

## 2. Pipeline architecture

The run executes seven phases (0–6). PubChem queries are rate-limited, retried
with backoff, and cached in `pubchem_cache.sqlite`, so re-runs are resumable and
offline once the cache is warm.

**Phase 0 — Self-consistency.** For each compound's SMILES: parse with RDKit;
compute the InChI and InChIKey if missing; verify the stored InChIKey matches
the one derived from the InChI; verify SMILES and InChI agree on connectivity;
and recompute formula and charge from the InChI, fixing them where they
disagree. R-group/polymer SMILES (`*`) and unparseable SMILES are skipped and
logged.

**Phase 1 — External-ID lookup.** ChEBI and KEGG aliases are resolved to PubChem
CIDs via the `xref/RegistryID` endpoint. When ChEBI and KEGG point a compound to
*different* CIDs (an xref conflict), the candidates are scored against the stored
InChIKey (MATCH > PROTONATION_DIFF > STEREO_DIFF) and the best is chosen;
unresolved conflicts are reported.

**Phase 2 — Name lookup.** Compounds not resolved by xref are queried by
compound name in batches. When a name maps to several CIDs, disambiguation tries
exact InChIKey match, then connectivity+stereo block match, then connectivity
block, then highest Tanimoto similarity to the stored structure.

**Phase 3 — Recovery.** Still-unresolved or mismatched compounds get a direct
InChIKey lookup, then a direct InChI lookup, to recover matches missed by the
xref/name routes.

**Phase 4 — Classification & reporting.** Stored vs PubChem InChIKeys are
compared into: `MATCH`, `PROTONATION_DIFF` (same connectivity+stereo block),
`STEREO_DIFF` (same connectivity block), `MISMATCH`, `NOT_FOUND`,
`XREF_CONFLICT`, `AMBIGUOUS`. Each `MISMATCH` is sub-classified using formula
comparison, InChI-layer comparison, O-ring (ring/chain) counting, tautomer
canonicalization, metal detection, and Morgan-fingerprint Tanimoto — yielding
tautomer / salt-form / metal-representation / wrong-mapping / positional-isomer /
different-compound labels. Some mismatches are reclassified to STEREO_DIFF or
PROTONATION_DIFF when the InChI layers show only stereo or protonation differ.

**Phase 5 — Corrections.** Three correction types are produced, each guarded:
- *STEREO_DIFF* — accepted only if PubChem's SMILES/InChIKey are mutually
  consistent, PubChem has **≥** the stored number of defined stereocenters
  (never fewer), and no shared stereocenter is inverted. The inversion check
  compares the InChI `/t` (relative configuration) **and** `/m` (absolute
  configuration) layers, so a full enantiomer — which shares `/t` and differs
  only in `/m` — is rejected, not adopted.
- *PROTONATION_DIFF / pKa* — protonation is adjusted toward the charge predicted
  from ChemAxon pKa at pH 7. Guards: borderline pKa 6–8 is skipped; corrections
  that would re-protonate a strong acid (carboxylic/sulfonic/sulfinic) are
  rejected; a pKa cross-validation keeps the stored charge when it fits the data
  better; and metal coordination complexes are skipped entirely (RDKit's
  uncharge/reprotonate cannot safely manipulate metal–ligand bonds).
- Rejected corrections are written to `pubchem_rejected_corrections.tsv` with the
  reason.

**Phase 6 — Apply & write output.** Accepted corrections are written into the
All file; formula and charge are recomputed from each corrected InChI; **all**
SMILES are re-canonicalized to RDKit canonical form. The `All_updated` file is
written atomically, then the `Unique_updated` file is rebuilt from it, and a
post-correction pass re-checks every corrected compound for SMILES/InChI/InChIKey
consistency.

## 3. Unique-file picker (matches dev-branch `List_ModelSEED_Structures.py`)
`build_unique_file` reproduces the dev-branch selection algorithm: reference
type/stage priority **InChI → SMILE**, **Charged → Original**; a compound passes
if it has one structure, or several structures sharing one formula; differing
formulas are reported as conflicts and no pick is made. Conflicts among
same-formula structures are resolved by **multi-database agreement → InChIKey
connectivity grouping → stereo preference (avoid `UHFFFAOYSA`) → source priority
MetaCyc > KEGG > ChEBI > Rhea**. Curated ignore lists are read from both the
legacy (`Biochemistry/Structures/Curation`) and dev
(`Biochemistry/Curation/ignores`) locations (identical 132-ID set).

## 4. Corrected compounds (68)

All 68 corrections are stereochemistry or protonation refinements on existing
structures — there are no connectivity rewrites. Breakdown:

| Change type | Count |
|---|---|
| Stereochemistry gained (PubChem adds detail to undefined centers) | 48 |
| Stereo/representation change (e.g. porphyrin double-bond/InChIKey layout) | 13 |
| Protonation + stereo gain (e.g. phenol → neutral at pH 7, plus stereo) | 6 |
| Protonation change | 1 |

Full before→after for all 68 is also written to
`Scripts/Structures/Corrected_Compounds_Before_After.tsv` (full strings,
machine-readable). Listing below is sorted by reaction usage (most-used first).

### Before → after (all 68)

**1. cpd01476 — 7,12-diethenyl-3,8,13,17-tetramethylporphyrin-2,18-dipropanoate**  (42 rxns; InChIKey ZCFFYALKHPIRKJ-UJJXFSCMSA-L → ZCFFYALKHPIRKJ-UHFFFAOYSA-L)
- before: `C=Cc1c(C)c2cc3[nH]c(cc4nc(cc5nc(cc1[nH]2)C(C)=C5CCC(=O)[O-])C(CCC(=O)[O-])=C4C)c(C)c3C=C`
- after:  `C=Cc1c(C)c2cc3[nH]c(cc4nc(cc5nc(cc1[nH]2)C(C)=C5CCC(=O)[O-])C(CCC(=O)[O-])=C4C)c(C)c3C=C`

**2. cpd39967 — 3,7,11,15-tetramethylhexadeca-2,6,10,14-tetraen-1-yl diphosphate**  (9 rxns; InChIKey OINNEUNVOZHBOX-UHFFFAOYSA-L → OINNEUNVOZHBOX-QIRCYJPOSA-L, stereo 0→3)
- before: `CC(C)=CCCC(C)=CCCC(C)=CCCC(C)=CCOP(=O)([O-])OP(=O)([O-])O`
- after:  `CC(C)=CCC/C(C)=C/CC/C(C)=C/CC/C(C)=C/COP(=O)([O-])OP(=O)([O-])O`

**3. cpd02168 — (1R,2R,3R)-prephytoene diphosphate**  (8 rxns; InChIKey RVCNKTPCHZNAAO-IMSLGMFESA-L → RVCNKTPCHZNAAO-UZDKSQMHSA-L)
- before: `CC(C)=CCC/C(C)=C/CC/C(C)=C/CC/C(C)=C/[C@H]1[C@H](COP(=O)([O-])OP(=O)([O-])O)[C@]1(C)CC/C=C(\C)CC/C=C(\C)CCC=C(C)C`
- after:  `CC(C)=CCC/C(C)=C/CC/C(C)=C/CC/C(C)=C/[C@@H]1[C@@H](COP(=O)([O-])OP(=O)([O-])O)[C@]1(C)CC/C=C(\C)CC/C=C(\C)CCC=C(C)C`

**4. cpd03913 — Hydrogenobyrinate a,c diamide**  (6 rxns; InChIKey ZGGWTIPDUOTHRA-FVLIVPSMSA-J → ZGGWTIPDUOTHRA-QQHNNQQLSA-J, stereo 10→11)
- before: `CC1=C2[NH+]=C(C=C3N/C(=C(/C)C4=N[C@@](C)([C@@H]5N=C1[C@](C)(CCC(=O)[O-])[C@H]5CC(=O)[O-])[C@@](C)(CC(N)=O)[C@@H]4CCC(=O)[O-])[C@@](C)(CC(N)=O)[C@@H]3CCC(=O)[O-])C(C)(C)[C@@H]2CCC(=O)[O-]`
- after:  `C/C1=C2/N=C(/C=C3\N/C(=C(/C)C4=N[C@@](C)(C5N=C1[C@](C)(CCC(=O)[O-])[C@H]5CC(=O)[O-])[C@@](C)(CC(N)=O)[C@@H]4CCC(=O)[O-])[C@@](C)(CC(N)=O)[C@@H]3CCC(=O)O)C(C)(C)[C@@H]2CCC(=O)[O-]`

**5. cpd23886 — neo-lambda-carratetraose**  (5 rxns; InChIKey SYQAURQZRLJKJT-GTBSAFNRSA-H → SYQAURQZRLJKJT-JYOYXRFGSA-H, stereo 19→20)
- before: `O=S(=O)([O-])OC[C@H]1O[C@H](O[C@H]2[C@@H](O)[C@@H](CO)O[C@@H](O[C@@H]3[C@H](O)[C@@H](OS(=O)(=O)[O-])[C@@H](O[C@H]4[C@@H](O)[C@@H](CO)OC(O)[C@@H]4OS(=O)(=O)[O-])O[C@@H]3COS(=O)(=O)[O-])[C@@H]2OS(=O)(=O)[O-])[C@H](OS(=O)(=O)[O-])[C@@H](O)[C@H]1O`
- after:  `O=S(=O)([O-])OC[C@H]1O[C@H](O[C@H]2[C@@H](O)[C@@H](CO)O[C@@H](O[C@@H]3[C@H](O)[C@@H](OS(=O)(=O)[O-])[C@@H](O[C@H]4[C@@H](O)[C@@H](CO)O[C@@H](O)[C@@H]4OS(=O)(=O)[O-])O[C@@H]3COS(=O)(=O)[O-])[C@@H]2OS(=O)(=O)[O-])[C@H](OS(=O)(=O)[O-])[C@@H](O)[C@H]1O`

**6. cpd23887 — neo-lambda-carrahexaose**  (5 rxns; InChIKey XCAVJOBVNSMNME-WRUHAMNASA-E → XCAVJOBVNSMNME-HOQVOYLOSA-E, stereo 29→30)
- before: `O=S(=O)([O-])OC[C@H]1O[C@H](O[C@H]2[C@@H](O)[C@@H](CO)O[C@@H](O[C@@H]3[C@H](O)[C@@H](OS(=O)(=O)[O-])[C@@H](O[C@H]4[C@@H](O)[C@@H](CO)O[C@@H](O[C@@H]5[C@H](O)[C@@H](OS(=O)(=O)[O-])[C@@H](O[C@H]6[C@@H](O)[C@@H](CO)OC(O)[C@@H]6OS(=O)(=O)[O-])O[C@@H]5COS(=O)(=O)[O-])[C@@H]4OS(=O)(=O)[O-])O[C@@H]3COS(=O)(=O)[O-])[C@@H]2OS(=O)(=O)[O-])[C@H](OS(=O)(=O)[O-])[C@@H](O)[C@H]1O`
- after:  `O=S(=O)([O-])OC[C@H]1O[C@H](O[C@H]2[C@@H](O)[C@@H](CO)O[C@@H](O[C@@H]3[C@H](O)[C@@H](OS(=O)(=O)[O-])[C@@H](O[C@H]4[C@@H](O)[C@@H](CO)O[C@@H](O[C@@H]5[C@H](O)[C@@H](OS(=O)(=O)[O-])[C@@H](O[C@H]6[C@@H](O)[C@@H](CO)O[C@@H](O)[C@@H]6OS(=O)(=O)[O-])O[C@@H]5COS(=O)(=O)[O-])[C@@H]4OS(=O)(=O)[O-])O[C@@H]3COS(=O)(=O)[O-])[C@@H]2OS(=O)(=O)[O-])[C@H](OS(=O)(=O)[O-])[C@@H](O)[C@H]1O`

**7. cpd03832 — Hydrogenobyrinate**  (4 rxns; InChIKey FJDBIDBCPIOVPG-FVLIVPSMSA-H → FJDBIDBCPIOVPG-QQHNNQQLSA-H, stereo 10→11)
- before: `CC1=C2[NH+]=C(C=C3N/C(=C(/C)C4=N[C@@](C)([C@@H]5N=C1[C@](C)(CCC(=O)[O-])[C@H]5CC(=O)[O-])[C@@](C)(CC(=O)[O-])[C@@H]4CCC(=O)[O-])[C@@](C)(CC(=O)[O-])[C@@H]3CCC(=O)[O-])C(C)(C)[C@@H]2CCC(=O)[O-]`
- after:  `C/C1=C2/N=C(/C=C3\N/C(=C(/C)C4=N[C@@](C)(C5N=C1[C@](C)(CCC(=O)[O-])[C@H]5CC(=O)[O-])[C@@](C)(CC(=O)[O-])[C@@H]4CCC(=O)[O-])[C@@](C)(CC(=O)O)[C@@H]3CCC(=O)[O-])C(C)(C)[C@@H]2CCC(=O)[O-]`

**8. cpd23840 — 2,2'-azino-bis-(3-ethylbenzothiazoline-6-sulfonate) radical cation**  (4 rxns; InChIKey HFNOFNDLBTVECO-YAFCTCPESA-L → HFNOFNDLBTVECO-UHFFFAOYSA-L)
- before: `CCn1/c(=N/N=C2/Sc3cc(S(=O)(=O)[O-])ccc3[N+]2CC)sc2cc(S(=O)(=O)[O-])ccc21`
- after:  `CCn1/c(=N/N=C2/Sc3cc(S(=O)(=O)[O-])ccc3[N+]2CC)sc2cc(S(=O)(=O)[O-])ccc21`

**9. cpd35693 — coniferyl alcohol radical**  (4 rxns; InChIKey ORAJWSYKRGVTDP-FPYGCLRLSA-N → ORAJWSYKRGVTDP-UHFFFAOYSA-N)
- before: `COC1=C/C(=C/[CH]CO)C=CC1=O`
- after:  `COC1=C/C(=C\[CH]CO)C=CC1=O`

**10. cpd03417 — Coproporphyrin I**  (3 rxns; InChIKey VCCUOZSDXVZCSK-JRHDEHKPSA-J → VCCUOZSDXVZCSK-UHFFFAOYSA-J)
- before: `CC1=C(CCC(=O)[O-])c2cc3nc(cc4[nH]c(cc5[nH]c(cc1n2)c(CCC(=O)[O-])c5C)c(CCC(=O)[O-])c4C)C(CCC(=O)[O-])=C3C`
- after:  `CC1=C(CCC(=O)[O-])c2cc3nc(cc4[nH]c(cc5[nH]c(cc1n2)c(CCC(=O)[O-])c5C)c(CCC(=O)[O-])c4C)C(CCC(=O)[O-])=C3C`

**11. cpd03617 — (1S,2R,3S)-2-Formyl-alpha-(hydroxymethylene)-3-methylcyclopentaneacetaldehyde**  (3 rxns; InChIKey PFGBAVLSGZLAAY-FXBDTBDDSA-N → PFGBAVLSGZLAAY-JYPFKNLRSA-N, stereo 3→4)
- before: `C[C@H]1CC[C@H](C(C=O)=CO)[C@@H]1C=O`
- after:  `C[C@H]1CC[C@H](/C(C=O)=C/O)[C@@H]1C=O`

**12. cpd23880 — 3,6-anhydro-alpha-D-galactopyranosyl-(1->3)-4-O-sulfo-D-galactose**  (3 rxns; InChIKey OGFGYTRGDDWUKC-PBFIDDPKSA-M → OGFGYTRGDDWUKC-YBFCQOOYSA-M, stereo 9→10)
- before: `O=S(=O)([O-])O[C@@H]1[C@H](O[C@H]2O[C@@H]3CO[C@@H]([C@H]3O)[C@H]2O)[C@@H](O)C(O)O[C@@H]1CO`
- after:  `O=S(=O)([O-])O[C@@H]1[C@H](O[C@H]2O[C@@H]3CO[C@@H]([C@H]3O)[C@H]2O)[C@@H](O)[C@H](O)O[C@@H]1CO`

**13. cpd23888 — neo-lambda-carrabiose**  (3 rxns; InChIKey RHDWPJKPWCDPBP-JZSVMVJISA-K → RHDWPJKPWCDPBP-WSWWMNSNSA-K, stereo 9→10)
- before: `O=S(=O)([O-])OC[C@H]1O[C@H](O[C@H]2[C@@H](O)[C@@H](CO)OC(O)[C@@H]2OS(=O)(=O)[O-])[C@H](OS(=O)(=O)[O-])[C@@H](O)[C@H]1O`
- after:  `O=S(=O)([O-])OC[C@H]1O[C@H](O[C@H]2[C@@H](O)[C@@H](CO)O[C@@H](O)[C@@H]2OS(=O)(=O)[O-])[C@H](OS(=O)(=O)[O-])[C@@H](O)[C@H]1O`

**14. cpd25687 — alpha-D-Galp2,6S2-(1->3)-beta-D-Galp2S-(1->4)-alpha-D-Galp2,6S2-(1->3)-D-Galp2S**  (3 rxns; InChIKey CFLYIPGVTFOQCI-SDNHUNQPSA-H → CFLYIPGVTFOQCI-SQESNNFJSA-H, stereo 19→20)
- before: `O=S(=O)([O-])OC[C@H]1O[C@H](O[C@H]2[C@@H](O)[C@@H](CO)O[C@@H](O[C@@H]3[C@H](O)[C@@H](OS(=O)(=O)[O-])C(O)O[C@@H]3COS(=O)(=O)[O-])[C@@H]2OS(=O)(=O)[O-])[C@H](OS(=O)(=O)[O-])[C@@H](O)[C@H]1O[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1OS(=O)(=O)[O-]`
- after:  `O=S(=O)([O-])OC[C@H]1O[C@H](O)[C@H](OS(=O)(=O)[O-])[C@@H](O)[C@H]1O[C@@H]1O[C@H](CO)[C@H](O)[C@H](O[C@H]2O[C@H](COS(=O)(=O)[O-])[C@H](O[C@@H]3O[C@H](CO)[C@H](O)[C@H](O)[C@H]3OS(=O)(=O)[O-])[C@H](O)[C@H]2OS(=O)(=O)[O-])[C@H]1OS(=O)(=O)[O-]`

**15. cpd31674 — quercetin-3-gentiobioside**  (3 rxns; InChIKey FDRQPMVGJOQVTL-BVYPKUFKSA-L → FDRQPMVGJOQVTL-DEFKTLOSSA-N, charge -2→0, stereo 9→10)
- before: `O=c1c(O[C@@H]2O[C@H](CO[C@@H]3O[C@H](CO)[C@@H](O)C(O)[C@H]3O)[C@@H](O)[C@H](O)[C@H]2O)c(-c2ccc(O)c(O)c2)oc2cc([O-])cc([O-])c12`
- after:  `O=c1c(O[C@@H]2O[C@H](CO[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)[C@@H](O)[C@H](O)[C@H]2O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12`

**16. cpd32789 — (S)-4-deoxygadusol**  (3 rxns; InChIKey ZONRIYAALKITGT-UHFFFAOYSA-M → ZONRIYAALKITGT-QMMMGPOBSA-N, charge -1→0, stereo 0→1)
- before: `COC1=C([O-])CC(O)(CO)CC1=O`
- after:  `COC1=C(O)C[C@@](O)(CO)CC1=O`

**17. cpd35178 — 4-acetylnivalenol**  (3 rxns; InChIKey XGCUCFKWVIWWNW-IZEYKREJSA-N → XGCUCFKWVIWWNW-CAYGJDLQSA-N, stereo 7→8)
- before: `CC(=O)O[C@@H]1[C@@H](O)[C@H]2O[C@@H]3C=C(C)C(=O)[C@@H](O)[C@]3(CO)[C@]1(C)C21CO1`
- after:  `CC(=O)O[C@@H]1[C@@H](O)[C@H]2O[C@@H]3C=C(C)C(=O)[C@@H](O)[C@]3(CO)[C@]1(C)[C@]21CO1`

**18. cpd23818 — (1,3)-beta-xylotriose**  (2 rxns; InChIKey WKSXOLNEYANUQB-RTAFPKLOSA-N → WKSXOLNEYANUQB-KWGKWNPMSA-N, stereo 11→12)
- before: `OC1OC[C@@H](O)[C@H](O[C@@H]2OC[C@@H](O)[C@H](O[C@@H]3OC[C@@H](O)[C@H](O)[C@H]3O)[C@H]2O)[C@H]1O`
- after:  `O[C@@H]1[C@@H](O)[C@H](O[C@@H]2[C@@H](O)[C@H](O[C@@H]3[C@@H](O)[C@H](O)OC[C@H]3O)OC[C@H]2O)OC[C@H]1O`

**19. cpd23819 — (1,3)-beta-xylotetraose**  (2 rxns; InChIKey YIMKISJWNJZSQI-GKDCLMHXSA-N → YIMKISJWNJZSQI-JTOFAPFQSA-N, stereo 15→16)
- before: `OC1OC[C@@H](O)[C@H](O[C@@H]2OC[C@@H](O)[C@H](O[C@@H]3OC[C@@H](O)[C@H](O[C@@H]4OC[C@@H](O)[C@H](O)[C@H]4O)[C@H]3O)[C@H]2O)[C@H]1O`
- after:  `O[C@@H]1[C@@H](O)[C@H](O[C@@H]2[C@@H](O)[C@H](O[C@@H]3[C@@H](O)[C@H](O[C@@H]4[C@@H](O)[C@H](O)OC[C@H]4O)OC[C@H]3O)OC[C@H]2O)OC[C@H]1O`

**20. cpd23881 — neocarratetraose 4-O-disulfate**  (2 rxns; InChIKey PPWBJIYNGRTPNY-SAEYTFRYSA-L → PPWBJIYNGRTPNY-OKPLKWSVSA-L, stereo 19→20)
- before: `O=S(=O)([O-])O[C@@H]1[C@H](O[C@H]2O[C@@H]3CO[C@@H]([C@H]3O)[C@H]2O)[C@@H](O)[C@H](O[C@@H]2[C@@H]3OC[C@H]2O[C@H](O[C@H]2[C@@H](OS(=O)(=O)[O-])[C@@H](CO)OC(O)[C@@H]2O)[C@@H]3O)O[C@@H]1CO`
- after:  `O=S(=O)([O-])O[C@@H]1[C@H](O[C@H]2O[C@@H]3CO[C@H]([C@H]2O)[C@H]3O[C@@H]2O[C@H](CO)[C@H](OS(=O)(=O)[O-])[C@H](O[C@H]3O[C@@H]4CO[C@@H]([C@H]4O)[C@H]3O)[C@H]2O)[C@@H](O)[C@H](O)O[C@@H]1CO`

**21. cpd25005 — quercetin 3,5-O-diglucoside**  (2 rxns; InChIKey YOXWSUCVDVXAMX-LOZUMUMPSA-M → YOXWSUCVDVXAMX-DEFKTLOSSA-N, charge -1→0, stereo 6→10)
- before: `O=c1c(OC2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)c(-c2ccc(O)c(O)c2)oc2cc([O-])cc(OC3O[C@H](CO)C(O)[C@H](O)C3O)c12`
- after:  `O=c1c(O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)c12`

**22. cpd32048 — 20-[alpha-L-arabinofuranosyl-(1->6)-beta-D-glucopyranosyloxy]dammar-24-ene-3beta,12-diol**  (2 rxns; InChIKey CJFGBCWGOQRURQ-NMUMGEHWSA-N → CJFGBCWGOQRURQ-JFJIKBJRSA-N, stereo 18→19)
- before: `CC(C)=CCCC(C)(O[C@@H]1O[C@H](CO[C@@H]2O[C@@H](CO)[C@H](O)[C@H]2O)[C@@H](O)[C@H](O)[C@H]1O)[C@H]1CC[C@]2(C)[C@@H]1[C@H](O)C[C@@H]1[C@@]3(C)CC[C@H](O)C(C)(C)[C@@H]3CC[C@]12C`
- after:  `CC(C)=CCC[C@](C)(O[C@@H]1O[C@H](CO[C@@H]2O[C@@H](CO)[C@H](O)[C@H]2O)[C@@H](O)[C@H](O)[C@H]1O)[C@H]1CC[C@]2(C)[C@@H]1[C@H](O)C[C@@H]1[C@@]3(C)CC[C@H](O)C(C)(C)[C@@H]3CC[C@]12C`

**23. cpd32552 — 4,5-dehydro-L-arginine**  (2 rxns; InChIKey FVSQXMGZDMFDEJ-BYPYZUCNSA-O → FVSQXMGZDMFDEJ-REJIZHKJSA-O, stereo 1→2)
- before: `NC(=[NH2+])NC=CC[C@H]([NH3+])C(=O)[O-]`
- after:  `NC(N)=[NH+]/C=C/C[C@H]([NH3+])C(=O)[O-]`

**24. cpd33745 — 7,8-dihydroxy-3,4,15-triacetoxyscirpenol**  (2 rxns; InChIKey URJKZBKWIWMEQI-LYTYWTHFSA-N → URJKZBKWIWMEQI-IGZVRCDWSA-N, stereo 8→9)
- before: `CC(=O)OC[C@]12[C@H](O)[C@H](O)C(C)=C[C@H]1O[C@@H]1[C@H](OC(C)=O)[C@@H](OC(C)=O)[C@@]2(C)C12CO2`
- after:  `CC(=O)OC[C@]12[C@H](O)[C@H](O)C(C)=C[C@H]1O[C@@H]1[C@H](OC(C)=O)[C@@H](OC(C)=O)[C@@]2(C)[C@]12CO2`

**25. cpd34103 — (20S)-ginsenoside C-Y**  (2 rxns; InChIKey YNBYFOIDLBTOMW-MOBYHOAZSA-N → YNBYFOIDLBTOMW-QHNUHGIDSA-N, stereo 18→19)
- before: `CC(C)=CCCC(C)(O[C@@H]1O[C@H](CO[C@@H]2OC[C@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@H](O)[C@H]1O)[C@H]1CC[C@]2(C)[C@@H]1[C@H](O)C[C@@H]1[C@@]3(C)CC[C@H](O)C(C)(C)[C@@H]3CC[C@]12C`
- after:  `CC(C)=CCC[C@](C)(O[C@@H]1O[C@H](CO[C@@H]2OC[C@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@H](O)[C@H]1O)[C@H]1CC[C@]2(C)[C@@H]1[C@H](O)C[C@@H]1[C@@]3(C)CC[C@H](O)C(C)(C)[C@@H]3CC[C@]12C`

**26. cpd34182 — (R)-4-deoxygadusol**  (2 rxns; InChIKey ZONRIYAALKITGT-UHFFFAOYSA-M → ZONRIYAALKITGT-MRVPVSSYSA-N, charge -1→0, stereo 0→1)
- before: `COC1=C([O-])CC(O)(CO)CC1=O`
- after:  `COC1=C(O)C[C@](O)(CO)CC1=O`

**27. cpd34487 — 12,13-epoxytrichothec-9-ene**  (2 rxns; InChIKey LZAJKCZTKKKZNT-HHHGZCDHSA-N → LZAJKCZTKKKZNT-QMIVOQANSA-N, stereo 4→5)
- before: `CC1=C[C@H]2O[C@@H]3CC[C@@](C)(C34CO4)[C@@]2(C)CC1`
- after:  `CC1=C[C@H]2O[C@@H]3CC[C@](C)([C@@]2(C)CC1)[C@]31CO1`

**28. cpd38198 — (20S)-20-(beta-D-glucopyranosyloxy)-dammara-24-ene-3beta-ol**  (2 rxns; InChIKey JQOUYGJYNQSCQP-SKWFVUTISA-N → JQOUYGJYNQSCQP-UZTNOJNJSA-N, stereo 13→14)
- before: `CC(C)=CCCC(C)(O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H]1CC[C@]2(C)[C@@H]1CC[C@@H]1[C@@]3(C)CC[C@H](O)C(C)(C)[C@@H]3CC[C@]12C`
- after:  `CC(C)=CCC[C@](C)(O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H]1CC[C@]2(C)[C@@H]1CC[C@@H]1[C@@]3(C)CC[C@H](O)C(C)(C)[C@@H]3CC[C@]12C`

**29. cpd39850 — (9Z)-11-{(3S)-3-[(2Z)-pent-2-en-1-yl]oxiran-2-ylidene}undec-9-enoate**  (2 rxns; InChIKey YZBZORUZOSCZRN-IRBKMFMFSA-M → YZBZORUZOSCZRN-YWHLHSFDSA-M, stereo 2→3)
- before: `CC/C=C\C[C@@H]1OC1=CC=CCCCCCCCC(=O)[O-]`
- after:  `CC/C=C\C[C@@H]1OC1=C/C=C\CCCCCCCC(=O)[O-]`

**30. cpd03386 — 2-Aminoethylphosphocholate**  (1 rxns; InChIKey BGUPNWPPECTFDP-RUUYYOGZSA-M → BGUPNWPPECTFDP-HZAMXZRMSA-M, stereo 10→11)
- before: `C[C@H](CCC(=O)NCCP(=O)([O-])O)[C@H]1CC[C@H]2[C@@H]3C(O)C[C@@H]4C[C@H](O)CC[C@]4(C)[C@H]3C[C@H](O)[C@]12C`
- after:  `C[C@H](CCC(=O)NCCP(=O)([O-])O)[C@H]1CC[C@H]2[C@@H]3[C@H](O)C[C@@H]4C[C@H](O)CC[C@]4(C)[C@H]3C[C@H](O)[C@]12C`

**31. cpd22207 — 7''-O-phosphohygromycin B**  (1 rxns; InChIKey DDJWTKQJOKVHBW-GNNWKVJVSA-O → DDJWTKQJOKVHBW-NZSRVPFOSA-O)
- before: `C[NH2+][C@H]1C[C@@H]([NH3+])[C@H](O)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@@H]3OC4(O[C@H]23)O[C@H]([C@H]([NH3+])COP(=O)([O-])[O-])[C@H](O)[C@H](O)[C@H]4O)[C@@H]1O`
- after:  `CN[C@H]1C[C@@H](N)[C@H](O)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@@H]3O[C@]4(O[C@H]23)O[C@H](C([NH3+])COP(=O)(O)O)[C@H](O)[C@H](O)[C@H]4O)[C@@H]1O`

**32. cpd23882 — iota-neocarratetraose sulfate**  (1 rxns; InChIKey QILNVJZKGZSQFJ-GCDDTEOQSA-J → QILNVJZKGZSQFJ-IXTYJRGDSA-J, stereo 19→20)
- before: `O=S(=O)([O-])O[C@@H]1[C@H](O[C@H]2O[C@@H]3CO[C@@H]([C@H]3O)[C@H]2OS(=O)(=O)[O-])[C@@H](O)[C@H](O[C@@H]2[C@@H]3OC[C@H]2O[C@H](O[C@H]2[C@@H](OS(=O)(=O)[O-])[C@@H](CO)OC(O)[C@@H]2O)[C@@H]3OS(=O)(=O)[O-])O[C@@H]1CO`
- after:  `O=S(=O)([O-])O[C@@H]1[C@H](O[C@H]2O[C@@H]3CO[C@@H]([C@H]3O[C@@H]3O[C@H](CO)[C@H](OS(=O)(=O)[O-])[C@H](O[C@H]4O[C@@H]5CO[C@@H]([C@H]5O)[C@H]4OS(=O)(=O)[O-])[C@H]3O)[C@H]2OS(=O)(=O)[O-])[C@@H](O)[C@H](O)O[C@@H]1CO`

**33. cpd23883 — iota-neocarrahexaose sulfate**  (1 rxns; InChIKey JKZMNXQJLUUIRT-ORGYOHQWSA-H → JKZMNXQJLUUIRT-YTZOPKLZSA-H, stereo 29→30)
- before: `O=S(=O)([O-])O[C@@H]1[C@H](O[C@H]2O[C@@H]3CO[C@@H]([C@H]3O[C@@H]3O[C@H](CO)[C@H](OS(=O)(=O)[O-])[C@H](O[C@H]4O[C@@H]5CO[C@@H]([C@H]5O)[C@H]4OS(=O)(=O)[O-])[C@H]3O)[C@H]2OS(=O)(=O)[O-])[C@@H](O)[C@H](O[C@@H]2[C@@H]3OC[C@H]2O[C@H](O[C@H]2[C@@H](OS(=O)(=O)[O-])[C@@H](CO)OC(O)[C@@H]2O)[C@@H]3OS(=O)(=O)[O-])O[C@@H]1CO`
- after:  `O=S(=O)([O-])O[C@@H]1[C@H](O[C@H]2O[C@@H]3CO[C@@H]([C@H]3O[C@@H]3O[C@H](CO)[C@H](OS(=O)(=O)[O-])[C@H](O[C@H]4O[C@@H]5CO[C@@H]([C@H]5O[C@@H]5O[C@H](CO)[C@H](OS(=O)(=O)[O-])[C@H](O[C@H]6O[C@@H]7CO[C@@H]([C@H]7O)[C@H]6OS(=O)(=O)[O-])[C@H]5O)[C@H]4OS(=O)(=O)[O-])[C@H]3O)[C@H]2OS(=O)(=O)[O-])[C@@H](O)[C@H](O)O[C@@H]1CO`

**34. cpd31237 — alpha-D-Gal-(1->3)-(alpha-L-Fuc-(1->2))-beta-D-Gal-(1->3)-alpha-D-GalNAc-(1->3)-alpha-D-GalNAc-diphospho-ditrans,octacis-undecaprenol**  (1 rxns; InChIKey UQYSJTDLRLYPNU-ABPCESQWSA-L → UQYSJTDLRLYPNU-OOXGFABVSA-L, stereo 34→35)
- before: `CC(=O)N[C@H]1[C@@H](OP(=O)([O-])OP(=O)([O-])OC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(\C)CC/C=C(\C)CCC=C(C)C)O[C@H](CO)[C@H](O)[C@@H]1O[C@H]1O[C@H](CO)[C@H](O)[C@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@H](OC3O[C@H](CO)[C@H](O)[C@H](O)[C@H]3O)[C@H]2O[C@@H]2O[C@@H](C)[C@@H](O)[C@@H](O)[C@@H]2O)[C@H]1NC(C)=O`
- after:  `CC(=O)N[C@H]1[C@@H](OP(=O)([O-])OP(=O)([O-])OC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(/C)CC/C=C(\C)CC/C=C(\C)CCC=C(C)C)O[C@H](CO)[C@H](O)[C@@H]1O[C@H]1O[C@H](CO)[C@H](O)[C@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@H](O[C@H]3O[C@H](CO)[C@H](O)[C@H](O)[C@H]3O)[C@H]2O[C@@H]2O[C@@H](C)[C@@H](O)[C@@H](O)[C@@H]2O)[C@H]1NC(C)=O`

**35. cpd37318 — (5S)-2,5-dihydroxy-5-(hydroxymethyl)-3-oxocyclohex-1-en-1-olate**  (1 rxns; InChIKey OWHGXOODGNBQRG-UHFFFAOYSA-M → OWHGXOODGNBQRG-SSDOTTSWSA-M, stereo 0→1)
- before: `O=C1CC(O)(CO)CC([O-])=C1O`
- after:  `O=C1C[C@@](O)(CO)CC(O)=C1[O-]`

**36. cpd38354 — (9S),10-epoxy-(10,12Z)-octadecadienoic acid(1-)**  (1 rxns; InChIKey LVVCDOSOKGYQFY-KRWDZBQOSA-M → LVVCDOSOKGYQFY-KKJXIILJSA-M, stereo 1→3)
- before: `CCCCCC=CC=C1O[C@H]1CCCCCCCC(=O)[O-]`
- after:  `CCCCC/C=C\C=C1/O[C@H]1CCCCCCCC(=O)[O-]`

**37. cpd39984 — 2-C-[(phosphonatooxy)methyl]-D-ribose**  (1 rxns; InChIKey RIOZVCDMYGAYCJ-QYRBDRAASA-L → RIOZVCDMYGAYCJ-HSUXUTPPSA-L, stereo 2→3)
- before: `O=CC(O)(COP(=O)([O-])[O-])[C@H](O)[C@H](O)CO`
- after:  `O=C[C@@](O)(COP(=O)([O-])[O-])[C@H](O)[C@H](O)CO`

**38. cpd46292 — 6'-sialyllactose**  (1 rxns; InChIKey TYALNJQZQRNQNQ-SKHIEIRWSA-M → TYALNJQZQRNQNQ-JLYOMPFMSA-M, stereo 14→16)
- before: `CC(=O)N[C@H]1[C@H]([C@H](O)[C@H](O)CO)OC(OC[C@H]2O[C@@H](O[C@H]3[C@H](O)[C@@H](O)C(O)O[C@@H]3CO)[C@H](O)[C@@H](O)[C@H]2O)(C(=O)[O-])C[C@@H]1O`
- after:  `CC(=O)N[C@H]1[C@H]([C@H](O)[C@H](O)CO)O[C@@](OC[C@H]2O[C@@H](O[C@H]3[C@H](O)[C@@H](O)[C@H](O)O[C@@H]3CO)[C@H](O)[C@@H](O)[C@H]2O)(C(=O)[O-])C[C@@H]1O`

**39. cpd03042 — Porphyrin**  (0 rxns; InChIKey JZRYQZJSTWVBBD-CEVVSZFKSA-N → JZRYQZJSTWVBBD-UHFFFAOYSA-N)
- before: `C1=Cc2cc3ccc(cc4ccc(cc5nc(cc1n2)C=C5)[nH]4)[nH]3`
- after:  `C1=Cc2cc3ccc(cc4ccc(cc5nc(cc1n2)C=C5)[nH]4)[nH]3`

**40. cpd04236 — Cefotetan**  (0 rxns; InChIKey SRZNHPXWXCNNDU-IXOPCIAXSA-L → SRZNHPXWXCNNDU-RHBCBLIFSA-L)
- before: `CO[C@@]1(NC(=O)C2SC(=C(C(N)=O)C(=O)[O-])S2)C(=O)N2C(C(=O)[O-])=C(CSc3nnnn3C)CS[C@@H]21`
- after:  `CO[C@@]1(NC(=O)C2SC(=C(C(N)=O)C(=O)[O-])S2)C(=O)N2C(C(=O)[O-])=C(CSc3nnnn3C)CS[C@@H]21`

**41. cpd05139 — Sermorelin**  (0 rxns; InChIKey WGWPRVFKDLAUQJ-UHFFFAOYSA-R → WGWPRVFKDLAUQJ-MITYVQBRSA-R, stereo 0→31)
- before: `CCC(C)C(NC(=O)C(C)NC(=O)C(CC(=O)[O-])NC(=O)C(C)NC(=O)C([NH3+])Cc1ccc(O)cc1)C(=O)NC(Cc1ccccc1)C(=O)NC(C(=O)NC(CC(N)=O)C(=O)NC(CO)C(=O)NC(Cc1ccc(O)cc1)C(=O)NC(CCCNC(N)=[NH2+])C(=O)NC(CCCC[NH3+])C(=O)NC(C(=O)NC(CC(C)C)C(=O)NCC(=O)NC(CCC(N)=O)C(=O)NC(CC(C)C)C(=O)NC(CO)C(=O)NC(C)C(=O)NC(CCCNC(N)=[NH2+])C(=O)NC(CCCC[NH3+])C(=O)NC(CC(C)C)C(=O)NC(CC(C)C)C(=O)NC(CCC(N)=O)C(=O)NC(CC(=O)[O-])C(=O)NC(C(=O)NC(CCSC)C(=O)NC(CO)C(=O)NC(CCCNC(N)=[NH2+])C(N)=O)C(C)CC)C(C)C)C(C)O`
- after:  `CC[C@H](C)[C@H](NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](C)NC(=O)[C@@H](N)Cc1ccc(O)cc1)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@H](C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H](CO)C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@@H](CCCNC(=N)[NH3+])C(=O)N[C@@H](CCCC[NH3+])C(=O)N[C@H](C(=O)N[C@@H](CC(C)C)C(=O)NCC(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CO)C(=O)N[C@@H](C)C(=O)N[C@@H](CCCNC(=N)[NH3+])C(=O)N[C@@H](CCCC[NH3+])C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H](C(=O)N[C@@H](CCSC)C(=O)N[C@@H](CO)C(=O)N[C@@H](CCCNC(=N)N)C(N)=O)[C@@H](C)CC)C(C)C)[C@@H](C)O`

**42. cpd08539 — Foscan**  (0 rxns; InChIKey LYPFDBRUNKHDGX-LWQDQPMZSA-N → LYPFDBRUNKHDGX-UHFFFAOYSA-N)
- before: `Oc1cccc(-c2c3nc(c(-c4cccc(O)c4)c4ccc([nH]4)c(-c4cccc(O)c4)c4nc(c(-c5cccc(O)c5)c5ccc2[nH]5)CC4)C=C3)c1`
- after:  `Oc1cccc(-c2c3nc(c(-c4cccc(O)c4)c4ccc([nH]4)c(-c4cccc(O)c4)c4nc(c(-c5cccc(O)c5)c5ccc2[nH]5)CC4)C=C3)c1`

**43. cpd09606 — 6,10-Diaza-3(1,3)8,(1,4)-dibenzena-1,5(1,4)- diquinolinacyclodecaphane**  (0 rxns; InChIKey HZWVJPDDZQOYGA-QUTRQNJUSA-P → HZWVJPDDZQOYGA-UHFFFAOYSA-P)
- before: `c1cc2cc(c1)C[n+]1ccc(c3ccccc31)NCc1ccc(cc1)CNc1cc[n+](c3ccccc13)C2`
- after:  `c1cc2cc(c1)C[n+]1ccc(c3ccccc31)NCc1ccc(cc1)CNc1cc[n+](c3ccccc13)C2`

**44. cpd09758 — Cefotetan disodium**  (0 rxns; InChIKey ZQQALMSFFARWPK-GLHLDKNHSA-L → ZQQALMSFFARWPK-ZTQQJVKJSA-L)
- before: `CO[C@@]1(NC(=O)C2SC(=C(C(N)=O)C(=O)[O-])S2)C(=O)N2C(C(=O)[O-])=C(CSc3nnnn3C)CS[C@@H]21.[Na+].[Na+]`
- after:  `CO[C@@]1(NC(=O)C2SC(=C(C(N)=O)C(=O)[O-])S2)C(=O)N2C(C(=O)[O-])=C(CSc3nnnn3C)CS[C@@H]21.[Na+].[Na+]`

**45. cpd17351 — Caulerpin**  (0 rxns; InChIKey PVWALMHWUOWPPA-HCOJNGAVSA-N → PVWALMHWUOWPPA-RPLHEXBESA-N)
- before: `COC(=O)/C1=C/c2c([nH]c3ccccc23)/C(C(=O)OC)=C\c2c1[nH]c1ccccc21`
- after:  `COC(=O)/C1=C/c2c([nH]c3ccccc23)/C(C(=O)OC)=C\c2c1[nH]c1ccccc21`

**46. cpd17664 — Chikusetsusaponin III**  (0 rxns; InChIKey ZICDJKZDHVLVOD-GREWDSJQSA-N → ZICDJKZDHVLVOD-HUGMCNGHSA-N, stereo 23→24)
- before: `CC(C)=CCCC(C)(O)[C@H]1CC[C@]2(C)[C@@H]1[C@H](O)C[C@@H]1[C@@]3(C)CC[C@H](O[C@@H]4O[C@H](CO[C@@H]5OC[C@@H](O)[C@H](O)[C@H]5O)[C@@H](O)[C@H](O)[C@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4O)C(C)(C)[C@@H]3CC[C@]12C`
- after:  `CC(C)=CCC[C@](C)(O)[C@H]1CC[C@]2(C)[C@@H]1[C@H](O)C[C@@H]1[C@@]3(C)CC[C@H](O[C@@H]4O[C@H](CO[C@@H]5OC[C@@H](O)[C@H](O)[C@H]5O)[C@@H](O)[C@H](O)[C@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4O)C(C)(C)[C@@H]3CC[C@]12C`

**47. cpd19283 — Protojujuboside A**  (0 rxns; InChIKey DAYKIQFQZVAKAA-PADDJQGGSA-N → DAYKIQFQZVAKAA-DJISZYPWSA-N, stereo 37→38)
- before: `CC(C)=C[C@@H](CC(C)(O)[C@H]1C(=O)C[C@]2(CO)[C@@H]1CC[C@@H]1[C@@]3(C)CC[C@H](O[C@@H]4OC[C@H](O)[C@H](O[C@@H]5O[C@H](CO[C@@H]6O[C@H](CO)[C@@H](O)[C@H](O)[C@H]6O)[C@@H](O)[C@H](O)[C@H]5O[C@@H]5OC[C@@H](O)[C@H](O)[C@H]5O)[C@H]4O[C@@H]4O[C@@H](C)[C@H](O)[C@@H](O)[C@H]4O)C(C)(C)[C@@H]3CC[C@]12C)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O`
- after:  `CC(C)=C[C@@H](C[C@](C)(O)[C@H]1C(=O)C[C@]2(CO)[C@@H]1CC[C@@H]1[C@@]3(C)CC[C@H](O[C@@H]4OC[C@H](O)[C@H](O[C@@H]5O[C@H](CO[C@@H]6O[C@H](CO)[C@@H](O)[C@H](O)[C@H]6O)[C@@H](O)[C@H](O)[C@H]5O[C@@H]5OC[C@@H](O)[C@H](O)[C@H]5O)[C@H]4O[C@@H]4O[C@@H](C)[C@H](O)[C@@H](O)[C@H]4O)C(C)(C)[C@@H]3CC[C@]12C)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O`

**48. cpd19284 — Protojujuboside B**  (0 rxns; InChIKey JWDBYPGADNNYPU-RNYOKPNMSA-N → JWDBYPGADNNYPU-KRACRRDZSA-N, stereo 32→33)
- before: `CC(C)=C[C@@H](CC(C)(O)[C@H]1C(=O)C[C@]2(CO)[C@@H]1CC[C@@H]1[C@@]3(C)CC[C@H](O[C@@H]4OC[C@H](O)[C@H](O[C@@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@H]5O[C@@H]5OC[C@@H](O)[C@H](O)[C@H]5O)[C@H]4O[C@@H]4O[C@@H](C)[C@H](O)[C@@H](O)[C@H]4O)C(C)(C)[C@@H]3CC[C@]12C)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O`
- after:  `CC(C)=C[C@@H](C[C@](C)(O)[C@H]1C(=O)C[C@]2(CO)[C@@H]1CC[C@@H]1[C@@]3(C)CC[C@H](O[C@@H]4OC[C@H](O)[C@H](O[C@@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@H]5O[C@@H]5OC[C@@H](O)[C@H](O)[C@H]5O)[C@H]4O[C@@H]4O[C@@H](C)[C@H](O)[C@@H](O)[C@H]4O)C(C)(C)[C@@H]3CC[C@]12C)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O`

**49. cpd19285 — Protojujuboside B1**  (0 rxns; InChIKey JWDBYPGADNNYPU-PNPNNWMTSA-N → JWDBYPGADNNYPU-HFNHGWQBSA-N, stereo 32→33)
- before: `CC(C)=C[C@@H](CC(C)(O)[C@H]1C(=O)C[C@]2(CO)[C@@H]1CC[C@@H]1[C@@]3(C)CC[C@H](O[C@@H]4OC[C@H](O)[C@H](O[C@@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@H]5O[C@@H]5OC[C@@H](O)[C@H](O)[C@H]5O)[C@H]4O[C@@H]4O[C@@H](C)[C@@H](O)[C@@H](O)[C@@H]4O)C(C)(C)[C@@H]3CC[C@]12C)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O`
- after:  `CC(C)=C[C@@H](C[C@](C)(O)[C@H]1C(=O)C[C@]2(CO)[C@@H]1CC[C@@H]1[C@@]3(C)CC[C@H](O[C@@H]4OC[C@H](O)[C@H](O[C@@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@H]5O[C@@H]5OC[C@@H](O)[C@H](O)[C@H]5O)[C@H]4O[C@@H]4O[C@@H](C)[C@@H](O)[C@@H](O)[C@@H]4O)C(C)(C)[C@@H]3CC[C@]12C)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O`

**50. cpd22704 — muricholate**  (0 rxns; InChIKey DKPMWHFRUGMUKF-QONOTNIVSA-M → DKPMWHFRUGMUKF-JDDNAIEOSA-M, stereo 4→9)
- before: `C[C@H](CCC(=O)[O-])C1CCC2C3C(O)C(O)C4C[C@H](O)CC[C@]4(C)C3CC[C@@]21C`
- after:  `C[C@H](CCC(=O)[O-])[C@H]1CC[C@H]2[C@@H]3C(O)C(O)[C@@H]4C[C@H](O)CC[C@]4(C)[C@H]3CC[C@]12C`

**51. cpd23878 — 3,6-anhydro-2-O-sulfo-alpha-D-galactopyranosyl-(1->3)-4-O-sulfo-D-galactose**  (0 rxns; InChIKey CTNFROXWAZXRKW-VDGMBKLFSA-L → CTNFROXWAZXRKW-VHBGUFLRSA-L, stereo 9→10)
- before: `O=S(=O)([O-])O[C@@H]1[C@H](O[C@H]2O[C@@H]3CO[C@@H]([C@H]3O)[C@H]2OS(=O)(=O)[O-])[C@@H](O)C(O)O[C@@H]1CO`
- after:  `O=S(=O)([O-])O[C@@H]1[C@H](O[C@H]2O[C@@H]3CO[C@@H]([C@H]3O)[C@H]2OS(=O)(=O)[O-])[C@@H](O)[C@H](O)O[C@@H]1CO`

**52. cpd23879 — 4-O-sulfo- beta-D-galactopyranosyl-(1->4)-3,6-anhydro-2-O-sulfo-D-galactose**  (0 rxns; InChIKey MWCCETPFHWJFDS-OSEFXCQPSA-L → MWCCETPFHWJFDS-XWTDTOCNSA-L, stereo 9→10)
- before: `O=S(=O)([O-])O[C@@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]2[C@@H]3OC[C@H]2OC(O)[C@@H]3OS(=O)(=O)[O-])O[C@@H]1CO`
- after:  `O=S(=O)([O-])O[C@@H]1[C@H]2OC[C@@H](O[C@@H]1O)[C@@H]2O[C@@H]1O[C@H](CO)[C@H](OS(=O)(=O)[O-])[C@H](O)[C@H]1O`

**53. cpd24275 — (11Z)-eicos-11-enoyl-CoA**  (0 rxns; InChIKey ZDRKXADSROCWCG-LURPFOQRSA-J → ZDRKXADSROCWCG-FVLDFCIYSA-J, stereo 5→6)
- before: `CCCCCCCC/C=C\CCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)([O-])[O-]`
- after:  `CCCCCCCC/C=C\CCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)([O-])[O-]`

**54. cpd24513 — sterigmatin**  (0 rxns; InChIKey CWYJYLXZMAUSNI-UHFFFAOYSA-M → CWYJYLXZMAUSNI-BWKAKNAASA-N, charge -1→0, stereo 0→2)
- before: `O=c1c2c(O)cccc2oc2cc3c(c([O-])c12)C1C=COC1O3`
- after:  `O=c1c2c(O)cccc2oc2cc3c(c(O)c12)[C@@H]1C=CO[C@@H]1O3`

**55. cpd26249 — S-2,3-dicarboxyaziridine**  (0 rxns; InChIKey IFCCPDAHQDGHMH-UHFFFAOYSA-M → IFCCPDAHQDGHMH-LWMBPPNESA-M, stereo 0→2)
- before: `O=C([O-])C1[NH2+]C1C(=O)[O-]`
- after:  `O=C([O-])[C@H]1N[C@@H]1C(=O)O`

**56. cpd26275 — (E)-alpha-monofluoromethyl-3,4-dehydroarginine**  (0 rxns; InChIKey UDWPEHOXJIIOEY-OWOJBTEDSA-O → UDWPEHOXJIIOEY-QOHHWTFISA-O, stereo 1→2)
- before: `NC(=[NH2+])NC/C=C/C([NH3+])(CF)C(=O)[O-]`
- after:  `NC(N)=[NH+]C/C=C/[C@@]([NH3+])(CF)C(=O)[O-]`

**57. cpd31660 — 2-O-(3'-hydroxy)phytanyl-3-O-phytanyl-sn-glycero-1-phosphoethanolamine**  (0 rxns; InChIKey VUKYXYWCFFKPMX-PLSZJGEGSA-N → VUKYXYWCFFKPMX-ARNVMZJBSA-N, stereo 6→7)
- before: `CC(C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCOC[C@@H](COP(=O)([O-])OCC[NH3+])OCCC(C)(O)CCC[C@H](C)CCC[C@H](C)CCCC(C)C`
- after:  `CC(C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCOC[C@@H](COP(=O)([O-])OCC[NH3+])OCC[C@](C)(O)CCC[C@H](C)CCC[C@H](C)CCCC(C)C`

**58. cpd34026 — manumycin B**  (0 rxns; InChIKey ZGICGDCGECBVTD-ULNOGCHZSA-M → ZGICGDCGECBVTD-LOWJZYKPSA-N, charge -1→0, stereo 7→8)
- before: `CCCCC(C)/C=C(\C)C(=O)NC1=C[C@@](O)(/C=C/C=C/C=C/C(=O)NC2=C([O-])CCC2=O)[C@@H]2O[C@@H]2C1=O`
- after:  `CCCC[C@@H](C)/C=C(\C)C(=O)NC1=C[C@@](O)(/C=C/C=C/C=C/C(=O)NC2=C(O)CCC2=O)[C@@H]2O[C@@H]2C1=O`

**59. cpd34426 — 2-O-(3'-hydroxy)phytanyl-3-O-phytanyl-sn-glycerol-1-phospho-3''-sn-glycerol**  (0 rxns; InChIKey BIVFMXPCWMONSY-XMRRTHDVSA-M → BIVFMXPCWMONSY-XHQCNPPWSA-M, stereo 7→8)
- before: `CC(C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCOC[C@@H](COP(=O)([O-])OC[C@H](O)CO)OCCC(C)(O)CCC[C@H](C)CCC[C@H](C)CCCC(C)C`
- after:  `CC(C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCOC[C@@H](COP(=O)([O-])OC[C@H](O)CO)OCC[C@](C)(O)CCC[C@H](C)CCC[C@H](C)CCCC(C)C`

**60. cpd36034 — 2-O-(3'-hydroxy)phytanyl-3-O-phytanyl-sn-glycero-1-phosphoserine**  (0 rxns; InChIKey YALRQOMGANIIQK-BOMFZTBMSA-M → YALRQOMGANIIQK-QQLPJPLTSA-M, stereo 7→8)
- before: `CC(C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCOC[C@@H](COP(=O)([O-])OC[C@H]([NH3+])C(=O)[O-])OCCC(C)(O)CCC[C@H](C)CCC[C@H](C)CCCC(C)C`
- after:  `CC(C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@@H](C)CCOC[C@@H](COP(=O)([O-])OC[C@H]([NH3+])C(=O)[O-])OCC[C@](C)(O)CCC[C@H](C)CCC[C@H](C)CCCC(C)C`

**61. cpd44412 — a chlorin**  (0 rxns; InChIKey UGADAJMDJZPKQX-CEVVSZFKSA-N → UGADAJMDJZPKQX-UHFFFAOYSA-N)
- before: `C1=Cc2cc3ccc(cc4nc(cc5ccc(cc1n2)[nH]5)CC4)[nH]3`
- after:  `C1=Cc2cc3ccc(cc4nc(cc5ccc(cc1n2)[nH]5)CC4)[nH]3`

**62. cpd44629 — normaritidine**  (0 rxns; InChIKey TZTBAJFJEZRQCV-UHFFFAOYSA-O → TZTBAJFJEZRQCV-RLCCDNCMSA-O, stereo 0→3)
- before: `COc1cc2c(cc1O)C[NH+]1CCC23C=CC(O)CC13`
- after:  `COc1cc2c(cc1O)C[NH+]1CC[C@@]23C=C[C@@H](O)C[C@H]13`

**63. cpd44677 — (5Z,8E,10E)-12-hydroxyheptadeca-5,8,10-trienoate**  (0 rxns; InChIKey KUKJHGXXZWHSBG-FQWDWNLFSA-M → KUKJHGXXZWHSBG-WBGSEQOASA-M, stereo 3→4)
- before: `CCCCCC(O)/C=C/C=C/C/C=C\CCCC(=O)[O-]`
- after:  `CCCCC[C@H](O)/C=C/C=C/C/C=C\CCCC(=O)[O-]`

**64. cpd46464 — ferrirhodin**  (0 rxns; InChIKey SEQVPEQKJMMKNO-SQLDAFDYSA-N → SEQVPEQKJMMKNO-RQCMKQRDSA-N, stereo 5→8)
- before: `CC(=CC1=[O+][Fe-3]2345O[C+](C=C(C)CCO)N(CCC[C@@H]6NC(=O)[C@H](CCCN(O2)C(C=C(C)CCO)=[O+]3)NC(=O)[C@H](CO)NC(=O)[C@H](CO)NC(=O)CNC(=O)[C@H](CCCN1O4)NC6=O)O5)CCO`
- after:  `C/C(=C/C(=O)N([O-])CCC[C@@H]1NC(=O)[C@H](CCCN([O-])C(=O)/C=C(/C)CCO)NC(=O)[C@H](CCCN([O-])C(=O)/C=C(/C)CCO)NC(=O)[C@H](CO)NC(=O)[C@H](CO)NC(=O)CNC1=O)CCO.[Fe+3]`

**65. cpd46646 — lomaiviticin B**  (0 rxns; InChIKey ZPEQGULUANUJNH-RJEFGVBDSA-N → ZPEQGULUANUJNH-RJEFGVBDSA-O, charge -1→0)
- before: `CC[C@]12O[C@]3(O)c4c5c(O)c6c(c([O-])c5c([N+]#N)[c-]4[C@H](O[C@H]4C[C@H](OC)[C@@H](O)[C@H](C)O4)[C@]4(CC)O[C@](O)(c5c7c(O)c8c(=O)ccc(=O)c=8c(=O)c-7c(N=N)c5[C@@H]1O[C@H]1C[C@H](OC)[C@@H](O)[C@H](C)O1)[C@H]2[C@@H]43)C(=O)C=CC6=O`
- after:  `CC[C@]12O[C@]3(O)c4c5c(O)c6c(c(O)c5c([N+]#N)[c-]4[C@H](O[C@H]4C[C@H](OC)[C@@H](O)[C@H](C)O4)[C@]4(CC)O[C@](O)(c5c7c(O)c8c(=O)ccc(=O)c=8c(=O)c-7c(N=N)c5[C@@H]1O[C@H]1C[C@H](OC)[C@@H](O)[C@H](C)O1)[C@H]2[C@@H]43)C(=O)C=CC6=O`

**66. cpd46783 — withanone**  (0 rxns; InChIKey FAZIYUIDUNHZRG-QUSHIGSFSA-N → FAZIYUIDUNHZRG-PCTWTJKKSA-N, stereo 10→11)
- before: `CC1=C(C)C(=O)O[C@@H]([C@@H](C)[C@@]2(O)CC[C@H]3[C@H]4C(CC[C@@]32C)[C@@]2(C)C(=O)C=CC[C@]2(O)[C@H]2O[C@@H]42)C1`
- after:  `CC1=C(C)C(=O)O[C@@H]([C@@H](C)[C@@]2(O)CC[C@H]3[C@@H]4[C@@H]5O[C@@H]5[C@@]5(O)CC=CC(=O)[C@]5(C)[C@H]4CC[C@@]32C)C1`

**67. cpd47191 — [(2-guanidinoethyl)sulfanyl]butanedioate**  (0 rxns; InChIKey VKVCLXDFOQQABP-UHFFFAOYSA-M → VKVCLXDFOQQABP-BYPYZUCNSA-M, stereo 0→1)
- before: `NC(=[NH2+])NCCSC(CC(=O)[O-])C(=O)[O-]`
- after:  `NC(N)=NCCS[C@@H](CC(=O)[O-])C(=O)O`

**68. cpd47672 — 2-deoxyparaherquamide A**  (0 rxns; InChIKey DYVLXWPZFQQUIU-BUWOMYMCSA-O → DYVLXWPZFQQUIU-WGNDVSEMSA-O, stereo 3→5)
- before: `CN1C(=O)C23C[C@@H]4C1(C[NH+]2CC[C@@]3(C)O)C[C@@]1(CNc2c1ccc1c2OC=CC(C)(C)O1)C4(C)C`
- after:  `CN1C(=O)[C@]23C[C@H]4C(C)(C)[C@]5(C[NH2+]c6c5ccc5c6OC=CC(C)(C)O5)C[C@@]41CN2CC[C@@]3(C)O`


### Usage context
The corrected compounds are peripheral/secondary-metabolism species, not core
metabolites; the heavily-used core compounds were already correct and untouched.
Reaction usage across the 56,012 reactions:

| Set | % orphan (0 rxns) | median rxns | mean | max |
|---|---|---|---|---|
| Corrected (68) | 44% | 1 | 2.1 | 42 |
| All structured compounds | 38% | 1 | 6.3 | 28,815 |
| Newly added (40) | 10% | 3 | 6.8 | 65 |

## 5. SMILES canonicalization (cosmetic, identity-preserving)
Beyond the 68 corrections, 415 compounds had their SMILES re-canonicalized
without any chemical change: all 415 had non-canonical input, 100% preserved
their InChIKey (no identity changed), and 407/415 are exact round-trip-stable
RDKit-canonical form (the 8 exceptions are metal/radical clusters, an RDKit
idempotency quirk, still valid).

## 6. Items left for manual review
- **cpd35694 (R)-2-(phosphomethyl)malate** — PubChem's (R) form matches the name
  while the stored value is (S); the automated pass conservatively kept the
  stored form because adopting it would invert a defined stereocenter. A manual
  flip to (R) is likely warranted.
- **cpd26222 DDF** — PubChem differs by a stereocenter inversion; the name has no
  R/S descriptor to adjudicate.
- **cpd23840 ABTS** — passes correction but raises one post-correction warning
  (its InChI fails RDKit re-parse due to radical mobile-H/charge representation);
  a pre-existing data quirk, not introduced by this run.

## 7. Reproducibility
`python Validate_PubChem_Structures.py --apply` (offline once the cache is warm,
~40 s). Outputs land in `Biochemistry/Structures/*_updated.txt` and
`Scripts/Structures/structure_review_output/`.
