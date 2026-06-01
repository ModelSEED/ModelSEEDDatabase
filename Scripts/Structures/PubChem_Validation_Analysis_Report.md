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

Full before→after SMILES for all 68 are in
`Scripts/Structures/Corrected_Compounds_Before_After.tsv`. Side-by-side 2D
diagrams (before vs after) are in `structure_review_output/corrected_2d/`
(one PNG per compound, plus `Corrected_Compounds_2D.pdf` covering all 68).
The table below is sorted by reaction usage (most-used first); the **2D** column
links to each before/after diagram.

### Before → after (all 68)

| # | cpd_id | name | rxns | change | InChIKey before → after | 2D |
|---|---|---|---|---|---|---|
| 1 | cpd01476 | 7,12-diethenyl-3,8,13,17-tetramethylporphyrin-2,18-dipr | 42 | repr. | `ZCFFYALKHPIRKJ-UJJXFSCMSA-L` → `ZCFFYALKHPIRKJ-UHFFFAOYSA-L` | [png](structure_review_output/corrected_2d/cpd01476.png) |
| 2 | cpd39967 | 3,7,11,15-tetramethylhexadeca-2,6,10,14-tetraen-1-yl di | 9 | stereo 0→3 | `OINNEUNVOZHBOX-UHFFFAOYSA-L` → `OINNEUNVOZHBOX-QIRCYJPOSA-L` | [png](structure_review_output/corrected_2d/cpd39967.png) |
| 3 | cpd02168 | (1R,2R,3R)-prephytoene diphosphate | 8 | repr. | `RVCNKTPCHZNAAO-IMSLGMFESA-L` → `RVCNKTPCHZNAAO-UZDKSQMHSA-L` | [png](structure_review_output/corrected_2d/cpd02168.png) |
| 4 | cpd03913 | Hydrogenobyrinate a,c diamide | 6 | stereo 10→11 | `ZGGWTIPDUOTHRA-FVLIVPSMSA-J` → `ZGGWTIPDUOTHRA-QQHNNQQLSA-J` | [png](structure_review_output/corrected_2d/cpd03913.png) |
| 5 | cpd23886 | neo-lambda-carratetraose | 5 | stereo 19→20 | `SYQAURQZRLJKJT-GTBSAFNRSA-H` → `SYQAURQZRLJKJT-JYOYXRFGSA-H` | [png](structure_review_output/corrected_2d/cpd23886.png) |
| 6 | cpd23887 | neo-lambda-carrahexaose | 5 | stereo 29→30 | `XCAVJOBVNSMNME-WRUHAMNASA-E` → `XCAVJOBVNSMNME-HOQVOYLOSA-E` | [png](structure_review_output/corrected_2d/cpd23887.png) |
| 7 | cpd03832 | Hydrogenobyrinate | 4 | stereo 10→11 | `FJDBIDBCPIOVPG-FVLIVPSMSA-H` → `FJDBIDBCPIOVPG-QQHNNQQLSA-H` | [png](structure_review_output/corrected_2d/cpd03832.png) |
| 8 | cpd23840 | 2,2'-azino-bis-(3-ethylbenzothiazoline-6-sulfonate) rad | 4 | repr. | `HFNOFNDLBTVECO-YAFCTCPESA-L` → `HFNOFNDLBTVECO-UHFFFAOYSA-L` | [png](structure_review_output/corrected_2d/cpd23840.png) |
| 9 | cpd35693 | coniferyl alcohol radical | 4 | repr. | `ORAJWSYKRGVTDP-FPYGCLRLSA-N` → `ORAJWSYKRGVTDP-UHFFFAOYSA-N` | [png](structure_review_output/corrected_2d/cpd35693.png) |
| 10 | cpd03417 | Coproporphyrin I | 3 | repr. | `VCCUOZSDXVZCSK-JRHDEHKPSA-J` → `VCCUOZSDXVZCSK-UHFFFAOYSA-J` | [png](structure_review_output/corrected_2d/cpd03417.png) |
| 11 | cpd03617 | (1S,2R,3S)-2-Formyl-alpha-(hydroxymethylene)-3-methylcy | 3 | stereo 3→4 | `PFGBAVLSGZLAAY-FXBDTBDDSA-N` → `PFGBAVLSGZLAAY-JYPFKNLRSA-N` | [png](structure_review_output/corrected_2d/cpd03617.png) |
| 12 | cpd23880 | 3,6-anhydro-alpha-D-galactopyranosyl-(1->3)-4-O-sulfo-D | 3 | stereo 9→10 | `OGFGYTRGDDWUKC-PBFIDDPKSA-M` → `OGFGYTRGDDWUKC-YBFCQOOYSA-M` | [png](structure_review_output/corrected_2d/cpd23880.png) |
| 13 | cpd23888 | neo-lambda-carrabiose | 3 | stereo 9→10 | `RHDWPJKPWCDPBP-JZSVMVJISA-K` → `RHDWPJKPWCDPBP-WSWWMNSNSA-K` | [png](structure_review_output/corrected_2d/cpd23888.png) |
| 14 | cpd25687 | alpha-D-Galp2,6S2-(1->3)-beta-D-Galp2S-(1->4)-alpha-D-G | 3 | stereo 19→20 | `CFLYIPGVTFOQCI-SDNHUNQPSA-H` → `CFLYIPGVTFOQCI-SQESNNFJSA-H` | [png](structure_review_output/corrected_2d/cpd25687.png) |
| 15 | cpd31674 | quercetin-3-gentiobioside | 3 | charge -2→0; stereo 9→10 | `FDRQPMVGJOQVTL-BVYPKUFKSA-L` → `FDRQPMVGJOQVTL-DEFKTLOSSA-N` | [png](structure_review_output/corrected_2d/cpd31674.png) |
| 16 | cpd32789 | (S)-4-deoxygadusol | 3 | charge -1→0; stereo 0→1 | `ZONRIYAALKITGT-UHFFFAOYSA-M` → `ZONRIYAALKITGT-QMMMGPOBSA-N` | [png](structure_review_output/corrected_2d/cpd32789.png) |
| 17 | cpd35178 | 4-acetylnivalenol | 3 | stereo 7→8 | `XGCUCFKWVIWWNW-IZEYKREJSA-N` → `XGCUCFKWVIWWNW-CAYGJDLQSA-N` | [png](structure_review_output/corrected_2d/cpd35178.png) |
| 18 | cpd23818 | (1,3)-beta-xylotriose | 2 | stereo 11→12 | `WKSXOLNEYANUQB-RTAFPKLOSA-N` → `WKSXOLNEYANUQB-KWGKWNPMSA-N` | [png](structure_review_output/corrected_2d/cpd23818.png) |
| 19 | cpd23819 | (1,3)-beta-xylotetraose | 2 | stereo 15→16 | `YIMKISJWNJZSQI-GKDCLMHXSA-N` → `YIMKISJWNJZSQI-JTOFAPFQSA-N` | [png](structure_review_output/corrected_2d/cpd23819.png) |
| 20 | cpd23881 | neocarratetraose 4-O-disulfate | 2 | stereo 19→20 | `PPWBJIYNGRTPNY-SAEYTFRYSA-L` → `PPWBJIYNGRTPNY-OKPLKWSVSA-L` | [png](structure_review_output/corrected_2d/cpd23881.png) |
| 21 | cpd25005 | quercetin 3,5-O-diglucoside | 2 | charge -1→0; stereo 6→10 | `YOXWSUCVDVXAMX-LOZUMUMPSA-M` → `YOXWSUCVDVXAMX-DEFKTLOSSA-N` | [png](structure_review_output/corrected_2d/cpd25005.png) |
| 22 | cpd32048 | 20-[alpha-L-arabinofuranosyl-(1->6)-beta-D-glucopyranos | 2 | stereo 18→19 | `CJFGBCWGOQRURQ-NMUMGEHWSA-N` → `CJFGBCWGOQRURQ-JFJIKBJRSA-N` | [png](structure_review_output/corrected_2d/cpd32048.png) |
| 23 | cpd32552 | 4,5-dehydro-L-arginine | 2 | stereo 1→2 | `FVSQXMGZDMFDEJ-BYPYZUCNSA-O` → `FVSQXMGZDMFDEJ-REJIZHKJSA-O` | [png](structure_review_output/corrected_2d/cpd32552.png) |
| 24 | cpd33745 | 7,8-dihydroxy-3,4,15-triacetoxyscirpenol | 2 | stereo 8→9 | `URJKZBKWIWMEQI-LYTYWTHFSA-N` → `URJKZBKWIWMEQI-IGZVRCDWSA-N` | [png](structure_review_output/corrected_2d/cpd33745.png) |
| 25 | cpd34103 | (20S)-ginsenoside C-Y | 2 | stereo 18→19 | `YNBYFOIDLBTOMW-MOBYHOAZSA-N` → `YNBYFOIDLBTOMW-QHNUHGIDSA-N` | [png](structure_review_output/corrected_2d/cpd34103.png) |
| 26 | cpd34182 | (R)-4-deoxygadusol | 2 | charge -1→0; stereo 0→1 | `ZONRIYAALKITGT-UHFFFAOYSA-M` → `ZONRIYAALKITGT-MRVPVSSYSA-N` | [png](structure_review_output/corrected_2d/cpd34182.png) |
| 27 | cpd34487 | 12,13-epoxytrichothec-9-ene | 2 | stereo 4→5 | `LZAJKCZTKKKZNT-HHHGZCDHSA-N` → `LZAJKCZTKKKZNT-QMIVOQANSA-N` | [png](structure_review_output/corrected_2d/cpd34487.png) |
| 28 | cpd38198 | (20S)-20-(beta-D-glucopyranosyloxy)-dammara-24-ene-3bet | 2 | stereo 13→14 | `JQOUYGJYNQSCQP-SKWFVUTISA-N` → `JQOUYGJYNQSCQP-UZTNOJNJSA-N` | [png](structure_review_output/corrected_2d/cpd38198.png) |
| 29 | cpd39850 | (9Z)-11-{(3S)-3-[(2Z)-pent-2-en-1-yl]oxiran-2-ylidene}u | 2 | stereo 2→3 | `YZBZORUZOSCZRN-IRBKMFMFSA-M` → `YZBZORUZOSCZRN-YWHLHSFDSA-M` | [png](structure_review_output/corrected_2d/cpd39850.png) |
| 30 | cpd03386 | 2-Aminoethylphosphocholate | 1 | stereo 10→11 | `BGUPNWPPECTFDP-RUUYYOGZSA-M` → `BGUPNWPPECTFDP-HZAMXZRMSA-M` | [png](structure_review_output/corrected_2d/cpd03386.png) |
| 31 | cpd22207 | 7''-O-phosphohygromycin B | 1 | repr. | `DDJWTKQJOKVHBW-GNNWKVJVSA-O` → `DDJWTKQJOKVHBW-NZSRVPFOSA-O` | [png](structure_review_output/corrected_2d/cpd22207.png) |
| 32 | cpd23882 | iota-neocarratetraose sulfate | 1 | stereo 19→20 | `QILNVJZKGZSQFJ-GCDDTEOQSA-J` → `QILNVJZKGZSQFJ-IXTYJRGDSA-J` | [png](structure_review_output/corrected_2d/cpd23882.png) |
| 33 | cpd23883 | iota-neocarrahexaose sulfate | 1 | stereo 29→30 | `JKZMNXQJLUUIRT-ORGYOHQWSA-H` → `JKZMNXQJLUUIRT-YTZOPKLZSA-H` | [png](structure_review_output/corrected_2d/cpd23883.png) |
| 34 | cpd31237 | alpha-D-Gal-(1->3)-(alpha-L-Fuc-(1->2))-beta-D-Gal-(1-> | 1 | stereo 34→35 | `UQYSJTDLRLYPNU-ABPCESQWSA-L` → `UQYSJTDLRLYPNU-OOXGFABVSA-L` | [png](structure_review_output/corrected_2d/cpd31237.png) |
| 35 | cpd37318 | (5S)-2,5-dihydroxy-5-(hydroxymethyl)-3-oxocyclohex-1-en | 1 | stereo 0→1 | `OWHGXOODGNBQRG-UHFFFAOYSA-M` → `OWHGXOODGNBQRG-SSDOTTSWSA-M` | [png](structure_review_output/corrected_2d/cpd37318.png) |
| 36 | cpd38354 | (9S),10-epoxy-(10,12Z)-octadecadienoic acid(1-) | 1 | stereo 1→3 | `LVVCDOSOKGYQFY-KRWDZBQOSA-M` → `LVVCDOSOKGYQFY-KKJXIILJSA-M` | [png](structure_review_output/corrected_2d/cpd38354.png) |
| 37 | cpd39984 | 2-C-[(phosphonatooxy)methyl]-D-ribose | 1 | stereo 2→3 | `RIOZVCDMYGAYCJ-QYRBDRAASA-L` → `RIOZVCDMYGAYCJ-HSUXUTPPSA-L` | [png](structure_review_output/corrected_2d/cpd39984.png) |
| 38 | cpd46292 | 6'-sialyllactose | 1 | stereo 14→16 | `TYALNJQZQRNQNQ-SKHIEIRWSA-M` → `TYALNJQZQRNQNQ-JLYOMPFMSA-M` | [png](structure_review_output/corrected_2d/cpd46292.png) |
| 39 | cpd03042 | Porphyrin | 0 | repr. | `JZRYQZJSTWVBBD-CEVVSZFKSA-N` → `JZRYQZJSTWVBBD-UHFFFAOYSA-N` | [png](structure_review_output/corrected_2d/cpd03042.png) |
| 40 | cpd04236 | Cefotetan | 0 | repr. | `SRZNHPXWXCNNDU-IXOPCIAXSA-L` → `SRZNHPXWXCNNDU-RHBCBLIFSA-L` | [png](structure_review_output/corrected_2d/cpd04236.png) |
| 41 | cpd05139 | Sermorelin | 0 | stereo 0→31 | `WGWPRVFKDLAUQJ-UHFFFAOYSA-R` → `WGWPRVFKDLAUQJ-MITYVQBRSA-R` | [png](structure_review_output/corrected_2d/cpd05139.png) |
| 42 | cpd08539 | Foscan | 0 | repr. | `LYPFDBRUNKHDGX-LWQDQPMZSA-N` → `LYPFDBRUNKHDGX-UHFFFAOYSA-N` | [png](structure_review_output/corrected_2d/cpd08539.png) |
| 43 | cpd09606 | 6,10-Diaza-3(1,3)8,(1,4)-dibenzena-1,5(1,4)- diquinolin | 0 | repr. | `HZWVJPDDZQOYGA-QUTRQNJUSA-P` → `HZWVJPDDZQOYGA-UHFFFAOYSA-P` | [png](structure_review_output/corrected_2d/cpd09606.png) |
| 44 | cpd09758 | Cefotetan disodium | 0 | repr. | `ZQQALMSFFARWPK-GLHLDKNHSA-L` → `ZQQALMSFFARWPK-ZTQQJVKJSA-L` | [png](structure_review_output/corrected_2d/cpd09758.png) |
| 45 | cpd17351 | Caulerpin | 0 | repr. | `PVWALMHWUOWPPA-HCOJNGAVSA-N` → `PVWALMHWUOWPPA-RPLHEXBESA-N` | [png](structure_review_output/corrected_2d/cpd17351.png) |
| 46 | cpd17664 | Chikusetsusaponin III | 0 | stereo 23→24 | `ZICDJKZDHVLVOD-GREWDSJQSA-N` → `ZICDJKZDHVLVOD-HUGMCNGHSA-N` | [png](structure_review_output/corrected_2d/cpd17664.png) |
| 47 | cpd19283 | Protojujuboside A | 0 | stereo 37→38 | `DAYKIQFQZVAKAA-PADDJQGGSA-N` → `DAYKIQFQZVAKAA-DJISZYPWSA-N` | [png](structure_review_output/corrected_2d/cpd19283.png) |
| 48 | cpd19284 | Protojujuboside B | 0 | stereo 32→33 | `JWDBYPGADNNYPU-RNYOKPNMSA-N` → `JWDBYPGADNNYPU-KRACRRDZSA-N` | [png](structure_review_output/corrected_2d/cpd19284.png) |
| 49 | cpd19285 | Protojujuboside B1 | 0 | stereo 32→33 | `JWDBYPGADNNYPU-PNPNNWMTSA-N` → `JWDBYPGADNNYPU-HFNHGWQBSA-N` | [png](structure_review_output/corrected_2d/cpd19285.png) |
| 50 | cpd22704 | muricholate | 0 | stereo 4→9 | `DKPMWHFRUGMUKF-QONOTNIVSA-M` → `DKPMWHFRUGMUKF-JDDNAIEOSA-M` | [png](structure_review_output/corrected_2d/cpd22704.png) |
| 51 | cpd23878 | 3,6-anhydro-2-O-sulfo-alpha-D-galactopyranosyl-(1->3)-4 | 0 | stereo 9→10 | `CTNFROXWAZXRKW-VDGMBKLFSA-L` → `CTNFROXWAZXRKW-VHBGUFLRSA-L` | [png](structure_review_output/corrected_2d/cpd23878.png) |
| 52 | cpd23879 | 4-O-sulfo- beta-D-galactopyranosyl-(1->4)-3,6-anhydro-2 | 0 | stereo 9→10 | `MWCCETPFHWJFDS-OSEFXCQPSA-L` → `MWCCETPFHWJFDS-XWTDTOCNSA-L` | [png](structure_review_output/corrected_2d/cpd23879.png) |
| 53 | cpd24275 | (11Z)-eicos-11-enoyl-CoA | 0 | stereo 5→6 | `ZDRKXADSROCWCG-LURPFOQRSA-J` → `ZDRKXADSROCWCG-FVLDFCIYSA-J` | [png](structure_review_output/corrected_2d/cpd24275.png) |
| 54 | cpd24513 | sterigmatin | 0 | charge -1→0; stereo 0→2 | `CWYJYLXZMAUSNI-UHFFFAOYSA-M` → `CWYJYLXZMAUSNI-BWKAKNAASA-N` | [png](structure_review_output/corrected_2d/cpd24513.png) |
| 55 | cpd26249 | S-2,3-dicarboxyaziridine | 0 | stereo 0→2 | `IFCCPDAHQDGHMH-UHFFFAOYSA-M` → `IFCCPDAHQDGHMH-LWMBPPNESA-M` | [png](structure_review_output/corrected_2d/cpd26249.png) |
| 56 | cpd26275 | (E)-alpha-monofluoromethyl-3,4-dehydroarginine | 0 | stereo 1→2 | `UDWPEHOXJIIOEY-OWOJBTEDSA-O` → `UDWPEHOXJIIOEY-QOHHWTFISA-O` | [png](structure_review_output/corrected_2d/cpd26275.png) |
| 57 | cpd31660 | 2-O-(3'-hydroxy)phytanyl-3-O-phytanyl-sn-glycero-1-phos | 0 | stereo 6→7 | `VUKYXYWCFFKPMX-PLSZJGEGSA-N` → `VUKYXYWCFFKPMX-ARNVMZJBSA-N` | [png](structure_review_output/corrected_2d/cpd31660.png) |
| 58 | cpd34026 | manumycin B | 0 | charge -1→0; stereo 7→8 | `ZGICGDCGECBVTD-ULNOGCHZSA-M` → `ZGICGDCGECBVTD-LOWJZYKPSA-N` | [png](structure_review_output/corrected_2d/cpd34026.png) |
| 59 | cpd34426 | 2-O-(3'-hydroxy)phytanyl-3-O-phytanyl-sn-glycerol-1-pho | 0 | stereo 7→8 | `BIVFMXPCWMONSY-XMRRTHDVSA-M` → `BIVFMXPCWMONSY-XHQCNPPWSA-M` | [png](structure_review_output/corrected_2d/cpd34426.png) |
| 60 | cpd36034 | 2-O-(3'-hydroxy)phytanyl-3-O-phytanyl-sn-glycero-1-phos | 0 | stereo 7→8 | `YALRQOMGANIIQK-BOMFZTBMSA-M` → `YALRQOMGANIIQK-QQLPJPLTSA-M` | [png](structure_review_output/corrected_2d/cpd36034.png) |
| 61 | cpd44412 | a chlorin | 0 | repr. | `UGADAJMDJZPKQX-CEVVSZFKSA-N` → `UGADAJMDJZPKQX-UHFFFAOYSA-N` | [png](structure_review_output/corrected_2d/cpd44412.png) |
| 62 | cpd44629 | normaritidine | 0 | stereo 0→3 | `TZTBAJFJEZRQCV-UHFFFAOYSA-O` → `TZTBAJFJEZRQCV-RLCCDNCMSA-O` | [png](structure_review_output/corrected_2d/cpd44629.png) |
| 63 | cpd44677 | (5Z,8E,10E)-12-hydroxyheptadeca-5,8,10-trienoate | 0 | stereo 3→4 | `KUKJHGXXZWHSBG-FQWDWNLFSA-M` → `KUKJHGXXZWHSBG-WBGSEQOASA-M` | [png](structure_review_output/corrected_2d/cpd44677.png) |
| 64 | cpd46464 | ferrirhodin | 0 | stereo 5→8 | `SEQVPEQKJMMKNO-SQLDAFDYSA-N` → `SEQVPEQKJMMKNO-RQCMKQRDSA-N` | [png](structure_review_output/corrected_2d/cpd46464.png) |
| 65 | cpd46646 | lomaiviticin B | 0 | charge -1→0 | `ZPEQGULUANUJNH-RJEFGVBDSA-N` → `ZPEQGULUANUJNH-RJEFGVBDSA-O` | [png](structure_review_output/corrected_2d/cpd46646.png) |
| 66 | cpd46783 | withanone | 0 | stereo 10→11 | `FAZIYUIDUNHZRG-QUSHIGSFSA-N` → `FAZIYUIDUNHZRG-PCTWTJKKSA-N` | [png](structure_review_output/corrected_2d/cpd46783.png) |
| 67 | cpd47191 | [(2-guanidinoethyl)sulfanyl]butanedioate | 0 | stereo 0→1 | `VKVCLXDFOQQABP-UHFFFAOYSA-M` → `VKVCLXDFOQQABP-BYPYZUCNSA-M` | [png](structure_review_output/corrected_2d/cpd47191.png) |
| 68 | cpd47672 | 2-deoxyparaherquamide A | 0 | stereo 3→5 | `DYVLXWPZFQQUIU-BUWOMYMCSA-O` → `DYVLXWPZFQQUIU-WGNDVSEMSA-O` | [png](structure_review_output/corrected_2d/cpd47672.png) |

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
