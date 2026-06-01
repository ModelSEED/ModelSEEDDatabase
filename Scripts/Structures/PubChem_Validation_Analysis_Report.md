# PubChem Structure Validation — Analysis Report

**Date:** 2026-06-01
**Script reviewed:** `Scripts/Structures/Validate_PubChem_Structures.py`
**Reference:** dev-branch `Scripts/Structures/List_ModelSEED_Structures.py`
**Outputs:** `Biochemistry/Structures/{All,Unique}_ModelSEED_Structures_updated.txt`

---

## 1. Scope & method
Full read of the 3,128-line validation script, line-by-line comparison of the
unique-structure picker against the dev-branch `List_ModelSEED_Structures.py`,
RDKit-based chemical inspection of every corrected compound, two bug fixes, and
a clean offline regeneration of the output files (~40 s; all PubChem lookups
served from cache).

## 2. Does the pipeline achieve its goal?
Yes, and conservatively. Only three change classes are ever auto-applied —
consistency fixes (formula/charge recomputed from InChI), STEREO_DIFF (PubChem
adds stereo detail), and pKa-validated protonation — each behind multiple
guards. `MISMATCH` is never auto-applied; it is classified and routed to review.

## 3. Consistency with the dev-branch picker
`build_unique_file` / `_resolve_struct_conflict` faithfully reproduce the dev
algorithm (InChI>SMILE reference priority; multi-DB agreement → connectivity
grouping → source priority MetaCyc>KEGG>ChEBI>Rhea). The dev diff is otherwise
additive (`Pick_Reasons.txt`, formula-loader refactor). One divergence was
fixed: the ignore-list directory was relocated on dev
(`Biochemistry/Curation/ignores/`); `_load_ignored_structures()` now reads both
locations (verified identical 132-ID set, byte-for-byte same as dev).

## 4. Bugs found and fixed
- **Bug A — enantiomer inversions accepted.** `_check_stereo_compatibility`
  compared InChI `/t` (relative config) and `/b` layers but not `/m` (absolute
  config). Enantiomers share `/t` and differ only in `/m`, so full inversions
  passed as "compatible." Sourced from name-based PubChem matches, which often
  return the wrong stereoisomer. **Fix:** added an `/m`-layer enantiomer guard.
- **Bug B — metal complexes mangled.** The pKa/protonation path runs RDKit's
  `Uncharger` + re-protonation, which cannot handle metal–ligand dative bonds
  and rewrote cobalt coordination shells. **Fix:** `_has_metal` guards on all
  three pKa/protonation paths.

These removed 13 wrong/unsafe corrections (5 enantiomer inversions, 6 cobalt
corrinoids, scopolamine 6→5 stereo loss, a porphyrin 4→0 stereo loss). All are
now logged in `pubchem_rejected_corrections.tsv`.

## 5. Canonicalization assessment (the cosmetic SMILES changes)
- 415 pure-canonicalization compounds; **100% identity-preserved** (0 InChIKeys
  broken).
- **415/415** had non-canonical input SMILES (genuinely normalized).
- **407/415** are now exact round-trip-stable RDKit-canonical form; the 8
  exceptions are metal/radical clusters (RDKit idempotency quirk), still valid.

**Verdict:** yes, better — standardized, deduplicated, valid notation.

## 6. How common are the corrected compounds?
Reaction usage across 56,012 reactions. Corrected compounds are peripheral,
low-traffic (secondary metabolism); core metabolites were already correct and
untouched.

| Set | % orphan (0 rxns) | median rxns | mean | max |
|---|---|---|---|---|
| Corrected (68) | 44% | 1 | 2.1 | 42 |
| All structured compounds | 38% | 1 | 6.3 | 28,815 |
| Newly added (40) | 10% | 3 | 6.8 | 65 |

## 7. Final results (before → after fixes)

| Metric | Pre-fix | Post-fix (final) |
|---|---|---|
| Real identity (InChIKey) corrections | 81 | **68** |
| Bad/unsafe corrections | 13 | **0** |
| New compounds added | 42 | 40 |
| Compounds removed | 3 | 3 |
| SMILES canonicalized (cosmetic) | — | 632 |
| Unique file rows | 97,608 | 97,601 |
| Corrections-log provenance | April (stale cache) | 2026-06-01 (clean) |

By change type: 48 stereo gained · 13 stereo/representation change ·
6 protonation + stereo gain · 1 protonation change.

## 8. Complete list of the 68 corrected compounds
Sorted by reaction usage (most-used first).

| # | cpd_id | rxns | change type | strategy | name |
|---|---|---|---|---|---|
| 1 | cpd01476 | 42 | stereo/representation change | STEREO_DIFF | 7,12-diethenyl-3,8,13,17-tetramethylporphyrin-2,18-dipropanoate |
| 2 | cpd39967 | 9 | stereo gained | STEREO_DIFF | 3,7,11,15-tetramethylhexadeca-2,6,10,14-tetraen-1-ol |
| 3 | cpd02168 | 8 | stereo/representation change | STEREO_DIFF | (1R,2R,3R)-prephytoene diphosphate |
| 4 | cpd03913 | 6 | stereo gained | STEREO_DIFF | Hydrogenobyrinate a,c diamide |
| 5 | cpd23886 | 5 | stereo gained | STEREO_DIFF | neo-lambda-carratetraose |
| 6 | cpd23887 | 5 | stereo gained | STEREO_DIFF | neo-lambda-carrahexaose |
| 7 | cpd03832 | 4 | stereo gained | STEREO_DIFF | Hydrogenobyrinate |
| 8 | cpd23840 | 4 | stereo/representation change | STEREO_DIFF | 2,2'-azino-bis-(3-ethylbenzothiazoline-6-sulfonate) (ABTS) |
| 9 | cpd35693 | 4 | stereo/representation change | STEREO_DIFF | coniferyl alcohol radical |
| 10 | cpd03417 | 3 | stereo/representation change | STEREO_DIFF | Coproporphyrin I |
| 11 | cpd03617 | 3 | stereo gained | STEREO_DIFF | (1S,2R,3S)-2-Formyl-alpha-(hydroxymethylene)-3-methylcyclopropaneacetate |
| 12 | cpd23880 | 3 | stereo gained | STEREO_DIFF | 3,6-anhydro-alpha-D-galactopyranosyl-(1->3)-4-O-... |
| 13 | cpd23888 | 3 | stereo gained | STEREO_DIFF | neo-lambda-carrabiose |
| 14 | cpd25687 | 3 | stereo gained | STEREO_DIFF | alpha-D-Galp2,6S2-(1->3)-beta-D-Galp2S-(1->4)-alpha-... |
| 15 | cpd31674 | 3 | protonation + stereo gain | STEREO_DIFF | quercetin-3-gentiobioside |
| 16 | cpd32789 | 3 | protonation + stereo gain | STEREO_DIFF | (S)-4-deoxygadusol |
| 17 | cpd35178 | 3 | stereo gained | STEREO_DIFF | 4-acetylnivalenol |
| 18 | cpd23818 | 2 | stereo gained | STEREO_DIFF | (1,3)-beta-xylotriose |
| 19 | cpd23819 | 2 | stereo gained | STEREO_DIFF | (1,3)-beta-xylotetraose |
| 20 | cpd23881 | 2 | stereo gained | STEREO_DIFF | neocarratetraose 4-O-disulfate |
| 21 | cpd25005 | 2 | protonation + stereo gain | STEREO_DIFF | quercetin 3,5-O-diglucoside |
| 22 | cpd32048 | 2 | stereo gained | STEREO_DIFF | 20-[alpha-L-arabinofuranosyl-(1->6)-beta-D-glucosyl]... |
| 23 | cpd32552 | 2 | stereo gained | STEREO_DIFF | 4,5-dehydro-L-arginine |
| 24 | cpd33745 | 2 | stereo gained | STEREO_DIFF | 7,8-dihydroxy-3,4,15-triacetoxyscirpenol |
| 25 | cpd34103 | 2 | stereo gained | STEREO_DIFF | (20S)-ginsenoside C-Y |
| 26 | cpd34182 | 2 | protonation + stereo gain | STEREO_DIFF | (R)-4-deoxygadusol |
| 27 | cpd34487 | 2 | stereo gained | STEREO_DIFF | 12,13-epoxytrichothec-9-ene |
| 28 | cpd38198 | 2 | stereo gained | STEREO_DIFF | (20S)-20-(beta-D-glucopyranosyloxy)-dammara-24-ene-... |
| 29 | cpd39850 | 2 | stereo gained | STEREO_DIFF | (9Z)-11-{(3S)-3-[(2Z)-pent-2-en-1-yl]oxiran-2-yl}... |
| 30 | cpd03386 | 1 | stereo gained | STEREO_DIFF | 2-Aminoethylphosphocholate |
| 31 | cpd22207 | 1 | stereo/representation change | STEREO_DIFF | 7''-O-phosphohygromycin B |
| 32 | cpd23882 | 1 | stereo gained | STEREO_DIFF | iota-neocarratetraose sulfate |
| 33 | cpd23883 | 1 | stereo gained | STEREO_DIFF | iota-neocarrahexaose sulfate |
| 34 | cpd31237 | 1 | stereo gained | STEREO_DIFF | alpha-D-Gal-(1->3)-(alpha-L-Fuc-(1->2))-beta-D-Gal-... |
| 35 | cpd37318 | 1 | stereo gained | STEREO_DIFF | (5S)-2,5-dihydroxy-5-(hydroxymethyl)-3-oxocyclohexa-... |
| 36 | cpd38354 | 1 | stereo gained | STEREO_DIFF | (9S),10-epoxy-(10,12Z)-octadecadienoic acid(1-) |
| 37 | cpd39984 | 1 | stereo gained | STEREO_DIFF | 2-C-[(phosphonatooxy)methyl]-D-ribose |
| 38 | cpd46292 | 1 | stereo gained | STEREO_DIFF | 6'-sialyllactose |
| 39 | cpd03042 | 0 | stereo/representation change | STEREO_DIFF | Porphyrin |
| 40 | cpd04236 | 0 | stereo/representation change | STEREO_DIFF | Cefotetan |
| 41 | cpd05139 | 0 | stereo gained | STEREO_DIFF | Sermorelin |
| 42 | cpd08539 | 0 | stereo/representation change | STEREO_DIFF | Foscan |
| 43 | cpd09606 | 0 | stereo/representation change | STEREO_DIFF | 6,10-Diaza-3(1,3)8,(1,4)-dibenzena-1,5(1,4)-diquinolinacyclodecaphane |
| 44 | cpd09758 | 0 | stereo/representation change | STEREO_DIFF | Cefotetan disodium |
| 45 | cpd17351 | 0 | stereo/representation change | STEREO_DIFF | Caulerpin |
| 46 | cpd17664 | 0 | stereo gained | STEREO_DIFF | Chikusetsusaponin III |
| 47 | cpd19283 | 0 | stereo gained | STEREO_DIFF | Protojujuboside A |
| 48 | cpd19284 | 0 | stereo gained | STEREO_DIFF | Protojujuboside B |
| 49 | cpd19285 | 0 | stereo gained | STEREO_DIFF | Protojujuboside B1 |
| 50 | cpd22704 | 0 | stereo gained | STEREO_DIFF | muricholate |
| 51 | cpd23878 | 0 | stereo gained | STEREO_DIFF | 3,6-anhydro-2-O-sulfo-alpha-D-galactopyranosyl-(1->3)-... |
| 52 | cpd23879 | 0 | stereo gained | STEREO_DIFF | 4-O-sulfo-beta-D-galactopyranosyl-(1->4)-3,6-anhydro-... |
| 53 | cpd24275 | 0 | stereo gained | STEREO_DIFF | (11Z)-eicos-11-enoyl-CoA |
| 54 | cpd24513 | 0 | protonation + stereo gain | STEREO_DIFF | sterigmatin |
| 55 | cpd26249 | 0 | stereo gained | STEREO_DIFF | S-2,3-dicarboxyaziridine |
| 56 | cpd26275 | 0 | stereo gained | STEREO_DIFF | (E)-alpha-monofluoromethyl-3,4-dehydroarginine |
| 57 | cpd31660 | 0 | stereo gained | STEREO_DIFF | 2-O-(3'-hydroxy)phytanyl-3-O-phytanyl-sn-glycerol-1-phospho-... |
| 58 | cpd34026 | 0 | protonation + stereo gain | STEREO_DIFF | manumycin B |
| 59 | cpd34426 | 0 | stereo gained | STEREO_DIFF | 2-O-(3'-hydroxy)phytanyl-3-O-phytanyl-sn-glycerol-1-phospho-... |
| 60 | cpd36034 | 0 | stereo gained | STEREO_DIFF | 2-O-(3'-hydroxy)phytanyl-3-O-phytanyl-sn-glycerol-1-phospho-... |
| 61 | cpd44412 | 0 | stereo/representation change | STEREO_DIFF | a chlorin |
| 62 | cpd44629 | 0 | stereo gained | STEREO_DIFF | normaritidine |
| 63 | cpd44677 | 0 | stereo gained | STEREO_DIFF | (5Z,8E,10E)-12-hydroxyheptadeca-5,8,10-trienoate |
| 64 | cpd46464 | 0 | stereo gained | STEREO_DIFF | ferrirhodin |
| 65 | cpd46646 | 0 | protonation change | PKA_CORRECTION | lomaiviticin B |
| 66 | cpd46783 | 0 | stereo gained | STEREO_DIFF | withanone |
| 67 | cpd47191 | 0 | stereo gained | STEREO_DIFF | [(2-guanidinoethyl)sulfanyl]butanedioate |
| 68 | cpd47672 | 0 | stereo gained | STEREO_DIFF | 2-deoxyparaherquamide A |

## 9. Remaining manual-review items
- **cpd35694 (R)-2-(phosphomethyl)malate** — PubChem (R) matches the name; the
  stored value (S) is likely wrong. Conservatively kept stored; consider a
  manual flip to R.
- **cpd26222 DDF** — enantiomer rejected; no R/S descriptor in the name to
  adjudicate.
- **cpd23840 ABTS** — 1 post-correction consistency warning (InChI fails RDKit
  parse: radical mobile-H/charge mismatch). Pre-existing data quirk.

## 10. Safety / reproducibility
- Pre-fix outputs snapshotted to `/tmp/preFix_snapshot/`; the run's own backups
  remain in `Biochemistry/Structures/0_run_backup/`.
- Code edits in `Scripts/Structures/Validate_PubChem_Structures.py`
  (not yet committed — review with `git diff`).
- Reproduce: `python Validate_PubChem_Structures.py --apply` (offline, ~40 s).
