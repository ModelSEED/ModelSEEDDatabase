# Results — Database connectivity + Ontological mappings (draft)

**Target sections:** Results, "Database connectivity" and "Ontological mappings for annotation comparison" (see `PAPER_2026_SKELETON.md`).
**Guide references:** `PAPER_2026_GUIDE.md` §10, §11.
**Status:** first pass. Both sections are refreshes of 2020 results with 2026 numbers.

---

## Database connectivity (whole-database FBA)

The 2020 release quantified the biochemistry database's own suitability for metabolic modeling by treating the entire database as a single reaction network and applying flux balance analysis. Two metrics were reported (Table 5 of the 2020 paper): (i) the number of reactions capable of carrying nonzero mass-balanced flux ("functional reactions") assuming every extracellular metabolite is available; (ii) the number of the 390 Biolog growth conditions under which the network is capable of producing biomass. The 2020 counts were 21,403 functional reactions and 355 / 390 Biolog conditions, comparable to the 2014 PlantSEED release but with substantially more mass-balanced reactions in the underlying database.

For the 2026 release we rerun the same analysis against the v2.0.0 tag. Refreshed numbers:

| | 2010 | 2014 | 2020 | 2026 [TBD] |
|---|---:|---:|---:|---:|
| Total reactions | 13,257 | 27,558 | 36,193 | 56,012 |
| Mass balanced | 10,263 | 17,264 | 25,457 | 34,370 |
| Reversible | 6,195 | 8,906 | 18,399 | [TBD] |
| Functional reactions | 8,505 | 21,917 | 21,403 | [TBD] |
| Functional Biolog conditions | 330 (85%) | 337 (86%) | 355 (91%) | [TBD] |

The absolute count of balanced reactions grew, but the functional-reaction fraction of the mass-balanced subset is the more informative metric [expected TBD range]. Under the Biolog conditions, [TBD interpretation whether the 355/390 figure holds, improves, or regresses].

## Ontological mappings for annotation comparison

The 2020 release described three ontology types built on top of the ModelSEED biochemistry: equivalent compound sets (mapping abstract or generic compounds to their structurally-specific representatives), lumped reaction sets (mapping a merged multi-step reaction to its component substeps), and context-specific reaction sets (mapping model-adapted reaction variants back to their standard representation). It also demonstrated the ontology's value by mapping the *E. coli* iJR904 model onto ModelSEED and showing that apparent-mismatch reactions dropped from 258 (unique to iJR904 by naive matching) to 159 (after ontology-guided reconciliation), primarily driven by isomer, phospholipid-abstraction, and fatty-acid-lumping equivalences.

For the 2026 release we refresh the same three ontology types and re-run the iJR904 comparison against the 2026 biochemistry.

Ontology growth:

| Ontology type | 2020 entries | 2026 entries [TBD] |
|---|---:|---:|
| Equivalent compound sets | [TBD 2020 count from source] | [TBD] |
| Lumped reaction sets | [TBD 2020 count] | [TBD] |
| Context-specific reaction sets | [TBD 2020 count] | [TBD] |

iJR904 reconciliation refresh: apparent-mismatch reactions dropped from [TBD 2026 baseline before ontology] to [TBD after ontology] via [TBD dominant equivalence categories]. The improvement over the 2020 result of 258 → 159 reflects both the added ontology entries and the improved compound-structure curation described earlier.

---

## Open loose ends flagged during drafting

- All numbers here are placeholders; the whole-DB FBA needs to be re-run against v2.0.0 before drafting can be finalized.
- The 2020 ontology-count numbers weren't given precisely in the 2020 paper's Table format; extract from the source manuscript's Results text or the ontology data files directly.
- The iJR904 comparison script may need locating / refactoring; find its 2020 counterpart before running for 2026.
