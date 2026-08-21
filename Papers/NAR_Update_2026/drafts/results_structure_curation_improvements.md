# Results — Structure-curation improvements (draft)

**Target section:** Results, "Structure-curation improvements" (see `PAPER_2026_SKELETON.md`).
**Guide references:** `PAPER_2026_GUIDE.md` §3, §4, §5.
**Status:** first draft. Numbers refresh once v2.0.0 is cut.

---

## Formula-conflict resolution

At the 2020 release, formula disagreement between source databases silently dropped the affected compound from `Unique_ModelSEED_Structures.txt`, cascading into a corresponding reaction becoming impossible to balance. The updated pipeline classifies each remaining conflict and routes it to one of three outcomes: auto-pick from the source variant that matches the historical `compound_*.json` formula (zero-disruption), curator-authored override with recorded rationale, or an entry in the mass-balance exclusion registry.

Applied to the priority-scope reaction set (~9,000 reactions used in v7.0 ModelSEEDTemplates and PlantSEED templates, comprising ~6,500 compounds), the pipeline resolved every previously-conflicting compound in three curation rounds during this release cycle:

- Round 1 (2026-07-03): 52 H-only conflicts auto-matched to the `compound_*.json` formula, plus two mass-balance-exclusion entries for cpd00766 (ascorbate radical) and cpd00849 (quercetin 3,3'-bissulfate) where no source formula matched.
- Round 2 (2026-07-04): 11 remaining priority-scope conflicts resolved — 9 auto-matched (spanning the alkyl-glycerol, sulfatide, Descarbamoylnovobiocin, and other classes), plus two additional mass-balance-exclusion entries for cpd27485 (Co(I)-corrinoid Fe-S protein — sole source structure includes the DMB nucleotide loop of vitamin B12, expanding formula from C₅₄H₈₀CoN₁₄O₉R₂ to C₆₇H₉₃CoN₁₆O₁₆PR₂) and cpd12396 (a pyruvate-kinase phosphoserine residue — closest source is one hydrogen off from any protonation state).
- Round 3 (2026-07-04): cpd00244 Ni²⁺ resolved via ChEBI 49786; the one remaining source variant (KEGG C00291 at Ni/0) was identified as elemental nickel mis-aliased to the ion, flagged for a separate alias cleanup.

After these rounds, the classifier reports **56 remaining formula-conflict compounds** across 143 source-variant rows (42 H-only, 12 heavy-atom skeleton, 1 element-set, 1 single-variant). Zero of these fall in the v7.0 priority template compound set. The priority-scope formula-conflict backlog is empty, and the balanced-reaction rate on the priority set held at 8,683 / 9,000 across all three PRs (no regressions).

## Mass-balance exclusion registry

The `mass_balance_excluded.tsv` registry — introduced in this release cycle — records four compounds whose stored formula in `compound_*.json` cannot be reproduced from any available source structure. For these compounds the picked structure's SMILES and InChIKey populate `compound_*.json` (so downstream atom-mapping tools have something to work with), while the historical formula and charge are preserved unchanged. As a consequence, the reactions in which these compounds appear balance against the historical formula rather than a formula computed from the picked structure. The registry is small and expected to remain so; each entry represents a species (a radical, an unusually deprotonated anion, or a protein-bound cofactor with database-inconsistent boundary conventions) whose canonical representation in ModelSEED predates its availability as a source structure.

## ACP formula override propagation

An earlier commit in this release cycle (2026-07-02, "Bulk-add 47 safe pantetheine-inclusive ACP formula overrides") added 47 rows to `acps_formula_charge.tsv` and regenerated `Unique_ModelSEED_Structures.txt`, but did not propagate the values into the `compound_*.json` shards. The Round 1 curation PR (2026-07-03) picked up those 47 propagations as a signed-off carryover, and the current release now carries a consistent pantetheine-inclusive ACP representation. This is the seed of the broader protein-carrier cofactor standardisation described below and detailed further in the Discussion.

## Placeholder for §5 and §8 sections

The following Results paragraphs will be added once the corresponding guide sections land:

- **Protein-carrier cofactor standardisation** (§5) — per-class reductions in false-positive mass-imbalance flags (acyl-ACP phosphopantetheine + biotinyl-BCCP + lipoyl-carrier + any extensions), broken out with before/after counts of reactions flagged as unbalanced. Expected to be a significant recovery of the 70% → 61% balanced-share regression noted in the growth statistics.
- **PubChem-alias mismap detection** — quantify the alias-mismap cases surfaced during this release cycle (MetaCyc CPD-15419 → clorobiocin, KEGG C00291 → elemental Ni, and others). Report as either a discovered-issues count or an actionable follow-up for the next release.
- **Atom mapping coverage** (§8) — fraction of priority-scope reactions with mappings delivered by the Nikoloski collaboration, with an illustrative use case in atom-tracking for C, P, S.

---

## Open loose ends flagged during drafting

- Numbers refresh once v2.0.0 is cut (per §1 decision).
- Confirm the alias mismap counts are worth surfacing as a Results sub-bullet vs a Discussion mention.
- The single-variant remaining-conflict entry (cpd where all source formulas are identical yet the row survived to `Formula_Conflicts.txt`) needs investigation before drafting the final Results text — likely an artefact of the picker's conflict-detection logic.
