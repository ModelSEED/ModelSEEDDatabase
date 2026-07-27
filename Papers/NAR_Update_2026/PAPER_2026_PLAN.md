# ModelSEED Biochemistry Database — 2026 Update Paper Plan

**Draft plan • 2026-07-27 • prepared for team review**

## Snapshot

- **Prior paper:** Seaver et al. 2020, *The ModelSEED Biochemistry Database…*, Nucleic Acids Research 49(D1):D575-D588. DOI [10.1093/nar/gkaa746](https://doi.org/10.1093/nar/gkaa746). (Note: the DOI [10.1093/nar/gkaa1143](https://doi.org/10.1093/nar/gkaa1143) is the erratum, not the main paper. Cite gkaa746.)
- **Target venue:** NAR Database Issue (same as 2020 — keeps citation lineage; ~4000-word cap).
- **Draft window:** TBD by team.
- **Novel empirical hook:** systematic study of reaction-direction heuristic sensitivity across the ModelSEED v2 draft-model corpus ([10.1101/2023.10.04.556561](https://doi.org/10.1101/2023.10.04.556561)).

## What the 2020 paper committed to (baseline for the update)

According to PubMed ([10.1093/nar/gkaa746](https://doi.org/10.1093/nar/gkaa746)), the 2020 paper's four positioning claims were:

1. **Rosetta-Stone integration** of KEGG + MetaCyc + BioCyc + 34 published models + BiGG/MetaNetX/Rhea aliases.
2. **Structure-first curation** — ChemAxon Marvin protonation @ pH 7; RDKit + OpenBabel for formula/charge derivation; mass-balance ledger via a `status` field.
3. **Thermodynamics** — eQuilibrator (primary) + Group Contribution (fallback), used to predict reaction reversibility.
4. **Community + Provenance** — GitHub PR workflow with Travis CI, quarterly releases.

**2020 baseline numbers to update:**

| | 2020 snapshot |
|---|---:|
| Compounds | 33,958 |
| Compounds with structures | 28,120 (83%) |
| Reactions | 36,193 |
| Mass+charge balanced | 25,457 (70%) |
| Functional in whole-DB FBA | 21,403 |
| Biolog conditions supported | 355 / 390 (91%) |
| Compounds with eQuilibrator ΔfG′ | 17,510 |
| Reactions with accepted ΔrG′ | 13,298 |

We'll refresh every row and add a `2026` column when the update numbers are in.

## Proposed section outline

### 1. Introduction

Frame: six years of community-driven biochemistry curation. Motivate: rise of atom mapping in metabolic-model applications, ML-based thermodynamics methods maturing, growing need for reaction-direction transparency and its downstream model impact.

### 2. Materials & Methods

- **Structure-curation pipeline** — PubChem validation stage (Ray's stereo-loss guard); RDKit canonicalization (one-time mass recanonicalization); curator override system (`Biochemistry/Curation/overrides/structure_picks/<curator>.tsv`); mass-balance exclusion mechanism (new file schema: `cpd_id, name, reason, date, curator`); per-source SMILES preserved; PubChem-alias mismap detection surfacing.
- **Formula-conflict resolver** — H-only vs skeleton conflicts; auto-match to `compound_*.json` for zero-disruption picks; explicit exclusion when no source matches (e.g., ascorbate radical, quercetin bissulfate).
- **ACP formula override framework** — pantetheine-inclusive formula overrides in `acps_formula_charge.tsv`.
- **Multi-source thermodynamics** — rebuild of Group Contribution from MFAToolkit; dGpredictor retrained on ModelSEED; eQuilibrator refresh; OpenTECR integration. Per-reaction-class heuristics applied globally.
- **Atom mapping** — collaboration with the Nikoloski lab; describe methodology and coverage-to-date.
- **PR-time validation** — GitHub Actions replacing Travis (Travis is retired). New curated entries validated on submission via a workflow that runs the curation self-check.
- **Packaging** — PyPi distribution with a documented API.

### 3. Results

- **Growth statistics 2020→2026** (mirror Tables 2-3 of the 2020 paper).
- **Structure-curation improvements** — number of compounds with newly-picked structures, mass-balance exclusions applied, PubChem pipeline yield.
- **Multi-source thermodynamics comparison** — heat map / correlation among the four ΔG sources; per-class heuristic direction assignments.
- **NEW: reaction-direction heuristic sensitivity study** — systematic sweep over the ModelSEED v2 draft-model corpus ([10.1101/2023.10.04.556561](https://doi.org/10.1101/2023.10.04.556561)). Measures per-model: predicted growth rate, essential-gene set overlap, mass-balance survival rate, feasibility under Biolog conditions, as we swap between direction sources (eQuilibrator vs GC vs dGpredictor vs heuristic overlay). This is the empirical hook that lifts the paper above a numbers-refresh.
- **Atom mapping coverage** — fraction of priority-scope reactions with mappings; example use cases in atom-tracking (C, P, S).

### 4. Discussion

- Where 6 years of community curation has taken the DB.
- Open problems: fragment-aware pKa/pKb (per Priorities list); EC number refresh beyond ExPASy name-matching (bound to have false / missing associations); obsolescence-leakage audit (obsolete compounds/reactions may be accidentally used); direction-field removal from the reaction schema now that thermoreversibility + template capture direction.

### 5. Data & software availability

- GitHub repo (dev/master branches, PR workflow, GitHub Actions CI).
- PyPi package + API.
- Updated website (with atom-mapping surfaced).
- Nikoloski-lab atom-mapping endpoint (as delivered).

## Empirical study — design sketch (reaction-direction sensitivity)

**Corpus:** ModelSEED v2 draft-model corpus published in [10.1101/2023.10.04.556561](https://doi.org/10.1101/2023.10.04.556561). Use the models as-published (no manual curation) to keep the study systematic.

**Direction sources to compare** (each is a full reversibility labeling over the ModelSEED reactions):

1. eQuilibrator-derived ΔrG′ + reversibility rules (this is the 2020 default).
2. Rebuilt Group Contribution ΔrG′ + reversibility rules.
3. dGpredictor (retrained on ModelSEED) ΔrG′ + reversibility rules.
4. Heuristic overlay on each of the above (per-reaction-class rules from the meeting notes — e.g., MetaCyc fatty-acyl ACP handling, sugar-isomer conventions).

**Per-model metrics:**

- Predicted growth rate under the model's original medium.
- Essential-gene set (single-gene knockout FBA) — Jaccard overlap across direction sources.
- Mass-balance survival rate (reactions removed as thermodynamically infeasible under each direction source).
- Feasibility of biomass production under the 390 Biolog conditions (mirroring the 2020 whole-DB FBA analysis).

**Reporting:** matrix of `(model × direction source × metric)` deltas; identify direction sources that systematically break vs preserve biology across the corpus.

**Explicitly out of scope:** ANIME-based heuristic-break detection. Even though ANIME appeared in the July 2024 meeting notes as a candidate, it is not part of this paper — save for a future methods paper if warranted.

## Open decisions to close with the team

Before the first full draft, the team should reach agreement on:

1. **Which direction sources make the final cut** for the empirical study — dropping any of the four saves complexity but weakens the systematic-sweep claim.
2. **Whether OpenTECR gets its own Methods paragraph or is folded into the eQuilibrator refresh.**
3. **How the Nikoloski atom-mapping work is credited** — co-authorship, collaboration acknowledgment, or joint methods citation depending on delivery timing.
4. **Timing of the direction-field removal** — do we describe it as done (schema-breaking release) or planned (Discussion)?
5. **Cutoff date** for the compound / reaction growth statistics (i.e. which release version becomes "2026 snapshot").

## Explicitly out of scope for this paper

- ANIME-based heuristic-break detection.
- Full PlantSEED refactor (three-way merge of PlantSEED-Curation / PlantSEED-Refactor / plantseed-v3) — separate concern.
- Reaction-similarity work (Laino paper embeddings, horizyn) — worth mentioning in Discussion as a future direction but not a Results section.
- Kinetic-constant integration (mentioned in 2020 Introduction as a related need; still not the focus here).

## Infrastructure note (Travis → GitHub Actions)

The 2020 paper's Methods described "We utilize Travis CI along with scripts for testing data immediately, and reporting whether or not data in the pull request is valid." Travis is no longer in use. The replacement is a **GitHub Actions workflow that runs at PR submission time** to validate any newly curated entries against the existing schema and consistency checks. This should be a short Methods paragraph in the update paper, framed as an infrastructure upgrade rather than a novel contribution.

## Relationship to prior in-house work

- **Structure curation stream** (`MSD_Structures/`) — the PubChem pipeline, stereo-loss guard, curator override system, and formula-conflict resolver are all recent products of this working directory. The 2020 paper described mass-balance in one paragraph; the update should devote several paragraphs to the *pipeline* that produces mass-balanced compounds today.
- **Priority scope** — the ~9K reactions / ~6.5K compounds from v7.0 ModelSEEDTemplates + PlantSEED (`/scratch/seaver/Claude_Projects/priority/`) is the working slice for structure curation; the empirical study extends to the full ModelSEED v2 model corpus.
- **PlantSEED** — currently spread across three sibling project directories; a merge is pending but out of scope for this paper.

---

*Regenerate this file whenever paper scope shifts materially. Attribution for source paper information: PubMed ([10.1093/nar/gkaa746](https://doi.org/10.1093/nar/gkaa746)).*
