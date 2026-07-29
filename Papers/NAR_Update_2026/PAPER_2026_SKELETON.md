# ModelSEED Biochemistry Database — 2026 Update Paper Skeleton

**Target venue:** Nucleic Acids Research, Database Issue.
**Reference outline:** matches the 2020 paper (Seaver et al., [10.1093/nar/gkaa746](https://doi.org/10.1093/nar/gkaa746)) section-for-section, updated with what's new since. Word cap ~4000. See `PAPER_2026_GUIDE.md` for the in-repository work needed to fill each section.

---

## Title (draft)

> The ModelSEED Biochemistry Database, 2026 update: curated molecular structures, multi-source thermodynamics, atom mapping, and the impact of reaction-direction assignments on metabolic-model output.

## Abstract (draft outline)

- Six years since the 2020 release; ModelSEED biochemistry remains the chemistry foundation for ModelSEED and KBase metabolic-model reconstruction.
- Four update themes: (i) improved molecular structures via a curator override + provenance pipeline with PubChem validation and RDKit canonicalization; (ii) systematic mass-balance completeness, including a database-wide standardization of protein-carrier cofactors (acyl-ACP phosphopantetheine, biotinyl-BCCP, lipoyl carriers); (iii) multi-source thermodynamics (eQuilibrator, rebuilt Group Contribution, dGpredictor retrained on ModelSEED, OpenTECR) and heuristic overlays for reaction direction; (iv) new capabilities — atom mapping via collaboration with the Nikoloski lab, a refreshed reaction-similarity matrix, and PyPi distribution with an API.
- Novel empirical result: systematic study of how the choice of reaction-direction source affects growth-rate, essential-gene, and mass-balance outcomes across the ModelSEED v2 draft-model corpus ([10.1101/2023.10.04.556561](https://doi.org/10.1101/2023.10.04.556561)).
- Updated numbers (compounds, reactions, structures, balanced, thermodynamics coverage) to close.

## Introduction

- Six-year retrospective on community curation and the growth of the field's demand on biochemistry databases.
- What is qualitatively new in the domain since 2020: **atom mapping** as a first-class capability, **ML-based thermodynamics** (dGpredictor and successors) reaching maturity, growing awareness that reaction-direction assignments materially affect model output.
- Position the update: same Rosetta-Stone integration mission, same community-contribution philosophy, expanded pipeline for structure curation, systematic thermodynamics, and a first empirical accounting of reaction-direction sensitivity.
- Section roadmap (three technical components — structures, similarities, dynamics — plus their integration into the repository, the website, and the manuscript).

## Materials and Methods

### Collation of biochemical data

Updated list of sources and their release versions (KEGG, MetaCyc, BioCyc databases, BiGG, MetaNetX, Rhea) as of the 2026 cutoff release. Note any newly added or dropped sources vs 2020.

### Biochemical integration

Largely as 2020 (staged primary / secondary / tertiary sources with structure-first matching, then identifier, then synonym). Note any refinements to matching heuristics.

### Provenance

Extends the 2020 provenance model (source-DB attribution on every compound and reaction) with **curator provenance**: per-curator override files (`Biochemistry/Curation/overrides/structure_picks/<curator>.tsv`) that record *what* was chosen, *why*, and *who* decided, for every reconciled structure. Reconciled entries are distinguishable from silently ingested ones.

### Transport

Unchanged from 2020 (compartment-generalized to 0 or 1 for reusability).

### Structure-curation pipeline

Substantially expanded from the 2020 "Protonation and conversion" and "Balancing" paragraphs. New pipeline consists of:

- **RDKit canonicalization** — migration from ChemAxon Marvin to RDKit-canonical SMILES via a one-time mass recanonicalization pass; per-source SMILES preserved.
- **PubChem validation stage** — for every candidate structure imported from a source, a validator (with a stereo-loss guard) verifies the structure against PubChem and rejects imports that would drop stereo information present in the canonical record.
- **Curator override system** — hand-picked replacements per curator, with rationale recorded.
- **Formula-conflict resolver** — classification of source disagreements (H-only vs skeleton vs element-set) and rules for auto-matching to the stored `compound_*.json` formula when possible; explicit exclusion when no source matches (via `mass_balance_excluded.tsv`).
- **MetaCyc→KEGG structure alignment policy** — MetaCyc-only structures adjusted to KEGG conventions when needed to preserve reaction balance in reactions already validated under KEGG-conforming compounds.
- **Protein-carrier cofactor standardization** — acyl-ACP compounds include the phosphopantetheine cofactor in their stored formula; the same pattern applied to biotinyl-BCCP carriers (biotin) and lipoyl carriers on H-protein / E2 subunits (lipoic acid). Motivated by cofactor biosynthesis modelability and by removal of false-positive mass-imbalance flags. Candidate extensions (covalent FMN/FAD, molybdopterin, heme-c) surveyed and either included or deferred with explanation.

### Balancing of reactions

Updated. Formula derivation still uses RDKit + OpenBabel; mass and charge balance still tracked via the `status` field. The novelty is the *pipeline that produces balanceable compounds* (above), not the balancing algorithm itself.

### Multi-source thermodynamics

Replaces the 2020 "eQuilibrator (primary) + Group Contribution (fallback)" pair with a **four-source integration**:

1. **eQuilibrator** — refreshed against the current compound set.
2. **Group Contribution** — rebuilt from MFAToolkit against the current compound set.
3. **dGpredictor** — retrained on ModelSEED.
4. **OpenTECR** — integrated where coverage exists.

Also introduces per-reaction-class heuristic overlays for reaction direction (e.g. hydrolases, cofactor-loading reactions, transport with proton gradients), applied globally on top of any of the four sources.

### Reaction similarity

New section. Pipeline for computing reaction similarity from updated compound structures; GPU regeneration completes in under ~10 minutes end-to-end, so the matrix is regenerated whenever compound curation triggers a material change. External reaction-embedding foundation-model selection is described: candidate models are evaluated for how discerning their embeddings are on ModelSEED reactions, and the paper defends the chosen model on that basis.

### Atom mapping

New section. Describes the collaboration with the Nikoloski lab and the mapping methodology being applied (RXNMapper-family or its successor, depending on delivery). Coverage-to-date on the priority-scope reaction set (v7.0 templates) is reported.

### Undetermined compounds

Updated. Retains the 2020 handling of R-groups and generic compounds. Adds a note on fragment-aware pKa/pKb (open work; currently the biochemistry object only stores pKa/pKb for atoms in the first fragment, and this needs a schema-plus-code update).

### Community contributions, branches, and CI

Retains the 2020 GitHub PR workflow (dev / master branches, external contributor PRs to dev). CI has migrated from **Travis to GitHub Actions**: the same intent as 2020 (validate any newly curated entries at PR time), implemented in the current CI tooling. Describe the actions workflow briefly.

### Release procedure

Unchanged from 2020 (quarterly release to `master`, tagged in GitHub, deployed to ModelSEED and KBase).

### Distribution: PyPi + API

New section. The database is now installable as a PyPi package with a documented programmatic API alongside the existing Solr REST endpoint. Brief examples of typical calls (fetch compound, walk reaction, query balance).

## Results

### Growth in compounds and reactions

Table mirroring the 2020 Table 2 (compounds, compounds with structures, compounds with generic groups, reactions, complete reactions, balanced reactions, reactions with generic groups) with the 2026 column added. Same for Table 3 (per-source counts).

### Structure-curation improvements

- Per-curator counts of overrides applied since the 2020 release.
- Number of compounds with mass-balance exclusions (via `mass_balance_excluded.tsv`).
- PubChem pipeline yield (compounds validated / rejected).
- Reduction in false-positive mass-imbalance flags attributable to protein-carrier-cofactor standardization, broken out per class (acyl-ACP phosphopantetheine + biotinyl-BCCP + lipoyl-carrier + any extensions that make the cut).
- Optionally: illustrative example figure of one reconciled compound with its provenance chain.

### Multi-source thermodynamics

- Per-source coverage counts (compounds with ΔfG′, reactions with accepted ΔrG′) across the four sources.
- Head-to-head comparison of ΔrG′ across sources (correlation, uncertainty distributions), mirroring the 2020 Figure 2.
- How the per-reaction-class heuristic overlays change direction assignments relative to the source values alone.

### Reaction similarity

- The refreshed similarity matrix statistics (density, clustering behavior on curated reaction classes).
- Head-to-head comparison of external reaction-embedding foundation models with the selection defended on discerning-power grounds.

### Empirical study — reaction-direction heuristic sensitivity on ModelSEED v2 models

The novel empirical hook. Systematic sweep over the ModelSEED v2 draft-model corpus ([10.1101/2023.10.04.556561](https://doi.org/10.1101/2023.10.04.556561)). Per-model metrics:

- Predicted growth rate under the model's original medium.
- Essential-gene set (single-gene knockout FBA) — Jaccard overlap across direction sources.
- Mass-balance survival rate (reactions removed as thermodynamically infeasible under each direction source).
- Feasibility of biomass production under the 390 Biolog conditions (mirroring the 2020 whole-DB FBA analysis).

Reported as a matrix of `(model × direction source × metric)` deltas, plus a summary figure identifying direction sources that systematically break vs preserve biology across the corpus.

### Atom mapping coverage

- Fraction of priority-scope reactions with atom mappings delivered.
- Illustrative use case: tracking C, P, S atoms across a metabolic pathway using the mappings.

### Database connectivity (whole-DB FBA)

Retained from 2020, refreshed against the 2026 database. Updated version of Table 5 (total reactions / mass-balanced / reversible / functional / functional growth conditions).

### Ontological mappings for annotation comparison

Retained from 2020. Report any growth in the ontology (equivalent compound sets, lumped reaction sets, context-specific reaction sets).

## Discussion

- Six years of community-driven curation: what worked, what accumulated as debt.
- The protein-carrier-cofactor standardization pattern as a general DB-cleanup approach (ACP → biotinyl → lipoyl → extensions).
- Open problems for the next release:
    - Fragment-aware pKa/pKb (biochemistry object + code refresh).
    - EC number refresh beyond ExPASy name matching (bound to have false or missing associations).
    - Obsolescence-leakage audit — obsolete compounds and reactions marked as such, but are they accidentally used elsewhere?
    - Direction-field removal from the reaction schema, now that thermoreversibility + template capture direction.
- Where atom mapping and multi-source thermodynamics take the community next.

## Data Availability

- GitHub repository, branches (dev / master), release tagging.
- GitHub Actions CI for PR-time validation of curated entries.
- PyPi package + documented API.
- Solr REST endpoint (retained from 2020).
- Updated ModelSEED website with atom-mapping surfaced.
- Nikoloski-lab atom-mapping endpoint (or dataset, depending on delivery format).
- KBase integration retained.

## Author contributions, funding, acknowledgments

To be composed at drafting time. Note: Nikoloski-lab atom-mapping credit format is an open decision (co-authorship vs collaboration acknowledgment vs joint methods citation) — see `PAPER_2026_PLAN.md`.
