# ModelSEED Biochemistry Database — 2026 Update Paper Guide

**A step-by-step, tackle-one-section-at-a-time companion to `PAPER_2026_SKELETON.md`.**

The skeleton describes what the paper will *look like*. This guide describes what needs to be *done in this repository* before each section can be written. Sections are ordered so that upstream work lands before the sections that depend on it. Each entry gives:

- **Repository work** — concrete tasks: analyses to run, scripts to write / refactor, data to publish, figures to generate, decisions to close.
- **Writing** — the prose deliverable.
- **Blocked by** — other guide sections that must be complete first.
- **Effort** — S (a few hours) / M (a day or two) / L (a week+).

Legend for status boxes: `[ ]` open · `[~]` in progress · `[x]` done.

---

## Order of attack

Preferred sequence (dependencies flow downward — anything below can be done in parallel with anything above it once its blockers are cleared):

1. Cutoff release decision
2. Source-database snapshot + growth statistics
3. Structure-curation pipeline (write-up of what is already built)
4. Formula-conflict resolver (write-up of what is already built)
5. Protein-carrier cofactor standardization (survey + execute + measure impact)
6. Multi-source thermodynamics (collect from four sources)
7. Reaction similarity (regenerate matrix + foundation-model selection)
8. Atom mapping (integrate Nikoloski deliverables + measure coverage)
9. Empirical study — reaction-direction sensitivity on ModelSEED v2 models (biggest lift)
10. Database connectivity (whole-DB FBA refresh)
11. Ontological mappings (mostly a refresh)
12. Distribution: PyPi + API
13. CI: GitHub Actions PR-time validation
14. Website + atom-mapping surfacing
15. Discussion + Data availability + Abstract + Introduction (last)

---

## §1 · Cutoff release decision

- [x] **Decision (2026-07-29):** the 2026 paper snapshot is tagged **`v2.0.0`**, cut *after* guide sections §5–§8 land (protein-carrier standardization + multi-source thermodynamics + reaction similarity refresh + atom mapping delivered). Until then, the current `dev` HEAD (**`3cc9a1d`**) is the provisional snapshot for tracking growth statistics as they move.
- [~] **Repository work**
    - [ ] After §5–§8 land, cut and push `v2.0.0` from `dev`; announce the cutoff so late-landing PRs are held for `v2.1.0`.
    - [x] Baseline decision recorded here so subsequent guide sections cite the right tag.
- [ ] **Writing** — every subsequent table caption cites `v2.0.0`; the Data Availability section links to the tagged release on GitHub.
- **Blocked by:** nothing (for the decision); tagging blocked by §5–§8.
- **Effort:** S.

## §2 · Source-database snapshot + growth statistics

- [~] **Repository work**
    - [x] Table 2 and Table 3 numbers computed against dev HEAD `1111754` and persisted at `Papers/NAR_Update_2026/data/snapshot_2026-07-29.md`. Provisional (will be re-run against the v2.0.0 tag when it's cut).
    - [x] Delta narrative captured — balanced-reaction share regressed 70%→61%, which motivates the §5 recovery story.
    - [x] ChEBI now appears in Table 3 alongside the sources named in the 2020 paper.
    - [x] Rhea count fixed — the earlier ad-hoc query used lowercase `rhea`; the alias key is capitalized `Rhea`. Actual: 1,931 compounds, **17,477 reactions** (roughly 2× the 2020 count of 8,786).
    - [ ] Promote the two ad-hoc snippets into a reproducible script under `Scripts/Statistics/Tables/` and refresh `Growth_Stats.tsv` when v2.0.0 is cut.
    - [ ] Record which KEGG / MetaCyc / BioCyc / BiGG / MetaNetX / Rhea version was ingested for the snapshot (Methods "Collation" paragraph needs source version numbers, not just names).
    - [ ] Decide: aggregate the 11 BioCyc-family sources in Table 3 (as in 2020) or itemize them (as they appear in the current alias table).
- [ ] **Writing** — Methods "Collation of biochemical data" paragraph updated with source versions; Results "Growth in compounds and reactions" section with the two updated tables and the balanced-share regression narrative.
- **Blocked by:** §1 (decided).
- **Effort:** S–M (mostly done; source-version capture + Rhea fix remaining).

## §3 · Structure-curation pipeline (write-up)

The pipeline is largely built already in this session's work. This step is *writing*, not building.

- [ ] **Repository work**
    - Confirm the following live in the repo and are pointer-linkable from the paper: PubChem validation stage with stereo-loss guard; RDKit canonicalization script and canonicalized `All_ModelSEED_Structures.txt`; curator override files at `Biochemistry/Curation/overrides/structure_picks/<curator>.tsv`; mass-balance exclusion file at `Biochemistry/Curation/exclusions/mass_balance_excluded.tsv`; formula-conflict report generator.
    - Draft one worked example per pipeline stage for the paper (e.g. cpd35693 for the stereo-loss guard; cpd00766 or cpd00849 for the mass-balance exclusion).
- [ ] **Writing** — Methods "Structure-curation pipeline" section (the biggest new Methods block). Provenance framing throughout (what / why / who).
- **Blocked by:** §1.
- **Effort:** M.

## §4 · Formula-conflict resolver (write-up)

Also mostly built.

- [ ] **Repository work**
    - Freeze the H-only vs skeleton vs element-set classifier logic in the repo (script that produces `Formula_Conflicts.txt`).
    - Report before/after counts of formula conflicts across the priority-scope set.
- [ ] **Writing** — Methods "Balancing of reactions" section refresh + Results "Structure-curation improvements" bullet on formula-conflict resolution.
- **Blocked by:** §1.
- **Effort:** S–M.

## §5 · Protein-carrier cofactor standardization

The load-bearing new methods contribution. Not fully built yet.

- [ ] **Repository work**
    - **Census** — for each candidate carrier class (acyl-ACP, biotinyl-BCCP, lipoyl-carrier + candidate extensions covalent-FMN/FAD, molybdopterin, heme-c), count how many compounds in the DB represent the loaded state and how many currently include vs omit the cofactor's atoms.
    - **Decide** which classes to include in this paper vs defer (open decision #7 in `PAPER_2026_PLAN.md`).
    - **Execute** the standardization pass per included class: extend the `acps_formula_charge.tsv` framework (or a sibling file per class) with the cofactor-inclusive formulas and charges.
    - **Measure** the before/after false-positive mass-imbalance flag counts.
    - **Sanity-check** with FBA on a representative model that the standardization doesn't introduce infeasibility (should only fix balance, not break flux).
- [ ] **Writing** — Methods "Protein-carrier cofactor standardization" subsection under Structure-curation pipeline; Results bullet on false-positive reduction broken out per class. Optionally a Discussion paragraph on this as a general DB-cleanup pattern.
- **Blocked by:** §3.
- **Effort:** L.

## §6 · Multi-source thermodynamics

- [ ] **Repository work**
    - **eQuilibrator refresh** — rerun eQuilibrator against the current compound set at pH 7.0, ionic strength 0.25 M, 298.15 K (matching the 2020 conditions).
    - **Group Contribution rebuild** — rebuild from MFAToolkit against the current compound set.
    - **dGpredictor retrain** — retrain on ModelSEED training data; commit the trained model artifact or a reproducible training script.
    - **OpenTECR integration** — pull OpenTECR values where coverage exists; document the mapping between OpenTECR reaction identifiers and ModelSEED reaction IDs.
    - **Per-reaction-class heuristic overlays** — codify the direction rules for hydrolases, cofactor-loading reactions, transport with proton gradients, and any others the team agrees on.
    - Produce a unified table: per compound / reaction, the ΔG′ and uncertainty from each of the four sources, plus the direction assignment from each source pre- and post-heuristic overlay.
- [ ] **Writing** — Methods "Multi-source thermodynamics" section; Results "Multi-source thermodynamics" section with a refreshed Figure 2 analogue (correlation + uncertainty distributions across sources).
- **Blocked by:** §3 (curated structures are the input); §5 (some ΔG values change when carrier cofactors are included).
- **Effort:** L.

## §7 · Reaction similarity

- [ ] **Repository work**
    - Regenerate the reaction-similarity matrix against the updated compound structures (Ray reports ~10 min on GPU end-to-end).
    - **Foundation-model bake-off** — evaluate the candidate external reaction-embedding models on ModelSEED reactions; measure discerning-power (silhouette on curated reaction classes, retrieval quality on held-out reaction pairs, or a similar metric).
    - Pick and defend the model.
- [ ] **Writing** — Methods "Reaction similarity" section; Results "Reaction similarity" section with the matrix statistics and the model-selection table.
- **Blocked by:** §3.
- **Effort:** M.

## §8 · Atom mapping

- [ ] **Repository work**
    - Land whatever the Nikoloski lab delivers into the repo (as a subdirectory under `Biochemistry/AtomMappings/` or similar, TBD).
    - Compute coverage: fraction of priority-scope reactions with a delivered mapping.
    - Build the illustrative use case (C/P/S atom tracking across a metabolic pathway).
    - Wire the atom mapping into the website (see §14).
- [ ] **Writing** — Methods "Atom mapping" section; Results "Atom mapping coverage" section; Discussion note.
- **Blocked by:** §3 (mappings depend on canonicalized structures); external delivery timeline.
- **Effort:** M–L (depends on Nikoloski delivery).

## §9 · Empirical study — reaction-direction sensitivity on ModelSEED v2 models

The novel empirical hook. This is the largest single unit of work.

- [ ] **Repository work**
    - Pull the ModelSEED v2 draft-model corpus from [10.1101/2023.10.04.556561](https://doi.org/10.1101/2023.10.04.556561) into the repo (or a sibling data folder) and pin the exact model set used.
    - For each `(model × direction source)` pair produce: predicted growth rate under original medium; essential-gene set from single-gene knockout FBA; number of reactions removed as thermodynamically infeasible; number of Biolog conditions under which biomass is feasible.
    - Aggregate into a `(model × direction source × metric)` matrix and identify direction sources that systematically break vs preserve biology.
    - Design the summary figure (options: heatmap of essential-gene Jaccard across sources; per-metric delta boxplots across models; a per-model "direction-sensitivity" score).
- [ ] **Writing** — Results "Empirical study" section — the paper's headline. Materials-and-Methods paragraph on the empirical design (corpus, metrics, statistical treatment).
- **Blocked by:** §6 (need all four direction sources computed) and §5 (ACP standardization should be applied before running FBA so the "mass-balance survival" metric is not dominated by carrier-cofactor artifacts).
- **Effort:** L.

## §10 · Database connectivity (whole-DB FBA refresh)

- [ ] **Repository work**
    - Rerun the whole-DB FBA analysis from the 2020 paper against the 2026 compound / reaction set: total reactions, mass-balanced, reversible, functional reactions, functional growth conditions across the 390 Biolog conditions.
    - Include a comparison row to 2020 for the Table 5 analogue.
- [ ] **Writing** — Results "Database connectivity" section refresh.
- **Blocked by:** §2, §5, §6.
- **Effort:** M.

## §11 · Ontological mappings

- [ ] **Repository work**
    - Report any growth in the three ontology types since 2020 (equivalent compound sets, lumped reaction sets, context-specific reaction sets).
    - Refresh the *E. coli* iJR904 vs ModelSEED comparison from the 2020 paper (Figure 4 analogue) — do the improved structures reduce the number of unique-in-iJR904 reactions further?
- [ ] **Writing** — Results "Ontological mappings" section — mostly a refresh with updated counts.
- **Blocked by:** §2, §3.
- **Effort:** S–M.

## §12 · Distribution: PyPi + API

- [ ] **Repository work**
    - Package the database + Python helpers as an installable PyPi package.
    - Document the programmatic API (fetch compound, walk reaction, query balance / thermodynamics / atom mapping).
    - Publish to test PyPi, then production PyPi.
    - Add a quickstart to the repo README pointing at the package.
- [ ] **Writing** — Methods "Distribution: PyPi + API" section; Data Availability entry.
- **Blocked by:** §1 (need a stable release to package).
- **Effort:** M.

## §13 · CI: GitHub Actions PR-time validation

- [ ] **Repository work**
    - Refactor the validation code (2020's Travis-run scripts) to run under GitHub Actions.
    - Wire it to run on every PR to `dev`.
    - Ensure it catches at least: mass-balance regression, formula-conflict introduction, invalid override file schema, missing provenance fields.
- [ ] **Writing** — Methods "Community contributions, branches, and CI" paragraph updated with the Actions workflow (short — this is an infrastructure paragraph, not a novel contribution).
- **Blocked by:** §3, §4 (the validation checks are defined by these pipelines).
- **Effort:** M.

## §14 · Website + atom-mapping surfacing

- [ ] **Repository work**
    - Update the ModelSEED website (`modelseed.org/biochem`) to reflect the current schema and to surface atom mappings on the reaction landing pages.
    - Coordinate with Alan on repository / website integration (per the 2026-07 working session).
- [ ] **Writing** — Data Availability entry updated.
- **Blocked by:** §8.
- **Effort:** M.

## §15 · Discussion + Data availability + Abstract + Introduction

Best written last so the specifics of every prior section are already fixed.

- [ ] **Repository work** — none (this is writing).
- [ ] **Writing**
    - **Introduction** — position the six-year update; motivate atom mapping, ML thermodynamics, reaction-direction transparency as drivers.
    - **Discussion** — six-year retrospective; general lessons from the protein-carrier cofactor standardization pattern; enumerate open problems for the next release (fragment-aware pKa/pKb, EC refresh beyond ExPASy name matching, obsolescence audit, direction-field removal).
    - **Data Availability** — GitHub, PyPi, Solr, website, atom-mapping endpoint, KBase integration.
    - **Abstract** — condense the above.
    - **Author contributions / funding / acknowledgments** — close per team decision on Nikoloski credit (open decision #3 in `PAPER_2026_PLAN.md`).
- **Blocked by:** everything above.
- **Effort:** M.

---

## Open decisions carried from `PAPER_2026_PLAN.md`

Repeated here so they are visible while tackling sections:

1. Which direction sources make the final cut for the empirical study (§9).
2. Whether OpenTECR gets its own Methods paragraph or is folded into eQuilibrator (§6).
3. How the Nikoloski atom-mapping work is credited (§15).
4. Timing of the direction-field removal — done or planned (§15 Discussion).
5. ~~Cutoff release for the 2026 snapshot~~ — **decided 2026-07-29:** `v2.0.0`, cut after §5–§8 land.
6. Which external reaction-embedding foundation model to use (§7).
7. Which protein-carrier-cofactor classes make the cut (§5).

---

*Update this file alongside `PAPER_2026_PLAN.md` and `PAPER_2026_SKELETON.md` whenever section scope or dependencies change. Attribution for source paper information: PubMed ([10.1093/nar/gkaa746](https://doi.org/10.1093/nar/gkaa746)).*
