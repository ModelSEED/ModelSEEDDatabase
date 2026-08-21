# The ModelSEED Biochemistry Database: 2026 update

**Reading-order assembly of the drafts under `drafts/`.**
This is a working manuscript, not a submission. Every section either (a) inlines the draft prose that already exists, or (b) marks a `[DRAFT PENDING]` placeholder pointing at what still needs writing. `[TBD]` markers within prose flag numbers or decisions that fill in as underlying work completes.

Prior paper: Seaver et al. 2020, Nucleic Acids Research 49(D1):D575–D588, DOI [10.1093/nar/gkaa746](https://doi.org/10.1093/nar/gkaa746). Cite this as the immediate predecessor throughout.

---

## Title

The ModelSEED Biochemistry Database, 2026 update: curated molecular structures, multi-source thermodynamics, atom mapping, and the impact of reaction-direction assignments on metabolic-model output.

## Authors

[Author list TBD — the 2020 paper's list plus this cycle's contributors (Ray, ctaylor, Nikoloski-lab collaborators per the crediting decision), minus anyone stepping off, in the agreed order.]

## Abstract

*(Source: drafts/discussion_and_availability.md — this is the outline; final abstract is drafted last so it can distill finalised results.)*

Six years since the 2020 release, the ModelSEED biochemistry database remains the chemistry foundation for ModelSEED and KBase metabolic-model reconstruction. Growth to [TBD] compounds, [TBD] reactions, [TBD %] balanced. Structure curation via PubChem validation, RDKit canonicalisation, per-curator provenance tracking, and a new mass-balance exclusion mechanism; extended to protein-carrier cofactor standardisation (acyl-ACP phosphopantetheine, biotinyl-BCCP, lipoyl carrier). Multi-source thermodynamics — eQuilibrator refresh, group-contribution rebuilt from MFAToolkit, dGpredictor retrained on ModelSEED, OpenTECR. New capabilities: atom mapping through collaboration with the [TBD] lab, a refreshed reaction-similarity matrix with defended foundation-model selection, PyPi distribution with a documented API. Novel empirical finding: reaction-direction source choice materially affects [TBD] on the ModelSEED v2 draft-model corpus published in [10.1101/2023.10.04.556561](https://doi.org/10.1101/2023.10.04.556561).

## Introduction

*(Source: drafts/discussion_and_availability.md.)*

The 2020 release of the ModelSEED biochemistry database ([10.1093/nar/gkaa746](https://doi.org/10.1093/nar/gkaa746)) established the database as a curated integration of KEGG, MetaCyc, BiGG, MetaNetX, Rhea, BioCyc-family databases, and dozens of published metabolic models — a "Rosetta Stone" for cross-database mapping backed by structure-first compound reconciliation, per-reaction thermodynamics, and a community-contribution workflow on GitHub. Six years later, we present an update that extends the database along four axes and provides the first systematic empirical accounting of a question the 2020 release could not: how much does the choice of a reaction-direction source actually matter for downstream metabolic-model output.

Three developments in the field motivate specific parts of this update. **Atom mapping** has matured from a specialty capability into a routine expectation for any biochemistry database supporting isotope-tracing or degradation-prediction workflows; we present a first-pass integration of atom mappings via a collaboration with the [TBD] group. **Machine-learning-based thermodynamics** methods (dGpredictor and successors) now cover a substantial fraction of the ModelSEED compound set and enable a multi-source approach that goes beyond the single-primary-plus-fallback strategy of 2020. **Reaction-direction transparency** has emerged as a growing concern in the metabolic-modeling community; each choice a curator makes about direction propagates through model output in ways that are rarely quantified, and we provide the first systematic quantification here.

The paper walks through three technical components — compound structures, reaction similarities, and dynamics (the reaction-direction empirical study) — and their integration into the repository, the website, and the manuscript. [TBD 1–2 sentence roadmap.]

## Materials and Methods

### Collation of biochemical data

**[DRAFT PENDING]** — needs a paragraph mirroring the 2020 paper's "Collation" subsection, updated with source versions of KEGG, MetaCyc, BioCyc-family databases, BiGG, MetaNetX, Rhea, and ChEBI ingested for the v2.0.0 snapshot. See `data/snapshot_2026-07-29.md` for the per-source counts; source versions are the remaining `[TBD]` piece.

### Biochemical integration

**[DRAFT PENDING]** — mostly a carryover from the 2020 paper's "Biochemical integration" subsection (staged primary / secondary / tertiary sources with structure-first matching, then identifier, then synonym). Note any refinements to matching heuristics introduced since 2020.

### Provenance

**[DRAFT PENDING for basics; extended discussion inlined in Structure-curation pipeline below.]** — the 2020 provenance model (source-DB attribution on every compound and reaction) is unchanged. The extension is per-curator provenance for reconciled structures, covered in the Structure-curation pipeline section.

### Transport

**[DRAFT PENDING]** — unchanged from 2020 (compartment-generalized to 0 or 1 for reusability across models).

### Structure-curation pipeline

*(Source: drafts/methods_structure_curation.md — inlined below.)*

The 2020 release described protonation and reaction balancing in two short paragraphs: ChemAxon Marvin was applied to protonate every source structure at pH 7, and a combination of RDKit and OpenBabel was used to derive the formula and charge from the resulting InChI and SMILES. Six years of curation experience have exposed the limits of a single-tool ingestion pass: source databases disagree on stereochemistry, hydrogen counts, and — in the case of protein-bound reactive groups — where to draw the boundary between the small-molecule substrate and its carrier. To address these limitations we have replaced the single-pass ingestion with a multi-stage structure-curation pipeline that (i) canonicalises structures deterministically, (ii) validates every source-provided structure against an external reference, (iii) allows explicit curator overrides with recorded provenance, (iv) resolves formula conflicts by class rather than case-by-case, and (v) marks compounds whose stored formula cannot be reproduced from any source as mass-balance-excluded rather than silently unbalancing the reactions that use them.

#### Canonicalisation

All ModelSEED SMILES are now stored in an RDKit-canonical form. The `Scripts/Structures/Recanonicalize_SMILES.py` script applies `Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)` and iterates until the output is a fixed point (in practice, at most three passes are needed). A one-time mass recanonicalisation converted every entry in `Biochemistry/Structures/Unique_ModelSEED_Structures.txt` and `All_ModelSEED_Structures.txt` to the canonical form; per-source files under `Biochemistry/Structures/{KEGG,MetaCyc,ChEBI,Rhea}/` are preserved unchanged so that source-provided SMILES remain auditable. Downstream tooling (`Scripts/Structures/List_ModelSEED_Structures.py`) is idempotent with respect to `Recanonicalize_SMILES.py`: running the picker after a recanonicalisation, or vice versa, produces no further changes.

#### PubChem validation stage

Every candidate structure imported from a source database is passed through a validator that queries PubChem for the reference structure and checks compatibility along four axes: InChI stereo layer presence, `/t` (tetrahedral) stereocenter configuration, `/b` (bond E/Z) configuration, and `/m` (mirror) enantiomer parity. When PubChem returns a structure that has strictly fewer stereo layers than the canonical ModelSEED record, the validator rejects the update rather than accepting it as a "correction". This stereo-loss guard was motivated by cpd35693 (coniferyl alcohol radical), whose canonical structure carried a `/b8-3+` marker specifying E geometry at the sinapyl double bond; PubChem returned an InChI with no `/b` layer at all, and the pre-guard validator accepted the update as compatible because the shared-center inversion count was zero on an empty shared-bond set. The guard is implemented as two additive checks in `Validate_PubChem_Structures._check_stereo_compatibility`: a layer-presence check that rejects when the canonical record has a `/t` or `/b` layer that PubChem drops entirely, and a specification-loss check that rejects when a shared `/t` centre transitions from a specified (`+` / `-`) configuration to an unspecified (`?`) one. Both rejections use the `stereo_loss:` prefix in the rejection reason so they group naturally in the validator's Phase 5 log.

#### Curator override system

Structures reconciled by hand rather than by the automated cascade are recorded in per-curator tab-separated files under `Biochemistry/Curation/overrides/structure_picks/<curator>.tsv`. Each row records the compound identifier, the picked source database and identifier (with the ingestion stage — Charged or Original), the format of the picked structure (InChI or SMILES), the structure string itself, the pick date, and a rationale. The rationale field is required: rather than allowing a curator to silently overwrite a compound, the system asks *why* the pick was made, and stores the answer. This distinguishes reconciled entries from those silently ingested from a source database. Two override files are currently populated — `samseaver.tsv` and `Ray16.tsv` — reflecting the primary and grad-student curation streams. The picker (`List_ModelSEED_Structures.py`) consults these files at ingestion time, applies the override in place of the cascade's structural pick, and records `manual_curation:<curator>` in the pick-reason ledger so that every override's provenance is recoverable.

#### Formula-conflict resolver

When multiple sources supply structures for the same ModelSEED compound and the resulting formulas disagree, the pipeline classifies the conflict as one of three types: **H-only**, where sources agree on the heavy-atom composition but disagree on hydrogen count; **heavy-atom** (skeleton), where sources disagree on heavy-atom counts within the same element set; and **element-set**, where sources disagree on the elements present (typically an accidental variant containing a halogen or metal). The classification is implemented in `Scripts/Structures/Classify_Formula_Conflicts.py`, which reads `Biochemistry/Structures/_reports/Formula_Conflicts.txt`, groups the source-variant rows by compound, and emits both a summary count and a per-compound classified TSV. An optional `--priority` flag adds a cross-tabulation against a compound-ID priority list (e.g. the v7.0 template compound set), so that curation work can be triaged by whether the conflict affects a template reaction. For H-only conflicts and for skeleton or element-set conflicts where at least one source's formula matches the value stored in `compound_*.json`, the resolver auto-picks the matching source as the least-disruptive choice; the compound acquires a structure without any reaction it participates in becoming newly imbalanced. For the residual cases where no source formula matches — a small class where the stored formula reflects a species (typically a radical or a doubly-deprotonated anion) that no source database catalogues — the compound is instead added to the mass-balance exclusion registry (below). Across three rounds of this session's curation, 64 formula-conflict compounds were resolved this way (52 H-only in the first round, 11 skeleton/element-set in the second, cpd00244 Ni²⁺ in the third), leaving 56 residual compounds — none of which fall in the v7.0 priority set.

#### Mass-balance exclusion mechanism

A new file, `Biochemistry/Curation/exclusions/mass_balance_excluded.tsv`, records compounds whose stored formula (in `compound_*.json`) cannot be reproduced from any available source structure. The schema is `cpd_id, name, reason, date, curator`. The compound is still assigned a source structure so downstream atom-mapping tools have something to work with, but `Update_Compound_Structures_Formulas_Charge.py` explicitly skips the formula and charge overwrite for excluded compounds; the picked structure's SMILES and InChIKey populate `compound_*.json` while the historical formula and charge are preserved unchanged. As a consequence, reactions in which an excluded compound appears remain balanced against the historical formula. The registry currently contains four compounds: cpd00766 (ascorbate radical, `C6H5O6/-2` — no source provides a matching one-electron oxidation intermediate; picked KEGG's `C6H7O6/0` for the atom-mapping structure); cpd00849 (quercetin 3,3'-bissulfate, `C15H8O13S2/-2` — sources are further-deprotonated); cpd27485 (Co(I)-corrinoid Fe–S protein, `C54H80CoN14O9R2/+1` — the sole available structure includes the DMB-nucleotide loop of vitamin B12, expanding the formula to `C67H93CoN16O16PR2`); and cpd12396 (a pyruvate-kinase phosphoserine residue, `C4H6N2O6PR2/-1` — off by one hydrogen from any source's protonation state).

#### MetaCyc→KEGG structure alignment policy

ModelSEED originated in a KEGG-first era, and many reactions were validated against KEGG-conforming compound formulas. When MetaCyc later supplied a structure absent from KEGG, adopting the MetaCyc structure directly could silently unbalance those pre-existing reactions. The pipeline now applies a MetaCyc-to-KEGG alignment policy for these cases: the MetaCyc-provided structure is adjusted so that its formula matches what a KEGG-conforming representation would produce. This policy is applied through the same override mechanism (a row in the curator TSV file with the alignment rationale in the `rationale` field) rather than through a separate file, so that every aligned structure carries its provenance forward.

#### Protein-carrier cofactor standardisation

Protein-bound intermediates that carry a reactive group between enzymes — most prominently the acyl-carrier-protein (ACP) family — are handled through a dedicated override file, `Biochemistry/Curation/overrides/acps_formula_charge.tsv`. Acyl-ACP compounds bind an acyl group to the phosphopantetheine cofactor of an acyl-carrier protein; unlike acyl-CoA compounds (where the full CoA structure is present in the compound record and no ambiguity arises), acyl-ACP compounds may or may not include the phosphopantetheine moiety in the stored formula. Where the phosphopantetheine is omitted, mass balance fails in any reaction that produces or consumes phosphopantetheine biosynthetically or that transfers acyl groups between the CoA-bound and ACP-bound pools. The override framework requires the phosphopantetheine formula to be included for every acyl-ACP compound; a bulk-add pass added 47 pantetheine-inclusive ACP formula rows in a single commit, verified as producing zero newly-imbalanced reactions across the whole database. The same standardisation pattern is being extended to biotinyl carriers (biotin covalently attached to a lysine ε-amine on the biotin-carboxyl-carrier proteins of the carboxylase family) and to lipoyl carriers (lipoic acid attached in the same way on the H-protein of the glycine cleavage system and on the E2 subunits of the α-ketoacid dehydrogenase complexes). Candidate extensions to survey but potentially defer for a future release include covalently-bound FMN and FAD (as in the histidyl-FAD of succinate dehydrogenase), molybdopterin cofactors, and heme-c attached via thioether bonds in c-type cytochromes. In every case the unifying claim is that the carrier cofactor is biosynthesised, its demand should be modellable, and reactions handling it should mass-balance; standardising the formula representation to always include the cofactor achieves all three goals in one pass.

### Balancing of reactions

**[DRAFT PENDING]** — carryover from 2020 (RDKit + OpenBabel derive formula and charge; balance status recorded in the `status` field). The novelty of the update is the pipeline that produces balanceable compounds, described above; the balancing algorithm itself is unchanged.

### Multi-source thermodynamics

*(Source: drafts/methods_multi_source_thermodynamics.md — inlined below.)*

The 2020 release combined two thermodynamics sources: eQuilibrator as the primary estimator of ΔfG′ and ΔrG′, and a group-contribution (GC) method as a fallback for compounds and reactions eQuilibrator could not estimate. For this update we integrate four sources, each providing an independent estimate wherever its coverage permits.

**eQuilibrator refresh.** [TBD version of eQuilibrator] was run against the 2026 compound set at pH 7.0, ionic strength 0.25 M and temperature 298.15 K, matching the 2020 conditions. Coverage: [N compounds with ΔfG′, N reactions with accepted ΔrG′; TBD].

**Group Contribution rebuilt from MFAToolkit.** The 2020 fallback GC method has been rebuilt from the MFAToolkit implementation against the current compound set, so the group definitions and parameter fits reflect the enlarged compound base rather than the 2020 snapshot. Coverage: [TBD].

**dGpredictor retrained on ModelSEED.** dGpredictor (an ML-based ΔG estimator) has been retrained using the ModelSEED-curated compound structures as training data. This produces a ModelSEED-native estimator rather than one trained on a competing biochemistry database, and lets us report ΔG estimates for compounds the other three methods cannot cover. Training procedure and hyperparameters: [TBD]. Coverage: [TBD].

**OpenTECR integration.** OpenTECR is a community effort to standardize experimentally-measured thermodynamic constants; we integrate its values where a mapping between OpenTECR reaction identifiers and ModelSEED reaction IDs exists. Coverage: [TBD].

**Direction assignment.** For each reaction with an accepted ΔrG′ (from any source, chosen with per-compound-consistency constraints described in the 2020 paper), we apply the same reversibility rules used in 2020 (based on Jankowski et al.) to assign one of `=` (reversible), `>` (irreversible left-to-right), or `<` (irreversible right-to-left). On top of this we add per-reaction-class heuristic overlays — hard rules for classes where the standard thermodynamic assignment is known to disagree with textbook biology. Rule set: [TBD; discussed in the working session but not yet codified in a repo script]. Overlays are applied after the base assignment and their scope is documented per class.

**Combined direction ledger.** For every reaction, the pipeline records: the ΔrG′ estimate from each of the four sources (with per-source uncertainty), the base direction assignment from each source, the heuristic overlay if applicable, and the final direction. This ledger is the input to the reaction-direction sensitivity study reported in Results.

### Reaction similarity

*(Source: drafts/methods_reaction_similarity.md.)*

Reaction similarity — a per-pair distance or similarity score across all mass-balanced ModelSEED reactions — is a new capability in the 2026 release. It supports two use cases: as a diagnostic for the reconciliation ontology described in the 2020 paper (reactions with high similarity that are not already ontology-linked are candidates for review), and as an input feature for downstream applications including annotation transfer, gap-filling, and pathway-inference.

**Pipeline.** Every ModelSEED reaction is represented by the concatenation of its balanced reactant and product structures. Structures are passed through an external reaction-embedding foundation model to produce a fixed-dimensional vector; the similarity between two reactions is the cosine similarity between their embeddings. GPU-accelerated batching allows the full pairwise similarity matrix to be regenerated in under approximately ten minutes on a single GPU, so the matrix is treated as a derived artifact and refreshed whenever compound curation materially changes the input structures.

**External reaction-embedding foundation-model selection.** Several external reaction-embedding models are available at the time of this release. Rather than adopting one by default, we evaluate the candidate models on their discerning power over ModelSEED reactions and defend the chosen model on that basis. The evaluation compares candidate models on: (i) intra-EC clustering — reactions sharing the same EC number should embed close together; (ii) intra-pathway clustering — reactions within a KEGG or MetaCyc pathway should embed close together; (iii) retrieval quality on a held-out set of known reaction-pair relationships (e.g. adjacent reactions in a pathway). Candidate models compared: [TBD list]. Chosen model: [TBD] on the basis of [TBD metric result].

**Publication of the matrix.** The similarity matrix is exported as a compressed sparse artifact under `Biochemistry/Reaction_Similarity/`. Users can regenerate it locally by rerunning the pipeline script against a specific ModelSEED release tag.

### Atom mapping

*(Source: drafts/methods_atom_mapping.md.)*

Atom mapping — the assignment of individual reactant atoms to their product-side counterparts across a mass-balanced reaction — is a new capability in the 2026 release. Mappings enable atom-level tracing of carbon, phosphorus, sulfur, and other atoms through a metabolic network, and are increasingly required for downstream applications including 13C metabolic flux analysis, pathway-degradation prediction, and cofactor-usage attribution.

Mappings for the ModelSEED biochemistry are being generated through a collaboration with the [TBD] group at the [TBD institution]. Mapping methodology: [TBD; will be described briefly with reference to the collaborators' methods paper]. Delivery format: [TBD — likely RXN files or a tabular format per reaction].

Mapped reactions are integrated into ModelSEED under `Biochemistry/AtomMappings/` (or equivalent, [TBD]). Each mapped reaction is stored alongside its ModelSEED reaction ID; downstream consumers can access mappings through the same Solr endpoint or the new PyPi API (see Distribution below).

Coverage as of the 2026 snapshot: [TBD fraction of priority-scope reactions with mappings delivered]. Coverage of the full ModelSEED reaction set: [TBD]. Reactions without a mapping fall into three classes: (i) reactions whose reactant compounds do not yet have a stored structure (necessarily unmappable); (ii) reactions with `R` groups whose mapping requires expansion to specific chain lengths; (iii) reactions whose complexity (e.g. very long acyl chains, PKS/NRPS iterative modules) exceeds the mapper's practical scope.

### Undetermined compounds

**[DRAFT PENDING]** — carryover from 2020 (R-groups treated as first-class atoms; generic compounds linked to their structurally-specific representatives). Add a paragraph on fragment-aware pKa/pKb as an open problem for the next release (Discussion also notes this).

### Community contributions, branches, and continuous integration

*(Source: drafts/methods_distribution_pypi_ci.md — CI section.)*

The 2020 release described a Travis CI setup that validated pull requests to `dev` at submission time. Travis is no longer used; the same intent — validating any newly curated entries at PR time — is now implemented via a GitHub Actions workflow. The workflow runs on every PR to `dev` and checks: (i) that the compound and reaction JSON shards parse and match their expected schema; (ii) that newly added or modified reactions remain mass- and charge-balanced against the (possibly modified) compound set; (iii) that newly added override or exclusion file rows conform to the appropriate schema; (iv) that every curated pick carries the required provenance fields. Failures block the PR from being merged until the contributor addresses them.

The workflow definition is under `.github/workflows/` in the repository; contributors can run the same checks locally before opening a PR by invoking the validation script directly (`Scripts/Validation/`, [TBD confirm script path]).

### Release procedure

**[DRAFT PENDING]** — carryover from 2020 (quarterly release to `master`, tagged in GitHub, deployed to ModelSEED and KBase). Note: the 2026 update ships as tag `v2.0.0` per the release-cutoff decision.

### Distribution: PyPi package + API

*(Source: drafts/methods_distribution_pypi_ci.md — PyPi section.)*

The 2020 release published the database as a Git repository and a Solr REST endpoint; both remain available. The 2026 release adds an installable PyPi package (`modelseed-biochem`, [TBD confirm exact package name]) that ships the database contents plus a documented Python API for programmatic access. Typical use cases include fetching a compound record by ID, walking a reaction's reagents, querying the thermodynamics ledger, and retrieving atom mappings. Installation and usage examples are documented in the repository README and mirror the shape of the existing Solr access patterns, so existing users can migrate incrementally.

Version pinning: the PyPi package version tracks the ModelSEED release tag (i.e. installing `modelseed-biochem==2.0.0` pulls the database contents as of the `v2.0.0` release tag). This lets downstream tools depend on a specific database state rather than the moving `master` head.

## Results

### Growth in compounds and reactions

*(Source: data/snapshot_2026-07-29.md — provisional numbers against dev HEAD `1111754`; refresh against v2.0.0 when tagged.)*

The database has grown substantially since the 2020 release, driven mainly by contributions from MetaCyc and the BioCyc-family databases (Table 1). Compound counts grew 34% (33,992 → 45,708), with structure coverage rising from 28,120 to 36,943 in absolute terms while the coverage rate slipped one percentage point (83% → 81%) because new compounds are ingested faster than the curation pipeline can add structures. Reaction counts grew 55% (36,193 → 56,012), with balanced-reaction counts growing 35% (25,457 → 34,370). The balanced-share of reactions regressed from 70% to 61%, reflecting the reaction-count growth outpacing balance-restoring curation; the structure-curation and protein-carrier-cofactor-standardisation pipelines described in Methods are designed to recover this gap, and the impact is quantified in the following section.

**Table 1 — Growth in compounds and reactions (2010–2026)**

| | 2010 | 2014 | 2020 | 2026 |
|---|---:|---:|---:|---:|
| Compounds | 16,275 | 27,694 | 33,992 | 45,708 |
| Compounds with structures | 13,821 (85%) | 19,605 (71%) | 28,120 (83%) | 36,943 (81%) |
| Compounds with generic groups | 1,261 (8%) | 1,402 (5%) | 4,416 (13%) | 7,071 (15%) |
| Reactions | 13,257 | 27,558 | 36,193 | 56,012 |
| Complete reactions | 11,338 (86%) | 7,898 (29%) | 27,991 (77%) | 44,877 (80%) |
| Balanced reactions | 10,263 (77%) | 17,264 (63%) | 25,457 (70%) | 34,370 (61%) |
| Reactions with generic groups | 1,988 (15%) | 2,939 (11%) | 9,772 (27%) | 15,357 (27%) |

**Table 2 — Per-source integration (2020 vs 2026)**

Per-source counts refreshed. Rhea integration doubled (8,786 → 17,477 reactions). ChEBI is added as a source column not present in the 2020 Table 3. [TBD decision: aggregate 11 BioCyc-family sources into one row as in 2020, or itemise separately as the current alias table does.]

| Source | 2020 compounds | 2026 compounds | 2020 reactions | 2026 reactions |
|---|---:|---:|---:|---:|
| ModelSEED | 33,958 | 45,708 | 36,193 | 56,012 |
| KEGG | 17,760 | 17,803 | 10,850 | 14,066 |
| MetaCyc | 19,138 | 25,740 | 21,830 | 31,805 |
| BiGG | 2,704 | 2,736 | 4,306 | 7,326 |
| MetaNetX | 30,858 | 30,896 | 23,758 | 29,323 |
| Rhea | — | 1,931 | 8,786 | 17,477 |
| ChEBI | — | 11,429 | — | — |

### Structure-curation improvements

*(Source: drafts/results_structure_curation_improvements.md — full text inlined below.)*

At the 2020 release, formula disagreement between source databases silently dropped the affected compound from `Unique_ModelSEED_Structures.txt`, cascading into a corresponding reaction becoming impossible to balance. The updated pipeline classifies each remaining conflict and routes it to one of three outcomes: auto-pick from the source variant that matches the historical `compound_*.json` formula (zero-disruption), curator-authored override with recorded rationale, or an entry in the mass-balance exclusion registry.

Applied to the priority-scope reaction set (~9,000 reactions used in v7.0 ModelSEEDTemplates and PlantSEED templates, comprising ~6,500 compounds), the pipeline resolved every previously-conflicting compound in three curation rounds during this release cycle. Round 1 (2026-07-03): 52 H-only conflicts auto-matched to the `compound_*.json` formula, plus two mass-balance-exclusion entries for cpd00766 (ascorbate radical) and cpd00849 (quercetin 3,3'-bissulfate) where no source formula matched. Round 2 (2026-07-04): 11 remaining priority-scope conflicts resolved — 9 auto-matched (spanning the alkyl-glycerol, sulfatide, Descarbamoylnovobiocin, and other classes), plus two additional mass-balance-exclusion entries for cpd27485 (Co(I)-corrinoid Fe-S protein) and cpd12396 (a pyruvate-kinase phosphoserine residue). Round 3 (2026-07-04): cpd00244 Ni²⁺ resolved via ChEBI 49786; the one remaining source variant (KEGG C00291 at Ni/0) was identified as elemental nickel mis-aliased to the ion, flagged for a separate alias cleanup.

After these rounds, the classifier reports 56 remaining formula-conflict compounds across 143 source-variant rows (42 H-only, 12 heavy-atom skeleton, 1 element-set, 1 single-variant). Zero of these fall in the v7.0 priority template compound set. The priority-scope formula-conflict backlog is empty, and the balanced-reaction rate on the priority set held at 8,683 / 9,000 across all three PRs (no regressions).

The `mass_balance_excluded.tsv` registry — introduced in this release cycle — records four compounds whose stored formula in `compound_*.json` cannot be reproduced from any available source structure. For these compounds the picked structure's SMILES and InChIKey populate `compound_*.json` (so downstream atom-mapping tools have something to work with), while the historical formula and charge are preserved unchanged. As a consequence, the reactions in which these compounds appear balance against the historical formula rather than a formula computed from the picked structure. The registry is small and expected to remain so.

For protein-carrier cofactor standardisation, an earlier commit in this release cycle added 47 pantetheine-inclusive acyl-ACP formula overrides in a single PR (verified as producing zero newly-imbalanced reactions). A follow-up analysis (Round 4, 2026-07-30) applied the same pipeline to the 260 open ACP name-matches: 179 survived a structural filter (`R` in formula or `*` in SMILES), 173 received a computed override proposal, and an iterative safe-filter converged at zero acceptable proposals — every candidate cascades into breakage through co-application effects. The zero-safe result reflects the tight interconnection of ACP-bound intermediates in specialty secondary-metabolism pathways (polyketide synthase chains, NRPS chains, actinorhodin, mupirocin, mithramycin), and zero of the 173 candidates touch any priority-scope reaction; systematic curation of these residuals via coupled-update reasoning is future work.

[TBD] Extensions to biotinyl-BCCP and lipoyl-carrier classes are pending as smaller PRs.

### Multi-source thermodynamics

*(Source: drafts/methods_multi_source_thermodynamics.md, Results portion — pending numbers.)*

**[NUMBERS PENDING]** Per-source coverage (compounds with ΔfG′, reactions with accepted ΔrG′) across eQuilibrator, rebuilt Group Contribution, dGpredictor retrained on ModelSEED, and OpenTECR. Head-to-head ΔrG′ comparison (correlation, uncertainty distributions) mirroring the 2020 Figure 2. How per-reaction-class heuristic overlays change direction assignments relative to the source values alone.

### Reaction similarity

*(Source: drafts/methods_reaction_similarity.md, Results portion — pending.)*

**[NUMBERS PENDING]** Refreshed similarity matrix statistics (density, clustering behaviour on curated reaction classes). Head-to-head comparison of external reaction-embedding foundation models with the selection defended on discerning-power grounds.

### Empirical study — reaction-direction heuristic sensitivity on ModelSEED v2 models

*(Source: drafts/results_reaction_direction_sensitivity.md — inlined below.)*

The most consequential open question for a mass-balance-first biochemistry database is not "is the reaction thermodynamically feasible" — that question is answered per reaction by the ΔrG′ ledger described in Methods — but rather "does the direction assignment materially change what a metabolic model built on top of this database predicts". The 2020 release described a single direction-assignment procedure (based on eQuilibrator-derived ΔrG′ with a group-contribution fallback and a reversibility rule set from Jankowski et al.). The 2026 release provides four direction sources (see Methods), and we present here the first systematic study of how the choice between them propagates through metabolic-model output.

**Model corpus.** We use the ModelSEED v2 draft-model corpus published in [10.1101/2023.10.04.556561](https://doi.org/10.1101/2023.10.04.556561) — [TBD N] draft genome-scale metabolic models spanning [TBD taxonomy summary]. Models are used as-published (no manual curation) so that the study measures direction-source sensitivity in isolation from other model-improvement effects.

**Direction sources compared.** Four sources, each producing a full reversibility labelling over the ModelSEED reactions used in the model corpus: (i) eQuilibrator-derived ΔrG′ with reversibility rules (the 2020 default); (ii) rebuilt group-contribution ΔrG′ with the same reversibility rules; (iii) dGpredictor retrained on ModelSEED with the same reversibility rules; (iv) [TBD — OpenTECR integrated as a fourth source, or heuristic-overlay variants of (i)-(iii)]. See Methods for the derivation of each.

**Per-model metrics.** For each `(model × direction source)` pair we compute: (a) predicted growth rate under the model's original medium via flux balance analysis; (b) essential-gene set derived from single-gene knockout FBA; (c) number of reactions marked infeasible under the direction source's constraints; (d) number of the 390 Biolog conditions from Ye et al. under which the model produces biomass; (e) [TBD any additional metric agreed with the co-authors].

**Aggregation and reporting.** We report: (1) a `(model × direction source × metric)` matrix as supplementary material; (2) per-metric distributional summaries (median + inter-quartile range across the model corpus, per direction source); (3) a summary figure with pairwise Jaccard similarity of essential-gene sets across direction sources, plus a per-model direction-sensitivity score defined as [TBD].

**Key findings [PLACEHOLDER — filled after the study runs]:**
- [TBD] direction source(s) produced growth-rate predictions within [TBD %] of the corpus median; [TBD] produced systematic deviations.
- Essential-gene Jaccard between the eQuilibrator and dGpredictor sources was [TBD] on median, indicating that [TBD interpretation].
- Under the Biolog condition set, [TBD] source predicted the largest number of feasible conditions, [TBD] the smallest; the difference of [TBD] conditions is (small / substantial) relative to the 390-condition baseline.

### Atom mapping coverage

*(Source: drafts/methods_atom_mapping.md, Results portion — pending.)*

**[NUMBERS PENDING]** Fraction of priority-scope reactions with mappings delivered by the Nikoloski collaboration, with an illustrative use case in atom-tracking for C, P, S.

### Database connectivity and ontological mappings

*(Source: drafts/results_db_connectivity_and_ontology.md — inlined below with placeholders.)*

The 2020 release quantified the biochemistry database's suitability for metabolic modelling by treating the entire database as a single reaction network and applying flux balance analysis. Two metrics were reported (Table 5 of the 2020 paper): the number of reactions capable of carrying nonzero mass-balanced flux ("functional reactions") assuming every extracellular metabolite is available, and the number of the 390 Biolog growth conditions under which the network is capable of producing biomass. The 2020 counts were 21,403 functional reactions and 355/390 Biolog conditions.

**Table 3 — Whole-database FBA refresh** (see Methods for procedure)

| | 2010 | 2014 | 2020 | 2026 [TBD] |
|---|---:|---:|---:|---:|
| Total reactions | 13,257 | 27,558 | 36,193 | 56,012 |
| Mass balanced | 10,263 | 17,264 | 25,457 | 34,370 |
| Reversible | 6,195 | 8,906 | 18,399 | [TBD] |
| Functional reactions | 8,505 | 21,917 | 21,403 | [TBD] |
| Functional Biolog conditions | 330 (85%) | 337 (86%) | 355 (91%) | [TBD] |

For the ontology mappings, the 2020 release described three ontology types (equivalent compound sets, lumped reaction sets, context-specific reaction sets) and demonstrated ontology-guided reconciliation reducing apparent-mismatch reactions between the *E. coli* iJR904 model and ModelSEED from 258 to 159. The 2026 refresh reports:

**Table 4 — Ontology growth 2020 → 2026 [numbers TBD]**

| Ontology type | 2020 entries | 2026 entries |
|---|---:|---:|
| Equivalent compound sets | [TBD] | [TBD] |
| Lumped reaction sets | [TBD] | [TBD] |
| Context-specific reaction sets | [TBD] | [TBD] |

iJR904 reconciliation refresh: apparent-mismatch reactions dropped from [TBD 2026 baseline] to [TBD after ontology] via [TBD dominant equivalence categories].

## Discussion

*(Source: drafts/discussion_and_availability.md.)*

Six years of community-driven curation have grown the ModelSEED biochemistry database by 34% in compounds and 55% in reactions. That growth has been driven mainly by contributions from MetaCyc and the BioCyc-family databases; KEGG's contribution is roughly flat over the interval. The balanced-reaction share of the database has regressed from 70% to 61% over the interval; the structure-curation pipeline described in Methods (particularly the protein-carrier cofactor standardisation) [TBD accounts for X percentage points of recovery on that regression once §5 completes].

**Broader lesson — protein-carrier cofactor standardisation as a general pattern.** The acyl-ACP work described in Methods is one instance of a general pattern: whenever a protein-covalent cofactor carries a reactive group between enzymes, some databases include the cofactor in the stored formula for the loaded state and others do not, and the inconsistency generates false-positive mass-imbalance flags in reactions handling the loaded and unloaded forms. We extend the same pattern to biotinyl-carrier (biotin-carboxyl-carrier protein) and lipoyl-carrier (H-protein of glycine cleavage, E2 subunits of α-ketoacid dehydrogenase complexes) classes; candidate extensions to covalent FMN/FAD, molybdopterin, and heme-c are surveyed and deferred to future work.

**Contrast — automated curation vs coupled-update residuals.** The bulk-add-safe pipeline that added 47 pantetheine-inclusive ACP overrides also identified 173 further candidates for which every proposal cascades into breakage via co-application effects. The dichotomy between the 47-safe and the 0-safe subsets — the low-hanging fruit vs the specialty-metabolism residuals whose reactions form tightly-coupled substrate chains — is itself a lesson: single-compound-at-a-time bulk curation reaches a ceiling that the next generation of tools will need coupled-update reasoning to break through.

**Open problems for the next release.**
- **Fragment-aware pKa/pKb.** The current biochemistry object stores pKa/pKb values only for atoms in the first InChI fragment; extending this to multi-fragment compounds requires both a schema change and a corresponding update to the protonation and thermodynamics pipelines.
- **EC number refresh.** The current ExPASy-name-based mapping between reactions and EC numbers is bound to have false and missing associations; a curation pass beyond name matching would improve annotation quality.
- **Obsolescence-leakage audit.** The database marks obsolete compounds and reactions as such, but some may still be inadvertently referenced by non-obsolete reactions or models; a systematic audit is warranted.
- **Direction-field removal.** With thermoreversibility computed per reaction and direction captured in the template, the reaction schema's `direction` field is arguably redundant; a removal pass has been considered but is [done / planned — TBD] as of the 2026 release.

**Where atom mapping and multi-source thermodynamics take the community next.** [TBD 2–3 forward-looking sentences.]

## Data Availability

- **GitHub repository:** `ModelSEED/ModelSEEDDatabase`, tagged as `v2.0.0` for the 2026 release.
- **PyPi package:** [TBD `modelseed-biochem==2.0.0`].
- **Solr REST endpoint:** retained from the 2020 release.
- **Web interface:** the ModelSEED website ([TBD final URL]) has been updated to include atom-mapping-aware reaction landing pages. KBase integration is retained.
- **Atom mapping endpoint:** [TBD delivery mechanism from the Nikoloski collaboration].
- **Continuous integration:** GitHub Actions workflow at `.github/workflows/` runs on every PR to `dev`.
- **License:** Creative Commons Attribution, matching the 2020 release; data derived from KEGG and MetaCyc remains subject to those databases' licenses.

## Author contributions, funding, acknowledgments

*(Deferred to submission time. Author list, Nikoloski credit format — see `PAPER_2026_PLAN.md` open decision #3.)*

---

## Manuscript-shape health check (2026-07-31)

| Section | Draft | Numbers | Status |
|---|---|---|---|
| Abstract | outline | [TBD] | needs full drafting after Results finalise |
| Introduction | ✓ prose | mostly done | reads well |
| M&M — Collation | [DRAFT PENDING] | numbers in `data/snapshot_2026-07-29.md` | paragraph needs writing |
| M&M — Integration | [DRAFT PENDING] | — | 2020 carryover |
| M&M — Provenance | [DRAFT PENDING] basics | — | extended piece already covered |
| M&M — Transport | [DRAFT PENDING] | — | 2020 carryover |
| M&M — Structure-curation pipeline | ✓ full prose | — | largest new Methods block; ready |
| M&M — Balancing | [DRAFT PENDING] | — | 2020 carryover |
| M&M — Multi-source thermodynamics | ✓ prose | [TBD numbers] | needs source-version capture |
| M&M — Reaction similarity | ✓ prose | [TBD list + results] | needs foundation-model bake-off |
| M&M — Atom mapping | ✓ prose | [TBD Nikoloski deliverables] | pending external delivery |
| M&M — Undetermined | [DRAFT PENDING] | — | 2020 carryover + fragment note |
| M&M — CI (GitHub Actions) | ✓ prose | — | needs workflow-file path |
| M&M — Release procedure | [DRAFT PENDING] | — | 2020 carryover + v2.0.0 note |
| M&M — Distribution (PyPi) | ✓ prose | [TBD package name] | needs PyPi publication |
| Results — Growth | ✓ full prose + tables | ✓ provisional numbers | ready pending v2.0.0 refresh |
| Results — Structure-curation improvements | ✓ full prose | ✓ | ready |
| Results — Multi-source thermodynamics | [NUMBERS PENDING] | [TBD] | numbers first |
| Results — Reaction similarity | [NUMBERS PENDING] | [TBD] | numbers first |
| Results — Reaction-direction empirical study | ✓ prose | [TBD numbers, figure, interpretation] | biggest gap; needs §6 done first |
| Results — Atom mapping coverage | [NUMBERS PENDING] | [TBD] | Nikoloski delivery |
| Results — Whole-DB FBA + Ontology | ✓ prose + table shells | [TBD numbers] | needs whole-DB FBA rerun |
| Discussion | ✓ prose | mostly done | some [TBD] specifics |
| Data Availability | ✓ prose | [TBD URLs] | ready pending publication |
| Author list / credit | deferred | — | submission-time |

**Word count (provisional first pass): ~[TBD final trim].** NAR Database Issue target ~4,000 words; this draft is over that (deliberately — trims happen after the second pass fills the [TBD] markers).
