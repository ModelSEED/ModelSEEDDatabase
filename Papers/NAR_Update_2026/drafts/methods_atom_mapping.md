# Methods — Atom mapping (draft)

**Target section:** Materials and Methods, "Atom mapping" (see `PAPER_2026_SKELETON.md`).
**Guide reference:** `PAPER_2026_GUIDE.md` §8.
**Status:** first pass filled in against the committed integration; illustrative-use-case figure and formal group/institution attribution still `[TBD]`.

---

Atom mapping — the assignment of individual reactant atoms to their product-side counterparts across a mass-balanced reaction — is a new capability in the 2026 release. Mappings enable atom-level tracing of carbon, phosphorus, sulfur, and other atoms through a metabolic network, and are increasingly required for downstream applications including 13C metabolic flux analysis, pathway-degradation prediction, and cofactor-usage attribution.

**Generation.** Mappings for the ModelSEED biochemistry are generated through a collaboration with the [Nikoloski lab, TBD confirm institutional attribution]. The underlying mapper is the **Reaction Decoder Tool** (RDT) v4.0.0 (Rahman et al., <a href="https://doi.org/10.1093/bioinformatics/btw096">10.1093/bioinformatics/btw096</a>), driven per-reaction from ModelSEED SMILES by the **UniversalRDT** wrapper maintained by Sebastian Huhn (<a href="https://github.com/sebahu/UniversalRDT">github.com/sebahu/UniversalRDT</a>, `ModelSEED/` subdirectory). The wrapper flattens each `reaction_*.tsv` equation to a SMILES-based reaction, invokes RDT with per-reaction timeout, and collects the atom-atom-pair output into a flat `all_mapping.txt` file. A subsequent row-level filter (`Scripts/Structures/Rebuild_AtomMappings_from_raw.py` in this repository, contributed upstream to `unite_and_filter_mappings.sh`) drops rows that cannot be salvaged and keeps every row that can — see below.

**Row-level filter.** RDT's raw output for a given reaction can contain rows that don't match the canonical single-pair form. Three common patterns require handling: (i) *run-on chains* like `A=B=C=D` that appear when RDT cannot uniquely resolve a symmetric group (the two O of CO₂ in a decarboxylation is the archetypal case); (ii) *dangling orphans* — a single atom reference with no `=` partner; (iii) *cross-element pairs* (e.g. `O=S`) from near-isomorphic subgraph confusion in the alignment. The upstream shell filter originally rejected the entire reaction whenever any single row was malformed; the ModelSEED filter drops only the bad rows at the row level and preserves the reaction's valid same-element pairs, expressing the element-pair constraint by construction (same-element check on both endpoints, 1-2 character element regex to accept two-letter symbols such as Cl, Fe, Mg, Zn, Hg). The row-level filter has been contributed back to `sebahu/UniversalRDT` upstream.

**Confidence tag.** To preserve information about which reactions the filter had to salvage from vs which came through clean, every mapped reaction carries an `atom_mapping_confidence` field alongside its `atom_mapping` array. A reaction is `clean` when every raw RDT row for it was already a canonical single-pair same-element row (no salvage was needed); `salvaged` when at least one row was a chain, orphan, cross-element pair, or malformed. The distinction matters: as one collaborator noted, *"as soon as there is one problem [in an RDT reaction], there are likely more (hidden) problems"* — so a salvaged mapping may carry correct-looking rows whose partner atoms are subtly misassigned elsewhere in the same reaction. Reachability / neighborhood-level applications are largely fine using both classes; mechanism-level tracing (13C flux, exact atom fate) should filter to `clean` only.

**Storage.** Mapped reactions are stored under `Biochemistry/Structures/AtomMappings/`. The flat pair file (`all_mapping_no_problem.txt`), the reaction-ID list (`rxns_no_problems.txt`), the confidence tags (`rxns_confidence.tsv`), the raw superset (`all_mapping.txt` for reproducibility and re-filtering), and the non-attempted / SMILES-gap / wildcard-SMILES metadata files are all committed. The corresponding `reaction_*.json` records carry `atom_mapping` (a list of `cpdA:E#N=cpdB:E#M` strings) and `atom_mapping_confidence` (`clean` | `salvaged`) fields. Downstream consumers can access mappings through the same Solr endpoint or the new PyPi API (see Methods "Distribution").

**Coverage.** As of the 2026 snapshot: **32,877 of 56,012 total ModelSEED reactions carry a clean atom mapping (59%)**, of which 25,058 (76%) are `clean` and 7,819 (24%) are `salvaged`. Coverage on the priority-scope reaction set (v7.0 ModelSEEDTemplates ∪ PlantSEED_v3 role assignments; 9,125 reactions total): **7,378 (80.9%)** mapped. Coverage on the plant reconstructions in `plantseed-v3` (the Arabidopsis TAIR10 model as the illustrative case; 782 unique base reaction IDs): **728 (93.1%)** mapped.

Reactions without a mapping fall into four classes:
- **Compound(s) lack a SMILES structure** (18,621 reactions; 1,611 in priority scope) — RDT cannot attempt these.
- **Wildcard SMILES** (1,725 reactions; 42 in priority scope) — the reaction's SMILES contains a `*` atom (unspecified R-group), fundamentally unmappable without picking a concrete placeholder.
- **RDT ran but no salvageable rows** (1,054 reactions; 57 in priority scope) — all raw rows were malformed or element-mismatched.
- **Pipeline timed out** (1,735 reactions; 37 in priority scope) — RDT did not complete even under an extended per-reaction budget.

**Illustrative use case: PlantSEED biomass reachability from CO₂.** As one concrete carbon-tracing demonstration, we build the undirected carbon-atom graph induced by the mapped reactions in the PlantSEED role set (1,058 reactions), seed it at CO₂ (cpd00011) and bicarbonate (cpd00242), and enumerate the carbon-containing biomass components reachable via carbon-atom transfers. Of the 76 carbon-containing biomass components, **74 are reachable** using the clean atom mapping alone. The two unreached components — glucotropaeolin (cpd28704) and sinalbin (cpd34716) — are a PlantSEED curation gap (no producing reaction is encoded for the benzenic glucosinolate branch), not an atom-mapping gap. Reference figure: [FIGURE TBD — the reachability graph, annotated by the six biomass components that specifically require the row-level filter's salvage recovery: Biotin, Leucine, Lysine, Phosphopantetheine, Thiamin diphosphate, UDP-Xylose].

---

## Open loose ends flagged during drafting

- Formal group / institution attribution (Nikoloski lab affiliation with Sebastian Huhn's UniversalRDT contribution) needs confirming before submission.
- Co-authorship for the atom-mapping collaborators is `PAPER_2026_PLAN.md` open decision #3.
- Illustrative-use-case figure needs generating from the reachability graph — recommended one page showing the 74 reachable biomass components colored by distance-from-CO₂ and the six specifically-filter-recovered ones highlighted.
- InChIKey-mismatch audit (Sebastian's `species_id_inchikey.diff.txt` identifies 129 compounds whose stored InChIKey doesn't match what OpenBabel regenerates from the same SMILES) is a companion structure-curation follow-up that should get a line in the Discussion.
