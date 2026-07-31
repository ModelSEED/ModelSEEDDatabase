# Methods — Atom mapping (draft)

**Target section:** Materials and Methods, "Atom mapping" (see `PAPER_2026_SKELETON.md`).
**Guide reference:** `PAPER_2026_GUIDE.md` §8.
**Status:** first pass, contingent on Nikoloski-lab delivery.

---

Atom mapping — the assignment of individual reactant atoms to their product-side counterparts across a mass-balanced reaction — is a new capability in the 2026 release. Mappings enable atom-level tracing of carbon, phosphorus, sulfur, and other atoms through a metabolic network, and are increasingly required for downstream applications including 13C metabolic flux analysis, pathway-degradation prediction, and cofactor-usage attribution.

Mappings for the ModelSEED biochemistry are being generated through a collaboration with the [TBD group name] group at the [TBD institution] (the Nikoloski lab). Mapping methodology: [TBD; will be described briefly with reference to the collaborators' methods paper]. Delivery format: [TBD — likely RXN files or a tabular format per reaction].

Mapped reactions are integrated into ModelSEED under `Biochemistry/AtomMappings/` (or equivalent, [TBD]). Each mapped reaction is stored alongside its ModelSEED reaction ID; downstream consumers can access mappings through the same Solr endpoint or the new PyPi API (see Methods "Distribution").

Coverage as of the 2026 snapshot: [TBD fraction of priority-scope reactions with mappings delivered]. Coverage of the full ModelSEED reaction set: [TBD]. Reactions without a mapping fall into three classes: (i) reactions whose reactant compounds do not yet have a stored structure (necessarily unmappable); (ii) reactions with `R` groups whose mapping requires expansion to specific chain lengths; (iii) reactions whose complexity (e.g. very long acyl chains, PKS/NRPS iterative modules) exceeds the mapper's practical scope.

**Illustrative use case.** [TBD — pick one central-metabolism pathway and demonstrate carbon-atom tracing from a labeled input through the mapped reactions]. Reference figure: [FIGURE TBD].

---

## Open loose ends flagged during drafting

- Everything about the Nikoloski-lab delivery is pending — group name / institution formal identification, exact mapping algorithm, delivery format, coverage statistics.
- Storage path (`Biochemistry/AtomMappings/`) is a placeholder; may end up under a different directory depending on file format.
- Attribution / co-authorship for the Nikoloski-lab collaborators is `PAPER_2026_PLAN.md` open decision #3.
- The illustrative use-case pathway needs picking; central carbon metabolism (glycolysis + TCA cycle) is an obvious default.
