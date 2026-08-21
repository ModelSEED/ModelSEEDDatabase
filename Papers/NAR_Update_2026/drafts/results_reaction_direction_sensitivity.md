# Results — Reaction-direction heuristic sensitivity (draft)

**Target section:** Results, "Empirical study — reaction-direction heuristic sensitivity on ModelSEED v2 models" (see `PAPER_2026_SKELETON.md`).
**Guide reference:** `PAPER_2026_GUIDE.md` §9.
**Status:** first pass — study design is fixed; numbers, figure, and conclusions are placeholders until the study runs.

---

The most consequential open question for a mass-balance-first biochemistry database is not "is the reaction thermodynamically feasible" — that question is answered per reaction by the ΔrG′ ledger described in the Methods section on thermodynamics — but rather "does the direction assignment materially change what a metabolic model built on top of this database predicts". The 2020 release described a single direction-assignment procedure (based on eQuilibrator-derived ΔrG′ with a group-contribution fallback and a reversibility rule set from Jankowski et al.). The 2026 release provides four direction sources (see Methods), and we present here the first systematic study of how the choice between them propagates through metabolic-model output.

**Model corpus.** We use the ModelSEED v2 draft-model corpus published in [10.1101/2023.10.04.556561](https://doi.org/10.1101/2023.10.04.556561) — [TBD N] draft genome-scale metabolic models spanning [TBD taxonomy summary]. Models are used as-published (no manual curation) so that the study measures direction-source sensitivity in isolation from other model-improvement effects.

**Direction sources compared.** Four sources, each producing a full reversibility labeling over the ModelSEED reactions used in the model corpus: (i) eQuilibrator-derived ΔrG′ with reversibility rules (the 2020 default); (ii) rebuilt group-contribution ΔrG′ with the same reversibility rules; (iii) dGpredictor retrained on ModelSEED with the same reversibility rules; (iv) [TBD — OpenTECR integrated as a fourth source, or heuristic-overlay variants of (i)-(iii)]. See Methods for the derivation of each.

**Per-model metrics.** For each `(model × direction source)` pair we compute: (a) predicted growth rate under the model's original medium via flux balance analysis; (b) essential-gene set derived from single-gene knockout FBA; (c) number of reactions marked infeasible under the direction source's constraints; (d) number of the 390 Biolog conditions from Ye et al. under which the model produces biomass; (e) [TBD any additional metric agreed with the co-authors — e.g. maximum-yield of a stated target metabolite].

**Aggregation and reporting.** We report: (1) a `(model × direction source × metric)` matrix as supplementary material; (2) per-metric distributional summaries (median + inter-quartile range across the model corpus, per direction source); (3) a summary figure with pairwise Jaccard similarity of essential-gene sets across direction sources, plus a per-model direction-sensitivity score defined as [TBD — e.g. one minus the mean pairwise Jaccard across all six source-pairs].

**Key findings** [PLACEHOLDER — filled after the study runs]:
- [TBD] direction source(s) produced growth-rate predictions within [TBD percent] of the corpus median; [TBD] produced systematic deviations.
- Essential-gene Jaccard between the eQuilibrator and dGpredictor sources was [TBD] on median, indicating that [TBD interpretation].
- Under the Biolog condition set, [TBD] source predicted the largest number of feasible conditions, [TBD] the smallest; the difference of [TBD] conditions is (small / substantial) relative to the 390-condition baseline.

**Interpretation.** [TBD paragraph tying the numerical findings to a takeaway about which direction source is preferable for which downstream use case, and what the community should watch for when choosing.]

---

## Open loose ends flagged during drafting

- Corpus size and taxonomy summary come directly from the bioRxiv paper — needs fetching and citing precisely.
- Whether the fourth direction source is OpenTECR-based or heuristic-overlay-based is an open scoping decision; if OpenTECR-based, coverage may be too sparse to be a full-corpus fourth column.
- The empirical study needs the multi-source thermodynamics work in §6 to be complete first (all four direction ledgers must exist before per-model FBA can be run).
- Whether direction-field removal from the reaction schema happens before or after this study runs is `PAPER_2026_PLAN.md` open decision #4.
