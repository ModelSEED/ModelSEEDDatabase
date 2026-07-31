# Methods — Multi-source thermodynamics (draft)

**Target section:** Materials and Methods, "Multi-source thermodynamics" (see `PAPER_2026_SKELETON.md`).
**Guide reference:** `PAPER_2026_GUIDE.md` §6.
**Status:** first pass with `[TBD]` markers for the numerical results and per-source configuration details.

---

The 2020 release combined two thermodynamics sources: eQuilibrator as the primary estimator of ΔfG′ and ΔrG′, and a group-contribution (GC) method as a fallback for compounds and reactions eQuilibrator could not estimate. For this update we integrate four sources, each providing an independent estimate wherever its coverage permits.

**eQuilibrator refresh.** [TBD version of eQuilibrator] was run against the 2026 compound set at pH 7.0, ionic strength 0.25 M and temperature 298.15 K, matching the 2020 conditions. Coverage: [N compounds with ΔfG′, N reactions with accepted ΔrG′; TBD].

**Group Contribution rebuilt from MFAToolkit.** The 2020 fallback GC method has been rebuilt from the MFAToolkit implementation against the current compound set, so the group definitions and parameter fits reflect the enlarged compound base rather than the 2020 snapshot. Coverage: [TBD].

**dGpredictor retrained on ModelSEED.** dGpredictor (an ML-based ΔG estimator) has been retrained using the ModelSEED-curated compound structures as training data. This produces a ModelSEED-native estimator rather than one trained on a competing biochemistry database, and lets us report ΔG estimates for compounds the other three methods cannot cover. Training procedure and hyperparameters: [TBD]. Coverage: [TBD].

**OpenTECR integration.** OpenTECR is a community effort to standardize experimentally-measured thermodynamic constants; we integrate its values where a mapping between OpenTECR reaction identifiers and ModelSEED reaction IDs exists. Coverage: [TBD].

**Direction assignment.** For each reaction with an accepted ΔrG′ (from any source, chosen with per-compound-consistency constraints described in the 2020 paper's Materials-and-Methods "Computation of thermodynamic properties" subsection), we apply the same reversibility rules used in 2020 (based on Jankowski et al.) to assign one of `=` (reversible), `>` (irreversible left-to-right), or `<` (irreversible right-to-left). On top of this we add **per-reaction-class heuristic overlays** — hard rules for classes where the standard thermodynamic assignment is known to disagree with textbook biology. Rule set: [TBD; discussed in the working session but not yet codified in a repo script]. Overlays are applied after the base assignment and their scope is documented per class.

**Combined direction ledger.** For every reaction, the pipeline records: the ΔrG′ estimate from each of the four sources (with per-source uncertainty), the base direction assignment from each source, the heuristic overlay if applicable, and the final direction. This ledger is the input to §9's empirical study of how the direction-source choice affects metabolic-model output.

---

## Open loose ends flagged during drafting

- Version numbers for eQuilibrator, dGpredictor, and OpenTECR need to be captured for reproducibility. Add to the Methods "Collation" subsection alongside the other source-database versions.
- The MFAToolkit-based GC rebuild's exact provenance (which fork / commit / parameter set) needs recording.
- The heuristic overlay rule set needs codifying as a repo script (probably a YAML or TSV file listing patterns and their direction assignments) before the paper's claim that overlays are "applied globally" is defensible.
- Whether OpenTECR gets its own Methods paragraph or is folded into the eQuilibrator refresh is an open decision (see `PAPER_2026_PLAN.md` open decision #2).
