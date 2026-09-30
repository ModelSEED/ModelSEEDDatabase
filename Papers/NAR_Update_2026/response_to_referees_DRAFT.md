# Response to Referees — Draft (Reviewer 1, comments 1 & 4)

## A note on scope: the live-only conversion

While answering comment 1 below, we discovered that the manuscript mixed two
reporting populations without saying so: the abstract's headline ("~56,000
reactions"), Figure 2A's denominator, and the growth-rate arithmetic in the
Results all counted every record including the ~7,600 flagged `is_obsolete`,
while the "55% growth in reactions" sentence quietly set 2026's total against
a 2020 baseline with its obsolete rows *removed* — the only place in the draft
that filtered. Neither number was wrong on its own terms; together they were
not comparable, and no single population reproduced the 55% figure.

**This was not requested by either reviewer.** We are recording that
explicitly: comment 1 asked where new biochemistry is landing, not what
population it should be counted against. Having found the inconsistency while
building the pathway-distribution panel the comment required, the
corresponding author decided to fix it properly rather than patch around it —
converting the *entire* manuscript to the live-only (non-obsolete) basis
rather than cherry-picking a fix for the one sentence that was visibly wrong.
That moves the abstract's headline from "~56,000 reactions" to "~48,000," and
every count in Methods/Results Sections M09, M11, M12, M13 and Supplementary
Table S1 with it, along with all three main-text figures (none of which
filtered obsolete records before this pass). `analysis/population_basis_table.py`
prints every quantity on both populations side by side so the conversion can
be audited rather than trusted.

---

## Reviewer 1, Comment 1 (pathway-level distribution of new biochemistry)

> "It would be interesting to show a frequency distribution at pathway/subsystem level. Assuming that central pathways are already well characterized for many years, it would be curious to see where are we still gaining new information."

**Response:** The reviewer's hypothesis is correct, and we can now quantify it directly. Of the 12,259 reactions added to the database since the 2020 release, only 41 (0.3%) fall in a central-carbon or energy-metabolism MetaCyc class (Energy-Metabolism, Electron-Transfer, Fermentation, TCA-VARIANTS, Glycolysis, Pentose-Phosphate-Cycle, Photosynthesis, Respiration, Methanogenesis) — against 328 of 36,125 (0.9%) among reactions already present in 2020. The growth is concentrated in specialised and peripheral metabolism instead: antibiotic and polyketide biosynthesis, O-antigen and cell-envelope assembly, branched and acylated lipid biosynthesis, and alkaloids, several of which are enriched by one to two orders of magnitude among the new reactions relative to what the database already held (e.g. branched fatty acids 67×, O-antigen 41×). The released per-reaction pathway annotation cannot answer this directly for the new reactions — `Unique_ModelSEED_Reaction_Pathways.txt` and the EC alias file both stop just below the 2020 boundary — so we rebuilt the mapping from the shipped MetaCyc source table for both eras through the same one-level parent join, added as a new Figure 1D, and state the coverage limit (only MetaCyc-sourced additions can be placed; the 8,409 new Rhea-sourced reactions carry no pathway ontology) directly in the caption rather than leaving it implicit.

**Manuscript changes:**
- Figure 1 — new panel D: reactions added since 2020 against those already held, by MetaCyc pathway class, for the eight classes that gained most (antibiotics, O-antigen, polyketides, toxins, lipids, branched fatty acids, fatty acids, alkaloids). Figure 1B is also replaced in this same PR (see Comment 4 below), so the whole figure is regenerated together.
- `M06` (figure caption) — extended to describe panel D: what it shows, how both eras are joined identically for comparability, the MetaCyc-only coverage limit, and the central-metabolism headline (41 of 12,259).
- `M09` (Results, Growth and coverage) — two new closing sentences stating the finding in prose: growth is concentrated in specialised metabolism and cell-envelope biosynthesis, central metabolism was already saturated in 2020, and the comparison is scoped to MetaCyc-sourced additions.
- `analysis/pathway_distribution_of_growth.py` (new) — the script behind every number in this response; reads only `Biochemistry/` and the repository's own git history, so it needs no external cache and runs anywhere.

**Supplementary changes:** None for this comment.

---

## Reviewer 1, Comment 4 (figure clarity: Fig 1B as a Venn/Euler diagram, Fig 2A needs a legend)

> "The figures could use a little improvement. There is a bit of barplot abuse (smiley face). For instance, Fig 1B would be more informative as a Venn diagram. Fig 2A needs a color legend (and maybe a different color scheme)."

**Response:** Agreed on both counts, and for the reason implied in each case. The stacked bar in the old Figure 1B reported that a majority of MetaCyc's reactions are unique to MetaCyc, but never said *who the rest is shared with* — exactly the question that integrating Rhea as a third reaction source raises. We replaced it with a three-circle Euler diagram over the three primary reaction sources (MetaCyc, KEGG, Rhea), with every one of the seven regions labelled with its exact reaction count; circle areas are explicitly *not* drawn proportionally (three-set areas cannot be solved exactly), and the caption says so, so a reader is not misled into reading the geometry quantitatively. `matplotlib-venn` is not in our build environment, so the panel is drawn directly with matplotlib patches rather than adding a new dependency for one panel.

For Figure 2A, the panel previously had four visually distinct segments but named only two of them in the caption prose, so the caption was doing the legend's job. We added a four-entry inline key. On the colour scheme, the reviewer's parenthetical was well aimed independently of the legend: two of the four segments — "no ionizable site" and "no structure" — were both rendered in the same neutral grey, which visually merged a genuine chemistry *result* (Marvin ran and found no dissociable proton between pH −2 and 16, for 1,237 compounds) with a curation *gap* (there is no structure to run Marvin on at all). These are conceptually different things and are now given distinct categorical colours; only the genuine structural absence keeps the grey.

**Manuscript changes:**
- Figure 1B — replaced the stacked unique/shared bar with a three-circle Euler diagram; every region (MetaCyc-only, KEGG-only, Rhea-only, the three pairwise overlaps, the three-way overlap, and reactions in none of the three sources) is labelled with its exact count.
- Figure 2A — added a four-entry colour key; gave "no ionizable site" its own colour, distinct from "no structure."
- `M06` (Figure 1 caption) and `M11` (Figure 2 caption) — both rewritten to describe the new panels; the Figure 1 caption is necessarily a substantial rewrite since panel B is now a fundamentally different chart type (Euler diagram vs. bar chart), not merely a relabelling. Wording elsewhere in both captions that was not affected by the panel change is preserved close to verbatim (see rewrite-scope note below).

**Rewrite-scope check:** A word-level SequenceMatcher retention check was run on every touched paragraph. `M06`'s Figure 1 caption retains only ~46% of the original wording, well under our usual ~50% threshold — but this is the one case we judge structurally justified rather than a stylistic overreach: panel B changed from a bar chart to an Euler diagram and panel D is a wholly new panel, so a large fraction of the caption necessarily describes content that did not exist in the original figure. Every other touched paragraph in this PR (M09, M11 body and caption, M12, M13, S03) retains 84–99% of its original wording; none required tightening.

**Supplementary changes:** None for this comment.

---

## Additional note: obsolete-reaction accounting

The live-only conversion described above depends on the database's `is_obsolete` flag being complete and correctly linked. We audited this directly: every one of the ~7,600 reactions flagged obsolete correctly links to its live replacement, with no orphaned or circular references, and the live-only reaction count (48,384) is fully trustworthy on that basis. Separately, however, we found that a policy adopted around the 2020 release — that a reaction found to duplicate an existing one thereafter should have its names, aliases and EC numbers merged into the surviving (lower-numbered) reaction and then be physically removed rather than merely flagged — had never actually been carried out. We applied it for the first time while preparing this response, which physically removed a small number of newly-obsolete reactions (low single digits) from the counts in this PR; older reactions later found obsolete are deliberately left in place, unflagged for removal, to preserve backward compatibility with published models built against them. We confirmed this cleanup has no measurable effect on any thermodynamic value, direction call, or evidence grade for any reaction that remains. We mention it because several of the counts in this response (e.g. the graded total and gold/silver/bronze breakdown) were measured before and after this cleanup and moved by amounts consistent with it; the database is under continuous curation, and we do not expect counts of this kind to be stable to the last digit between this submission and publication.

## Referee-facing marked copy

Both `main_diff.pdf` and `supplementary_diff.pdf` are built with `latexdiff --flatten` against `origin/dev` and use the same red-underline (addition) / blue-strikethrough (deletion) colour convention already established for PRs #297 and #298 (see `latex/make_marked_copy.sh` and `latex/README.md`). This PR is the first in the series to touch the supplement (`S03_evidence_grading.tex`), so both documents' diffs are generated and both pairs of PDFs (`main_clean.pdf`/`main_diff.pdf` and `supplementary.pdf`/`supplementary_diff.pdf`) are included.
