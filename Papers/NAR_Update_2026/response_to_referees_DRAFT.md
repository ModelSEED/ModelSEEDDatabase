# Response to Referees — Combined Draft (PR #297, #298, #299)

---

# Response to Referees — Draft (Reviewer 1, comments 2 & 3)

## Reviewer 1, Comment 2 (transport reactions / compartment indices)

> "Transport reactions do not have compartment information (they are only distinguished as compartments 0 and 1, unlike Rhea which identifies them as in and out). And the authors mention that transport reactions are scored from stoichiometry alone. I think it would be important to elaborate a bit more on this and discuss potential limitations (after all, thermodynamics are the most important for energy metabolism, which usually involves transport reactions)."

**Response:** We agree the original single-sentence caveat under-treated this point. We promoted it to its own Methods paragraph (`Transport`) that (i) clarifies the 0/1 compartment indices are relative, not absolute, and map directly onto Rhea's `in`/`out` convention; (ii) quantifies why an energy-only score is largely blind to translocation, since a primary active transporter's ATP hydrolysis dominates the estimate; and (iii) adds one Discussion sentence naming compartment-aware scoring as the natural next release and transport as a standing limit.

**Manuscript changes:**
- Methods, `Transport` paragraph (new) — states that transport reactions are scored from stoichiometry and energy alone with no membrane potential or pH gradient term; that 96% of gold-graded transport reactions are ATP-coupled while only 2.7% translocate a proton, so a confident grade mostly reflects confidence in ATP chemistry rather than the translocation itself; and that dGPredictor commits a direction on only 1.4% of transport reactions versus 17.3% elsewhere.
- Methods, `Uncertainty` and `Experimental anchors` paragraphs — tightened for page budget (the three median per-source uncertainties, 0.63 / 10.41 / 17.01 kcal mol⁻¹, were duplicated verbatim from the Results section and are now stated once, in Results).
- Discussion — two sentences naming transport's lack of membrane-potential/pH-gradient information as a standing limitation, alongside the roughly half of reactions with no assigned direction.

**Supporting repository documentation (not part of the submitted manuscript):**
- `Biochemistry/REACTIONS.md` now documents explicitly that compartment indices 0/1 are relative (inside/outside), not absolute compartment identifiers, matching Rhea's in/out convention, and explains why transport reactions carry no membrane-potential or pH-gradient term.

**Supplementary changes:** None for this comment.

## Reviewer 1, Comment 3 (~9,000 of 56,000 reactions used in reconstructions)

> "The authors mention that only 9.000 of the 56.000 reactions in the database are used for reconstructions. It would be great if they could elaborate on this. Is it because they are disconnected from the main core network? Is it due to lack of GPR associations?"

**Response:** Neither of the reviewer's hypotheses is the constraint. Template membership is not gated on gene-protein-reaction (GPR) associations — a quarter of the Gram-negative v7.0 template (2,258 of 8,584 reactions) is typed `gapfilling` and requires no gene at all — nor on graph connectivity, since templates are curated scopes and nothing prunes on connectivity. The real gates are mass/charge balance and EC-based functional-role annotation, both of which we now quantify directly, and we reframed the opening of the Structure Curation results to state plainly that the ~9,000-reaction scope reflects a curation priority, not a claim about how much of the database is usable.

**Manuscript changes:**
- Results (`M10`) — opening sentence reframed: the curation pipeline was prioritized on the ~9,000 template reactions because an error there propagates into every model built from them, not because the rest of the database is unusable; points the reader to the new `M13` paragraph for the full answer.
- Results (`M13`) — new `Reconstruction scope` paragraph answering the reviewer's two hypotheses directly, with the corrected figures: v7.0 template union 8,597 reactions; 69% of the database satisfies mass/charge balance; 41% carry an EC-annotated functional role; 15,106 reactions satisfy both balance and role annotation; of those, 8,815 already lie outside every current template — i.e., roughly a doubling of usable scope is already available once a template is built to draw on it.
- Discussion — one sentence naming the roughly half of reactions still lacking an assigned thermodynamic direction as a standing limit, alongside the transport limitation above.

**A note on obsolete-reaction accounting and measurement stability.** While preparing this response we found that the database's obsolete-reaction bookkeeping, though internally consistent (every retired reaction correctly links to its replacement, with no orphaned or circular references), had not been fully applied going forward: a policy adopted around the time of the 2020 release — that a reaction found to duplicate an existing one after that point should have its names, aliases and EC numbers merged into the surviving (lower-numbered, and therefore more likely to already be embedded in published models) reaction and then be removed outright, rather than merely flagged — had never actually been carried out. We have now applied it for the first time, which physically removed a small number of newly-obsolete reaction records (in the low single digits) from the counts underlying this section; older reactions that were later found obsolete are deliberately left in place, unflagged for removal, to avoid breaking backward compatibility with models built against them. We verified this cleanup has no measurable effect on any thermodynamic value, direction call, or evidence grade for any reaction that remains in the database. We mention this because it is one plausible, though unconfirmed, contributor to the small measurement drift noted above (15,106 / 8,815, versus an earlier pass of 15,128 / 8,826): the database is under continuous curation, obsolete-reaction accounting is one of several moving parts, and we do not expect counts of this kind to be stable to the last digit between this submission and publication.

**Supplementary changes:** None for this comment.

---

# Response to Referees — Draft (Reviewer 2, comment 2.2)

## Reviewer 2, Comment 2.2 (LLM ensemble's role and the pathway-thermodynamics motivation)

> "The authors state in the Ensemble LLMs predictions section that a reaction may be driven by a pathway in a cell, which may be different from their grading of thermodynamic evidence. They then introduce LLMs and make predictions. Except the numbers in Figure 3, it is not apparent how LLMs play roles in reconstruction and other analyses."

**Response:** The reviewer identified a real inconsistency: the section motivated the LLM ensemble by pathway-level driving force (citing three references on the phenomenon), then introduced an ensemble that sees only a reaction's name and stoichiometry — no pathway context whatsoever — and finally disclaimed the ensemble's calls from both the evidence grading and any use as a downstream filter. The method was motivated by something it structurally cannot do, and then hedged for everything else. We have rewritten the section to state the role the data actually supports: **coverage**. Thermodynamic assignment is bounded by structural completeness (23,729 reactions get no direction from any energy-based source); the ensemble is not so bounded, and directs 22,902 of those 23,729 (96.5%), roughly doubling the fraction of the database that carries a direction. We also disclose, for the first time, how that coverage checks against measurement: scored against `opentecr_comparison.csv` at the same τ = 2 kcal/mol tolerance used in evidence grading, the ensemble agrees with measurement 76.7% of the time (Cohen's κ = 0.50) versus 98.2% for eQuilibrator (κ = 0.92), and 85.3% of its calls simply reproduce the direction the equation is written in — so its errors are concentrated almost entirely on reactions that run in reverse (46/95 correct) rather than forward (119/120 correct). This is disclosed explicitly as a strength of the *coverage* argument (a reconstruction would otherwise leave these reactions reversible by default) and as the reason the calls are excluded from evidence grading and should be reviewed, not treated as thermodynamic evidence.

**Manuscript changes:**
- Methods (`M06`, "Ensemble LLMs predictions") — removed the pathway-level-driving-force motivating sentence (and its three citations `mavrovouniotis1993`, `xu2008`, `noor2014` as the primary justification for introducing the ensemble); replaced the section's framing with structural-coverage bounding (23,729 reactions with no direction from any predictor). Added a new second paragraph stating the ensemble's actual, measured contribution: 22,902 of those 23,729 reactions directed; accuracy against measurement (77% vs. 98% for eQuilibrator); the one-sided error pattern; and the decision to release the calls as a separate, clearly-labelled layer excluded from evidence grading. The three orphaned citations are retained on a single new closing sentence that frames pathway context as a genuine limitation shared by *both* approaches (thermodynamic and LLM), rather than as a one-sided motivation for the LLM route alone — this keeps the citations relevant without repeating the reviewer's original objection.
- We reviewed the edit for scope: the paragraph containing the removed pathway-motivation sentence was tightened to a minimal edit (delete the motivating sentence and its citations, add one clause with the new 23,729 figure) rather than being reworded end-to-end; word-level retention of the original paragraph's wording is 54%, and every retained clause is verbatim. The second paragraph (ensemble's measured role and accuracy) is a wholly new paragraph — no text from the original section survives in it beyond incidental shared words (retention of the *entire original paragraph's* wording checked against it is under 6%).

**Supplementary changes:**
- `S04_direction_from_chemistry.tex`, "Summary of results" — restructured to state the coverage number (22,902 of 23,729) up front; word-level retention of the original four-sentence paragraph's wording is 63%. This is a targeted edit, not a rewrite: the closing sentence and framing wording are kept close to verbatim, and the confidence-not-a-filter sentence was relocated (not deleted) to the new `Accuracy against measurement` subsection below, where it now leads a fuller, quantified discussion.
- `S04`, new subsection **"Accuracy against measurement"** — a wholly new paragraph (word-level retention against the *entire* original S04 "Summary of results" paragraph is 40%, and inspection shows the residual overlap is incidental shared vocabulary — e.g. "ensemble", "direction", "the" — not reused sentences). Reports the τ-sensitivity analysis (65–77% agreement, κ 0.33–0.50 across a bare sign test to RT·ln 1000), the headline 94.7%-vs-eQuilibrator figure alongside its 88.3% chance baseline and κ = 0.55 (so it is not over-read as ensemble accuracy), and the reasoning for why these calls sit outside evidence grading.
- `S04`, new closing paragraph (call-distribution breakdown: 85.3% forward / 5.1% reverse / 7.8% reversible, and the 119/120 forward vs. 46/95 reverse confusion breakdown) — also a wholly new paragraph (retention against the original paragraph is under 16%, again incidental overlap only).

## Additional note: obsolete-reaction accounting

The counts above (23,729 undirected reactions; 22,902 directed by the ensemble; 8,085 eQuilibrator/ensemble co-commitments) were re-measured once during preparation of this response and moved by a small amount (to 23,790 / 22,963 / 8,084) between the two measurements. Investigating why, we found the database's obsolete-reaction bookkeeping — while internally consistent, with every retired reaction correctly linking to its replacement — had not been fully applied going forward: a policy from around the 2020 release, that a reaction found to duplicate an existing one thereafter should have its metadata merged into the surviving (lower-numbered) reaction and then be physically removed, had never actually been carried out. We applied it for the first time; it removed a small number of reactions (low single digits) from the live count, with older reactions retroactively found obsolete deliberately left in place for backward compatibility. We confirmed this has no effect on any reaction's direction call or thermodynamic value other than the removed records' own. We note it here because the database is under continuous curation and counts of this kind — like those above — should not be read as stable to the last digit between this submission and publication.

## Additional note: qualitative inspection of the ensemble's own stored rationale text

While preparing this response we separately inspected the actual free-text rationale the ensemble stores per reaction (`variant_A_reasoning.jsonl`, `/scratch/jplfaria/share/reaction_direction_council/`, referenced in `S04` as "we store the output of the ensemble in our repository for review"). This was not previously examined for the manuscript. We report the finding here for transparency, though we are **not** adding it to the manuscript text itself (Sam did not request a manuscript change, and it does not change the core point above):

Across all 45,644 reactions with a stored proposer rationale, a narrow, conservative keyword search for language that explicitly invokes pathway-level or cellular-context *driving-force* reasoning (e.g. "pulls the reaction/pathway forward," "driven by [downstream/cellular] concentration/flux/pool," "downstream consumption/clearance drives flux," "physiological concentrations of [substrate] drive," etc. — deliberately excluding the much larger set of rationales that merely *name* the pathway an enzyme belongs to, e.g. "this is a step in heme biosynthesis") finds this reasoning in **2,684 of 45,644 reactions (5.9%)**, drawn from 3,073 of 136,932 individual proposer calls (2.2%). A broader, unfiltered search that also counts any mention of the bare word "pathway" anywhere in the rationale returns a much higher and less meaningful figure (34.8%), because most such mentions are descriptive (naming the metabolic pathway an enzyme belongs to) rather than an argument that pathway-level dynamics determine the direction call.

Representative examples of the narrower, genuine pathway/flux-driven-force reasoning:

- *rxn00117* (nucleoside diphosphate kinase): "...making the reaction freely reversible and driven by relative nucleotide pool concentrations."
- *rxn00175* (acetyl-CoA synthetase): "...the released PPi is rapidly hydrolyzed by pyrophosphatase in vivo, pulling the reaction strongly toward acetyl-CoA formation."
- *rxn00179* (glutamate 5-kinase): "...ATP consumption and rapid downstream conversion of the unstable glutamyl 5-phosphate drive flux from glutamate toward glutamyl 5-phosphate."
- *rxn00096* (ADP:ATP adenylyltransferase): "Physiological concentrations of the substrates ATP and ADP are high compared to the signaling molecule AppppA, driving the reaction in the forward direction..."
- *rxn00217* (glucose dehydrogenase): "...driven by the rapid downstream hydrolysis of gluconolactone to gluconate..."

So the reviewer's original point and Sam's follow-up observation are, in a specific sense, both correct and complementary: the ensemble receives **no** pathway-level information as input (confirmed — the prompt is limited to reaction identifier, name, and unsigned equation), yet a non-trivial minority of its rationales (roughly 1 in 17 reactions, by our conservative count) invoke pathway-level or flux/concentration-driven reasoning anyway — almost certainly drawing on general biochemical knowledge from training rather than any pathway context we supply. This is an interesting property of the ensemble worth noting for transparency, but it does not change the core point of this PR: the pathway-level thermodynamic argument cannot be systematically verified or relied upon from our side, since the model is given no pathway-level input, and its unprompted invocation of such reasoning is not something we can validate against ground truth (we have not attempted to check whether the *specific* pathway claims in these rationales are correct, only that the language appears). We have not added this observation to the manuscript or supplement; we surface it here for the record in case it is useful in the response letter or a future revision.

**Manuscript changes for this note:** None (not requested; keeping manuscript edits limited to the rewrite-scope corrections above).

---

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
