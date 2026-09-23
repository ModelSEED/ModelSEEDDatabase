# Reviewer response plan — NAR 2026 update

Five comments (four from reviewer 1, one from reviewer 2). For each: where it
comes from in the manuscript, what the data actually says, and the proposed
change. Numbers below are measured on the manuscript's own population (all
records, obsolete included) against `upstream/dev` at `41b20c21`, with two
scripts added for this purpose:

    ~/Documents/py_venv/bin/python analysis/pathway_distribution_of_growth.py
    ~/Documents/py_venv/bin/python analysis/review_transport_and_llm.py

The changes described here are implemented across four pull requests against
`dev`; each PR body names the comment it answers. Draft text in this document
is the intent — the PR is the authority where they differ.

---

## Read this first: one problem the reviewers did not raise

### The population question, and a false alarm I raised first

**Correction, recorded because it shaped an earlier draft of this document.**
I first reported that the shipped tree no longer reproduced the manuscript's
grade counts — 27,355 graded against the published 33,099, and 365 anchors
against 806 — and attributed it to the 09-11 and 09-15 rebuilds. That was
wrong. It was my own filter: I had excluded `is_obsolete` records and the
manuscript does not. Counting every record, the tree reproduces the paper
**exactly**: 33,099 graded, 3,434 / 18,388 / 11,277, and 806 anchors, matching
the 806 `stereo_exact` rows in `opentecr_comparison.csv`. There is no rebuild
regression and nothing to re-measure. M12's own arithmetic confirms the basis —
30,157 directed + 25,855 undirected = 56,012, the all-records total.

Every number in this document is now on the manuscript's basis. Both analysis
scripts default to it and take `--live` for the other.

### What *is* wrong: "55% growth in reactions"

`M09` opens: *"The database has grown 34% in compounds and 55% in reactions
since 2020."* Against the 2020 baseline the figure scripts already use
(`fd6c7849`, 2020-11-10):

| basis | compounds | reactions |
|---|---:|---:|
| **all records — the paper's basis elsewhere** | 33,992 → 45,708 = **+34.5%** | 43,774 → 56,012 = **+28.0%** |
| live only (obsolete dropped both sides) | 33,958 → 45,662 = +34.5% | 36,197 → 48,403 = +33.7% |
| 2026 all records vs 2020 live | — | 36,197 → 56,012 = **+54.7%** |

The compound figure is the first row and is correct. The 55% is the third row:
2026's totals set against a 2020 baseline with its 7,577 obsolete rows removed.
It is not reproducible on any single consistent basis.

**Fix: quote 28%.** That keeps the whole paper on one population — the
abstract's "~56,000 reactions", Figure 2A's 56,002, M12's 25,855 and
Supplementary Table S1 are all all-records figures, and Figure 1's source
counts (which apply no obsolete filter — verified in `_reaction_sources()`)
become consistent rather than wrong. The live-only basis is equally defensible
but would require converting every other number in the paper, including the
abstract's headline. Recommend 28%; this is the corresponding author's call and
is raised as such in the PR.

## R1 #1 — pathway/subsystem distribution of the new biochemistry

> *"It would be interesting to show a frequency distribution at pathway/subsystem
> level. Assuming that central pathways are already well characterized for many
> years, it would be curious to see where are we still gaining new information."*

**Origin.** `M02_introduction.tex` ("This update reports the growth of the
database") and `M09_results_growth.tex`, which reports growth only by source
database and by structural completeness — never by biology.

**This is a fair and answerable request, and the answer is a good one for us.**
It is also currently blocked by a data gap that we should fix rather than
explain away.

### The blocker

Of the 12,261 reactions added since 2020 (`rxn48576`–`rxn60859`), **not one
carries a pathway or an EC annotation** — neither in the per-reaction record
nor in the released alias files:

| | pre-2020 reactions (43,751) | new since 2020 (12,261) |
|---|---:|---:|
| `pathways` populated | 16,938 (38.7%) | **0 (0.0%)** |
| `ec_numbers` populated | 20,019 (45.8%) | **0 (0.0%)** |

`Unique_ModelSEED_Reaction_Pathways.txt` stops at `rxn48568` and
`Unique_ModelSEED_Reaction_ECs.txt` at `rxn48573`, both immediately below the
2020 boundary at `rxn48575`. `Unique_ModelSEED_Reaction_Aliases.txt` and
`..._Names.txt` do run to `rxn60859`, so this is the pathway/EC annotation step
never having been re-run over the MetaCyc and Rhea intake — not a source-data
problem. A user asking "what pathway is this new reaction in?" gets nothing
today, which is worth fixing on its own merits.

### It is reconstructable, and the answer is clear

Joining ModelSEED → MetaCyc alias → `Scripts/Provenance/MetaCyc/MetaCyc_pathways.tsv`
reaches 2,282 of the 4,206 new MetaCyc-sourced reactions without any new
downloads, and 16,499 of the carried-over reactions through the same join, so
the two columns are comparable. Distinct reactions, not reaction–pathway pairs,
all records:

| MetaCyc class | new | % of new | carried over | % of old | enrichment |
|---|---:|---:|---:|---:|---:|
| Antibiotic-Biosynthesis | 265 | 11.6% | 824 | 5.0% | 2.3× |
| O-Antigen-Biosynthesis | 218 | 9.6% | 29 | 0.2% | **54×** |
| POLYKETIDE-SYN | 169 | 7.4% | 224 | 1.4% | 5.5× |
| Toxin-Biosynthesis | 142 | 6.2% | 223 | 1.4% | 4.6× |
| Lipid-Biosynthesis | 140 | 6.1% | 438 | 2.7% | 2.3× |
| Branched-Fatty-Acids-Biosynthesis | 98 | 4.3% | 22 | 0.1% | **32×** |
| Fatty-acid-biosynthesis | 77 | 3.4% | 637 | 3.9% | 0.9× |
| ALKALOIDS-SYN | 53 | 2.3% | 223 | 1.4% | 1.7× |
| Sterol-Biosynthesis | 50 | 2.2% | 249 | 1.5% | 1.5× |
| Lipid-IV-A-Biosynthesis | 37 | 1.6% | 9 | 0.1% | **30×** |
| *Super-Pathways* | *299* | *13.1%* | *5,706* | *34.6%* | *0.38×* |

**The reviewer's hypothesis is right, and we can quantify it.** Only **41 of
the 12,261 new reactions (0.3%)** fall in a central-carbon or energy-metabolism
class (Energy-Metabolism, Electron-Transfer, Fermentation, TCA-VARIANTS,
Glycolysis, Pentose-Phosphate-Cycle, Photosynthesis, Respiration,
Methanogenesis), against 770 of 43,751 (1.8%) among the reactions already
there. The gain is in specialised metabolism — antibiotics, polyketides,
toxins, alkaloids — and in cell-envelope and lipid biosynthesis, where
O-antigen, lipid IV-A and branched fatty acids are enriched by one to two
orders of magnitude. That is a more interesting growth story than "MetaCyc got
bigger", and it is the story the Discussion already wants when it talks about
serving pathway design.

Named individual pathways make the point even more concretely — the fifteen
that gained most are mycolate biosynthesis (57 reactions), arachidonate
metabolites (37), colibactin (34), the iso- and anteiso-branched-chain fatty
acids (32 each), icosapentaenoate and docosahexaenoate metabolites, apratoxin A,
platensimycin, acylsucrose, and *H. pylori* O-antigen. Mycobacterial cell wall,
a bacterial genotoxin, eicosanoid signalling and two antibiotics: none of this
was reachable in 2020, and none of it is central metabolism.

**Two method points the panel must get right**, both of which caught this
analysis first:

- `MetaCyc_pathways.tsv` mixes individual pathways and ontology classes in one
  id column — `Degradation`, `Biosynthesis` and `Energy-Metabolism` all appear
  as ids. Treating every id as a pathway finds zero classes.
- The two columns must be walked to the **same depth**. The pre-2020 side in
  `Unique_ModelSEED_Reaction_Pathways.txt` carries a full ancestor closure,
  while a direct-`parent` join gives one level; comparing them would be
  meaningless. The table above puts both sides through the direct-parent join.
  `Super-Pathways` is an organisational class rather than a biological one and
  should be dropped from the panel (italicised above); `Lipid-Biosynthesis`,
  `Fatty-acid-biosynthesis` and `Branched-Fatty-Acids-Biosynthesis` are three
  depths of one lineage and should be collapsed or shown nested.

### Proposed change

1. **Re-run the pathway and EC annotation over the post-2020 intake** and ship
   it. This is the substantive fix; the figure is a by-product. For the 8,409
   Rhea-sourced reactions, Rhea carries EC cross-references directly and ChEBI
   carries no pathway concept, so expect EC coverage to be good and MetaCyc-class
   coverage to stay poor there — report that gap rather than hiding it.
2. **Add one panel** showing new-vs-carried-over reaction counts per MetaCyc
   class, sorted by the new count. Cheapest home is Figure 1 (see comment 4 —
   this replaces a bar panel rather than adding a fourth figure, which also
   answers the "barplot abuse" complaint).
3. **Two or three sentences in `M09`** after the completeness sentence. Draft:

   > The new biochemistry is not distributed like the old. Mapping reactions
   > onto MetaCyc's pathway ontology (Figure 1D), the classes that grew are
   > those of specialised metabolism — antibiotic and polyketide biosynthesis,
   > O-antigen and cell-envelope assembly, branched and acylated lipids,
   > alkaloids — while central carbon and energy metabolism account for 41 of
   > the 12,261 reactions added. Central metabolism was already saturated in
   > 2020; what a 2026 reconstruction gains is the periphery.

**What to run.** `grep -rl` finds the pathway and EC alias files written by
`Scripts/Biochemistry/Reset_Biochemistry_in_Git.sh` and, historically, by
`Scripts/Archived_Perl_Scripts/Compile_External_Pathways.pl` and
`Find_Unique_ModelSEED_Reaction_ECs.pl` — both archived Perl, which is
consistent with the annotation step having been dropped when the pipeline moved
to Python. There is no live Python equivalent, so this is a small port, not a
re-run. `Scripts/Provenance/{MetaCyc,KEGG}/Refactor_*_Pathway_Table.py` are the
modern parsers for the source tables and are the natural place to hang it.

**Author decision.** Whether to hold the revision for the annotation rebuild,
or submit the figure built from the MetaCyc join alone (2,282 reactions) and
state the coverage limit. Recommend the rebuild — the gap is a real defect a
later user would hit, and it is a port of two archived scripts rather than new
work.

**Budget.** +1 panel, +3 sentences. Offset: `M09`'s per-source reaction
sentence is long and partly duplicates Figure 1B; it can lose a clause.

---

## R1 #2 — transport reactions, compartments, and the limits of stoichiometry-only scoring

> *"Transport reactions do not have compartment information (they are only
> distinguished as compartments 0 and 1, unlike Rhea which identifies them as in
> and out). And the authors mention that transport reactions are scored from
> stoichiometry alone. I think it would be important to elaborate a bit more on
> this and discuss potential limitations (after all, thermodynamics are the most
> important for energy metabolism, which usually involves transport reactions)."*

**Origin.** Last sentence of `M05_methods_thermodynamics.tex`:

> *"A key caveat: as we publish a generalized database that can be applied across
> the kingdoms of life, transport reactions are scored from stoichiometry and
> energy alone, and neither an estimate of membrane potential or pH gradients
> are used."*

One sentence, at the end of a paragraph about reporting conditions, with no
follow-up anywhere in Results or Discussion. The reviewer is right that this is
under-treated.

**The compartment half of the comment is a documentation failure, not a data
failure.** ModelSEED indices are *relative* — 0 is intracellular and 1 is
extracellular relative to the reaction, resolved to named compartments
(`c0`, `e0`, and organelles) when a template instantiates the reaction into a
model. This is what lets one reaction serve a bacterial inner membrane, a
mitochondrial membrane and a plastid envelope. `Biochemistry/REACTIONS.md` calls
`m` a "compartment index number" and never explains this. It should, and the
paper should say it in one clause.

### What it costs — measured

6,313 transport reactions (11.3% of the database).

| | transport | non-transport |
|---|---:|---:|
| eQuilibrator commits to a direction | 24.4% | 17.4% |
| dGPredictor commits to a direction | **1.2%** | 17.7% |
| group contribution commits | 17.6% | 21.8% |
| graded at all | 76.5% | 56.9% |
| **gold / silver / bronze** | **27.7 / 14.6 / 57.7%** | 7.4 / 62.6 / 30.0% |

dGPredictor effectively refuses transport altogether — 1.2% against 17.7% — so
the transport grades rest on eQuilibrator and group contribution alone, and the
cross-source corroboration the grading scheme depends on is thinner here than
anywhere else in the database.

Transport reactions are graded gold nearly four times as often as the rest of
the database — which looks alarming until you see what the gold rests on. Of
the 1,338 gold transport reactions, **95.3% carry ATP and ADP**, and only 2.5%
translocate a proton at all. The deciding source is eQuilibrator for 1,231 of
them.

**That is the limitation, stated precisely:** for a primary active transporter
the estimate is dominated by ATP hydrolysis, which every predictor estimates
tightly, so the reaction earns a confident grade on the strength of its
chemistry while the translocation — the part that needs a membrane potential
and a pH gradient — contributes nothing to the score. The grade is trustworthy
as a statement about the hydrolysis and should not be read as a statement about
whether the transporter runs in that direction *in vivo*.

Worth saying in our favour: the scheme does **not** over-claim on the easy
cases. 4,056 transport reactions (64.2%) are uniport-like, with no net chemical
change, and **none of them is gold** — 2,647 bronze, 975 ungraded, and 434
silver. The 434 are worth understanding rather than glossing: with no net
chemistry every source returns ΔG ≈ 0, the sources therefore scatter very
little, and rule 3 (*corroborated*, $R \le 2$ and $z \le 2$) promotes them out
of bronze. They agree because there is nothing to disagree about. That is a
narrow, defensible artefact of the corroboration rule and it is better to name
it than to let a reader find it. The failure mode that matters is confined to
coupled transport.

One more gap to disclose: `S04` excludes symmetrical transport from the LLM
ensemble, so only 32.0% of transport reactions carry an LLM call against 81.5%
of the database. Figure 3A does not show this, and its caption implies uniform
coverage.

### Proposed change

1. **Expand the `M05` caveat** from one sentence to three, and move it out of
   the conditions paragraph into its own `\paragraph{Transport.}`. Draft:

   > \paragraph{Transport.} Compartments are recorded as relative indices —
   > 0 intracellular, 1 extracellular — rather than named compartments, so that
   > one reaction can be instantiated against a bacterial membrane, a
   > mitochondrion or a plastid when a template builds a model. Because the
   > database is not committed to a cell type, transport reactions are scored
   > from stoichiometry and energy alone: no membrane potential and no
   > transmembrane pH gradient enter the estimate. The consequence is specific.
   > A primary active transporter carries ATP hydrolysis in its stoichiometry,
   > every predictor estimates that tightly, and the reaction earns a confident
   > grade from chemistry that is not the part in question — 96% of the
   > gold-graded transport reactions are ATP-coupled, and only 2.5% translocate a
   > proton. Those grades should be read as confidence in the coupled chemistry,
   > not as evidence that the transporter runs in that direction in a given
   > organism. Reactions whose only change is a compartment change carry no net
   > chemistry, and none of them is graded gold: they are bronze or ungraded,
   > apart from a small number the corroboration rule promotes to silver because
   > every source agrees on an energy of approximately zero.

2. **One sentence in `M14_discussion.tex`** naming this as future work — a
   compartment-aware scoring layer needs organism context (Δψ, ΔpH,
   stoichiometry of the coupling ion) that belongs to a model, not to a
   reference database, and would be the natural next release.

3. **Flag transport in the released grades.** Cheap and it forecloses the
   objection: since `is_transport` already ships, either cap transport grades at
   silver, or add an `evidence_caveat: transport` key to `thermo-evidence`.
   Recommend the flag over the cap — capping destroys a real signal about the
   hydrolysis, whereas a flag lets the user decide.

4. **Fix `Biochemistry/REACTIONS.md`** to define the 0/1 convention as relative
   indices, and note that Rhea's in/out maps onto it directly.

**Author decision.** Item 3 changes shipped data, so it needs Chris/Sam. Items
1, 2 and 4 are text-only.

**Budget.** +5 sentences in Methods, +1 in Discussion. Offset: the magnesium
detail in `M05` duplicates `S02` and can be cut to a pointer.

---

## R1 #3 — why only ~9,000 reactions are used for reconstructions

> *"The authors mention that only 9,000 of the 56,000 reactions in the database
> are used for reconstructions. It would be great if they could elaborate on
> this. Is it because they are disconnected from the main core network? Is it
> due to lack of GPR associations?"*

**Origin.** First sentence of `M10_results_structure_curation.tex`:

> *"We built and ran our curation pipeline on all structures that would impact
> the set of roughly 9,000 reactions used in the ModelSEED and PlantSEED
> reconstruction templates."*

**The reviewer has misread this, and the sentence invites it.** That clause is
scoping the *structure-curation effort* — we prioritised curation on the
reactions the templates touch — not asserting that the other 47,000 reactions
are unusable. Nothing else in the paper says how much of the database reaches a
model. Given two reviewers read the paper and one drew this conclusion, the
fix is to rewrite the sentence *and* add the explanation, not to correct the
reviewer.

### What is actually true

The v7.0 templates (`ModelSEEDTemplates`, Gram-negative, Gram-positive,
Archaea, Core) hold **8,597 unique base reactions** in union; 8,594 are live in
the database. With PlantSEED_v3 roles that is the ~9,000 in the sentence.

Both of the reviewer's hypotheses are worth answering directly, and neither is
the main reason:

- **Not GPR.** The Gram-negative template types its 8,584 reactions explicitly:
  6,284 `conditional`, **2,258 `gapfilling`**, 31 `spontaneous`, 11 `universal`.
  A quarter of the template is therefore present with no gene requirement at
  all, as gapfilling candidates. (Complex references are sparser still — only
  2,590 of the 6,284 `conditional` reactions carry one in the shipped JSON, the
  rest presumably resolving their complexes through the role mapping at build
  time; the `type` field is the authoritative classification and is what the
  draft below cites.) A template is not a set of reactions with genes attached,
  so missing GPRs cannot be what keeps reactions out.
- **Not connectivity.** The templates are curated scopes, not the largest
  connected component; nothing prunes on graph connectivity.

The real reasons, in order:

1. **A template is a curated prokaryotic scope, not a usability verdict.** The
   database is a reference namespace spanning all kingdoms plus plants, plus
   reactions from secondary databases and published models (8,107 reactions
   carry only secondary-database or published-model aliases). Any one template
   deliberately covers a fraction.
2. **Mass and charge balance.** Only 28,096 of 48,403 live reactions (58%) are
   balanced; an unbalanced reaction cannot enter a template. This is the same
   structure-coverage limit the paper already reports, in a different guise.
3. **Functional-role association.** Gene-driven reconstruction needs a mapping
   to an annotated role. 20,019 live reactions (41.4%) carry an EC number.

And the headroom is real and quantifiable: **14,041 live reactions are both
balanced and EC-annotated, of which 8,089 are not in any v7.0 template.** That
is roughly a doubling of template scope available without any new curation —
a better answer than a defensive one.

### Proposed change

1. **Rewrite the `M10` opening sentence** so the scope claim cannot be read as
   a coverage claim. Draft:

   > We prioritised the curation pipeline on the compounds whose structures
   > affect the ~9,000 reactions carried by the ModelSEED and PlantSEED
   > reconstruction templates, since an error there propagates into every model
   > built from them.

2. **Add a short paragraph to `M13` or the Discussion** answering the question
   the reviewer actually asked. Draft:

   > A reconstruction template carries far fewer reactions than the database
   > holds — the v7.0 ModelSEED templates span 8,597 — and the gap is not a
   > statement that the remainder is unusable. The database is a reference
   > namespace across all kingdoms and includes biochemistry from secondary
   > databases and published models that no single template targets. What
   > limits promotion into a template is mass and charge balance, which 58% of
   > reactions currently achieve, and an association to an annotated functional
   > role, which 41% carry; gene associations are not the constraint, since a
   > quarter of the reactions already in a template are typed as gapfilling and
   > carry no gene requirement at all. 14,041 reactions
   > now satisfy both balance and role annotation, of which 8,089 sit outside
   > every current template — the structure curation reported here is what makes
   > that pool available.

**Author decision.** Confirm the mechanism by which reactions enter a template
(curated role mapping vs gapfill promotion) — that is a ModelSEEDTemplates
process question I cannot settle from this repository, and the draft above is
written to stay on the safe side of it.

**Budget.** +1 paragraph. Offset: `M13` atom-mapping is two sentences and could
absorb this without a new subsection.

---

## R1 #4 — the figures

> *"The figures could use a little improvement. There is a bit of barplot abuse
> (smiley face). For instance, Fig 1B would be more informative as a Venn
> diagram. Fig 2A needs a color legend (and maybe a different color scheme)."*

**Origin.** All three figures. The count is nine bar panels and one histogram
panel across three figures — the complaint is fair on its face.

| | current | verdict |
|---|---|---|
| 1A | grouped barh, 2020 vs 2026 per source | keep — magnitude over time, correct form |
| 1B | stacked barh, unique vs shared | **replace** — see below |
| 1C | stacked barh, complete vs incomplete | keep, or merge into 1A as a second encoding |
| 2A | stacked barh, pKa provenance | **needs a legend** — see below |
| 2B | plain barh, coverage per source | **weakest panel** — three numbers |
| 2C | three histograms | keep — this one earns its space |
| 3A | stacked barh, direction per source | keep — small-multiple comparison, correct form |
| 3B/3C | stacked barh, grade composition | keep, but merge is possible |

### 1B → Euler diagram

The reviewer is right that the stacked bar throws away the interesting part:
it says 67% of MetaCyc's reactions are unique but never says *who the other
33% are shared with*. The regions exist and are well separated (live reactions,
obsolete dropped):

| region | reactions |
|---|---:|
| MetaCyc only | 18,896 |
| Rhea only | 8,210 |
| KEGG only | 5,294 |
| all three | 3,465 |
| MetaCyc ∩ Rhea | 2,267 |
| MetaCyc ∩ KEGG | 1,515 |
| KEGG ∩ Rhea | 569 |
| in none of the three | 8,187 |

Note these differ from the counts in `M09` (MetaCyc 31,805/21,258 etc.), which
include obsolete reactions — see problem B above.

A three-set Euler reads well at these proportions. `matplotlib-venn` is **not
installed** in `~/Documents/py_venv`; either add it
(`uv pip install matplotlib-venn`) or draw three circles as matplotlib patches,
which is about fifteen lines and avoids a dependency for one panel. The 8,187
reactions in none of the three should sit as a labelled box outside the
circles, not be silently dropped.

### 2A needs a legend

The panel has four segments and labels two. Solid blue ("Marvin 29.5k") and
hatched blue ("6.9k") are direct-labelled; the narrow dark-grey segment and the
light-grey "8.2k" are not identified in the figure at all. The caption explains
them in prose — *"Marvin reported no dissociable proton ... for 1,242 compounds.
The remainder carry no structure to compute from"* — which is exactly the work
a legend should do. Add a four-entry key: **Marvin (query)**, **Marvin (from
SMILES)**, **no dissociable proton**, **no structure**.

On the colour scheme: the panel currently spends its only hue on the two Marvin
routes and renders both non-Marvin outcomes in grey, which visually merges a
chemistry result ("no dissociable proton") with a curation gap ("no structure").
Give the two grey segments distinct, non-sequential colours — they are not two
degrees of the same thing.

### Reduce the bar count

Two changes remove two bar panels and add the panel comment 1 needs, leaving
the figure count unchanged:

- **Drop 2B**, or fold it into 2C as a coverage annotation on each histogram.
  It carries three numbers (39%, 53%, 53%) and the text states all three.
- **Replace 1B with the Euler**, and **add the pathway-class panel** from
  comment 1 as 1D.

### Also: Figure 3 has no panel D

`M12_results_direction_sensitivity.tex` cites `Figure~\ref{fig:direction}D` for
the grade counts. Figure 3 has panels A, B and C only; the caption defines no D.
The reference should be to B and C, or the sentence should drop the pointer.
This is a live cross-reference error and needs fixing regardless of the redesign.

The same check across all three figures found this is the only broken pointer —
but also that **panels 3B and 3C are never cited by letter anywhere in the
text**, the only reference into Figure 3 being `A` and the non-existent `D`.
Since 3B and 3C are the two panels that show what the grading scheme produces,
retargeting the grading sentence at them fixes both problems at once.

**Budget.** Neutral to negative — one panel removed, one added, one replaced.

---

## R2 #2 — what the LLM ensemble is actually for

> *"The authors state in the Ensemble LLMs predictions section that a reaction may
> be driven by a pathway in a cell, which may be different from their grading of
> thermodynamic evidence. They then introduce LLMs and make predictions. Except
> the numbers in Figure 3, it is not apparent how LLMs play roles in
> reconstruction and other analyses."*

**Origin.** `M06_methods_council_direction.tex`, `\subsection{Ensemble LLMs
predictions}`, and `S04_direction_from_chemistry.tex`.

**The reviewer has found a real gap between what the section promises and what
it delivers.** The motivation opens on pathway-level driving force — *"every
reaction is a member of a pathway within a cell, and the overarching drive of
the pathway ... means a reaction may be driven in a direction on a level that
supersedes our evaluation"* — cites three papers on exactly that, and then
introduces an ensemble that sees only *"its name and stoichiometry alone"* and
therefore has no pathway context at all. The section then withdraws twice: the
calls *"do not [enter] our grading of the thermodynamic evidence"* (M06), and
their confidence *"does not appear to distinguish correct calls ... and is not
appropriate to use as a downstream filter"* (S04). A reader is left with a
method that is motivated by something it does not do and disclaimed for
everything else.

**The role we can defend is coverage — but the accuracy has to be stated with
it, and it is weaker than the headline agreement suggests.**

Coverage, which is the good news:

| | |
|---|---:|
| reactions carrying an LLM call | 45,644 (81.5%) |
| reactions **no** thermodynamic source directs | 25,855 |
| … of which the ensemble does direct | **22,902 (88.6%)** |

Accuracy. `Biochemistry/Thermodynamics/SourceGrading/opentecr_comparison.csv`
ships a measured `opentecr_dG_kJ` for 1,365 reactions, so the ensemble can be
scored against **measurement** rather than against another prediction. Taking
the measured direction from the sign of that energy:

| source | margin | n (both commit) | agree | rate | chance | kappa |
|---|---|---:|---:|---:|---:|---:|
| LLM ensemble | sign only | 361 | 246 | 68.1% | 52.8% | 0.33 |
| LLM ensemble | abs(dG) > 5 kJ | 264 | 200 | **75.8%** | 53.7% | 0.48 |
| LLM ensemble | abs(dG) > 11.7 kJ | 158 | 115 | **72.8%** | 48.0% | 0.48 |
| *eQuilibrator* | *sign only* | *127* | *126* | *99.2%* | *76.2%* | *0.97* |
| *eQuilibrator* | *abs(dG) > 11.7 kJ* | *82* | *81* | *98.8%* | *66.4%* | *0.96* |

*(Orientation check: the shipped eQuilibrator energy matches the sign of the
measured energy on 183 of 185 anchored reactions, which is what establishes
that `opentecr_dG_kJ` is already written in the ModelSEED equation's
orientation. Applying the `ms_orientation_vs_canonical` column as a flip drops
eQuilibrator to 45% and is the wrong reading — worth recording because it is an
easy mistake to repeat.)*

So the ensemble is right about three times in four where it commits and a
measurement exists, against eQuilibrator's ~99%. The 94.7% agreement with
eQuilibrator across the whole database (8,085 co-commitments, 7,659 agreeing)
reads far better than that, but its **chance baseline is 88.3%** and Cohen's
kappa is only **0.55**, because both sources overwhelmingly say "forward".

The mechanism is worth reporting because it is so clean. **85.3% of the
ensemble's calls are "forward as written"** (`>`), against 5.1% reverse and 7.8%
reversible, and the errors sit entirely on one side. At the 11.7 kJ margin:

| | LLM says `>` | LLM says `<` |
|---|---:|---:|
| **measured `>`** | 72 | 0 |
| **measured `<`** | 43 | 43 |

The ensemble is perfect when the answer is "forward" and a coin flip when it is
not. It has learned the convention that reaction equations are written in the
physiological direction — a genuine and useful prior — but it is not reading
thermodynamics. **This is the most important thing to disclose**: it is
recoverable from shipped data in an afternoon, so it is far better coming from
us, and it independently justifies the decision the paper already made to keep
these calls out of the evidence grading.

### Proposed change

1. **Rewrite the `M06` motivation** to stop promising pathway context. Draft
   replacement for the subsection's first two sentences:

   > Thermodynamic assignment is bounded by structural coverage: an energy
   > cannot be computed for a reaction whose participants lack complete
   > structures, however well understood the reaction is. That bound leaves
   > 23,729 reactions with no direction from any predictor. We therefore run a
   > complementary approach that does not depend on structures at all, using an
   > ensemble — a "council" — of large language models to interpret a reaction
   > from its name and stoichiometry.

2. **State the role explicitly**, replacing the bare *"We integrate these
   predictions transparently in our database, but do not include them in our
   grading"*:

   > The ensemble directs 22,902 of the reactions no thermodynamic source will
   > commit on, roughly doubling the fraction of the database that carries a
   > direction. Scored against the measured reactions where both commit, it
   > recovers the correct direction 76% of the time, against 99% for
   > eQuilibrator over the same anchors, and its errors fall almost entirely on
   > reactions that run in reverse: 85% of its calls are the direction the
   > equation is written in. We therefore release these calls as a separate,
   > clearly-labelled layer and exclude them from the evidence grading. They are
   > intended for reactions a reconstruction would otherwise have to leave
   > reversible by default, and should be reviewed before use rather than
   > trusted as evidence.

3. **Add the forward-bias caveat to `S04`**, in the results subsection: the
   85.3 / 5.1 / 7.8 call split, the 72/0 vs 43/43 confusion matrix, and the
   point that agreement with eQuilibrator (94.7%) sits against an 88.3% chance
   baseline, kappa 0.55. `S04` currently says only that the ensemble's
   *confidence* does not discriminate; the sharper and more useful finding is
   that its *direction* does not discriminate either, once the forward prior is
   accounted for.

4. **Fix the motivation's citations.** `mavrovouniotis1993`, `xu2008` and
   `noor2014` support the pathway-driving-force argument being removed. If that
   argument goes, either drop them or keep one sentence acknowledging that
   pathway context can override a per-reaction thermodynamic call and that
   neither approach here models it — which is true, and closes the loop with
   R1's transport comment.

**Author decision — the one real fork.** Either (a) keep the ensemble as a
released dataset with a coverage role and a stated accuracy, as drafted above,
or (b) give it a measured role in reconstruction by testing it on something.
(b) needs work we have not done: take a draft model, compare the reactions left
reversible by default against the ensemble's calls, and report how many
directions change and whether the model's predictions improve. That is the
analysis that would make the section unarguable, and it is the one the *title*
used to promise before it was cut on 2026-09-15.

The accuracy measurement above shifts the recommendation toward (a), and
firmly. A 76%-accurate direction call whose errors concentrate on exactly the
reactions a curator most needs help with is not something to build a
reconstruction claim on; it is a coverage layer with a known bias, which is
what the paper already treats it as. Recommend (a), with the bias disclosed.

There is a third option worth naming, (c): **keep only the reverse and
reversible calls as the useful signal.** The ensemble's `>` calls are nearly
uninformative given the 85% prior, but a `<` call is the model going against
that prior, and those are 2,350 reactions a curator could review in a
tractable pass. If anyone wants the ensemble to do real work in this release,
that framing costs nothing extra to compute and is defensible.

**Budget.** Roughly neutral — the rewritten motivation is shorter than the one
it replaces; the role sentences add about three lines.

---

## Summary of proposed actions

| # | change | kind | blocked on |
|---|---|---|---|
| B | fix the 55% growth figure to 28% (one consistent population) | text | author sign-off |
| 1a | re-run pathway/EC annotation over post-2020 reactions | **data/pipeline** | author sign-off |
| 1b | new Figure 1 panel: MetaCyc class distribution, new vs old | figure | 1a |
| 1c | 3 sentences in M09 | text | 1a |
| 2a | expand M05 transport caveat to its own paragraph | text | — |
| 2b | Discussion sentence on compartment-aware scoring | text | — |
| 2c | flag transport reactions in released grades | **data** | Chris/Sam |
| 2d | document the 0/1 convention in REACTIONS.md | docs | — |
| 3a | rewrite M10 opening sentence | text | — |
| 3b | paragraph on template scope vs database scope | text | template process Q |
| 4a | Fig 1B → Euler (+`matplotlib-venn` or hand-drawn) | figure | — |
| 4b | Fig 2A legend + recoloured non-Marvin segments | figure | — |
| 4c | drop or fold Fig 2B | figure | — |
| 4d | **fix `\ref{fig:direction}D` — panel D does not exist** | text | — |
| 5a | rewrite M06 motivation, state coverage role + 76% accuracy | text | — |
| 5b | forward-bias caveat + measured accuracy in S04 | text | — |
| 5c | resolve the three orphaned citations | text | — |

Net page budget: roughly neutral. The additions are concentrated in Methods
(transport) and Results (pathways, templates); the offsets available are the
magnesium detail in `M05`, the per-source reaction sentence in `M09`, and
Figure 2B.
