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
- Results (`M13`) — new `Reconstruction scope` paragraph answering the reviewer's two hypotheses directly, with the corrected figures: v7.0 template union 8,597 reactions; 69% of the database satisfies mass/charge balance; 41% carry an EC-annotated functional role; 15,128 reactions satisfy both balance and role annotation; of those, 8,826 already lie outside every current template — i.e., roughly a doubling of usable scope is already available once a template is built to draw on it.
- Discussion — one sentence naming the roughly half of reactions still lacking an assigned thermodynamic direction as a standing limit, alongside the transport limitation above.

**Supplementary changes:** None for this comment.
