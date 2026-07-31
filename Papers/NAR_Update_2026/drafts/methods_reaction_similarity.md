# Methods — Reaction similarity (draft)

**Target section:** Materials and Methods, "Reaction similarity" (see `PAPER_2026_SKELETON.md`).
**Guide reference:** `PAPER_2026_GUIDE.md` §7.
**Status:** first pass.

---

Reaction similarity — a per-pair distance or similarity score across all mass-balanced ModelSEED reactions — is a new capability in the 2026 release. It supports two use cases: as a diagnostic for the reconciliation ontology described in the 2020 paper (reactions with high similarity that are not already ontology-linked are candidates for review), and as an input feature for downstream applications including annotation transfer, gap-filling, and pathway-inference.

**Pipeline.** Every ModelSEED reaction is represented by the concatenation of its balanced reactant and product structures. Structures are passed through an external reaction-embedding foundation model to produce a fixed-dimensional vector; the similarity between two reactions is the cosine similarity between their embeddings. GPU-accelerated batching allows the full pairwise similarity matrix to be regenerated in under approximately ten minutes on a single GPU, so the matrix is treated as a derived artifact and refreshed whenever compound curation materially changes the input structures.

**External reaction-embedding foundation-model selection.** Several external reaction-embedding models are available at the time of this release. Rather than adopting one by default, we evaluate the candidate models on their discerning power over ModelSEED reactions and defend the chosen model on that basis. The evaluation compares candidate models on: (i) intra-EC clustering — reactions sharing the same EC number should embed close together; (ii) intra-pathway clustering — reactions within a KEGG or MetaCyc pathway should embed close together; (iii) retrieval quality on a held-out set of known reaction-pair relationships (e.g. adjacent reactions in a pathway). Candidate models compared: [TBD list]. Chosen model: [TBD] on the basis of [TBD metric result].

**Publication of the matrix.** The similarity matrix is exported as a compressed sparse artifact under `Biochemistry/Reaction_Similarity/`. Users can regenerate it locally by rerunning the pipeline script (`Scripts/Reaction_Similarity/regenerate.py`, [TBD]) against a specific ModelSEED release tag.

---

## Open loose ends flagged during drafting

- The candidate-model list and the evaluation results are the main gap here; they'll be filled once the foundation-model bake-off in §7 runs.
- The regeneration script does not yet exist under `Scripts/Reaction_Similarity/`; needs to be authored and committed before the Methods paragraph's pointer is real.
- The similarity matrix isn't yet published under `Biochemistry/Reaction_Similarity/`; also needs to land.
- Related tools not currently in the paper but worth mentioning in Discussion: the Laino et al. paper on reaction similarity (referenced in the 2026-07 working-session notes); the `horizyn` project (restrictive license — evaluate whether we can cite / benchmark against without licensing complications).
