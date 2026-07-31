# Discussion + Data Availability + Introduction + Abstract (drafts)

**Target sections:** Introduction, Discussion, Data Availability, Abstract (see `PAPER_2026_SKELETON.md`).
**Guide reference:** `PAPER_2026_GUIDE.md` §15.
**Status:** first pass. Best written last; these drafts are placeholders to be revisited after the earlier sections are finalized.

---

## Introduction (draft)

The 2020 release of the ModelSEED biochemistry database (Seaver et al., [10.1093/nar/gkaa746](https://doi.org/10.1093/nar/gkaa746)) established the database as a curated integration of KEGG, MetaCyc, BiGG, MetaNetX, Rhea, BioCyc-family databases, and dozens of published metabolic models — a "Rosetta Stone" for cross-database mapping backed by structure-first compound reconciliation, per-reaction thermodynamics, and a community-contribution workflow on GitHub. Six years later, we present an update that extends the database along four axes and provides the first systematic empirical accounting of a question the 2020 release could not: how much does the choice of a reaction-direction source actually matter for downstream metabolic-model output.

Three developments in the field motivate specific parts of this update. **Atom mapping** has matured from a specialty capability into a routine expectation for any biochemistry database supporting isotope-tracing or degradation-prediction workflows; we present a first-pass integration of atom mappings via a collaboration with the [TBD] group. **Machine-learning-based thermodynamics** methods (dGpredictor and successors) now cover a substantial fraction of the ModelSEED compound set and enable a multi-source approach that goes beyond the single-primary-plus-fallback strategy of 2020. **Reaction-direction transparency** has emerged as a growing concern in the metabolic-modeling community; each choice a curator makes about direction propagates through model output in ways that are rarely quantified, and we provide the first systematic quantification here.

We describe [TBD 3-4 sentence summary of the paper's structure, dropping in section markers].

## Discussion (draft)

Six years of community-driven curation have grown the ModelSEED biochemistry database by [TBD %] in compounds and [TBD %] in reactions ([TBD delta of balanced-reaction share]). That growth has been driven mainly by contributions from MetaCyc and the BioCyc-family databases; KEGG's contribution is roughly flat over the interval. The balanced-reaction share of the database has [regressed / held / recovered — TBD after final numbers] over the interval; the structure-curation pipeline described in Methods (particularly the protein-carrier cofactor standardisation) accounts for most of the [regression-recovery / continued-improvement] we report.

**Broader lesson — protein-carrier cofactor standardisation as a general pattern.** The acyl-ACP work described in Methods is one instance of a general pattern: whenever a protein-covalent cofactor carries a reactive group between enzymes, some databases include the cofactor in the stored formula for the loaded state and others do not, and the inconsistency generates false-positive mass-imbalance flags in reactions handling the loaded and unloaded forms. We extend the same pattern to biotinyl-carrier (biotin-carboxyl-carrier protein) and lipoyl-carrier (H-protein of glycine cleavage, E2 subunits of α-ketoacid dehydrogenase complexes) classes; candidate extensions to covalent FMN/FAD, molybdopterin, and heme-c are surveyed and deferred to future work.

**Open problems for the next release.** Several curation frontiers remain open:
- **Fragment-aware pKa/pKb.** The current biochemistry object stores pKa/pKb values only for atoms in the first InChI fragment; extending this to multi-fragment compounds requires both a schema change and a corresponding update to the protonation and thermodynamics pipelines.
- **EC number refresh.** The current ExPASy-name-based mapping between reactions and EC numbers is bound to have false and missing associations; a curation pass beyond name matching would improve annotation quality.
- **Obsolescence-leakage audit.** The database marks obsolete compounds and reactions as such, but some may still be inadvertently referenced by non-obsolete reactions or models; a systematic audit is warranted.
- **Direction-field removal.** With thermoreversibility computed per reaction and direction captured in the template, the reaction schema's `direction` field is arguably redundant; a removal pass has been considered but is [done / planned — TBD] as of the 2026 release.

**Where atom mapping and multi-source thermodynamics take the community next.** [TBD 2-3 forward-looking sentences.]

## Data Availability (draft)

- **GitHub repository:** `ModelSEED/ModelSEEDDatabase`, tagged as `v2.0.0` for the 2026 release. All curation history, override files, exclusion registry, thermodynamics ledger, and script sources are in this repository.
- **PyPi package:** [TBD `modelseed-biochem==2.0.0`, https://pypi.org/project/modelseed-biochem/].
- **Solr REST endpoint:** retained from the 2020 release; documentation in the repository's `Solr/` folder.
- **Web interface:** the ModelSEED website ([TBD final URL, e.g. `https://modelseed.org/biochem`]) has been updated to include atom-mapping-aware reaction landing pages. KBase integration is retained.
- **Atom mapping endpoint:** [TBD delivery mechanism from the Nikoloski collaboration — dataset in `Biochemistry/AtomMappings/`, external API, or both].
- **Continuous integration:** GitHub Actions workflow at `.github/workflows/` runs on every PR to `dev`.
- **License:** Creative Commons Attribution, matching the 2020 release; data derived from KEGG and MetaCyc remains subject to those databases' licenses.

## Abstract (draft — target ~250 words)

[TBD, drafted last so it can distill the finalised results. Working outline:]

1. Six years since 2020; ModelSEED biochemistry remains the chemistry foundation for ModelSEED and KBase metabolic-model reconstruction.
2. Growth to [TBD] compounds, [TBD] reactions, [TBD %] balanced.
3. Structure curation via PubChem validation + RDKit canonicalisation + per-curator provenance + mass-balance exclusion; extended to protein-carrier cofactor standardisation (acyl-ACP phosphopantetheine, biotinyl-BCCP, lipoyl carrier).
4. Multi-source thermodynamics (eQuilibrator refresh, rebuilt group-contribution, dGpredictor retrained on ModelSEED, OpenTECR).
5. New capabilities: atom mapping (Nikoloski collaboration), reaction similarity, PyPi distribution.
6. Novel empirical finding: reaction-direction source choice materially affects [TBD] on the ModelSEED v2 draft-model corpus.

---

## Open loose ends flagged during drafting

- Every `[TBD]` marker needs filling. Discussion and Introduction are best drafted after all Results sections are finalized.
- Author-contribution paragraph, funding, acknowledgments not drafted here — separate section, composed near submission.
- "How the Nikoloski atom-mapping work is credited" is `PAPER_2026_PLAN.md` open decision #3.
