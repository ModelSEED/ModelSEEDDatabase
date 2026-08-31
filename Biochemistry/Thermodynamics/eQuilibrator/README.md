# eQuilibrator thermodynamics inputs

## Live inputs

| file | rows | feeds |
|---|---|---|
| `ModelSEED_Reaction_Energies.tsv` | 56,012 | `Scripts/Thermodynamics/Update_Reaction_eQuilibrator_Energies.py` |
| `ModelSEED_Compound_Energies.tsv` | 45,708 | `Scripts/Thermodynamics/Update_Compound_eQuilibrator_Energies.py` |

Both carry a `# cache=... params=... p_h=... ionic_strength=... p_mg=... T=...`
provenance line. **Do not strip it** — it is the only record of which compound
cache and which component-contribution parameter set produced the numbers.

Generated 2026-08 from:

- **cache** — eQuilibrator's Zenodo cache with every `seed:` accession
  repointed at the structure ModelSEED holds (Path A), then training `kegg:`
  accessions repointed so `kegg:X` and `seed:cpd#####` resolve to the same
  compound (Path B). 97.1% of the 28,075 compounds where ModelSEED has a
  structure resolve to that structure; the remaining 2.9% lacked pKas, so the
  prior mapping was kept and is flagged unverified.
- **parameters** — component-contribution retrained on that cache with the
  TECRDB de-duplicated (78 measurement groups eQuilibrator counted 2–6 times;
  88 redundant rows removed, improving cross-validated prediction at
  p = 1.7e-05).
- **conditions** — pH 7.0, ionic strength 0.25 M, pMg 3.0, 298.15 K.

A row exists for every ModelSEED reaction and compound. The `status` column
says why a row has no energy rather than omitting it — that column is the point
of the format.

### Status vocabulary, and the split that matters

`GibbsEnergyPredictor.standard_dg` has two distinct outcomes and only one is a
failure:

```python
mu, sigma_fin, sigma_inf, residual = get_reaction_prediction(reaction)
if residual:
    return Q_(0, "kJ/mol").plus_minus(RMSE_inf)      # mu DISCARDED
else:
    std = |sigma_fin| + RMSE_inf * |sigma_inf|
    return Q_(mu, "kJ/mol").plus_minus(max(RMSE_eps, std))
```

When `residual` is non-empty, some reactant cannot be expressed in the model's
basis at all, `mu` is thrown away, and the return is literally zero plus the
pH/ionic-strength transform. **There is no estimate.** Status `undecomposable`,
with the offending compounds named so curation has a work queue.

When `residual` is empty the prediction is real, even if σ is enormous. Status
`ok`, and the true σ is recorded.

An earlier version keyed on `err > 1e4`, which cannot tell these apart — both
return a vast σ — and so filed 3,146 reactions and 9,081 compounds as
"outside CC span" without distinguishing them. That label is retired.

**Note the input.** `ComponentContribution.standard_dg_prime` calls
`reaction.separate_stored_dg_prime(...)` first and passes the *residual*
reaction to the predictor, so compounds carrying a measured dG never reach it.
Asking `get_reaction_prediction` about the full reaction reports those measured
compounds as undecomposable and declines reactions the library would answer.
`tools/verify_regen.py` exists to catch exactly that class of mistake: it fails
on any reaction that was `ok` and no longer is.

### σ = RMSE_inf is a sentinel, not an error bar

`RMSE_inf` is 100,000 kJ/mol = 23,901 kcal/mol. A reaction with an
unconstrained direction comes back at or beyond that scale. Of the 3,146
reactions that gained a value in the split, **none** is usable:

| σ band, kcal/mol | reactions | compounds |
|---|---|---|
| usable (< 10) | 0 | 0 |
| weak (10–100) | 0 | 0 |
| very weak (100–1,000) | 0 | 0 |
| unconstrained (> 1,000) | 3,146 | 9,081 |

σ min 11,065, median 23,902 (= RMSE_inf), max 861,957 kcal/mol.

They are stored anyway, with their true σ, so the decline is explicit and
quantified rather than an absent key — a consumer can see the method was asked
and answered "I cannot constrain this". Any σ filter removes them instantly,
and the `EQ` direction cascade already returns `?` for them via
`eq_undecomposable_heuristic`. **Do not quote the raw record count as
coverage:** 24,999 reactions carry a value, 21,853 carry a usable one.

### A note on pMg 3.0

eQuilibrator ships two default sets: "standard" (pMg 10, effectively no
magnesium) and "physiological" (pMg 3.0, 1 mM free Mg²⁺).
`ComponentContribution.__init__` uses the physiological one. The pre-2026
ModelSEED scripts set pH 7.0 explicitly but never set pMg, so they silently ran
at **pH 7.0 with physiological pMg 3.0** — a mix of the two convention sets.
This regeneration reproduces that mix deliberately, so old and new numbers stay
comparable. It is not internally consistent and is worth deciding on its own
terms. See "Open review items" below.

## Archived

| file | what it is |
|---|---|
| `eQuilibrator-2020_{reactions,compounds}.tsv` | the values stored in ModelSEED *before* this regeneration (25,028 / 30,607 entries, kcal/mol). The paper's baseline. |
| `MetaNetX_{Reaction,Compound}_Energies.tbl` | inputs of the superseded MetaNetX-mediated retrieval. No longer feed the database, but still read by `Check_eQuilibrator_Energy_Errors.py`, `Scripts/Tests/test_eq_heuristics.py`, `Scripts/Tests/test_reversibility_index.py`, `Preview_Structure_Update.py` and `Build_Compound_Field_Provenance.py`. |
| `InChIKey_Fallback_Compound_Energies.tbl` | the InChIKey-ladder fallback extract. Superseded. |

The MetaNetX retrieval matched ModelSEED compounds to MetaNetX ids down a
progressively loosening InChIKey ladder — full key, then key minus protonation,
then connectivity alone — so a stored energy could rest on a different
stereoisomer or protonation state than the structure ModelSEED holds. Removing
that is what this regeneration is for.

## What changed, precisely

Counting real energies (not keys — the old compound script wrote the
`10000000` sentinel for structureless compounds, so 8,765 of its 30,607 entries
were never energies):

| | before | after (usable) | after (all stored) | lost | gained | net usable |
|---|---|---|---|---|---|---|
| reactions | 25,028 | 21,853 | 24,999 | 6,238 | 3,063 | −3,175 |
| compounds | 21,842 | 16,372 | 25,453 | 9,223 | 3,753 | −5,470 |

"All stored" adds the unconstrained rows described above. They are records, not
coverage.

The new ingests write no sentinels at all: an entry either carries an energy or
carries no `eQuilibrator` key.

Most of the reaction loss is deliberate. Of the ~6,200 reactions that had a
value and no longer do, ~4,000 were outside the component-contribution span and
~1,900 do not balance — the old retrieval reported a number for both without
checking the infinite-uncertainty term. Roughly 270 are a genuine cache gap.

## Open review items

**These are not settled. They are recorded here so they do not get lost.**

### 1. The component-contribution span

The largest single reason for having no energy, on both sides:

| | undecomposable | unconstrained | combined share |
|---|---|---|---|
| reactions | 7,547 | 3,146 | 19.1% |
| compounds | 10,903 | 9,081 | 43.7% |

(These are the two halves of what the old `outside CC span` label lumped
together: 7,547 + 3,146 = 10,693 and 10,903 + 9,081 = 19,984.)

"Outside the span" means the molecule (or some participant) cannot be expressed
in the reactant- plus group-contribution basis the model was fitted in, so
`standard_dg` short-circuits and returns no estimate. **Nobody has yet audited what those
10,903 undecomposable compounds actually are.** The questions worth answering:

- How much of it is one structural class? Tetrapyrroles, metal complexes and
  polyprenyl chains are the known group-decomposer failures — if the 43.7% is
  mostly those, it is a bounded chemistry problem rather than a broad one.
- How much is reachable by curation? A compound outside the span because its
  structure is wrong is fixable; one outside because the basis genuinely does
  not cover its chemistry is not.
- How much of it matters? Restrict the count to the priority scope (the ~9,125
  reactions / ~6,500 compounds in v7.0 ModelSEEDTemplates + PlantSEED) before
  deciding it is a big number.

### 2. Compounds absent from the cache

| | not in cache | share |
|---|---|---|
| reactions | 8,367 | 14.9% |
| compounds | 9,352 | 20.5% |

`missing_compounds.tsv` ranks the ~4,650 compounds that block a reaction by how
many reactions each blocks — that ranking is the work queue. Known already:
roughly 1,100 are proteins or protein-bound cofactors carrying no formula at
all. Those can never be fixed by structure curation and would need redox-couple
handling instead, so they should be separated out before anyone treats the rest
as a curation backlog.

### 3. Reactions that stopped balancing

`lost_reactions_review.tsv` — 15 reactions that balanced under eQuilibrator's
structures and do not under ModelSEED's, implicating 9 compounds. Five are
polyprenyl / repeat-unit chain-length differences (cpd11713 dolichyl-PP is C80
in ModelSEED, C25 in eQuilibrator) — a representation convention clash rather
than an error in either. Four are small molecules 1–2 atoms apart (cpd22805,
cpd22806, cpd25746, cpd25506), where at least one of the two databases is
simply wrong. The small-molecule four are the tractable ones.

### 4. Structures we serve confidently and may have wrong

`modelseed_structure_checklist.tsv` (in the eQuilibrator working repo, 2,456
rows) ranks ModelSEED structures by how much external support they have from
KEGG / ChEBI / MetaCyc / Rhea. Highest-yield subset: the 13 "majority
disagrees" and 17 "only KEGG agrees". Every wrong structure there is now a
confidently-served wrong energy — the regeneration raised the stakes on this,
because the energy is computed from our structure with no fuzzy matching to
soften a mistake.

### 5. The pMg convention

See the note above. The pipeline runs pH 7.0 against physiological pMg 3.0,
inherited rather than chosen. Deciding it deliberately means re-running the
generation, so it was deliberately not bundled with this pass — changing it
here would have confounded the old-vs-new comparison.
