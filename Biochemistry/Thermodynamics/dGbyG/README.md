# dGbyG thermodynamics inputs

dGbyG (Fan et al. 2025, *Cell Systems*,
[doi:10.1016/j.cels.2025.101393](https://doi.org/10.1016/j.cels.2025.101393);
code [f-wc/dGbyG](https://github.com/f-wc/dGbyG), MIT) is a graph neural network
that predicts standard transformed formation energies, ΔfG′°. It is an ensemble
of 100 message-passing networks trained on TECRDB-derived reaction and formation
energies; its uncertainty is the spread of the 100 heads. The shipped model
(`models/mpnn_A139_B23_E300_L2_v2`) is run here, unmodified, over the structures
ModelSEED holds, and installed **additively** as `thermodynamics['dGbyG']`, next
to (never replacing) group contribution, eQuilibrator and dGPredictor.

## Live inputs

| file | rows | carry an energy | feeds |
|---|---|---|---|
| `ModelSEED_Compound_Energies.tsv` | 45,708 | 30,422 | `Scripts/Thermodynamics/Update_Compound_dGbyG_Energies.py` |
| `ModelSEED_Reaction_Energies.tsv` | 56,006 | 27,066 | `Scripts/Thermodynamics/Update_Reaction_dGbyG_Energies.py` |

One row per compound / reaction record, live and obsolete, in the layout of
`../eQuilibrator/ModelSEED_*_Energies.tsv` (without the reaction table's trailing
`formula` column): a `status` column says why a record has no energy, and only
`ok` rows carry one. The first line records the dGbyG commit, the model dGbyG
actually loaded, its head count and the conditions, and must not be dropped:

```
# dgbyg=https://github.com/f-wc/dGbyG@2202606a88f09423bcbf9e0989b7e12d79f32cc0 model=mpnn_A139_B23_E300_L2_v2 heads=100 p_h=7.0 ionic_strength=0.25M p_mg=14 T=298.15K
```

Energies are in kJ/mol to 3 decimals (dGbyG's native unit) and in kcal/mol
**already rounded to the 2 decimals the database stores**, from the unrounded
kJ value, so the updaters' 2-decimal formatting is the identity.

Stored record shape:

| kind | record | note |
|---|---|---|
| reaction | `[dg, err, operator]` | kcal/mol; operator from the DGB rule set (below) |
| compound | `[dg, err]` | kcal/mol; no operator — a formation energy has no direction |

Installed: 27,066 reactions (22,713 live) and 30,422 compounds (30,381 live).

## Provenance and regeneration

Both tables are produced by `Scripts/Thermodynamics/Generate_dGbyG_Energies.py`
from `Biochemistry/compound_*.json` and `reaction_*.json`; the run that wrote
them used dGbyG at the commit in the header line, on CPU (24 workers, about 4
minutes). The dGbyG checkout and its weights are not redistributed here (see
`Scripts/Thermodynamics/DATA_DEPENDENCIES.md`). To regenerate, then install:

```bash
cd Scripts/Thermodynamics
./Generate_dGbyG_Energies.py --dgbyg /path/to/dGbyG --workers 24   # dGbyG's environment.yml dependencies
./Update_Compound_dGbyG_Energies.py
./Update_Reaction_dGbyG_Energies.py
```

Keep inference on CPU: GPU float32 arithmetic can move a value across a
2-decimal rounding boundary. The values are static — like dGPredictor's, they
are not refreshed by the structure-update cascade; rerun the generator after
structure changes.

## Conditions and convention

- **pH 7.0, ionic strength 0.25 M, pMg 14, 298.15 K** — dGbyG's native
  condition, the one its training data were transformed to. The eQuilibrator
  tables are at pMg 3.0. dGbyG's Legendre transform has no Mg-binding constants
  (it moves only explicit Mg atoms), so re-transforming would not model Mg
  binding anyway. Supplementary S2 of the 2026 update paper puts the pMg 14 → 3
  shift at a median 0.044 kcal/mol over a 400-reaction sample of Mg-binding
  reactions — a median, not a bound.
- **Convention B** (Alberty, transformed): H⁺ is 0 (cpd00067 is stored as
  `[0.0, 0.0]`), water −37.36 kcal/mol — the same convention as eQuilibrator
  and dGPredictor, not group contribution's Convention A (see `../README.md`).

## How each record is computed

**Compounds.** Each ModelSEED SMILES goes through dGbyG's own normalisation
(RDKit `Normalize` + `Uncharger`); the value is the mean over the 100 heads, the
uncertainty their standard deviation. Multi-fragment SMILES (salts,
counter-ions) are predicted as given: 467 installed compounds, used by 60
installed reactions (45 live).

**Reactions.** Stoichiometry is netted by compound across compartments — the
compartment-collapsed reaction eQuilibrator and dGPredictor also score — and H⁺
is dropped: its transformed formation energy is 0 at fixed pH, and dGbyG's
normalisation changes the protonation state of many species from ModelSEED's
pH-7 form. The value is the mean over heads of the **summed** per-head vector
and the uncertainty is the standard deviation of that summed vector — dGbyG's
own `Reaction.standard_dGr_prime`. Compound uncertainties are not added in
quadrature: the heads' errors are strongly correlated (ATP's compound σ is
5.7 kcal/mol, ATP hydrolysis's 0.8).

**Direction.** `Update_Reaction_dGbyG_Energies.py` computes the operator at
write time with the rule set registered for `'dGbyG'` in
`reversibility_heuristics.SOURCE_HEURISTIC_SET` (`DGB`): the Noor 2012
reversibility index at one sigma, built by the same factory and arguments as the
eQuilibrator and dGPredictor sets. Live calls: `=` 8,035, `>` 7,292, `?` 6,486,
`<` 900 — 71.4% resolved, 36.1% irreversible. Median σ 2.87 kcal/mol on live
records (5th–95th percentile 0.26–11.73, max 341.70; the 2,500 gate never fires).

## Why a record has no energy

Reactions, in the generator's order of precedence (live / all records):

| status | live | all | meaning |
|---|---:|---:|---|
| `no stoichiometry` | 1 | 10 | — |
| `translocation only` | 2,752 | 4,120 | every participant cancels across compartments |
| `missing structure` | 10,314 | 10,900 | a participant has no SMILES |
| `compound failed` | 71 | 89 | a participant could not be parsed or predicted |
| `R group` | 9,965 | 11,099 | a participant carries an attachment point (`*`); dGbyG would featurise the dummy atom and return a number, but it is an extrapolation outside anything it was trained on. (12,506 live reactions have an R-group participant; the earlier statuses claim the rest.) |
| `unbalanced` | 2,140 | 2,213 | fails dGbyG's balance check on dGbyG-normalised structures; 2,095 of the live ones are flagged `MI` by ModelSEED itself |
| `degenerate` | 428 | 509 | every head gives the same value, σ < 0.005 (see below) |
| `ok` | 22,713 | 27,066 | installed |

Compounds: `ok` 30,381 live (30,422 all); `no structure` 8,755 (8,756);
`R group` 6,484 (6,488); `not parsed` 42 (42).

**The balance check enforces heavy atoms and water, not electrons.** It is
dGbyG's own `Reaction.is_balanced` (`ignore_H`, `ignore_H_ion`): hydrogen is
ignored and a charge difference that matches a hydrogen difference is absorbed
as H⁺, so electron-unbalanced half-reactions pass. 1,646 live installed
reactions (1,855 all) carry ModelSEED's own `CI` flag; the eQuilibrator tables
also ship such reactions.

**`degenerate` reactions are withheld on purpose.** dGbyG is a 2-layer
message-passing network with no chirality or E/Z feature, so a reaction whose two
sides are the same graph up to stereochemistry (racemases, epimerases), or that
leaves every atom's 2-bond neighbourhood unchanged (NAD + NADPH ⇌ NADH + NADP,
nucleotide exchanges such as GTP:GTP guanylyltransferase), is predicted as
exactly 0 by every head. Shipping `[0.0, 0.0, "="]` would tell every consumer
that reads σ that the model is certain, when the zero is a blind spot.
eQuilibrator covers 413 of the 428 live ones and puts 91 beyond ±1 kcal/mol.

## What to expect from the numbers

**Measured on this ingest** (live reactions unless stated; reproducible from the
database):

- **Agreement with the other sources.** Where both commit to an irreversible
  call (`>` or `<`), dGbyG agrees with eQuilibrator 99.4% (n = 4,885), group
  contribution 99.1% (4,746) and dGPredictor 95.9% (4,287). These pairs are
  dominated by large energies from sources that share TECRDB training data, so
  they show consistency, not accuracy. Against eQuilibrator energies: Pearson r
  0.923, median |Δ| 2.11 kcal/mol (n = 17,269), with a wide tail (27.1% differ
  by more than 5).
- **Coverage it adds.** 140 reactions get their first usable energy from dGbyG
  and 2,227 go from one source to two. Of the 23,790 reactions no other
  predictor gives a call, dGbyG gives one (`>`, `<` or `=`) for 1,284, 707 of
  them irreversible.
- **Formation energies sit on dGbyG's own reference.** Against eQuilibrator
  (the 16,651 live compounds where both carry a usable value) they correlate at
  r 0.995 but differ by a median 8.0 kcal/mol, an offset hydrogen count and
  charge do not explain (R² 0.002) and that largely cancels in balanced
  reactions. Compare sources at the reaction level.

**From the comparison analysis behind this ingest** (run in the dGbyG
checkout, including a grouped 10-fold retraining of dGbyG; not reproducible
from this repository alone):

- **Against measurement** (the 348 non-degenerate of the 363 live openTECR
  anchors): MAE 1.67 kcal/mol, 69.8% within 2 (eQuilibrator 0.85 / 89.1%,
  dGPredictor 0.62 / 93.1%). All of these are in-sample — every anchor is a
  dGbyG training reaction, and eQuilibrator and dGPredictor are fitted to the
  same TECRDB measurements. Retrained with each anchor withheld, dGbyG's MAE is
  2.12.
- **Structure representation costs dGbyG accuracy.** ModelSEED stores many
  compounds in charged pH-7 forms that survive dGbyG's uncharger (it adds
  explicit hydrogens first, which blocks deprotonation of a positive N–H: NH₄⁺,
  amino-acid zwitterions, `[nH+]` adenine), and draws amides as imidates,
  `N=C([O-])`, which the uncharger turns into imidic acids, `N=C(O)`. Neither
  form occurs in dGbyG's training structures. On dGbyG's own KEGG structures the
  held-out anchor MAE is 1.55 instead of 2.12. NH₄⁺ alone carries a compound σ
  of 13.6 kcal/mol, inflating σ on ammonia reactions. The values here use
  ModelSEED's structures as they are.
- **The ensemble spread understates the error on familiar inputs.** On its own
  training reactions, held out by reaction-grouped cross-validation, dGbyG's
  error exceeds its σ by a median factor of 2.3 and 26% fall within 1σ (68% if
  calibrated; σ = 0 rows excluded) — about eQuilibrator's position in its own
  cross-validation and far from dGPredictor's, although those two fold over
  measurement rows, a less strict split. On ModelSEED's structures the stored σ
  is wider, inflated by the unfamiliar forms above: on the held-out anchors 58%
  fall within 1σ (median |z| 0.79).

## Known limitations of dGbyG as applied here

- No stereochemistry (see `degenerate`); a 2-bond receptive field.
- `Compound(smiles)` cannot fully uncharge ModelSEED's pH-7 SMILES (above).
- dGbyG's `Reaction.balance()` is meant to add water for an oxygen deficit but
  never does (`src/dGbyG/api.py` tests `is_balanced` on the unmodified
  reaction); 760 records that would balance with water are `unbalanced` here.
- Electron balance is not checked (above).
- Transport is compartment-collapsed with no membrane potential or pH gradient,
  as for every source.

## Not in the evidence grading

dGbyG is an additional per-source record only. It is not one of the sources
`Scripts/Thermodynamics/SourceGrading/` grades, so `thermo-evidence`, the
recommended source and canonical `reversibility` are unchanged. Adding it to the
grading takes more than an entry in `SOURCES`
(`optimize_thermo_source_assignment.py`): the `proxy`/`trusted` tables there and
in `grade_thermo_sources.py`, `K`/`LABEL`/`PRECEDENCE` in
`grade_thermo_sources.py` and `recommend_thermo_source.py`, `FULL` in
`build_grade_table.py`, and the `KEY`/`PRECEDENCE` maps in
`Apply_Evidence_Grades_And_Recommendation.py` and
`Promote_Graded_Direction_to_Canonical.py` — and it changes the shipped grades.
