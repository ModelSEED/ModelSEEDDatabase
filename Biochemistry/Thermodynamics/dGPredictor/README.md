# dGPredictor thermodynamics inputs

dGPredictor (Wang et al. 2021) estimates reaction ΔG by decomposing molecules
into radius-1 and radius-2 atom-centred fragments and fitting a BayesianRidge
over the fragment counts. The model here is **retrained**, not the Wang lab's
shipped one, and it now produces per-compound ΔfG as well as per-reaction ΔrG.

## Live inputs

| file | rows | feeds |
|---|---|---|
| `retrained_dG.json` | 29,617 | `Scripts/Thermodynamics/Update_Reaction_dGPredictor_Energies.py` |
| `retrained_dG_compounds.json` | 30,185 | `Scripts/Thermodynamics/Update_Compound_dGPredictor_Energies.py` |

Schema: `{id: {dG_mean, dG_uncer, coverage, frags_used, frags_oov}}`, kJ/mol,
converted to kcal on ingest by 4.184.

**Every** prediction is installed, each carrying its own uncertainty —
29,617 reactions and 30,185 compounds. There is no coverage floor; see
"Why there is no coverage floor" below.

Stored record shape:

| kind | record | note |
|---|---|---|
| reaction | `[dg, err, operator, coverage]` | kcal/mol; operator from the dGPredictor rule set |
| compound | `[dg, err, coverage]` | kcal/mol; no operator — a formation energy has no direction |

Coverage is appended last on reactions, so positional readers of `[0]`/`[1]`/`[2]`
are unaffected. **Filter on `err`, not on `coverage`** — see below for why.

## Why the model was retrained

The predecessor `dGPredictor-ModelSEED` record, added Aug 2026, was
systematically broken: median dG_uncer 88.6 kJ/mol overall and 259.5 for
quinones, 44.5% of predictions past 100 kJ/mol, worst error 8,531 kJ/mol.

The cause was **RDKit canonicalization drift**. The Wang lab shipped a fragment
vocabulary written by an older RDKit (`CC(C)=O`); current RDKit's
`MolFragmentToSmiles(canonical=True)` emits `CC(=O)C` for the same fragment.
Reindexing the vocabulary against the shipped strings silently dropped most
fragment counts, which starved BayesianRidge on near-empty input. It then
extrapolated from its prior, and the prior covariance is what those enormous
error bars were.

The fix is not a patch to the vocabulary but a guarantee: the vocabulary is
rebuilt from current-RDKit output, and training and prediction call the same
`count_substructures`, so the two can no longer disagree about what a fragment
is called.

Separately — and this is why the archived baseline below must be read with
care — the *original* KEGG-mediated integration carried a stale-variable defect
in its 2023 staging notebook: the KEGG-alias extraction loop never reset its
holding variable, so every KEGG-less reaction inherited the id of the nearest
preceding KEGG-bearing one. 17,271 of the 27,715 archived records (62%) hold an
energy computed for a KEGG reaction that is not among that reaction's aliases.
The retrained model is structurally immune — it is keyed directly by ModelSEED
reaction id, with no KEGG detour.

## Training set: de-duplicated TECRDB

Fitted to `TECRDB_dedup.csv` (4,456 rows) rather than the raw Zenodo
`TECRDB.tsv` (4,544), **the same table component-contribution's
`cc_params_dedup.npz` was fitted to**. The 88 removed rows are redundant copies
of 78 measurement groups that eQuilibrator's table counts 2–6 times; one copy
of each is kept, so no measurement is lost.

This matters because the paper compares the two methods head to head. Fitting
them to different row sets of the same measurements is not defensible.

Held-out check that de-duplication is right here and not just consistent —
10 folds, both arms predicting an identical test set, every copy of a tested
measurement withheld from both arms so the raw arm cannot win on leakage
(`crossval_training_sets.tsv`, 3,633 held-out measurements):

| | raw TECRDB | TECRDB_dedup |
|---|---|---|
| RMSE kJ/mol | 14.5646 | 14.5588 |
| median \|error\| | 1.7125 | **1.6952** |
| mean \|error\| | 4.4274 | **4.4174** |
| 90th pct \|error\| | **7.3310** | 7.4279 |

De-duplicated is closer on 53.3%, paired Wilcoxon p = 1.5e-07.

Read that honestly: **highly significant and small.** Median error improves by
about 1%, RMSE is flat to four figures, and the tail is marginally worse; the
significance comes from the improvement being consistent across thousands of
rows rather than from its size. Same shape and direction as eQuilibrator's own
result for component-contribution (p = 1.7e-05). The conclusion to draw is that
consistency with eQuilibrator costs nothing measurable, not that
de-duplication is a large win.

Training-set RMSE *rises*, 5.4 → 5.85 kJ/mol. That is expected and is not
evidence of harm: removing 88 rows the fit was getting for free must raise
training error whether or not prediction improves. It is the reason the
held-out comparison above exists.

Fit quality: R² 0.9985, RMSE 5.85 kJ/mol, MAE 2.89, over 3,861 training rows ×
1,582 active features.

## Why there is no coverage floor

`coverage` is the fraction of a reaction's or compound's r=2 fragments present
in the training vocabulary. An earlier version of this work installed only
`coverage >= 0.9` — 8,279 reactions, 2,410 compounds — on the reasoning that
below that, a prediction is dominated by BayesianRidge's prior (variance
1/λ ≈ 1836 per uncovered fragment count) and is really "we don't know" wearing
a wide error bar.

That was tested and does not hold.

### The error bars are honest at every coverage level

20,691 reactions carry both a dGPredictor and an eQuilibrator estimate. Taking
`z = (dGP − eQ) / sqrt(σ_dGP² + σ_eQ²)`, a calibrated predictor puts ~68% of
cases within `|z| < 1`:

| stratum | n | median \|deviation\| | median σ | median \|z\| | frac \|z\|<1 |
|---|---|---|---|---|---|
| coverage = 1.00 | 3,011 | 7.7 | 9.6 | 0.49 | 74.3% |
| 0.90 ≤ cov < 1.00 | 4,236 | 16.8 | 50.9 | 0.43 | 75.4% |
| 0.75 ≤ cov < 0.90 | 7,316 | 28.2 | 65.2 | 0.47 | 76.1% |
| 0.50 ≤ cov < 0.75 | 4,484 | 35.1 | 81.0 | 0.46 | 78.7% |
| coverage < 0.50 | 1,644 | 87.1 | 110.7 | 0.58 | 67.4% |

kJ/mol. The calibration holds all the way down — at coverage < 0.5 the model
says ±111 and is off by 87, which is the error bar doing its job. `|z|` medians
of 0.43–0.58 against a nominal 0.67 mean the bars are if anything mildly
conservative. They never understate error, which is the only failure that would
justify withholding the prediction.

Caveat: eQuilibrator is the reference here and is itself an estimate, not
ground truth. Its σ is folded into `z`, which is the right treatment, but a
systematic eQuilibrator error would read as dGPredictor error.

### Coverage is a poor proxy for accuracy

Spearman correlation against |deviation from eQuilibrator|:

| | ρ |
|---|---|
| coverage | −0.349 |
| **σ** | **+0.650** |

The two correlate only −0.498 with each other, so they are not interchangeable.
Splitting at each metric's own median:

| | n | median \|deviation\| |
|---|---|---|
| coverage ≥ 0.9 but **σ high** | 1,947 | **43.3** |
| coverage < 0.9 but **σ low** | 5,038 | **12.3** |

The 0.9 floor was admitting ~1,950 predictions that are badly wrong while
rejecting ~5,000 that are nearly as good as the fully-covered set. It was
filtering on the wrong axis.

### What to filter on instead

σ. As a reference point, thresholds measured on the same both-methods set:

| σ cut | kept | median \|deviation\| |
|---|---|---|
| ≤ 2 kcal/mol | 3,682 | 2.2 |
| ≤ 5 | 5,266 | 2.6 |
| ≤ 10 | 7,650 | 5.2 |
| ≤ 15 | 12,501 | 9.7 |
| ≤ 20 | 18,210 | 13.9 |
| (old coverage ≥ 0.9) | 8,279 | 11.7 |

σ ≤ 10 kcal/mol keeps slightly fewer records than the old floor at less than
half the deviation. Pick a cut to suit the application rather than inheriting
one.

Coverage is still recorded, because it explains *why* σ is large and is useful
provenance when auditing a specific prediction. It is not the filter.

### Consequences for direction calls

The dGPredictor rule set (`DGP_HEURISTICS`, the Noor 2012 reversibility index
at one sigma) reads those honest error bars and declines to call a direction it
cannot support: of 29,617 records, 20,723 are `=`, 6,870 `>` and 2,024 `<`.
That is the intended behaviour — a wide bar produces an ambiguous call rather
than a confident one, so shipping low-coverage predictions does not smuggle
false confidence into the direction field.

Reported uncertainty by stratum, reactions:

| stratum | n | median \|ΔG\| | median dG_uncer | err > 100 |
|---|---|---|---|---|
| coverage 100% | 3,390 | 28.1 | 9.22 | 3.2% |
| coverage ≥ 90% | 8,279 | 29.6 | 36.00 | 9.1% |
| coverage ≥ 75% | 17,774 | 37.5 | 57.21 | 17.6% |
| coverage < 50% | 4,165 | 94.6 | 109.61 | 56.0% |
| (previous broken model) | 31,924 | — | 88.6 | 44.5% |

`DGP_COVERAGE_FLOOR` in the environment reinstates a floor if one is ever
wanted; it defaults to 0.0.

## Archived

`dGPredictor-2020.tsv` — 27,715 rows, the KEGG-trained predictions that used to
live in the reaction JSONs. Columns: `rxn_id`, `kegg_rxn_id`,
`dG_mean_kcal_per_mol`, `dG_uncer_kcal_per_mol`, `operator`. This is the
paper's "dGPredictor-2020" baseline.

**It is a baseline for the method, not a set of values that were ever correct
per-reaction** — see the carry-over defect above, which affects 62% of these
rows. Any head-to-head comparison should either restrict to the 10,444
correctly-mapped rows or state the contamination explicitly.

Removed at the same time, as superseded staging: `json_files/` (61 per-shard
Wang-lab prediction files), `pickle_files/` (61 pickles no current script ever
read), and `modelseed_retrained_dG.json` (the broken retrained payload).

## Regenerating

The pipeline is a self-contained project at
`MSD_Structures/dgpredictor_retrain/`, deliberately **not** committed here —
it is the build system, not the artefact. Seven phases, ~4 min total:

```bash
cd MSD_Structures/dgpredictor_retrain
for i in 01 02 03 04 05 06 06b; do
  micromamba run -n equilibrator python3 scripts/${i}_*.py
done
MSD_ROOT=/path/to/ModelSEEDDatabase \
  micromamba run -n equilibrator python3 scripts/07_atomic_swap.py
```

`01` takes `--tecrdb` / `--raw` to select the training table. `08` reruns the
raw-vs-dedup cross-validation. The raw-TECRDB run is preserved whole under
`dgpredictor_retrain/run_raw_tecrdb/` as the comparison arm.

Re-running phases 01–06b then the two `Update_*_dGPredictor_Energies.py`
scripts is idempotent.

## Known limitations

1. **No Legendre transform.** Wang lab's `get_ddG0` pH/ionic-strength
   correction is not implemented. The training data spans mixed pH and ionic
   strength and the model averages over them. Strict pH-7 predictions would
   need the transform reimplemented against ModelSEED's own pKa data.
2. **Abstract compounds are unreachable.** 22,642 reactions are skipped for
   want of a decomposition — Donor/Acceptor, ACP, ferredoxins, hemoproteins
   have no structure to fragment. This is a limit of any fragment-based
   method, not of this model.
3. **The vocabulary is capped by the training set.** ~600 training compounds
   is what sets the 3,390 fully-covered reactions. Only new experimental ΔG
   data on ModelSEED-specific chemistry moves that number.
