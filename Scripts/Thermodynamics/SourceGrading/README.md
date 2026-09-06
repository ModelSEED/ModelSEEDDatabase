# Thermodynamic source grading and recommendation

Relocated here 2026-09-03 from Cooper Taylor's `core_models_analysis`
(branch `template-direction-heuristics`, `cd2b93e`). The method is documented
in full in that repository at `reports/thermoSourceMethod/THERMO_SOURCE_METHOD.md`.

## Why it belongs in this repository

It answers a database question, not an analysis question: *given several
disagreeing ΔG′° values for a reaction, which should a consumer use?* It also
already depended on this repository — `recommend_thermo_source.py` imported
`reversibility_heuristics` out of `MSDB_CODE/Scripts/Thermodynamics`. That is
now a sibling import.

## What moved, and what did not

Moved — stages 1 to 3, calibration through recommendation:

| file | role |
|---|---|
| `optimize_thermo_source_assignment.py` | `SOURCES`, `load_db`, error-model fitting. The base module the other two import. |
| `grade_thermo_sources.py` | per-(reaction × source) grade, gold / silver / bronze |
| `recommend_thermo_source.py` | picks one source per reaction |

Left behind — stages 4 and 5, direction maps and core-model FBA
(`build_graded_direction_maps.py`, `run_graded_fba_all_models.py`,
`analyze_graded_fba.py`, `plot_graded_fba.py`). Those need 5,683 metabolic
models and cobra/GLPK. That is analysis, and it should stay with the models.

## The headline result, so nobody re-derives it

Direction is **not** chosen by uncertainty. Held-out accuracy over 20 random
70/30 splits of the 802-reaction anchor:

| strategy | accuracy | coverage |
|---|---:|---:|
| **priority EQ > DG > GC** | **95.9% ± 1.1** | 100% |
| eQuilibrator only | 95.9% ± 1.2 | 98.3% |
| argmin calibrated τ | 93.7% ± 1.1 | 100% |
| argmin expected magnitude error | 93.7% ± 1.1 | 100% |
| dGPredictor-ModelSEED only | 91.6% ± 1.1 | 100% |

Every uncertainty-based selector loses to a fixed priority order, because
direction errors concentrate where magnitude error is *smallest* — at
ΔG′° ≈ 0, which is exactly where the cascade's ±2 band decides.

## Paths

All defaults now resolve inside this repository. Override by environment:

| variable | default |
|---|---|
| `MSDB_ROOT` | repository root (reads `Biochemistry/reaction_*.json`) |
| `MSDB_CODE` | repository root |
| `THERMO_GRADING_OUT` | `Biochemistry/Thermodynamics/SourceGrading/` |
| `TECRDB_COMPARISON` | `Biochemistry/Thermodynamics/SourceGrading/tecrdb_vs_dgpredictor_modelseed.csv` |

## MISSING INPUT — this will not run yet

`tecrdb_vs_dgpredictor_modelseed.csv` is **not in either repository**. Three
scripts read it and none produce it; it was generated in a third place
(`/scratch/ctaylor/dgpredictor_tecrdb/results/`) during the dGPredictor
retraining work. It is the ground truth — 1,550 TECRDB matches, 802 of which
form the anchor every accuracy number is measured against.

Ask Cooper for it, or regenerate it by matching TECRDB measurements to
reactions and joining the dGPredictor-ModelSEED predictions. Until then the
harness is here and repathed but not runnable.

`recommend_thermo_source.py` also reads `results/eq_vs_dgpms/` from the
analysis repository. Check whether that is a hard dependency before running.

## Adding a source

`SOURCES` in `optimize_thermo_source_assignment.py`:

```python
SOURCES = {"Group contribution": "GC", "eQuilibrator": "EQ",
           "dGPredictor-ModelSEED": "DGPMS", ...}
```

It reads `thermodynamics[label]` from each reaction record, so a new source is
one entry. This is the intended route for grading `eQuilibrator-Marvin` against
`eQuilibrator-MolGpKa` once the split lands — see the two-energy-set plan.

## Caveats worth carrying

From the report's own limitations, and they bear on how far the numbers travel:

* the 802-reaction anchor is central metabolism, so accuracies are conditioned
  on the well-measured part of the database
* eQuilibrator is *fitted* on TECRDB and the anchor *is* TECRDB, so its margin
  is partly in-sample; only Group Contribution is genuinely out-of-sample
* results are snapshot-specific. The report was built against `dev @ 49563c6f`
  on 2026-08-12. The 2026-09-03 pKa rebuild moved 42.4% of eQuilibrator
  reactions by more than 1 kcal/mol with 1,580 direction changes, so
  "EQ first at 98.9%" is verified for the Marvin-era energies only and needs
  re-running before it can be claimed for the MolGpKa-era set.
